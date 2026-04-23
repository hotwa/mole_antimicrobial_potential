"""Runtime planning and multiprocessing topology helpers for process screening."""

from __future__ import annotations

import asyncio
import argparse
import json
import multiprocessing
import os
import queue
import tarfile
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Mapping

import pandas as pd

from src.batch_screening import iter_frame_batches, screen_frame_sync
from src.preprocess_screening_input import preprocess_many_to_parquet
from src.screening_planner import PlannerConfig, ScreeningPlanner, WorkUnit
from src.screening_sources import (
    DEFAULT_ARCHIVE_CHEM_ID_COL,
    DEFAULT_ARCHIVE_SMILES_COL,
    DEFAULT_INPUT_CHUNK_SIZE,
    process_work_unit,
)
from src.service import create_scheduler


ScreenInputPath = str | Path


@dataclass(frozen=True)
class ProcessScreenConfig:
    input_paths: list[ScreenInputPath]
    output_dir: str
    execution_mode: str = "thread"
    producer_processes: int | str = "auto"
    predict_queue_max_batches: int | str = "auto"
    result_queue_max_batches: int | str = "auto"
    batch_checkpoint_size: int = 2048
    rows_per_shard: int = 100_000
    row_group_size: int = 4096


@dataclass(frozen=True)
class RuntimePlan:
    predictor_processes: int
    producer_processes: int
    work_queue_max_batches: int
    predict_queue_max_batches: int
    result_queue_max_batches: int
    batch_checkpoint_size: int
    graph_workers: int


@dataclass(frozen=True)
class ProcessWorkItem:
    work_unit: WorkUnit
    smiles_colname: str = "smiles"
    chem_id_colname: str = "chem_id"
    archive_smiles_colname: str = DEFAULT_ARCHIVE_SMILES_COL
    archive_chem_id_colname: str = DEFAULT_ARCHIVE_CHEM_ID_COL
    chunk_size: int = DEFAULT_INPUT_CHUNK_SIZE


@dataclass
class ProducerBatch:
    producer_id: int
    work_unit: WorkUnit
    batch_index: int
    frame: pd.DataFrame


@dataclass
class PredictedBatch:
    producer_id: int
    work_unit: WorkUnit
    batch_index: int
    input_frame: pd.DataFrame
    predicted_frame: pd.DataFrame


@dataclass(frozen=True)
class ProcessErrorRecord:
    stage: str
    worker_name: str
    message: str
    traceback_text: str


@dataclass
class ProcessTopology:
    context: Any
    work_queue: Any
    predict_queue: Any
    result_queue: Any
    error_queue: Any

    def close(self) -> None:
        for queue in (self.work_queue, self.predict_queue, self.result_queue, self.error_queue):
            cancel_join_thread = getattr(queue, "cancel_join_thread", None)
            if callable(cancel_join_thread):
                cancel_join_thread()
            close = getattr(queue, "close", None)
            if callable(close):
                close()


@dataclass(frozen=True)
class PreparedInput:
    input_path: str
    prepared_path: str
    source_group: str
    manifest: dict[str, Any]


def _normalize_input_paths(
    input_path: ScreenInputPath | list[ScreenInputPath] | tuple[ScreenInputPath, ...] | None,
) -> list[ScreenInputPath]:
    if input_path is None:
        return []
    if isinstance(input_path, (str, Path)):
        return [input_path]
    if isinstance(input_path, (list, tuple)):
        return list(input_path)
    return [input_path]


def _resolve_existing_input_paths(
    input_paths: list[ScreenInputPath] | tuple[ScreenInputPath, ...] | None,
) -> list[Path]:
    resolved_paths: list[Path] = []
    for raw_path in _normalize_input_paths(input_paths):
        path = Path(raw_path).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"screen input path does not exist: {path}")
        resolved_paths.append(path)
    return resolved_paths


def _is_supported_process_passthrough(path: Path) -> bool:
    if path.is_file() and path.suffix.lower() == ".parquet":
        return True
    if path.is_dir() and any(candidate.is_file() for candidate in path.rglob("*.parquet")):
        return True
    return False


def _is_supported_process_tabular(path: Path) -> bool:
    if not path.is_file():
        return False
    if path.suffix.lower() == ".parquet":
        return False
    if path.suffix.lower() in {".sqlite", ".sqlite3", ".db", ".db3"}:
        return False
    if tarfile.is_tarfile(path):
        return False
    return True


def _build_passthrough_manifest(input_path: Path) -> PreparedInput:
    if input_path.is_dir():
        source_group = input_path.name
    else:
        source_group = input_path.stem

    manifest = {
        "input_path": str(input_path),
        "prepared_path": str(input_path),
        "source_group": source_group,
        "mode": "passthrough",
    }
    return PreparedInput(
        input_path=str(input_path),
        prepared_path=str(input_path),
        source_group=source_group,
        manifest=manifest,
    )


def _write_prepared_manifests(
    *,
    output_dir: Path,
    prepared_inputs: list[PreparedInput],
) -> tuple[Path, list[str]]:
    manifest_path = output_dir / "prepared_manifest.json"
    manifest_dir = output_dir / "prepared_manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "prepared_root": str(output_dir / "prepared"),
        "source_groups": sorted(item.source_group for item in prepared_inputs),
        "inputs": [item.manifest for item in prepared_inputs],
    }
    _write_json_atomic(manifest_path, payload)

    per_input_paths: list[str] = []
    for item in prepared_inputs:
        item_path = manifest_dir / f"{_sanitize_path_segment(item.source_group)}.json"
        _write_json_atomic(item_path, item.manifest)
        per_input_paths.append(str(item_path))

    return manifest_path, per_input_paths


def prepare_process_inputs(
    *,
    input_paths: list[ScreenInputPath] | tuple[ScreenInputPath, ...],
    output_dir: str | Path,
    smiles_colname: str,
    chem_id_colname: str,
    rows_per_shard: int,
    row_group_size: int,
) -> tuple[list[PreparedInput], Path, list[str]]:
    resolved_inputs = _resolve_existing_input_paths(input_paths)
    prepared_root = Path(output_dir).expanduser().resolve() / "prepared"

    tabular_inputs: list[Path] = []
    passthrough_inputs: list[PreparedInput] = []
    for input_path in resolved_inputs:
        if _is_supported_process_passthrough(input_path):
            passthrough_inputs.append(_build_passthrough_manifest(input_path))
            continue
        if _is_supported_process_tabular(input_path):
            tabular_inputs.append(input_path)
            continue
        raise ValueError(
            "process mode currently supports CSV/TSV-style tabular inputs or prepared Parquet files/directories"
        )

    prepared_inputs: list[PreparedInput] = []
    if tabular_inputs:
        source_groups = [path.stem for path in tabular_inputs]
        if len(source_groups) != len(set(source_groups)):
            raise ValueError(
                "process mode requires unique input basenames for tabular inputs because source_group defaults to Path.stem"
            )

        manifest = preprocess_many_to_parquet(
            input_paths=tabular_inputs,
            output_dir=prepared_root,
            smiles_colname=smiles_colname,
            chem_id_colname=chem_id_colname,
            rows_per_shard=rows_per_shard,
            row_group_size=row_group_size,
        )
        for item in manifest["inputs"]:
            source_group = str(item["source_group"])
            prepared_path = prepared_root / source_group
            item_manifest = dict(item)
            item_manifest["prepared_path"] = str(prepared_path)
            item_manifest["mode"] = "preprocessed"
            prepared_inputs.append(
                PreparedInput(
                    input_path=str(item["input_path"]),
                    prepared_path=str(prepared_path),
                    source_group=source_group,
                    manifest=item_manifest,
                )
            )

    prepared_inputs.extend(passthrough_inputs)
    manifest_path, per_input_manifest_paths = _write_prepared_manifests(
        output_dir=Path(output_dir).expanduser().resolve(),
        prepared_inputs=prepared_inputs,
    )
    return prepared_inputs, manifest_path, per_input_manifest_paths


def build_process_work_items(
    *,
    prepared_inputs: list[PreparedInput],
    grouping_mode: str,
    cpu_workers: int | str,
    target_rows_per_group: int | str,
    target_bytes_per_group: int | str,
    smiles_colname: str,
    chem_id_colname: str,
    archive_smiles_colname: str,
    archive_chem_id_colname: str,
    chunk_size: int,
) -> list[ProcessWorkItem]:
    planner = ScreeningPlanner(
        PlannerConfig(
            grouping_mode=grouping_mode,
            cpu_workers=cpu_workers,
            target_rows_per_group=target_rows_per_group,
            target_bytes_per_group=target_bytes_per_group,
            chunk_size=chunk_size,
        )
    )

    work_items: list[ProcessWorkItem] = []
    for prepared in prepared_inputs:
        units = planner.plan(
            prepared.prepared_path,
            archive_pattern="*",
            sqlite_table=None,
            sqlite_query=None,
            smiles_colname="smiles",
        )
        for unit in units:
            work_items.append(
                ProcessWorkItem(
                    work_unit=unit,
                    smiles_colname=smiles_colname,
                    chem_id_colname=chem_id_colname,
                    archive_smiles_colname=archive_smiles_colname,
                    archive_chem_id_colname=archive_chem_id_colname,
                    chunk_size=chunk_size,
                )
            )
    return work_items


def process_screen_config_from_args(args: argparse.Namespace) -> ProcessScreenConfig:
    return ProcessScreenConfig(
        input_paths=_normalize_input_paths(getattr(args, "input_path", None)),
        output_dir=str(getattr(args, "output_dir")),
        execution_mode=str(getattr(args, "execution_mode", "thread")),
        producer_processes=getattr(args, "producer_processes", "auto"),
        predict_queue_max_batches=getattr(args, "predict_queue_max_batches", "auto"),
        result_queue_max_batches=getattr(args, "result_queue_max_batches", "auto"),
        batch_checkpoint_size=int(getattr(args, "batch_checkpoint_size", 2048)),
        rows_per_shard=int(getattr(args, "rows_per_shard", 100_000)),
        row_group_size=int(getattr(args, "row_group_size", 4096)),
    )


def _coerce_positive_int(value: int | str, *, field_name: str) -> int:
    try:
        resolved = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be an integer or 'auto'") from exc

    if resolved <= 0:
        raise ValueError(f"{field_name} must be positive")
    return resolved


def resolve_graph_workers(value: int | str = "auto", *, cpu_count: int) -> int:
    if str(value).lower() != "auto":
        try:
            resolved = int(value)
        except (TypeError, ValueError) as exc:
            raise ValueError("graph_workers must be an integer or 'auto'") from exc
        if resolved < 0:
            raise ValueError("graph_workers must be non-negative")
        return resolved

    if cpu_count <= 2:
        return 0
    return max(1, min(4, cpu_count // 8))


def resolve_producer_processes(value: int | str, *, cpu_count: int, graph_workers: int) -> int:
    if str(value).lower() != "auto":
        return _coerce_positive_int(value, field_name="producer_processes")

    reserved = max(2, cpu_count // 12)
    available = max(1, cpu_count - reserved - max(0, graph_workers))
    scaled = max(4, min(32, available // 2))
    return min(max(1, cpu_count), scaled)


def resolve_queue_depth(value: int | str, *, producer_processes: int) -> int:
    if str(value).lower() != "auto":
        return _coerce_positive_int(value, field_name="queue_depth")
    return max(8, producer_processes * 2)


def plan_runtime(
    config: ProcessScreenConfig,
    *,
    cpu_count: int | None = None,
    gpu_count: int = 1,
    graph_workers: int | str = "auto",
) -> RuntimePlan:
    if config.execution_mode != "process":
        raise ValueError("plan_runtime requires execution_mode='process'")

    resolved_cpu_count = max(1, int(cpu_count or os.cpu_count() or 1))
    resolved_graph_workers = resolve_graph_workers(graph_workers, cpu_count=resolved_cpu_count)
    producer_processes = resolve_producer_processes(
        config.producer_processes,
        cpu_count=resolved_cpu_count,
        graph_workers=resolved_graph_workers,
    )

    return RuntimePlan(
        predictor_processes=max(1, int(gpu_count)),
        producer_processes=producer_processes,
        work_queue_max_batches=max(2, producer_processes * 2),
        predict_queue_max_batches=resolve_queue_depth(
            config.predict_queue_max_batches,
            producer_processes=producer_processes,
        ),
        result_queue_max_batches=resolve_queue_depth(
            config.result_queue_max_batches,
            producer_processes=producer_processes,
        ),
        batch_checkpoint_size=max(1, int(config.batch_checkpoint_size)),
        graph_workers=resolved_graph_workers,
    )


def create_process_topology(runtime: RuntimePlan, *, context: Any | None = None) -> ProcessTopology:
    ctx = context or multiprocessing.get_context("spawn")
    if ctx.get_start_method() != "spawn":
        raise ValueError("process screening requires a spawn multiprocessing context")

    return ProcessTopology(
        context=ctx,
        work_queue=ctx.JoinableQueue(maxsize=runtime.work_queue_max_batches),
        predict_queue=ctx.Queue(maxsize=runtime.predict_queue_max_batches),
        result_queue=ctx.Queue(maxsize=runtime.result_queue_max_batches),
        error_queue=ctx.Queue(),
    )


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _sanitize_path_segment(value: str | None) -> str:
    text = str(value or "").strip() or "input"
    separators = [os.sep]
    if os.altsep:
        separators.append(os.altsep)
    for separator in separators:
        text = text.replace(separator, "_")
    return text


def batch_manifest_path(output_dir: str | Path) -> Path:
    return Path(output_dir).expanduser().resolve() / "batch_manifest.jsonl"


def writer_run_state_path(output_dir: str | Path) -> Path:
    return Path(output_dir).expanduser().resolve() / "run_state.json"


def predicted_batch_shard_id(work_unit: WorkUnit) -> str:
    return _sanitize_path_segment(work_unit.group_id or work_unit.source_group or Path(work_unit.source_path).stem)


def predicted_batch_name(batch_index: int) -> str:
    return f"batch_{int(batch_index):06d}"


def predicted_batch_manifest_id(batch: PredictedBatch) -> str:
    return f"{predicted_batch_shard_id(batch.work_unit)}:{predicted_batch_name(batch.batch_index)}"


def batch_hits_output_path(output_dir: str | Path, shard_id: str, batch_id: str) -> Path:
    return Path(output_dir).expanduser().resolve() / "hits" / shard_id / f"{batch_id}.parquet"


def load_batch_manifest(path: str | Path) -> dict[str, dict[str, Any]]:
    manifest_path = Path(path).expanduser().resolve()
    if not manifest_path.exists():
        return {}

    records: dict[str, dict[str, Any]] = {}
    for line in manifest_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        records[str(record["batch_id"])] = record
    return records


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f"{path.name}.tmp")
    tmp_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp_path.replace(path)


def _write_batch_manifest_atomic(path: Path, records: Mapping[str, Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f"{path.name}.tmp")
    rendered = "".join(json.dumps(records[key], ensure_ascii=False) + "\n" for key in sorted(records))
    tmp_path.write_text(rendered, encoding="utf-8")
    tmp_path.replace(path)


def _drain_error_queue(error_queue: Any) -> ProcessErrorRecord | None:
    try:
        return error_queue.get_nowait()
    except queue.Empty:
        return None


def _terminate_processes(processes: list[Any]) -> None:
    for process in processes:
        if process.is_alive():
            process.terminate()
    for process in processes:
        process.join(timeout=5)


def _raise_process_error(error: ProcessErrorRecord) -> None:
    detail = f"{error.stage} worker {error.worker_name} failed: {error.message}"
    if error.traceback_text:
        detail = f"{detail}\n{error.traceback_text}"
    raise RuntimeError(detail)


def _derive_batch_row_range(batch: PredictedBatch) -> tuple[int | None, int | None]:
    if "source_row" in batch.input_frame.columns and not batch.input_frame.empty:
        source_rows = batch.input_frame["source_row"].dropna()
        if not source_rows.empty:
            return int(source_rows.min()), int(source_rows.max())

    if batch.work_unit.start_row is None:
        return None, None

    start_row = int(batch.work_unit.start_row)
    end_row = start_row + max(0, len(batch.input_frame) - 1)
    return start_row, end_row


def _build_batch_manifest_record(
    *,
    batch_id: str,
    shard_id: str,
    source_group: str,
    start_row: int | None,
    end_row: int | None,
    attempted_count: int,
    hit_count: int,
    output_path: str,
    status: str,
    error: str | None,
) -> dict[str, Any]:
    return {
        "batch_id": batch_id,
        "shard_id": shard_id,
        "source_group": source_group,
        "start_row": start_row,
        "end_row": end_row,
        "status": status,
        "attempted_count": int(attempted_count),
        "hit_count": int(hit_count),
        "output_path": output_path,
        "updated_at": _utc_now(),
        "error": error,
    }


def _write_writer_run_state(
    *,
    path: Path,
    manifest_records: Mapping[str, Mapping[str, Any]],
    status: str,
    error: str | None = None,
) -> None:
    committed = [record for record in manifest_records.values() if record.get("status") == "committed"]
    attempted_count = sum(int(record.get("attempted_count", 0)) for record in committed)
    hit_count = sum(int(record.get("hit_count", 0)) for record in committed)
    latest_output_path = ""
    if committed:
        latest_record = max(
            committed,
            key=lambda record: (str(record.get("updated_at", "")), str(record.get("batch_id", ""))),
        )
        latest_output_path = str(latest_record.get("output_path") or "")

    payload = {
        "status": status,
        "updated_at": _utc_now(),
        "attempted_count": attempted_count,
        "hit_count": hit_count,
        "output_path": latest_output_path,
        "progress": {
            "committed_batches": len(committed),
            "attempted_count": attempted_count,
            "hit_count": hit_count,
        },
        "error": error,
    }
    _write_json_atomic(path, payload)


def should_skip_batch(batch_id: str, manifest_records: Mapping[str, Mapping[str, Any]]) -> bool:
    record = manifest_records.get(batch_id)
    if not record or record.get("status") != "committed":
        return False

    output_path = str(record.get("output_path") or "").strip()
    return bool(output_path) and Path(output_path).exists()


def commit_prediction_batch(
    *,
    output_dir: str | Path,
    shard_id: str,
    batch_id: str,
    predicted_frame: pd.DataFrame,
    manifest_records: dict[str, dict[str, Any]],
    manifest_path: str | Path | None = None,
    source_group: str = "",
    start_row: int | None = None,
    end_row: int | None = None,
) -> dict[str, Any]:
    if "broad_spectrum" not in predicted_frame.columns:
        raise KeyError("predicted_frame must include a broad_spectrum column")

    hits = predicted_frame.loc[predicted_frame["broad_spectrum"] == 1].copy()
    target_path = batch_hits_output_path(output_dir, shard_id, batch_id)
    tmp_path = target_path.with_name(f"{target_path.name}.tmp")
    tmp_path.parent.mkdir(parents=True, exist_ok=True)
    hits.to_parquet(tmp_path, index=False)
    tmp_path.replace(target_path)
    output_path = str(target_path)

    record = _build_batch_manifest_record(
        batch_id=batch_id,
        shard_id=shard_id,
        source_group=source_group,
        start_row=start_row,
        end_row=end_row,
        attempted_count=len(predicted_frame),
        hit_count=len(hits),
        output_path=output_path,
        status="committed",
        error=None,
    )
    manifest_records[batch_id] = record

    if manifest_path is not None:
        _write_batch_manifest_atomic(Path(manifest_path).expanduser().resolve(), manifest_records)

    return record


def _build_error_record(stage: str, worker_name: str, exc: BaseException) -> ProcessErrorRecord:
    return ProcessErrorRecord(
        stage=stage,
        worker_name=worker_name,
        message=str(exc),
        traceback_text=traceback.format_exc(),
    )


def producer_main(
    *,
    producer_id: int,
    work_queue: Any,
    predict_queue: Any,
    error_queue: Any,
    batch_checkpoint_size: int,
) -> None:
    worker_name = f"producer-{producer_id}"

    try:
        while True:
            item = work_queue.get()
            try:
                if item is None:
                    return

                if not isinstance(item, ProcessWorkItem):
                    raise TypeError(f"{worker_name} received unexpected work item {type(item)!r}")

                batch_index = 0
                for frame in process_work_unit(
                    item.work_unit,
                    item.smiles_colname,
                    item.chem_id_colname,
                    item.archive_smiles_colname,
                    item.archive_chem_id_colname,
                    chunk_size=item.chunk_size,
                ):
                    for batch in iter_frame_batches(frame, batch_checkpoint_size):
                        predict_queue.put(
                            ProducerBatch(
                                producer_id=producer_id,
                                work_unit=item.work_unit,
                                batch_index=batch_index,
                                frame=batch,
                            )
                        )
                        batch_index += 1
            finally:
                task_done = getattr(work_queue, "task_done", None)
                if callable(task_done):
                    task_done()
    except Exception as exc:
        error_queue.put(_build_error_record("producer", worker_name, exc))
    finally:
        predict_queue.put(None)


def predictor_main(
    *,
    predict_queue: Any,
    result_queue: Any,
    error_queue: Any,
    producer_processes: int,
    aggregate_scores: bool,
    app_threshold: float,
    min_nkill: int,
    num_graph_workers: int | str = "auto",
    graph_batch_size: int = 1024,
    prefetch_batches: int = 2,
    classifier_workers: int | str = "auto",
    classifier_inflight_batches: int | str = "auto",
    deterministic_representation: bool = False,
    enable_profiling: bool = False,
) -> None:
    worker_name = "predictor-0"
    completed_producers = 0

    try:
        scheduler = create_scheduler(
            num_graph_workers=num_graph_workers,
            graph_batch_size=graph_batch_size,
            prefetch_batches=prefetch_batches,
            classifier_workers=classifier_workers,
            classifier_inflight_batches=classifier_inflight_batches,
            deterministic_representation=deterministic_representation,
        )

        while completed_producers < producer_processes:
            item = predict_queue.get()
            if item is None:
                completed_producers += 1
                continue

            if not isinstance(item, ProducerBatch):
                raise TypeError(f"{worker_name} received unexpected predict item {type(item)!r}")

            predicted = screen_frame_sync(
                frame=item.frame,
                aggregate_scores=aggregate_scores,
                app_threshold=app_threshold,
                min_nkill=min_nkill,
                scheduler=scheduler,
                num_graph_workers=num_graph_workers,
                graph_batch_size=graph_batch_size,
                prefetch_batches=prefetch_batches,
                classifier_workers=classifier_workers,
                classifier_inflight_batches=classifier_inflight_batches,
                enable_profiling=enable_profiling,
            )
            result_queue.put(
                PredictedBatch(
                    producer_id=item.producer_id,
                    work_unit=item.work_unit,
                    batch_index=item.batch_index,
                    input_frame=item.frame,
                    predicted_frame=predicted,
                )
            )
    except Exception as exc:
        error_queue.put(_build_error_record("predictor", worker_name, exc))
    finally:
        result_queue.put(None)


def writer_main(
    *,
    result_queue: Any,
    error_queue: Any,
    write_batch: Callable[[PredictedBatch], Any] | None = None,
    output_dir: str | Path | None = None,
    manifest_path: str | Path | None = None,
    run_state_path: str | Path | None = None,
) -> None:
    worker_name = "writer-0"
    sink = write_batch or _discard_batch
    resolved_manifest_path = None
    resolved_run_state_path = None
    manifest_records: dict[str, dict[str, Any]] = {}
    current_item: PredictedBatch | None = None
    current_item_committed = False

    if output_dir is not None:
        resolved_manifest_path = Path(manifest_path).expanduser().resolve() if manifest_path else batch_manifest_path(output_dir)
        resolved_run_state_path = Path(run_state_path).expanduser().resolve() if run_state_path else writer_run_state_path(output_dir)
        manifest_records = load_batch_manifest(resolved_manifest_path)
        _write_writer_run_state(path=resolved_run_state_path, manifest_records=manifest_records, status="running")

    try:
        while True:
            item = result_queue.get()
            if item is None:
                if resolved_run_state_path is not None:
                    _write_writer_run_state(
                        path=resolved_run_state_path,
                        manifest_records=manifest_records,
                        status="completed",
                    )
                return

            if not isinstance(item, PredictedBatch):
                raise TypeError(f"{worker_name} received unexpected result item {type(item)!r}")

            current_item = item
            current_item_committed = False
            if resolved_manifest_path is not None:
                batch_id = predicted_batch_manifest_id(item)
                if should_skip_batch(batch_id, manifest_records):
                    current_item = None
                    current_item_committed = False
                    continue

                shard_id = predicted_batch_shard_id(item.work_unit)
                start_row, end_row = _derive_batch_row_range(item)
                commit_prediction_batch(
                    output_dir=output_dir,
                    shard_id=shard_id,
                    batch_id=batch_id,
                    predicted_frame=item.predicted_frame,
                    manifest_records=manifest_records,
                    manifest_path=resolved_manifest_path,
                    source_group=item.work_unit.source_group,
                    start_row=start_row,
                    end_row=end_row,
                )
                current_item_committed = True
                if resolved_run_state_path is not None:
                    _write_writer_run_state(
                        path=resolved_run_state_path,
                        manifest_records=manifest_records,
                        status="running",
                    )

            sink(item)
            current_item = None
            current_item_committed = False
    except Exception as exc:
        if resolved_manifest_path is not None and current_item is not None and not current_item_committed:
            batch_id = predicted_batch_manifest_id(current_item)
            shard_id = predicted_batch_shard_id(current_item.work_unit)
            start_row, end_row = _derive_batch_row_range(current_item)
            manifest_records[batch_id] = _build_batch_manifest_record(
                batch_id=batch_id,
                shard_id=shard_id,
                source_group=current_item.work_unit.source_group,
                start_row=start_row,
                end_row=end_row,
                attempted_count=len(current_item.predicted_frame),
                hit_count=0,
                output_path="",
                status="failed",
                error=str(exc),
            )
            _write_batch_manifest_atomic(resolved_manifest_path, manifest_records)
        if resolved_run_state_path is not None:
            _write_writer_run_state(
                path=resolved_run_state_path,
                manifest_records=manifest_records,
                status="failed",
                error=str(exc),
            )
        error_queue.put(_build_error_record("writer", worker_name, exc))


def _discard_batch(_: PredictedBatch) -> None:
    return None


def run_process_pipeline(
    *,
    work_items: list[ProcessWorkItem],
    runtime: RuntimePlan,
    output_dir: str | Path,
    aggregate_scores: bool,
    app_threshold: float,
    min_nkill: int,
    num_graph_workers: int | str,
    graph_batch_size: int,
    prefetch_batches: int,
    classifier_workers: int | str,
    classifier_inflight_batches: int | str,
    deterministic_representation: bool,
    enable_profiling: bool,
) -> dict[str, Any]:
    resolved_output_dir = Path(output_dir).expanduser().resolve()
    batch_manifest = batch_manifest_path(resolved_output_dir)
    run_state = writer_run_state_path(resolved_output_dir)
    hits_dir = resolved_output_dir / "hits"
    hits_dir.mkdir(parents=True, exist_ok=True)

    topology = create_process_topology(runtime)
    processes: list[Any] = []
    try:
        writer_process = topology.context.Process(
            target=writer_main,
            kwargs={
                "result_queue": topology.result_queue,
                "error_queue": topology.error_queue,
                "output_dir": resolved_output_dir,
                "manifest_path": batch_manifest,
                "run_state_path": run_state,
            },
            name="writer-0",
        )
        processes.append(writer_process)

        predictor_process = topology.context.Process(
            target=predictor_main,
            kwargs={
                "predict_queue": topology.predict_queue,
                "result_queue": topology.result_queue,
                "error_queue": topology.error_queue,
                "producer_processes": runtime.producer_processes,
                "aggregate_scores": aggregate_scores,
                "app_threshold": app_threshold,
                "min_nkill": min_nkill,
                "num_graph_workers": runtime.graph_workers if str(num_graph_workers).lower() == "auto" else num_graph_workers,
                "graph_batch_size": graph_batch_size,
                "prefetch_batches": prefetch_batches,
                "classifier_workers": classifier_workers,
                "classifier_inflight_batches": classifier_inflight_batches,
                "deterministic_representation": deterministic_representation,
                "enable_profiling": enable_profiling,
            },
            name="predictor-0",
        )
        processes.append(predictor_process)

        for producer_index in range(runtime.producer_processes):
            processes.append(
                topology.context.Process(
                    target=producer_main,
                    kwargs={
                        "producer_id": producer_index,
                        "work_queue": topology.work_queue,
                        "predict_queue": topology.predict_queue,
                        "error_queue": topology.error_queue,
                        "batch_checkpoint_size": runtime.batch_checkpoint_size,
                    },
                    name=f"producer-{producer_index}",
                )
            )

        for process in processes:
            process.start()

        for item in work_items:
            while True:
                error = _drain_error_queue(topology.error_queue)
                if error is not None:
                    _terminate_processes(processes)
                    _raise_process_error(error)
                try:
                    topology.work_queue.put(item, timeout=0.1)
                    break
                except queue.Full:
                    continue

        for _ in range(runtime.producer_processes):
            while True:
                error = _drain_error_queue(topology.error_queue)
                if error is not None:
                    _terminate_processes(processes)
                    _raise_process_error(error)
                try:
                    topology.work_queue.put(None, timeout=0.1)
                    break
                except queue.Full:
                    continue

        while True:
            error = _drain_error_queue(topology.error_queue)
            if error is not None:
                _terminate_processes(processes)
                _raise_process_error(error)

            alive = False
            for process in processes:
                process.join(timeout=0.1)
                alive = alive or process.is_alive()
            if not alive:
                break

        error = _drain_error_queue(topology.error_queue)
        if error is not None:
            _raise_process_error(error)

        for process in processes:
            if process.exitcode not in (0, None):
                raise RuntimeError(f"{process.name} exited with code {process.exitcode}")

        return {
            "batch_manifest_path": str(batch_manifest),
            "run_state_path": str(run_state),
            "hits_dir": str(hits_dir),
        }
    finally:
        topology.close()


async def screen_paths_multiprocess(
    *,
    input_paths: list[ScreenInputPath] | tuple[ScreenInputPath, ...],
    output_dir: str | Path,
    execution_mode: str = "process",
    producer_processes: int | str = "auto",
    predict_queue_max_batches: int | str = "auto",
    result_queue_max_batches: int | str = "auto",
    batch_checkpoint_size: int = 2048,
    rows_per_shard: int = 100_000,
    row_group_size: int = 4096,
    smiles_colname: str = "smiles",
    chem_id_colname: str = "chem_id",
    archive_smiles_colname: str = DEFAULT_ARCHIVE_SMILES_COL,
    archive_chem_id_colname: str = DEFAULT_ARCHIVE_CHEM_ID_COL,
    grouping_mode: str = "source",
    cpu_workers: int | str = "auto",
    target_rows_per_group: int | str = "auto",
    target_bytes_per_group: int | str = "auto",
    chunk_size: int = DEFAULT_INPUT_CHUNK_SIZE,
    aggregate_scores: bool = True,
    app_threshold: float = 0.04374140128493309,
    min_nkill: int = 10,
    num_graph_workers: int | str = "auto",
    graph_batch_size: int = 1024,
    prefetch_batches: int = 2,
    classifier_workers: int | str = "auto",
    classifier_inflight_batches: int | str = "auto",
    deterministic_representation: bool = False,
    enable_profiling: bool = False,
) -> dict[str, Any]:
    if execution_mode != "process":
        raise ValueError("screen_paths_multiprocess requires execution_mode='process'")
    if not aggregate_scores:
        raise ValueError("process mode only supports aggregate score output because it persists hit-only batches")

    resolved_output_dir = Path(output_dir).expanduser().resolve()
    resolved_output_dir.mkdir(parents=True, exist_ok=True)

    config = ProcessScreenConfig(
        input_paths=[str(path) for path in _normalize_input_paths(input_paths)],
        output_dir=str(resolved_output_dir),
        execution_mode=execution_mode,
        producer_processes=producer_processes,
        predict_queue_max_batches=predict_queue_max_batches,
        result_queue_max_batches=result_queue_max_batches,
        batch_checkpoint_size=batch_checkpoint_size,
        rows_per_shard=rows_per_shard,
        row_group_size=row_group_size,
    )

    prepared_inputs, prepared_manifest, prepared_manifest_paths = prepare_process_inputs(
        input_paths=config.input_paths,
        output_dir=resolved_output_dir,
        smiles_colname=smiles_colname,
        chem_id_colname=chem_id_colname,
        rows_per_shard=config.rows_per_shard,
        row_group_size=config.row_group_size,
    )
    work_items = build_process_work_items(
        prepared_inputs=prepared_inputs,
        grouping_mode=grouping_mode,
        cpu_workers=cpu_workers,
        target_rows_per_group=target_rows_per_group,
        target_bytes_per_group=target_bytes_per_group,
        smiles_colname="smiles",
        chem_id_colname="chem_id",
        archive_smiles_colname=archive_smiles_colname,
        archive_chem_id_colname=archive_chem_id_colname,
        chunk_size=chunk_size,
    )
    runtime = plan_runtime(
        config,
        graph_workers=num_graph_workers,
    )
    pipeline_summary = run_process_pipeline(
        work_items=work_items,
        runtime=runtime,
        output_dir=resolved_output_dir,
        aggregate_scores=aggregate_scores,
        app_threshold=app_threshold,
        min_nkill=min_nkill,
        num_graph_workers=num_graph_workers,
        graph_batch_size=graph_batch_size,
        prefetch_batches=prefetch_batches,
        classifier_workers=classifier_workers,
        classifier_inflight_batches=classifier_inflight_batches,
        deterministic_representation=deterministic_representation,
        enable_profiling=enable_profiling,
    )

    return {
        "execution_mode": execution_mode,
        "input_paths": [str(path) for path in _resolve_existing_input_paths(config.input_paths)],
        "prepared_input_paths": [item.prepared_path for item in prepared_inputs],
        "source_groups": sorted(item.source_group for item in prepared_inputs),
        "prepared_manifest_path": str(prepared_manifest),
        "prepared_manifest_paths": prepared_manifest_paths,
        "batch_manifest_path": pipeline_summary["batch_manifest_path"],
        "run_state_path": pipeline_summary["run_state_path"],
        "hits_dir": pipeline_summary["hits_dir"],
        "work_unit_count": len(work_items),
        "runtime": {
            "predictor_processes": runtime.predictor_processes,
            "producer_processes": runtime.producer_processes,
            "work_queue_max_batches": runtime.work_queue_max_batches,
            "predict_queue_max_batches": runtime.predict_queue_max_batches,
            "result_queue_max_batches": runtime.result_queue_max_batches,
            "batch_checkpoint_size": runtime.batch_checkpoint_size,
            "graph_workers": runtime.graph_workers,
        },
        "prediction_runtime": None,
    }
