"""Batch screening helpers for TSV/CSV inputs and archive bundles."""

from __future__ import annotations

import asyncio
import csv
import heapq
import itertools
import shutil
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any, TYPE_CHECKING

import pandas as pd

from src.models import MoleculeInfo
from src.service import get_scheduler
from src.screening_planner import ScreeningPlanner, PlannerConfig
from src.screening_sources import (
    process_work_unit,
    _uniquify_chem_ids,
    DEFAULT_ARCHIVE_PATTERN,
    DEFAULT_ARCHIVE_SMILES_COL,
    DEFAULT_ARCHIVE_CHEM_ID_COL,
    DEFAULT_INPUT_CHUNK_SIZE,
)

if TYPE_CHECKING:
    from src.prediction_scheduler import PredictionScheduler


def _build_molecules(frame: pd.DataFrame) -> list[MoleculeInfo]:
    return [
        MoleculeInfo(smiles=row.smiles, chem_id=row.chem_id)
        for row in frame[["smiles", "chem_id"]].itertuples(index=False)
    ]


def iter_frame_batches(frame: pd.DataFrame, batch_size: int):
    """Yield stable row-order batches for process-mode producer checkpoints."""

    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

    for start in range(0, len(frame), batch_size):
        yield frame.iloc[start : start + batch_size].reset_index(drop=True).copy()


def screen_frame_sync(**kwargs) -> pd.DataFrame:
    """Run ``screen_frame`` from a synchronous multiprocessing worker entrypoint."""

    return asyncio.run(screen_frame(**kwargs))


async def screen_frame(
    frame: pd.DataFrame,
    aggregate_scores: bool,
    app_threshold: float,
    min_nkill: int,
    scheduler: "PredictionScheduler | None" = None,
    num_graph_workers: int | str = "auto",
    graph_batch_size: int = 1024,
    prefetch_batches: int = 2,
    classifier_workers: int | str = "auto",
    classifier_inflight_batches: int | str = "auto",
    enable_profiling: bool = False,
) -> pd.DataFrame:
    """Run batch prediction on a normalized frame and attach metadata."""

    if frame.empty:
        raise ValueError("No screening rows available")

    scheduler = scheduler or get_scheduler()

    items = await scheduler.predict_molecules(
        molecules=_build_molecules(frame),
        aggregate_scores=aggregate_scores,
        app_threshold=app_threshold,
        min_nkill=min_nkill,
        num_graph_workers=num_graph_workers,
        graph_batch_size=graph_batch_size,
        prefetch_batches=prefetch_batches,
        classifier_workers=classifier_workers,
        classifier_inflight_batches=classifier_inflight_batches,
        enable_profiling=enable_profiling,
    )

    chunk_result = pd.DataFrame(items)
    if aggregate_scores:
        chunk_result = chunk_result.merge(frame, on="chem_id", how="left")
    else:
        split = chunk_result["pred_id"].astype(str).str.rsplit(":", n=1, expand=True)
        chunk_result["chem_id"] = split[0]
        chunk_result["strain_name"] = split[1]
        chunk_result = chunk_result.merge(frame, on="chem_id", how="left")

    if "source_group" not in chunk_result.columns:
        chunk_result["source_group"] = "input"
    if "input_order" in chunk_result.columns:
        sort_columns = ["input_order"]
        if "strain_name" in chunk_result.columns:
            sort_columns.append("strain_name")
        chunk_result = chunk_result.sort_values(sort_columns, kind="stable")
    return chunk_result.reset_index(drop=True)


@dataclass
class ScreeningSummary:
    normalized_rows: int
    predicted_rows: int
    normalized_input_path: Path
    predictions_all_path: Path
    grouped_outputs: List[Dict[str, str]]
    grouping_mode: str = "auto"
    cpu_workers_selected: int = 1
    prefetch_queue_size_selected: int = 4
    work_unit_count: int = 0
    target_rows_per_group: int = 0
    target_bytes_per_group: int = 0
    profiling: Dict[str, Any] | None = None


def _new_screening_profile() -> dict[str, Any]:
    return {
        "queue_wait_seconds": 0.0,
        "screen_frame_seconds": 0.0,
        "chunks_predicted": 0,
        "prediction_seconds": 0.0,
        "representation_seconds": 0.0,
        "strain_expand_seconds": 0.0,
        "xgboost_seconds": 0.0,
        "prediction_frame_seconds": 0.0,
        "growth_inhibition_seconds": 0.0,
        "aggregate_scores_seconds": 0.0,
        "result_records_seconds": 0.0,
        "graph_build_seconds": 0.0,
        "graph_items": 0,
    }


def _merge_prediction_profile(
    aggregate: dict[str, Any],
    runtime_snapshot: dict[str, Any] | None,
) -> None:
    if not runtime_snapshot:
        return
    last_profile = runtime_snapshot.get("last_profile")
    if not isinstance(last_profile, dict):
        return

    for key in (
        "representation_seconds",
        "strain_expand_seconds",
        "xgboost_seconds",
        "prediction_frame_seconds",
        "growth_inhibition_seconds",
        "aggregate_scores_seconds",
        "result_records_seconds",
    ):
        value = float(last_profile.get(key, 0.0))
        aggregate[key] += value
        aggregate["prediction_seconds"] += value

    graph_build = last_profile.get("graph_build")
    if isinstance(graph_build, dict):
        aggregate["graph_build_seconds"] += float(graph_build.get("graph_total_seconds", 0.0))
        aggregate["graph_items"] += int(graph_build.get("graph_items", 0))


def _append_tsv(df: pd.DataFrame, path: Path | str, write_header: bool) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if write_header else "a"
    df.to_csv(p, sep="\t", index=False, header=write_header, mode=mode)


def _aggregate_sort_key_from_row(row: dict[str, str]) -> tuple[int, int, float, int]:
    def _as_int(value: str | None, default: int = 0) -> int:
        try:
            return int(float(value)) if value is not None else default
        except (TypeError, ValueError):
            return default

    def _as_float(value: str | None, default: float = float("inf")) -> float:
        try:
            return float(value) if value is not None else default
        except (TypeError, ValueError):
            return default

    return (
        -_as_int(row.get("broad_spectrum")),
        -_as_int(row.get("ginhib_total")),
        _as_float(row.get("apscore_total")),
        _as_int(row.get("input_order")),
    )


def _rewrite_sorted_aggregate_outputs(
    chunk_paths: list[Path],
    predictions_all_path: Path,
    by_source_root: Path,
) -> list[dict[str, str]]:
    if not chunk_paths:
        return []

    merged_output_path = predictions_all_path.with_suffix(".sorted.tmp")
    sorted_by_source_root = by_source_root.with_name(f"{by_source_root.name}.sorted_tmp")
    if sorted_by_source_root.exists():
        shutil.rmtree(sorted_by_source_root)
    sorted_by_source_root.mkdir(parents=True, exist_ok=True)

    readers: list[csv.DictReader] = []
    handles = []
    heap: list[tuple[tuple[int, int, float, int], int, int, dict[str, str]]] = []
    tie_breaker = itertools.count()
    grouped_outputs: list[dict[str, str]] = []
    source_handles: dict[str, Any] = {}
    source_writers: dict[str, csv.DictWriter] = {}

    try:
        fieldnames: list[str] | None = None
        for chunk_index, chunk_path in enumerate(chunk_paths):
            handle = chunk_path.open("r", encoding="utf-8", newline="")
            reader = csv.DictReader(handle, delimiter="\t")
            row = next(reader, None)
            if row is None:
                handle.close()
                continue
            handles.append(handle)
            readers.append(reader)
            if fieldnames is None:
                fieldnames = list(reader.fieldnames or row.keys())
            heapq.heappush(
                heap,
                (_aggregate_sort_key_from_row(row), next(tie_breaker), len(readers) - 1, row),
            )

        if fieldnames is None:
            return []

        with merged_output_path.open("w", encoding="utf-8", newline="") as merged_handle:
            merged_writer = csv.DictWriter(merged_handle, fieldnames=fieldnames, delimiter="\t")
            merged_writer.writeheader()

            while heap:
                _, _, reader_index, row = heapq.heappop(heap)
                merged_writer.writerow(row)

                source_group = str(row.get("source_group", "input"))
                if source_group not in source_writers:
                    source_dir = sorted_by_source_root / source_group
                    source_dir.mkdir(parents=True, exist_ok=True)
                    source_path = source_dir / "predictions.tsv"
                    source_handle = source_path.open("w", encoding="utf-8", newline="")
                    source_writer = csv.DictWriter(source_handle, fieldnames=fieldnames, delimiter="\t")
                    source_writer.writeheader()
                    source_handles[source_group] = source_handle
                    source_writers[source_group] = source_writer
                    grouped_outputs.append(
                        {
                            "source_group": source_group,
                            "path": str(by_source_root / source_group / "predictions.tsv"),
                        }
                    )

                source_writers[source_group].writerow(row)

                next_row = next(readers[reader_index], None)
                if next_row is not None:
                    heapq.heappush(
                        heap,
                        (_aggregate_sort_key_from_row(next_row), next(tie_breaker), reader_index, next_row),
                    )
    finally:
        for handle in source_handles.values():
            handle.close()
        for handle in handles:
            handle.close()

    merged_output_path.replace(predictions_all_path)
    if by_source_root.exists():
        shutil.rmtree(by_source_root)
    sorted_by_source_root.rename(by_source_root)
    return grouped_outputs


def _finalize_frame_in_source_order(
    frame: pd.DataFrame,
    *,
    dedupe_smiles: bool,
    seen_smiles: set[str],
    seen_ids: set[str],
    global_order_state: list[int],
) -> pd.DataFrame:
    frame = frame.dropna(subset=["smiles"]).copy()
    frame["smiles"] = frame["smiles"].astype(str)
    frame["chem_id"] = frame["chem_id"].astype(str)

    if dedupe_smiles:
        frame = frame.drop_duplicates(subset=["smiles"], keep="first")
        mask = ~frame["smiles"].isin(seen_smiles)
        frame = frame[mask].copy()
        seen_smiles.update(frame["smiles"])

    if frame.empty:
        return frame.reset_index(drop=True)

    frame = _uniquify_chem_ids(frame, seen_ids)
    frame["input_order"] = range(global_order_state[0], global_order_state[0] + len(frame))
    global_order_state[0] += len(frame)
    return frame.reset_index(drop=True)


async def screen_path(
    input_path: str | Path,
    output_dir: str | Path,
    smiles_colname: str = "smiles",
    chem_id_colname: str = "chem_id",
    archive_pattern: str = DEFAULT_ARCHIVE_PATTERN,
    archive_smiles_colname: str = DEFAULT_ARCHIVE_SMILES_COL,
    archive_chem_id_colname: str = DEFAULT_ARCHIVE_CHEM_ID_COL,
    sqlite_table: str | None = None,
    sqlite_query: str | None = None,
    dedupe_smiles: bool = True,
    aggregate_scores: bool = True,
    app_threshold: float = 0.04374140128493309,
    min_nkill: int = 10,
    chunk_size: int = DEFAULT_INPUT_CHUNK_SIZE,
    grouping_mode: str = "auto",
    cpu_workers: int | str = "auto",
    prefetch_queue_size: int | str = "auto",
    target_rows_per_group: int | str = "auto",
    target_bytes_per_group: int | str = "auto",
    scheduler: "PredictionScheduler | None" = None,
    num_graph_workers: int | str = "auto",
    graph_batch_size: int = 1024,
    prefetch_batches: int = 2,
    classifier_workers: int | str = "auto",
    classifier_inflight_batches: int | str = "auto",
    enable_profiling: bool = False,
    prediction_row_budget: int | None = None,
) -> ScreeningSummary:

    out_dir = Path(output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    scheduler = scheduler or get_scheduler()

    normalized_input_path = out_dir / "normalized_input.tsv"
    predictions_all_path = out_dir / "predictions_all.tsv"

    config = PlannerConfig(
        grouping_mode=grouping_mode,
        cpu_workers=cpu_workers,
        prefetch_queue_size=prefetch_queue_size,
        target_rows_per_group=target_rows_per_group,
        target_bytes_per_group=target_bytes_per_group,
        chunk_size=chunk_size
    )
    planner = ScreeningPlanner(config)
    work_units = planner.plan(
        input_path,
        archive_pattern=archive_pattern,
        sqlite_table=sqlite_table,
        sqlite_query=sqlite_query,
        smiles_colname=smiles_colname
    )

    seen_sources: set[str] = set()
    grouped_outputs: list[dict[str, str]] = []

    normalized_rows_count = 0
    predicted_rows_count = 0
    first_chunk = True

    queue = asyncio.Queue(maxsize=max(1, planner.prefetch_queue_size))
    stop_event = threading.Event()
    global_order_state = [0]

    seen_smiles: set[str] = set()
    seen_ids: set[str] = set()
    profiling = _new_screening_profile() if enable_profiling else None
    buffered_frames: dict[int, dict[int, pd.DataFrame]] = {}
    completed_units: set[int] = set()
    next_unit_index = 0
    next_chunk_index = 0
    producer_finished = False
    aggregate_chunk_dir = out_dir / ".aggregate_sort_chunks"
    aggregate_chunk_paths: list[Path] = []
    aggregate_chunk_counter = 0

    def _producer_task_func(unit_index, unit):
        if stop_event.is_set():
            return
        try:
            generator = process_work_unit(
                unit,
                smiles_colname,
                chem_id_colname,
                archive_smiles_colname,
                archive_chem_id_colname,
                chunk_size=chunk_size
            )
            for chunk_index, frame in enumerate(generator):
                if stop_event.is_set():
                    break

                asyncio.run_coroutine_threadsafe(
                    queue.put(("frame", unit_index, chunk_index, frame.reset_index(drop=True))),
                    loop,
                ).result()

            if not stop_event.is_set():
                asyncio.run_coroutine_threadsafe(queue.put(("done", unit_index, None, None)), loop).result()
        except Exception as e:
            if not stop_event.is_set():
                stop_event.set()
                try:
                    asyncio.run_coroutine_threadsafe(queue.put(e), loop).result()
                except Exception:
                    pass

    def _produce_all(loop, q, stop):
        try:
            with ThreadPoolExecutor(max_workers=planner.cpu_workers) as executor:
                from concurrent.futures import wait
                futures = [executor.submit(_producer_task_func, unit_index, unit) for unit_index, unit in enumerate(work_units)]
                wait(futures)

            if not stop.is_set():
                asyncio.run_coroutine_threadsafe(q.put(None), loop).result()
        except Exception as e:
            if not stop.is_set():
                stop.set()
                try:
                    asyncio.run_coroutine_threadsafe(q.put(e), loop).result()
                except Exception:
                    pass

    loop = asyncio.get_running_loop()
    producer_task = asyncio.create_task(asyncio.to_thread(_produce_all, loop, queue, stop_event))

    def _buffer_item(item: Any) -> None:
        nonlocal producer_finished
        if item is None:
            producer_finished = True
            return
        if isinstance(item, Exception):
            raise item
        kind, unit_index, chunk_index, payload = item
        if kind == "frame":
            buffered_frames.setdefault(unit_index, {})[int(chunk_index)] = payload
        elif kind == "done":
            completed_units.add(unit_index)

    def _pop_next_ready_frame() -> pd.DataFrame | None:
        nonlocal next_unit_index, next_chunk_index
        while True:
            unit_frames = buffered_frames.get(next_unit_index)
            if unit_frames and next_chunk_index in unit_frames:
                raw_frame = unit_frames.pop(next_chunk_index)
                if not unit_frames:
                    buffered_frames.pop(next_unit_index, None)
                next_chunk_index += 1
                return _finalize_frame_in_source_order(
                    raw_frame,
                    dedupe_smiles=dedupe_smiles,
                    seen_smiles=seen_smiles,
                    seen_ids=seen_ids,
                    global_order_state=global_order_state,
                )

            if next_unit_index in completed_units:
                completed_units.remove(next_unit_index)
                next_unit_index += 1
                next_chunk_index = 0
                continue

            return None

    async def _consume_ready_frames() -> None:
        nonlocal normalized_rows_count, predicted_rows_count, first_chunk
        nonlocal aggregate_chunk_counter

        while True:
            frames_to_predict: list[pd.DataFrame] = []
            rows_to_predict = 0
            while True:
                frame = _pop_next_ready_frame()
                if frame is None:
                    break
                if frame.empty:
                    continue
                frames_to_predict.append(frame)
                rows_to_predict += len(frame)
                if prediction_row_budget is None:
                    break
                if prediction_row_budget is not None and rows_to_predict >= prediction_row_budget:
                    break
                try:
                    queue_wait_start = time.perf_counter() if enable_profiling else None
                    item = await asyncio.wait_for(queue.get(), timeout=0.01)
                    if enable_profiling and profiling is not None:
                        profiling["queue_wait_seconds"] += time.perf_counter() - queue_wait_start
                    _buffer_item(item)
                except asyncio.TimeoutError:
                    break

            if frames_to_predict:
                frame = pd.concat(frames_to_predict, ignore_index=True)

                screen_frame_start = time.perf_counter() if enable_profiling else None
                prediction_task = asyncio.create_task(
                    screen_frame(
                        frame=frame,
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
                )
                try:
                    while True:
                        if prediction_task.done():
                            predicted = await prediction_task
                            break

                        queue_wait_start = time.perf_counter() if enable_profiling else None
                        queue_get_task = asyncio.create_task(queue.get())
                        done, _ = await asyncio.wait(
                            {prediction_task, queue_get_task},
                            return_when=asyncio.FIRST_COMPLETED,
                        )

                        if queue_get_task in done:
                            if enable_profiling and profiling is not None and queue_wait_start is not None:
                                profiling["queue_wait_seconds"] += time.perf_counter() - queue_wait_start
                            item = queue_get_task.result()
                            _buffer_item(item)
                            if prediction_task in done:
                                predicted = await prediction_task
                                break
                            continue

                        queue_get_task.cancel()
                        try:
                            await queue_get_task
                        except asyncio.CancelledError:
                            pass
                except Exception:
                    if not prediction_task.done():
                        prediction_task.cancel()
                        try:
                            await prediction_task
                        except asyncio.CancelledError:
                            pass
                    raise
                if enable_profiling and profiling is not None:
                    profiling["screen_frame_seconds"] += time.perf_counter() - screen_frame_start
                    profiling["chunks_predicted"] += 1
                    _merge_prediction_profile(profiling, scheduler.runtime_snapshot())

                if aggregate_scores:
                    rank_columns: list[str] = []
                    ascending: list[bool] = []
                    if "broad_spectrum" in predicted.columns:
                        rank_columns.append("broad_spectrum")
                        ascending.append(False)
                    if "ginhib_total" in predicted.columns:
                        rank_columns.append("ginhib_total")
                        ascending.append(False)
                    if "apscore_total" in predicted.columns:
                        rank_columns.append("apscore_total")
                        ascending.append(True)
                    if "input_order" in predicted.columns:
                        rank_columns.append("input_order")
                        ascending.append(True)
                    if rank_columns:
                        predicted = predicted.sort_values(rank_columns, ascending=ascending, kind="stable")
                else:
                    rank_columns = [column for column in ["input_order", "strain_name"] if column in predicted.columns]
                    if rank_columns:
                        predicted = predicted.sort_values(rank_columns, kind="stable")

                _append_tsv(frame, normalized_input_path, write_header=first_chunk)
                _append_tsv(predicted, predictions_all_path, write_header=first_chunk)

                if aggregate_scores:
                    aggregate_chunk_dir.mkdir(parents=True, exist_ok=True)
                    chunk_path = aggregate_chunk_dir / f"predicted_chunk_{aggregate_chunk_counter:06d}.tsv"
                    predicted.to_csv(chunk_path, sep="\t", index=False)
                    aggregate_chunk_paths.append(chunk_path)
                    aggregate_chunk_counter += 1

                if "source_group" in predicted.columns:
                    for source_group, group_frame in predicted.groupby("source_group", sort=True):
                        group_dir = out_dir / "by_source" / str(source_group)
                        group_path = group_dir / "predictions.tsv"

                        is_first_for_source = source_group not in seen_sources
                        _append_tsv(group_frame, group_path, write_header=is_first_for_source)

                        if is_first_for_source:
                            seen_sources.add(source_group)
                            grouped_outputs.append({"source_group": str(source_group), "path": str(group_path)})

                normalized_rows_count += len(frame)
                predicted_rows_count += len(predicted)
                first_chunk = False
                continue

            break

    try:
        while True:
            await _consume_ready_frames()
            if producer_finished:
                await _consume_ready_frames()
                break
            queue_wait_start = time.perf_counter() if enable_profiling else None
            item = await queue.get()
            if enable_profiling and profiling is not None:
                profiling["queue_wait_seconds"] += time.perf_counter() - queue_wait_start
            _buffer_item(item)

        await producer_task

        if aggregate_scores and aggregate_chunk_paths:
            grouped_outputs = _rewrite_sorted_aggregate_outputs(
                chunk_paths=aggregate_chunk_paths,
                predictions_all_path=predictions_all_path,
                by_source_root=out_dir / "by_source",
            )

    finally:
        stop_event.set()
        if not producer_task.done():
            producer_task.cancel()
            while not producer_task.done():
                try:
                    queue.get_nowait()
                except asyncio.QueueEmpty:
                    pass
                await asyncio.sleep(0.01)

        try:
            await producer_task
        except asyncio.CancelledError:
            pass
        except Exception:
            pass
        if aggregate_chunk_dir.exists():
            shutil.rmtree(aggregate_chunk_dir, ignore_errors=True)

    if normalized_rows_count == 0:
        raise ValueError(f"No valid screening rows could be loaded from {input_path}")

    return ScreeningSummary(
        normalized_rows=normalized_rows_count,
        predicted_rows=predicted_rows_count,
        normalized_input_path=normalized_input_path,
        predictions_all_path=predictions_all_path,
        grouped_outputs=grouped_outputs,
        grouping_mode=planner.grouping_mode,
        cpu_workers_selected=planner.cpu_workers,
        prefetch_queue_size_selected=planner.prefetch_queue_size,
        work_unit_count=len(work_units),
        target_rows_per_group=planner.target_rows_per_group,
        target_bytes_per_group=planner.target_bytes_per_group,
        profiling=profiling,
    )
