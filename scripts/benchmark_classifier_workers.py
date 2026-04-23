from __future__ import annotations

import argparse
import asyncio
import json
import os
import shutil
import statistics
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.models import MoleculeInfo
from src.prediction_scheduler import PredictionScheduler
from src.predictor import AntimicrobialPredictor
from src.stream_enumeration_screen import (
    DEFAULT_SCAFFOLD_FILE,
    IndexSpace,
    _load_fragment_smiles,
    _load_scaffolds,
    _materialize_batch,
    stream_enumeration_screen,
)

DEFAULT_WORKER_CANDIDATES = ["1", "2", "4", "6", "8", "12"]
DEFAULT_SAMPLE_SIZE = 8192
DEFAULT_REPEATS = 3


def _parse_worker_value(raw: str) -> int | str:
    token = str(raw).strip().lower()
    if token == "auto":
        return "auto"
    value = int(token)
    if value <= 0:
        raise ValueError("worker candidates must be positive integers or 'auto'")
    return value


def _worker_label(worker: int | str) -> str:
    return str(worker)


def _summarize_candidate_runs(
    *,
    mode: str,
    worker: int | str,
    sample_size: int,
    runs: list[dict[str, Any]],
) -> dict[str, Any]:
    wall_seconds = [float(run["wall_seconds"]) for run in runs]
    median_wall_seconds = float(statistics.median(wall_seconds))
    return {
        "mode": mode,
        "worker": worker,
        "sample_size": int(sample_size),
        "repeat_count": len(runs),
        "median_wall_seconds": median_wall_seconds,
        "median_molecules_per_second": float(sample_size) / median_wall_seconds if median_wall_seconds else 0.0,
        "runs": runs,
    }


def _select_best_candidate(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    if not candidates:
        raise ValueError("at least one candidate summary is required")
    return max(
        candidates,
        key=lambda item: (
            float(item["median_molecules_per_second"]),
            -float(item["median_wall_seconds"]),
        ),
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Benchmark classifier_workers across predictor and stream screening modes.",
    )
    parser.add_argument("--mode", choices=["predictor", "stream"], required=True)
    parser.add_argument(
        "--workers",
        nargs="+",
        default=DEFAULT_WORKER_CANDIDATES,
        help="Candidate classifier_workers values to benchmark. Accepts integers or 'auto'.",
    )
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    parser.add_argument("--sample-size", type=int, default=DEFAULT_SAMPLE_SIZE)
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument("--output", help="Optional JSON output path")
    parser.add_argument(
        "--scratch-root",
        default="data/05.stream_tasks/benchmarks/classifier_workers/_scratch",
        help="Temporary run directory root for stream benchmarks",
    )
    parser.add_argument("--keep-runs", action="store_true", help="Keep stream benchmark run directories")

    parser.add_argument(
        "--scaffold-file",
        default=str(DEFAULT_SCAFFOLD_FILE),
        help="Scaffold file used for predictor or stream benchmark materialization",
    )
    parser.add_argument("--ordinary-library", required=True, help="Ordinary fragment library CSV")
    parser.add_argument("--pos13-library", required=True, help="Position-13 fragment library CSV")
    parser.add_argument("--run-state", dest="run_state_source", help="Optional upstream run_state.json provenance")
    parser.add_argument("--chunk-manifest", dest="chunk_manifest_source", help="Optional upstream chunk manifest provenance")

    parser.add_argument("--classifier-backend", choices=["auto", "pickle", "timber"], default="pickle")
    parser.add_argument("--app-threshold", type=float, default=0.04374140128493309)
    parser.add_argument("--min-nkill", type=int, default=10)
    parser.add_argument("--num-graph-workers", default="auto")
    parser.add_argument("--graph-batch-size", type=int, default=1024)
    parser.add_argument("--prefetch-batches", type=int, default=2)
    parser.add_argument("--classifier-inflight-batches", default="auto")
    parser.add_argument("--deterministic-representation", action="store_true")
    parser.add_argument("--shard-size", type=int, default=4096)
    parser.add_argument("--prediction-batch-size", type=int, default=512)
    parser.add_argument("--enumeration-workers", default="auto")
    parser.add_argument("--enumeration-prefetch-batches", default="auto")
    return parser


def _build_index_space(
    *,
    ordinary_library: str | Path,
    pos13_library: str | Path,
    scaffold_file: str | Path | None,
) -> tuple[list[Any], list[str], list[str], IndexSpace]:
    scaffolds = _load_scaffolds(
        scaffold_file=scaffold_file,
        scaffold_dir=None,
        scaffold_catalog=None,
    )
    ordinary_fragments = _load_fragment_smiles(ordinary_library)
    pos13_fragments = _load_fragment_smiles(pos13_library)
    space = IndexSpace(
        scaffold_count=len(scaffolds),
        pos3_count=len(ordinary_fragments),
        pos4_count=len(ordinary_fragments),
        pos12_count=len(ordinary_fragments),
        pos13_count=len(pos13_fragments),
    )
    return scaffolds, ordinary_fragments, pos13_fragments, space


def _materialize_molecules(
    *,
    scaffold_file: str | Path,
    ordinary_library: str | Path,
    pos13_library: str | Path,
    start_index: int,
    sample_size: int,
) -> list[MoleculeInfo]:
    scaffolds, ordinary_fragments, pos13_fragments, space = _build_index_space(
        ordinary_library=ordinary_library,
        pos13_library=pos13_library,
        scaffold_file=scaffold_file,
    )
    stop_index = min(space.total_combinations, int(start_index) + int(sample_size))
    rows = _materialize_batch(
        start_idx=int(start_index),
        end_idx=stop_index,
        space=space,
        scaffolds=scaffolds,
        ordinary_fragments=ordinary_fragments,
        pos13_fragments=pos13_fragments,
    )
    return [MoleculeInfo(smiles=row["smiles"], chem_id=row["chem_id"]) for row in rows]


async def _benchmark_predictor_once(
    *,
    predictor: AntimicrobialPredictor,
    molecules: list[MoleculeInfo],
    worker: int | str,
    args: argparse.Namespace,
) -> dict[str, Any]:
    scheduler = PredictionScheduler(
        predictor=predictor,
        num_graph_workers=args.num_graph_workers,
        graph_batch_size=args.graph_batch_size,
        prefetch_batches=args.prefetch_batches,
        classifier_workers=worker,
        classifier_inflight_batches=args.classifier_inflight_batches,
        deterministic_representation=args.deterministic_representation,
    )
    start = time.perf_counter()
    records = await scheduler.predict_molecules(
        molecules=molecules,
        aggregate_scores=True,
        app_threshold=args.app_threshold,
        min_nkill=args.min_nkill,
        enable_profiling=True,
    )
    wall_seconds = time.perf_counter() - start
    return {
        "worker": worker,
        "attempted_count": len(records),
        "wall_seconds": wall_seconds,
        "molecules_per_second": float(len(records)) / wall_seconds if wall_seconds else 0.0,
        "prediction_runtime": scheduler.runtime_snapshot(),
    }


async def _benchmark_stream_once(
    *,
    predictor: AntimicrobialPredictor,
    worker: int | str,
    args: argparse.Namespace,
    repeat_index: int,
) -> dict[str, Any]:
    scheduler = PredictionScheduler(
        predictor=predictor,
        num_graph_workers=args.num_graph_workers,
        graph_batch_size=args.graph_batch_size,
        prefetch_batches=args.prefetch_batches,
        classifier_workers=worker,
        classifier_inflight_batches=args.classifier_inflight_batches,
        deterministic_representation=args.deterministic_representation,
    )
    scratch_root = Path(args.scratch_root).expanduser().resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)
    output_dir = Path(
        tempfile.mkdtemp(
            prefix=f"stream_{_worker_label(worker)}_r{repeat_index}_",
            dir=scratch_root,
        )
    )
    try:
        stop_index = int(args.start_index) + int(args.sample_size)
        start = time.perf_counter()
        summary = await stream_enumeration_screen(
            output_dir=output_dir,
            scaffold_file=args.scaffold_file,
            scaffold_dir=None,
            scaffold_catalog=None,
            ordinary_library=args.ordinary_library,
            pos13_library=args.pos13_library,
            run_state_source=args.run_state_source,
            chunk_manifest_source=args.chunk_manifest_source,
            start_index=args.start_index,
            stop_index=stop_index,
            shard_size=args.shard_size,
            prediction_batch_size=args.prediction_batch_size,
            scheduler=scheduler,
            app_threshold=args.app_threshold,
            min_nkill=args.min_nkill,
            classifier_backend=args.classifier_backend,
            num_graph_workers=args.num_graph_workers,
            graph_batch_size=args.graph_batch_size,
            prefetch_batches=args.prefetch_batches,
            classifier_workers=worker,
            classifier_inflight_batches=args.classifier_inflight_batches,
            enumeration_workers=args.enumeration_workers,
            enumeration_prefetch_batches=args.enumeration_prefetch_batches,
            deterministic_representation=args.deterministic_representation,
            enable_profiling=True,
        )
        wall_seconds = time.perf_counter() - start
        return {
            "worker": worker,
            "attempted_count": summary.attempted_count,
            "hit_count": summary.hit_count,
            "wall_seconds": wall_seconds,
            "molecules_per_second": float(summary.attempted_count) / wall_seconds if wall_seconds else 0.0,
            "prediction_runtime": scheduler.runtime_snapshot(),
            "output_dir": str(output_dir),
        }
    finally:
        if not args.keep_runs and output_dir.exists():
            shutil.rmtree(output_dir, ignore_errors=True)


async def _run_predictor_benchmark(args: argparse.Namespace) -> dict[str, Any]:
    os.environ["MOLE_CLASSIFIER_BACKEND"] = args.classifier_backend
    predictor = AntimicrobialPredictor()
    await predictor.ensure_loaded()
    molecules = _materialize_molecules(
        scaffold_file=args.scaffold_file,
        ordinary_library=args.ordinary_library,
        pos13_library=args.pos13_library,
        start_index=args.start_index,
        sample_size=args.sample_size,
    )
    worker_candidates = [_parse_worker_value(worker) for worker in args.workers]
    candidate_summaries: list[dict[str, Any]] = []
    for worker in worker_candidates:
        runs = []
        for repeat_index in range(args.repeats):
            run = await _benchmark_predictor_once(
                predictor=predictor,
                molecules=molecules,
                worker=worker,
                args=args,
            )
            run["repeat_index"] = repeat_index
            runs.append(run)
        candidate_summaries.append(
            _summarize_candidate_runs(
                mode="predictor",
                worker=worker,
                sample_size=len(molecules),
                runs=runs,
            )
        )
    return {
        "mode": "predictor",
        "classifier_backend": args.classifier_backend,
        "sample_size": len(molecules),
        "repeats": args.repeats,
        "worker_candidates": worker_candidates,
        "candidates": candidate_summaries,
        "best_candidate": _select_best_candidate(candidate_summaries),
    }


async def _run_stream_benchmark(args: argparse.Namespace) -> dict[str, Any]:
    os.environ["MOLE_CLASSIFIER_BACKEND"] = args.classifier_backend
    predictor = AntimicrobialPredictor()
    await predictor.ensure_loaded()
    worker_candidates = [_parse_worker_value(worker) for worker in args.workers]
    candidate_summaries: list[dict[str, Any]] = []
    for worker in worker_candidates:
        runs = []
        for repeat_index in range(args.repeats):
            run = await _benchmark_stream_once(
                predictor=predictor,
                worker=worker,
                args=args,
                repeat_index=repeat_index,
            )
            run["repeat_index"] = repeat_index
            runs.append(run)
        candidate_summaries.append(
            _summarize_candidate_runs(
                mode="stream",
                worker=worker,
                sample_size=args.sample_size,
                runs=runs,
            )
        )
    return {
        "mode": "stream",
        "classifier_backend": args.classifier_backend,
        "sample_size": args.sample_size,
        "repeats": args.repeats,
        "worker_candidates": worker_candidates,
        "candidates": candidate_summaries,
        "best_candidate": _select_best_candidate(candidate_summaries),
    }


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.repeats <= 0:
        raise ValueError("repeats must be positive")
    if args.sample_size <= 0:
        raise ValueError("sample-size must be positive")

    if args.mode == "predictor":
        payload = asyncio.run(_run_predictor_benchmark(args))
    else:
        payload = asyncio.run(_run_stream_benchmark(args))

    rendered = json.dumps(payload, indent=2, ensure_ascii=False)
    if args.output:
        output_path = Path(args.output).expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
