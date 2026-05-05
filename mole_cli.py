#!/usr/bin/env python3
"""Unified command line interface for MolE / Timber / REINVENT4."""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Sequence

import pandas as pd
import torch

from src.batch_screening import screen_path
from src.mole_representation import process_representation
from src.panel_scoring import load_panel_config, compute_panel_scores_from_dataframe
from src.preprocess_screening_input import preprocess_to_parquet
from src.screening_process_pipeline import process_screen_config_from_args, screen_paths_multiprocess
from src.stream_enumeration_screen import (
    DEFAULT_SCAFFOLD_FILE as DEFAULT_STREAM_SCAFFOLD_FILE,
    stream_enumeration_screen,
)
from scripts.benchmark_screening_inputs import benchmark_paths
from src.models import MoleculeInput, ReinventScoreRequest, ScoreObjective
from src.reinvent4_workflow import load_objective_spec, resolve_path
from src.reinvent_scoring import score_reinvent_predictions
from src.service import get_predictor, get_scheduler, create_scheduler
from src.classifier_backend import inspect_classifier_backends
from src.site_reward import scaffold_from_file


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_MOLE_MODEL = os.environ.get(
    "MOLE_MOLE_MODEL_PATH",
    str(REPO_ROOT / "pretrained_model" / "model_ginconcat_btwin_100k_d8000_l0.0001"),
)
DEFAULT_OPTIMIZE_SCRIPT = REPO_ROOT / "workflows" / "reinvent4" / "scripts" / "run_long_rl.py"


def _resolve_device(requested: str) -> str:
    requested = requested.strip().lower()
    if requested not in {"auto", "cpu", "cuda"} and not requested.startswith("cuda:"):
        raise ValueError("device must be one of auto, cpu, cuda, or cuda:<index>")
    if requested == "cpu":
        return "cpu"
    if requested.startswith("cuda:"):
        return requested
    if requested == "cuda":
        return "cuda:0"
    if torch.cuda.is_available():
        cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
        if cuda_visible:
            gpu_id = cuda_visible.split(",")[0].strip()
            return f"cuda:{gpu_id}"
        return "cuda:0"
    return "cpu"


def _dump_json(payload: Any, output: str | None) -> None:
    rendered = json.dumps(payload, indent=2, ensure_ascii=False)
    if output:
        Path(output).expanduser().resolve().write_text(rendered + "\n", encoding="utf-8")
    else:
        sys.stdout.write(rendered + "\n")


def _write_embedding_output(df: pd.DataFrame, output: str | None, format_name: str) -> None:
    if format_name == "tsv":
        text = df.to_csv(sep="\t", index=False)
        if output:
            Path(output).expanduser().resolve().write_text(text, encoding="utf-8")
        else:
            sys.stdout.write(text)
        return
    _dump_json(df.to_dict(orient="records"), output)


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(df.to_csv(sep="\t", index=False), encoding="utf-8")


def _apply_classifier_backend_arg(args: argparse.Namespace) -> None:
    backend = getattr(args, "classifier_backend", None)
    if backend:
        os.environ["MOLE_CLASSIFIER_BACKEND"] = backend


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mole", description="Unified MolE / Timber / REINVENT4 CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    doctor = subparsers.add_parser("doctor", help="Check the local environment and model availability")
    doctor.add_argument(
        "--strict-gpu",
        action="store_true",
        help="Fail when CUDA is not available instead of only warning",
    )
    doctor.add_argument(
        "--env-file",
        default=str(REPO_ROOT / "workflows" / "reinvent4" / "configs" / "local.env.recommended.example"),
        help="Optional REINVENT4 env file to validate",
    )
    doctor.add_argument(
        "--scaffold-file",
        default=str(REPO_ROOT / "workflows" / "reinvent4" / "inputs" / "scaffolds" / "mother_scaffold.template.smi"),
        help="Optional scaffold file to validate",
    )
    doctor.add_argument(
        "--objective-file",
        default=str(REPO_ROOT / "workflows" / "reinvent4" / "inputs" / "objectives" / "pathogen_group_a.site_reward.prototype.json"),
        help="Optional objective file to validate",
    )

    embed = subparsers.add_parser("embed", help="Generate MolE embeddings from SMILES")
    embed.add_argument("--smiles", nargs="+", action="extend", required=True, help="One or more SMILES strings")
    embed.add_argument("--chem-id", nargs="+", action="extend", dest="chem_ids", help="Optional chemical identifiers")
    embed.add_argument("--mole-model", default=DEFAULT_MOLE_MODEL, help="MolE checkpoint directory")
    embed.add_argument("--device", default="auto", help="MolE device: auto, cpu, cuda, or cuda:<index>")
    embed.add_argument(
        "--deterministic-representation",
        action="store_true",
        help="Force deterministic CUDA MolE forward passes for reproducibility checks.",
    )
    embed.add_argument("--output", help="Optional output file path")
    embed.add_argument("--format", choices=["json", "tsv"], default="json", help="Output format")

    predict = subparsers.add_parser("predict", help="Predict strain-level antimicrobial activity")
    predict.add_argument("--smiles", nargs="+", action="extend", required=True, help="One or more SMILES strings")
    predict.add_argument("--chem-id", nargs="+", action="extend", dest="chem_ids", help="Optional chemical identifiers")
    predict.add_argument("--aggregate-scores", action="store_true", help="Return aggregate antimicrobial scores")
    predict.add_argument(
        "--classifier-backend",
        choices=["auto", "timber", "pickle"],
        help="Classifier backend for strain probabilities. Use pickle for high-throughput screening.",
    )
    predict.add_argument("--app-threshold", type=float, default=0.04374140128493309, help="Growth inhibition threshold")
    predict.add_argument("--min-nkill", type=int, default=10, help="Broad-spectrum threshold")
    predict.add_argument(
        "--num-graph-workers",
        default="auto",
        help="CPU workers used to build MolE graph mini-batches",
    )
    predict.add_argument(
        "--graph-batch-size",
        type=int,
        default=1024,
        help="Mini-batch size for MolE graph construction and forward passes",
    )
    predict.add_argument(
        "--prefetch-batches",
        type=int,
        default=2,
        help="Prefetched graph mini-batches per worker",
    )
    predict.add_argument(
        "--classifier-workers",
        default="auto",
        help="CPU workers used for post-MolE classifier/aggregation overlap",
    )
    predict.add_argument(
        "--classifier-inflight-batches",
        default="auto",
        help="Maximum aggregate-classifier batches kept in flight ahead of ordered merge",
    )
    predict.add_argument(
        "--profiling",
        action="store_true",
        help="Emit stage profiling for MolE graph build and prediction",
    )
    predict.add_argument(
        "--deterministic-representation",
        action="store_true",
        help="Force deterministic CUDA MolE forward passes for reproducibility checks.",
    )
    predict.add_argument("--output", help="Optional output file path")

    preprocess = subparsers.add_parser(
        "preprocess-screening-input",
        help="Convert large CSV/TSV screening tables into Parquet shards",
    )
    preprocess.add_argument("--input-path", required=True, help="Path to a CSV/TSV input file")
    preprocess.add_argument("--output-dir", required=True, help="Directory where Parquet shards will be written")
    preprocess.add_argument("--smiles-colname", default="smiles", help="SMILES column name in the input file")
    preprocess.add_argument("--chem-id-colname", default="chem_id", help="chem_id column name in the input file")
    preprocess.add_argument("--source-group", default="input", help="Logical source group label for the output shards")
    preprocess.add_argument("--rows-per-shard", type=int, default=100000, help="Rows per Parquet shard")
    preprocess.add_argument("--row-group-size", type=int, default=4000, help="Parquet row group size")
    preprocess.add_argument("--output", help="Optional path for writing the preprocessing manifest JSON")

    benchmark = subparsers.add_parser(
        "benchmark-screening-inputs",
        help="Benchmark screening input formats",
    )
    benchmark.add_argument("--input-path", action="append", required=True, help="One or more input paths to benchmark")
    benchmark.add_argument("--output", help="Optional path for writing benchmark JSON")

    screen = subparsers.add_parser(
        "screen",
        help="Batch screen CSV/TSV/archive/SQLite inputs for antimicrobial potential",
    )
    screen.add_argument(
        "--input-path",
        action="append",
        required=True,
        help="One or more CSV/TSV/Parquet/SQLite/tar inputs. Repeat the flag for multiple sources.",
    )
    screen.add_argument("--output-dir", required=True, help="Directory where screening outputs will be written")
    screen.add_argument(
        "--execution-mode",
        choices=["thread", "process"],
        default="thread",
        help="Use thread mode for current behavior or process mode for CPU-heavy producer pools.",
    )
    screen.add_argument("--smiles-colname", default="smiles", help="SMILES column name for tabular inputs")
    screen.add_argument("--chem-id-colname", default="chem_id", help="chem_id column name for tabular inputs")
    screen.add_argument(
        "--archive-pattern",
        default="*_scheme_b_unique_products.csv",
        help="Filename pattern for CSV members inside tar archives",
    )
    screen.add_argument(
        "--archive-smiles-colname",
        default="product_smiles_canonical",
        help="SMILES column name inside archive CSV members",
    )
    screen.add_argument(
        "--archive-chem-id-colname",
        default="example_combo_id",
        help="chem_id column name inside archive CSV members",
    )
    screen.add_argument("--sqlite-table", help="SQLite table name to read when input is a database file")
    screen.add_argument("--sqlite-query", help="SQLite query to read when input is a database file")
    screen.add_argument(
        "--classifier-backend",
        choices=["auto", "timber", "pickle"],
        help="Classifier backend for strain probabilities. Use pickle for high-throughput screening.",
    )
    screen.add_argument(
        "--no-dedupe-smiles",
        dest="dedupe_smiles",
        action="store_false",
        help="Keep duplicate SMILES instead of collapsing to the first occurrence",
    )
    screen.set_defaults(dedupe_smiles=True)
    score_group = screen.add_mutually_exclusive_group()
    score_group.add_argument(
        "--aggregate-scores",
        dest="aggregate_scores",
        action="store_true",
        help="Return aggregate screening scores (default)",
    )
    score_group.add_argument(
        "--per-strain",
        dest="aggregate_scores",
        action="store_false",
        help="Return per-strain probabilities instead of aggregate scores",
    )
    screen.set_defaults(aggregate_scores=True)
    screen.add_argument("--app-threshold", type=float, default=0.04374140128493309, help="Growth inhibition threshold")
    screen.add_argument("--min-nkill", type=int, default=10, help="Broad-spectrum threshold")

    panel_group = screen.add_argument_group("Pathogen panel scoring")
    panel_group.add_argument(
        "--panel-file",
        default=None,
        help="Path to pathogen-selective panel JSON config. Enables panel scoring on aggregate output.",
    )
    panel_group.add_argument(
        "--panel-lambda",
        type=float,
        default=None,
        help="Override lambda (commensal penalty weight) from panel config.",
    )
    panel_group.add_argument(
        "--panel-min-pathogen-hard",
        type=int,
        default=None,
        help="Minimum pathogen_hard count to retain a hit when --panel-file is set.",
    )
    panel_group.add_argument(
        "--panel-min-selectivity",
        type=float,
        default=None,
        help="Minimum selectivity_score to retain a hit when --panel-file is set.",
    )

    tune_group = screen.add_argument_group("Tuning")
    tune_group.add_argument(
        "--grouping-mode", default="auto", choices=["auto", "source", "chunk", "none"],
        help="How to group input records (default: auto)"
    )
    tune_group.add_argument(
        "--cpu-workers", default="auto",
        help="Number of CPU workers for preprocessing (default: auto)"
    )
    tune_group.add_argument(
        "--producer-processes",
        default="auto",
        help="Number of producer processes to use in process execution mode",
    )
    tune_group.add_argument(
        "--target-rows-per-group", default="auto",
        help="Target number of rows per group for chunking (default: auto)"
    )
    tune_group.add_argument(
        "--target-bytes-per-group", default="auto",
        help="Target bytes per group for chunking in auto mode (default: auto)"
    )
    tune_group.add_argument(
        "--predict-queue-max-batches",
        default="auto",
        help="Bound on ready-to-predict batches queued in process execution mode",
    )
    tune_group.add_argument(
        "--result-queue-max-batches",
        default="auto",
        help="Bound on predicted result batches queued in process execution mode",
    )
    tune_group.add_argument(
        "--batch-checkpoint-size",
        type=int,
        default=2048,
        help="Rows per checkpoint commit for resumable process execution mode",
    )
    tune_group.add_argument(
        "--rows-per-shard",
        type=int,
        default=100000,
        help="Rows per prepared Parquet shard in process execution mode",
    )
    tune_group.add_argument(
        "--row-group-size",
        type=int,
        default=4096,
        help="Parquet row-group size for process-mode preprocessing",
    )
    tune_group.add_argument(
        "--input-chunk-size", "--chunk-size", type=int, default=10000, dest="input_chunk_size",
        help="Rows generated per chunk while reading inputs (default 10000)"
    )
    tune_group.add_argument(
        "--max-batch-size", type=int, default=16384,
        help="Hard cap on GPU batch size"
    )
    tune_group.add_argument(
        "--target-gpu-memory-fraction", type=float, default=0.8,
        help="Fraction of free VRAM to fill"
    )
    tune_group.add_argument(
        "--prefetch-queue-size", default="auto",
        help="Number of pre-normalized batches to queue for GPU"
    )
    tune_group.add_argument(
        "--num-graph-workers",
        default="auto",
        help="CPU workers used to build MolE graph mini-batches",
    )
    tune_group.add_argument(
        "--graph-batch-size",
        type=int,
        default=1024,
        help="Mini-batch size for MolE graph construction and forward passes",
    )
    tune_group.add_argument(
        "--prefetch-batches",
        type=int,
        default=2,
        help="Prefetched graph mini-batches per worker",
    )
    tune_group.add_argument(
        "--classifier-workers",
        default="auto",
        help="CPU workers used for post-MolE classifier/aggregation overlap",
    )
    tune_group.add_argument(
        "--classifier-inflight-batches",
        default="auto",
        help="Maximum aggregate-classifier batches kept in flight ahead of ordered merge",
    )
    tune_group.add_argument(
        "--deterministic-representation",
        action="store_true",
        help="Force deterministic CUDA MolE forward passes for reproducibility checks.",
    )
    tune_group.add_argument(
        "--prediction-row-budget",
        type=int,
        default=8192,
        help="Combine consecutive ready screening chunks into prediction calls up to this many rows",
    )
    tune_group.add_argument(
        "--profiling",
        action="store_true",
        help="Write stage profiling into manifest.json for screening runs",
    )

    stream = subparsers.add_parser(
        "stream-enumeration-screen",
        help="Enumerate scaffold/fragment combinations on the fly and persist only hit shards",
    )
    stream.add_argument("--output-dir", required=True, help="Directory where run_state, manifest, and hit shards are written")
    scaffold_group = stream.add_mutually_exclusive_group()
    scaffold_group.add_argument(
        "--scaffold-file",
        help=f"Single scaffold .smi file to enumerate (default: {DEFAULT_STREAM_SCAFFOLD_FILE})",
    )
    scaffold_group.add_argument("--scaffold-dir", help="Directory of scaffold .smi files for multi-scaffold runs")
    scaffold_group.add_argument("--scaffold-catalog", help="CSV/TSV catalog with scaffold_slug and scaffold_smiles/scaffold_file")
    stream.add_argument("--ordinary-library", required=True, help="CSV containing the shared ordinary fragment pool")
    stream.add_argument("--pos13-library", required=True, help="CSV containing the position-13 sugar fragment pool")
    stream.add_argument("--run-state", dest="run_state_source", help="Optional upstream run_state.json for provenance")
    stream.add_argument("--chunk-manifest", dest="chunk_manifest_source", help="Optional upstream chunk manifest for provenance")
    stream.add_argument("--start-index", type=int, default=0, help="Inclusive global combination start index")
    stream.add_argument("--stop-index", type=int, help="Exclusive global combination stop index")
    stream.add_argument("--shard-size", type=int, default=100000, help="Combinations per resumable shard")
    stream.add_argument("--prediction-batch-size", type=int, default=1024, help="Combinations per prediction call inside a shard")
    stream.add_argument(
        "--enumeration-workers",
        default="auto",
        help="CPU workers used to materialize enumeration batches ahead of prediction",
    )
    stream.add_argument(
        "--enumeration-prefetch-batches",
        default="auto",
        help="Bounded number of enumeration batches kept in flight ahead of ordered prediction",
    )
    stream.add_argument(
        "--classifier-backend",
        choices=["auto", "timber", "pickle"],
        help="Classifier backend for aggregate screening predictions.",
    )
    stream.add_argument("--app-threshold", type=float, default=0.04374140128493309, help="Growth inhibition threshold")
    stream.add_argument("--min-nkill", type=int, default=10, help="Broad-spectrum threshold")
    stream.add_argument(
        "--num-graph-workers",
        default="auto",
        help="CPU workers used to build MolE graph mini-batches",
    )
    stream.add_argument(
        "--graph-batch-size",
        type=int,
        default=1024,
        help="Mini-batch size for MolE graph construction and forward passes",
    )
    stream.add_argument(
        "--prefetch-batches",
        type=int,
        default=2,
        help="Prefetched graph mini-batches per worker",
    )
    stream.add_argument(
        "--classifier-workers",
        default="auto",
        help="CPU workers used for post-MolE classifier/aggregation overlap",
    )
    stream.add_argument(
        "--classifier-inflight-batches",
        default="auto",
        help="Maximum aggregate-classifier batches kept in flight ahead of ordered merge",
    )
    stream.add_argument(
        "--deterministic-representation",
        action="store_true",
        help="Force deterministic CUDA MolE forward passes for reproducibility checks.",
    )
    stream.add_argument(
        "--profiling",
        action="store_true",
        help="Enable prediction profiling while screening enumerated combinations",
    )

    score = subparsers.add_parser("score", help="Compute REINVENT4 rewards from MolE probabilities")
    score.add_argument("--objective-file", required=True, help="Path to a normalized objective JSON file")
    score.add_argument("--smiles", nargs="+", action="extend", required=True, help="One or more SMILES strings")
    score.add_argument("--chem-id", nargs="+", action="extend", dest="chem_ids", help="Optional chemical identifiers")
    score.add_argument(
        "--scaffold-file",
        default=str(REPO_ROOT / "workflows" / "reinvent4" / "inputs" / "scaffolds" / "mother_scaffold.template.smi"),
        help="Optional scaffold file used when site_reward is enabled and scaffold_smiles is absent",
    )
    score.add_argument("--app-threshold", type=float, default=0.04374140128493309, help="Shared inhibition threshold")
    score.add_argument("--min-nkill", type=int, default=10, help="Broad-spectrum reference threshold")
    score.add_argument("--output", help="Optional output file path")

    optimize = subparsers.add_parser("optimize", help="Run the REINVENT4 optimization workflow")
    optimize.add_argument("--template", required=True, help="Path to a REINVENT4 TOML template")
    optimize.add_argument("--env-file", required=True, help="Path to a local.env style runtime file")
    optimize.add_argument("--experiment-spec", required=True, help="Path to a long-run RL spec JSON file")
    optimize.add_argument("--scaffold-file", required=True, help="Path to the LibInvent scaffold .smi file")
    optimize.add_argument(
        "--output-root",
        default=str(REPO_ROOT / "workflows" / "reinvent4" / "results" / "runs"),
        help="Root directory for long-run outputs",
    )
    optimize.add_argument("--run-id", default="", help="Optional deterministic run identifier")
    optimize.add_argument("--dry-run", action="store_true", help="Render config and stop before launching REINVENT4")
    optimize.add_argument(
        "--script",
        default=str(DEFAULT_OPTIMIZE_SCRIPT),
        help="Workflow script to invoke (defaults to chunked long-run runner)",
    )

    return parser


def _command_doctor(args: argparse.Namespace) -> int:
    probe = inspect_classifier_backends()
    outputs: list[str] = []
    ok = True

    outputs.append(f"Python: {sys.version.split()[0]}")
    outputs.append(f"PyTorch: {torch.__version__}")
    outputs.append(f"CUDA available: {torch.cuda.is_available()}")
    if not torch.cuda.is_available():
        message = "WARNING: CUDA is not available; MolE will run on CPU."
        outputs.append(message)
        if args.strict_gpu:
            ok = False

    mol_e_model = Path(os.environ.get("MOLE_MOLE_MODEL_PATH", DEFAULT_MOLE_MODEL)).expanduser().resolve()
    mol_e_config = mol_e_model / "config.yaml"
    mol_e_weights = mol_e_model / "model.pth"
    mol_e_ok = mol_e_model.is_dir() and mol_e_config.is_file() and mol_e_weights.is_file()
    outputs.append(
        "MolE model directory: "
        f"{mol_e_model} [{'ok' if mol_e_ok else 'missing'}] "
        f"(config.yaml={'ok' if mol_e_config.is_file() else 'missing'}, "
        f"model.pth={'ok' if mol_e_weights.is_file() else 'missing'})"
    )
    if not mol_e_ok:
        ok = False

    strain_panel_primary = REPO_ROOT / "data" / "01.prepare_training_data" / "maier_screening_results.tsv.gz"
    strain_panel_fallback = REPO_ROOT / "workflows" / "reinvent4" / "inputs" / "strain_index.tsv"
    if strain_panel_primary.is_file():
        outputs.append(f"Strain panel: primary ok -> {strain_panel_primary}")
    elif strain_panel_fallback.is_file():
        outputs.append(
            f"Strain panel: fallback ok -> {strain_panel_fallback} (primary missing: {strain_panel_primary})"
        )
    else:
        outputs.append(
            f"Strain panel: missing primary and fallback ({strain_panel_primary} / {strain_panel_fallback})"
        )
        ok = False

    outputs.append(
        f"Pickle backend: {'ok' if probe.pickle_available else 'missing'} -> {probe.pickle_path}"
    )
    outputs.append(
        f"Timber backend: {'ok' if probe.timber_available else 'missing'} -> {probe.timber_library_path or probe.timber_model_dir}"
    )
    outputs.append(f"Selected classifier backend: {probe.selected_backend}")

    env_file = Path(args.env_file).expanduser().resolve()
    outputs.append(f"REINVENT4 env file: {env_file} [{'ok' if env_file.is_file() else 'missing'}]")
    if not env_file.is_file():
        ok = False

    scaffold_file = Path(args.scaffold_file).expanduser().resolve()
    outputs.append(f"REINVENT4 scaffold file: {scaffold_file} [{'ok' if scaffold_file.is_file() else 'missing'}]")
    if not scaffold_file.is_file():
        ok = False

    objective_file = Path(args.objective_file).expanduser().resolve()
    outputs.append(f"REINVENT4 objective file: {objective_file} [{'ok' if objective_file.is_file() else 'missing'}]")
    if not objective_file.is_file():
        ok = False

    for line in outputs:
        print(line)
    return 0 if ok else 1


def _command_embed(args: argparse.Namespace) -> int:
    device = _resolve_device(args.device)
    smiles = list(args.smiles)
    chem_ids = list(args.chem_ids) if args.chem_ids is not None else None
    if chem_ids is not None and len(chem_ids) != len(smiles):
        raise ValueError("--chem-id length must match --smiles length")

    frame = pd.DataFrame(
        {
            "smiles": smiles,
            "chem_id": chem_ids if chem_ids is not None else [f"mol{i + 1}" for i in range(len(smiles))],
        }
    )
    representation = process_representation(
        dataset_path=frame,
        smile_column_str="smiles",
        id_column_str="chem_id",
        pretrained_dir=args.mole_model,
        device=device,
        deterministic_representation=bool(getattr(args, "deterministic_representation", False)),
    ).reset_index(drop=False)
    _write_embedding_output(representation, args.output, args.format)
    return 0


async def _predict_async(args: argparse.Namespace) -> dict[str, Any]:
    _apply_classifier_backend_arg(args)
    scheduler = get_scheduler(
        num_graph_workers=getattr(args, "num_graph_workers", "auto"),
        graph_batch_size=int(getattr(args, "graph_batch_size", 1024)),
        prefetch_batches=int(getattr(args, "prefetch_batches", 2)),
        classifier_workers=getattr(args, "classifier_workers", "auto"),
        classifier_inflight_batches=getattr(args, "classifier_inflight_batches", "auto"),
        deterministic_representation=bool(getattr(args, "deterministic_representation", False)),
    )
    normalized = MoleculeInput(
        smiles=list(args.smiles),
        chem_id=list(args.chem_ids) if args.chem_ids is not None else None,
        aggregate_scores=bool(args.aggregate_scores),
        app_threshold=float(args.app_threshold),
        min_nkill=int(args.min_nkill),
    ).normalize()
    items = await scheduler.predict_molecules(
        molecules=normalized.molecules,
        aggregate_scores=normalized.aggregate_scores,
        app_threshold=normalized.app_threshold,
        min_nkill=normalized.min_nkill,
        num_graph_workers=getattr(args, "num_graph_workers", "auto"),
        graph_batch_size=int(getattr(args, "graph_batch_size", 1024)),
        prefetch_batches=int(getattr(args, "prefetch_batches", 2)),
        classifier_workers=getattr(args, "classifier_workers", "auto"),
        classifier_inflight_batches=getattr(args, "classifier_inflight_batches", "auto"),
        enable_profiling=bool(getattr(args, "profiling", False)),
        deterministic_representation=bool(getattr(args, "deterministic_representation", False)),
    )
    return {
        "mode": "aggregate" if normalized.aggregate_scores else "per_strain",
        "items": items,
        "prediction_runtime": scheduler.runtime_snapshot(),
    }


def _command_predict(args: argparse.Namespace) -> int:
    payload = asyncio_run(_predict_async(args))
    _dump_json(payload, args.output)
    return 0


def _command_preprocess_screening_input(args: argparse.Namespace) -> int:
    manifest = preprocess_to_parquet(
        input_path=args.input_path,
        output_dir=args.output_dir,
        smiles_colname=args.smiles_colname,
        chem_id_colname=args.chem_id_colname,
        source_group=args.source_group,
        rows_per_shard=int(args.rows_per_shard),
        row_group_size=int(args.row_group_size),
    )
    _dump_json(manifest, getattr(args, "output", None))
    return 0


def _command_benchmark_screening_inputs(args: argparse.Namespace) -> int:
    result = benchmark_paths(list(args.input_path))
    _dump_json(result, getattr(args, "output", None))
    return 0


def _apply_panel_scoring_to_predictions(
    *,
    predictions_path: Path,
    panel_config,
    lambda_override: float | None,
    min_pathogen_hard: int | None,
    min_selectivity: float | None,
) -> dict[str, Any]:
    """Compute panel scores on aggregate predictions and optionally filter rows."""
    from src.panel_scoring import compute_panel_scores_from_dataframe

    df = pd.read_csv(predictions_path, sep="\t")
    if "pred_id" in df.columns:
        per_strain_df = df
    else:
        split = df["chem_id"].astype(str) + ":" + df["strain_name"].astype(str)
        per_strain_df = df.copy()
        per_strain_df["pred_id"] = split

    scores = compute_panel_scores_from_dataframe(
        per_strain_df,
        panel_config,
        pred_id_col="pred_id",
        probability_col="1",
        lambda_=lambda_override,
    )
    scores = scores.reset_index()
    df = df.merge(scores, on="chem_id", how="left")

    rows_before = len(df)
    if min_pathogen_hard is not None or min_selectivity is not None:
        mask = pd.Series(True, index=df.index)
        if min_pathogen_hard is not None:
            mask &= df["pathogen_hard"] >= min_pathogen_hard
        if min_selectivity is not None:
            mask &= df["selectivity_score"] >= min_selectivity
        if "broad_spectrum" in df.columns:
            mask |= df["broad_spectrum"] == 1
        df = df[mask].copy()

    df.to_csv(predictions_path, sep="\t", index=False)
    return {
        "panel_scoring": {
            "panel_label": panel_config.label,
            "panel_mode": panel_config.mode,
            "pathogen_count": len(panel_config.pathogen_panel),
            "commensal_count": len(panel_config.commensal_sparing_panel),
            "lambda": float(lambda_override if lambda_override is not None else panel_config.lambda_),
            "threshold": panel_config.app_threshold,
            "tau": panel_config.tau,
            "min_pathogen_hard": min_pathogen_hard,
            "min_selectivity": min_selectivity,
            "rows_before_filter": rows_before,
            "rows_after_filter": len(df),
        }
    }


def _command_screen(args: argparse.Namespace) -> int:
    async def _run() -> dict[str, Any]:
        _apply_classifier_backend_arg(args)
        process_config = process_screen_config_from_args(args)
        output_dir = Path(process_config.output_dir).expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        resolved_input_paths = [
            str(Path(input_path).expanduser().resolve())
            for input_path in process_config.input_paths
        ]

        if process_config.execution_mode == "thread":
            if len(process_config.input_paths) != 1:
                raise ValueError(
                    "screen thread mode currently supports exactly one --input-path; "
                    "use --execution-mode process for multiple inputs"
                )
            screen_input_path = process_config.input_paths[0]
            scheduler = create_scheduler(
                max_batch_size=getattr(args, "max_batch_size", 2048),
                target_memory_fraction=getattr(args, "target_gpu_memory_fraction", 0.8),
                num_graph_workers=getattr(args, "num_graph_workers", "auto"),
                graph_batch_size=int(getattr(args, "graph_batch_size", 1024)),
                prefetch_batches=int(getattr(args, "prefetch_batches", 2)),
                classifier_workers=getattr(args, "classifier_workers", "auto"),
                classifier_inflight_batches=getattr(args, "classifier_inflight_batches", "auto"),
                deterministic_representation=bool(getattr(args, "deterministic_representation", False)),
            )

            summary = await screen_path(
                input_path=screen_input_path,
                output_dir=output_dir,
                smiles_colname=getattr(args, "smiles_colname", "smiles"),
                chem_id_colname=getattr(args, "chem_id_colname", "chem_id"),
                archive_pattern=getattr(args, "archive_pattern", "*_scheme_b_unique_products.csv"),
                archive_smiles_colname=getattr(args, "archive_smiles_colname", "product_smiles_canonical"),
                archive_chem_id_colname=getattr(args, "archive_chem_id_colname", "example_combo_id"),
                sqlite_table=getattr(args, "sqlite_table", None),
                sqlite_query=getattr(args, "sqlite_query", None),
                dedupe_smiles=bool(getattr(args, "dedupe_smiles", True)),
                aggregate_scores=bool(getattr(args, "aggregate_scores", True)),
                app_threshold=float(getattr(args, "app_threshold", 0.04374140128493309)),
                min_nkill=int(getattr(args, "min_nkill", 10)),
                chunk_size=int(getattr(args, "input_chunk_size", 10000)),
                prefetch_queue_size=getattr(args, "prefetch_queue_size", "auto"),
                grouping_mode=getattr(args, "grouping_mode", "auto"),
                cpu_workers=getattr(args, "cpu_workers", "auto"),
                target_rows_per_group=getattr(args, "target_rows_per_group", "auto"),
                target_bytes_per_group=getattr(args, "target_bytes_per_group", "auto"),
                scheduler=scheduler,
                num_graph_workers=getattr(args, "num_graph_workers", "auto"),
                graph_batch_size=int(getattr(args, "graph_batch_size", 1024)),
                prefetch_batches=int(getattr(args, "prefetch_batches", 2)),
                classifier_workers=getattr(args, "classifier_workers", "auto"),
                classifier_inflight_batches=getattr(args, "classifier_inflight_batches", "auto"),
                enable_profiling=bool(getattr(args, "profiling", False)),
                prediction_row_budget=int(getattr(args, "prediction_row_budget", 8192)),
            )
            manifest = {
                "input_path": resolved_input_paths[0],
                "input_paths": resolved_input_paths,
                "output_dir": str(output_dir),
                "execution_mode": process_config.execution_mode,
                "producer_processes": process_config.producer_processes,
                "predict_queue_max_batches": process_config.predict_queue_max_batches,
                "result_queue_max_batches": process_config.result_queue_max_batches,
                "batch_checkpoint_size": process_config.batch_checkpoint_size,
                "dedupe_smiles": bool(args.dedupe_smiles),
                "aggregate_scores": bool(args.aggregate_scores),
                "app_threshold": float(args.app_threshold),
                "min_nkill": int(args.min_nkill),
                "input_chunk_size": int(getattr(args, "input_chunk_size", 10000)),
                "normalized_rows": summary.normalized_rows,
                "predicted_rows": summary.predicted_rows,
                "normalized_input": str(summary.normalized_input_path),
                "predictions_all": str(summary.predictions_all_path),
                "grouped_outputs": summary.grouped_outputs,
                "grouping_mode": summary.grouping_mode,
                "cpu_workers_selected": summary.cpu_workers_selected,
                "prefetch_queue_size_selected": summary.prefetch_queue_size_selected,
                "work_unit_count": summary.work_unit_count,
                "target_rows_per_group": summary.target_rows_per_group,
                "target_bytes_per_group": summary.target_bytes_per_group,
                "prediction_row_budget": int(getattr(args, "prediction_row_budget", 8192)),
                "deterministic_representation": bool(getattr(args, "deterministic_representation", False)),
                "classifier_backend": getattr(args, "classifier_backend", None) or os.environ.get("MOLE_CLASSIFIER_BACKEND", "auto"),
                "profiling": summary.profiling,
                "prediction_runtime": scheduler.runtime_snapshot(),
            }
        else:
            process_summary = await screen_paths_multiprocess(
                input_paths=process_config.input_paths,
                output_dir=str(output_dir),
                execution_mode=process_config.execution_mode,
                producer_processes=process_config.producer_processes,
                predict_queue_max_batches=process_config.predict_queue_max_batches,
                result_queue_max_batches=process_config.result_queue_max_batches,
                batch_checkpoint_size=process_config.batch_checkpoint_size,
                rows_per_shard=process_config.rows_per_shard,
                row_group_size=process_config.row_group_size,
                smiles_colname=getattr(args, "smiles_colname", "smiles"),
                chem_id_colname=getattr(args, "chem_id_colname", "chem_id"),
                archive_smiles_colname=getattr(args, "archive_smiles_colname", "product_smiles_canonical"),
                archive_chem_id_colname=getattr(args, "archive_chem_id_colname", "example_combo_id"),
                grouping_mode=getattr(args, "grouping_mode", "source"),
                cpu_workers=getattr(args, "cpu_workers", "auto"),
                target_rows_per_group=getattr(args, "target_rows_per_group", "auto"),
                target_bytes_per_group=getattr(args, "target_bytes_per_group", "auto"),
                chunk_size=int(getattr(args, "input_chunk_size", 10000)),
                aggregate_scores=bool(getattr(args, "aggregate_scores", True)),
                app_threshold=float(getattr(args, "app_threshold", 0.04374140128493309)),
                min_nkill=int(getattr(args, "min_nkill", 10)),
                num_graph_workers=getattr(args, "num_graph_workers", "auto"),
                graph_batch_size=int(getattr(args, "graph_batch_size", 1024)),
                prefetch_batches=int(getattr(args, "prefetch_batches", 2)),
                classifier_workers=getattr(args, "classifier_workers", "auto"),
                classifier_inflight_batches=getattr(args, "classifier_inflight_batches", "auto"),
                deterministic_representation=bool(getattr(args, "deterministic_representation", False)),
                enable_profiling=bool(getattr(args, "profiling", False)),
            )
            manifest = {
                "input_path": resolved_input_paths[0],
                "input_paths": resolved_input_paths,
                "output_dir": str(output_dir),
                "execution_mode": process_config.execution_mode,
                "producer_processes": process_config.producer_processes,
                "predict_queue_max_batches": process_config.predict_queue_max_batches,
                "result_queue_max_batches": process_config.result_queue_max_batches,
                "batch_checkpoint_size": process_config.batch_checkpoint_size,
                "rows_per_shard": process_config.rows_per_shard,
                "row_group_size": process_config.row_group_size,
                "dedupe_smiles": bool(args.dedupe_smiles),
                "aggregate_scores": bool(args.aggregate_scores),
                "app_threshold": float(args.app_threshold),
                "min_nkill": int(args.min_nkill),
                "input_chunk_size": int(getattr(args, "input_chunk_size", 10000)),
                "prepared_manifest_path": process_summary["prepared_manifest_path"],
                "prepared_manifest_paths": process_summary["prepared_manifest_paths"],
                "prepared_input_paths": process_summary["prepared_input_paths"],
                "batch_manifest_path": process_summary["batch_manifest_path"],
                "run_state_path": process_summary["run_state_path"],
                "hits_dir": process_summary["hits_dir"],
                "source_groups": process_summary["source_groups"],
                "work_unit_count": process_summary["work_unit_count"],
                "prediction_row_budget": None,
                "normalized_rows": None,
                "predicted_rows": None,
                "normalized_input": None,
                "predictions_all": None,
                "grouped_outputs": [],
                "grouping_mode": getattr(args, "grouping_mode", "source"),
                "cpu_workers_selected": None,
                "prefetch_queue_size_selected": None,
                "target_rows_per_group": getattr(args, "target_rows_per_group", "auto"),
                "target_bytes_per_group": getattr(args, "target_bytes_per_group", "auto"),
                "deterministic_representation": bool(getattr(args, "deterministic_representation", False)),
                "classifier_backend": getattr(args, "classifier_backend", None) or os.environ.get("MOLE_CLASSIFIER_BACKEND", "auto"),
                "profiling": None,
                "prediction_runtime": process_summary.get("prediction_runtime"),
                "runtime": process_summary["runtime"],
            }

        panel_file = getattr(args, "panel_file", None)
        if panel_file:
            panel_config = load_panel_config(panel_file)
            panel_info = {
                "panel_file": str(Path(panel_file).expanduser().resolve()),
                "panel_label": panel_config.label,
                "panel_mode": panel_config.mode,
                "pathogen_count": len(panel_config.pathogen_panel),
                "commensal_count": len(panel_config.commensal_sparing_panel),
                "lambda": float(getattr(args, "panel_lambda", None) or panel_config.lambda_),
                "tau": panel_config.tau,
                "threshold": panel_config.app_threshold,
            }
            if process_config.execution_mode == "thread" and args.aggregate_scores:
                panel_result = _apply_panel_scoring_to_predictions(
                    predictions_path=summary.predictions_all_path,
                    panel_config=panel_config,
                    lambda_override=getattr(args, "panel_lambda", None),
                    min_pathogen_hard=getattr(args, "panel_min_pathogen_hard", None),
                    min_selectivity=getattr(args, "panel_min_selectivity", None),
                )
                panel_info.update(panel_result["panel_scoring"])
            manifest["panel"] = panel_info

        (output_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        return manifest

    payload = asyncio_run(_run())
    _dump_json(payload, None)
    return 0


def _command_stream_enumeration_screen(args: argparse.Namespace) -> int:
    async def _run() -> dict[str, Any]:
        _apply_classifier_backend_arg(args)
        scheduler = create_scheduler(
            num_graph_workers=getattr(args, "num_graph_workers", "auto"),
            graph_batch_size=int(getattr(args, "graph_batch_size", 1024)),
            prefetch_batches=int(getattr(args, "prefetch_batches", 2)),
            classifier_workers=getattr(args, "classifier_workers", "auto"),
            classifier_inflight_batches=getattr(args, "classifier_inflight_batches", "auto"),
            deterministic_representation=bool(getattr(args, "deterministic_representation", False)),
        )
        summary = await stream_enumeration_screen(
            output_dir=getattr(args, "output_dir"),
            scaffold_file=getattr(args, "scaffold_file", None),
            scaffold_dir=getattr(args, "scaffold_dir", None),
            scaffold_catalog=getattr(args, "scaffold_catalog", None),
            ordinary_library=getattr(args, "ordinary_library"),
            pos13_library=getattr(args, "pos13_library"),
            run_state_source=getattr(args, "run_state_source", None),
            chunk_manifest_source=getattr(args, "chunk_manifest_source", None),
            start_index=int(getattr(args, "start_index", 0)),
            stop_index=getattr(args, "stop_index", None),
            shard_size=int(getattr(args, "shard_size", 100000)),
            prediction_batch_size=int(getattr(args, "prediction_batch_size", 1024)),
            scheduler=scheduler,
            app_threshold=float(getattr(args, "app_threshold", 0.04374140128493309)),
            min_nkill=int(getattr(args, "min_nkill", 10)),
            classifier_backend=getattr(args, "classifier_backend", None) or os.environ.get("MOLE_CLASSIFIER_BACKEND", "auto"),
            num_graph_workers=getattr(args, "num_graph_workers", "auto"),
            graph_batch_size=int(getattr(args, "graph_batch_size", 1024)),
            prefetch_batches=int(getattr(args, "prefetch_batches", 2)),
            classifier_workers=getattr(args, "classifier_workers", "auto"),
            classifier_inflight_batches=getattr(args, "classifier_inflight_batches", "auto"),
            enumeration_workers=getattr(args, "enumeration_workers", "auto"),
            enumeration_prefetch_batches=getattr(args, "enumeration_prefetch_batches", "auto"),
            deterministic_representation=bool(getattr(args, "deterministic_representation", False)),
            enable_profiling=bool(getattr(args, "profiling", False)),
        )
        return {
            "output_dir": str(summary.output_dir),
            "run_state_path": str(summary.run_state_path),
            "shard_manifest_path": str(summary.shard_manifest_path),
            "attempted_count": summary.attempted_count,
            "hit_count": summary.hit_count,
            "completed_shards": summary.completed_shards,
            "start_index": summary.start_index,
            "stop_index": summary.stop_index,
            "total_combinations": summary.total_combinations,
            "prediction_runtime": scheduler.runtime_snapshot(),
        }

    payload = asyncio_run(_run())
    _dump_json(payload, None)
    return 0


def _command_score(args: argparse.Namespace) -> int:
    objective_dict = load_objective_spec(args.objective_file)
    site_reward = objective_dict.get("site_reward")
    if site_reward and site_reward.get("enabled") and not site_reward.get("scaffold_smiles"):
        scaffold_file = Path(args.scaffold_file).expanduser().resolve()
        if not scaffold_file.is_file():
            raise ValueError(
                "site_reward.enabled=true requires --scaffold-file or site_reward.scaffold_smiles"
            )
        objective_dict = {
            **objective_dict,
            "site_reward": {
                **site_reward,
                "scaffold_smiles": scaffold_from_file(scaffold_file),
            },
        }
    objective = ScoreObjective.model_validate(objective_dict)
    request = ReinventScoreRequest(
        smiles=list(args.smiles),
        chem_id=list(args.chem_ids) if args.chem_ids is not None else None,
        objective=objective,
        app_threshold=float(args.app_threshold),
        min_nkill=int(args.min_nkill),
    )

    async def _run() -> dict[str, Any]:
        scheduler = get_scheduler()
        items = await scheduler.predict_molecules(
            molecules=request.to_molecule_input().molecules,
            aggregate_scores=False,
            app_threshold=request.app_threshold,
            min_nkill=request.min_nkill,
        )
        scored_items = score_reinvent_predictions(items, request)
        return {
            "mode": "reinvent_score",
            "objective": {
                "mode": request.objective.mode,
                "strains": request.objective.resolved_strains(),
                "app_threshold": request.app_threshold,
                "min_nkill": request.min_nkill,
                "tau": request.objective.tau,
                "site_reward": (
                    request.objective.site_reward.as_payload()
                    if request.objective.site_reward is not None
                    else None
                ),
            },
            "items": scored_items,
        }

    payload = asyncio_run(_run())
    _dump_json(payload, args.output)
    return 0


def _command_optimize(args: argparse.Namespace) -> int:
    command = [
        sys.executable,
        str(Path(args.script).expanduser().resolve()),
        "--template",
        str(resolve_path(args.template)),
        "--env-file",
        str(resolve_path(args.env_file)),
        "--experiment-spec",
        str(resolve_path(args.experiment_spec)),
        "--scaffold-file",
        str(resolve_path(args.scaffold_file)),
        "--output-root",
        str(resolve_path(args.output_root)),
    ]
    if args.run_id:
        command += ["--run-id", args.run_id]
    if args.dry_run:
        command.append("--dry-run")

    completed = subprocess.run(command, check=False)
    return int(completed.returncode)


def asyncio_run(coro):
    return asyncio.run(coro)


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "doctor":
        return _command_doctor(args)
    if args.command == "embed":
        return _command_embed(args)
    if args.command == "predict":
        return _command_predict(args)
    if args.command == "preprocess-screening-input":
        return _command_preprocess_screening_input(args)
    if args.command == "benchmark-screening-inputs":
        return _command_benchmark_screening_inputs(args)
    if args.command == "screen":
        return _command_screen(args)
    if args.command == "stream-enumeration-screen":
        return _command_stream_enumeration_screen(args)
    if args.command == "score":
        return _command_score(args)
    if args.command == "optimize":
        return _command_optimize(args)
    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
