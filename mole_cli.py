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

from mole_representation import process_representation
from src.models import MoleculeInput, ReinventScoreRequest, ScoreObjective
from src.reinvent4_workflow import load_objective_spec, resolve_path
from src.reinvent_scoring import score_reinvent_predictions
from src.service import get_predictor
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
    embed.add_argument("--output", help="Optional output file path")
    embed.add_argument("--format", choices=["json", "tsv"], default="json", help="Output format")

    predict = subparsers.add_parser("predict", help="Predict strain-level antimicrobial activity")
    predict.add_argument("--smiles", nargs="+", action="extend", required=True, help="One or more SMILES strings")
    predict.add_argument("--chem-id", nargs="+", action="extend", dest="chem_ids", help="Optional chemical identifiers")
    predict.add_argument("--aggregate-scores", action="store_true", help="Return aggregate antimicrobial scores")
    predict.add_argument("--app-threshold", type=float, default=0.04374140128493309, help="Growth inhibition threshold")
    predict.add_argument("--min-nkill", type=int, default=10, help="Broad-spectrum threshold")
    predict.add_argument("--output", help="Optional output file path")

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
    ).reset_index(drop=False)
    _write_embedding_output(representation, args.output, args.format)
    return 0


async def _predict_async(args: argparse.Namespace) -> dict[str, Any]:
    predictor = get_predictor()
    await predictor.ensure_loaded()
    normalized = MoleculeInput(
        smiles=list(args.smiles),
        chem_id=list(args.chem_ids) if args.chem_ids is not None else None,
        aggregate_scores=bool(args.aggregate_scores),
        app_threshold=float(args.app_threshold),
        min_nkill=int(args.min_nkill),
    ).normalize()
    items = await predictor.predict(normalized)
    return {
        "mode": "aggregate" if normalized.aggregate_scores else "per_strain",
        "items": items,
    }


def _command_predict(args: argparse.Namespace) -> int:
    payload = asyncio_run(_predict_async(args))
    _dump_json(payload, args.output)
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
        predictor = get_predictor()
        await predictor.ensure_loaded()
        raw_items = await predictor.predict(request.to_molecule_input())
        scored_items = score_reinvent_predictions(raw_items, request)
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
    if args.command == "score":
        return _command_score(args)
    if args.command == "optimize":
        return _command_optimize(args)
    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
