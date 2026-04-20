#!/usr/bin/env python3
"""Run chunked long-horizon REINVENT4 RL with external early-stop monitoring."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.reinvent4_workflow import (  # noqa: E402
    build_reinvent_command,
    build_reinvent_context,
    build_score_request,
    chunk_reward_summary,
    create_run_directory,
    extract_top_unique_smiles,
    http_health_check,
    http_json_request,
    load_long_run_spec,
    load_runtime_settings,
    read_reinvent_csv,
    render_template_file,
    resolve_path,
    run_reinvent_command,
    unique_murcko_count,
    validate_scaffold_file,
    write_json_file,
)

REINVENT_SIMPLE_MAX_SCORE = 0.999
WARNING_WINDOW = 3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True, help="Path to the REINVENT4 TOML template")
    parser.add_argument("--env-file", required=True, help="Path to local.env style runtime settings")
    parser.add_argument("--experiment-spec", required=True, help="Path to a long-run RL spec JSON file")
    parser.add_argument("--scaffold-file", required=True, help="Path to the LibInvent scaffold .smi file")
    parser.add_argument(
        "--output-root",
        default=str(REPO_ROOT / "workflows" / "reinvent4" / "results" / "runs"),
        help="Root directory for long-run outputs",
    )
    parser.add_argument("--run-id", default="", help="Optional deterministic run identifier")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Render only the first chunk and manifest without launching REINVENT4",
    )
    return parser.parse_args()


def _recent_chunks(chunks: Sequence[Dict[str, Any]], window: int) -> List[Dict[str, Any]]:
    if len(chunks) < window:
        return []
    return list(chunks[-window:])


def _consecutive_success(
    chunks: Sequence[Dict[str, Any]],
    min_chunks: int,
    success_consecutive_chunks: int,
) -> bool:
    if len(chunks) < max(min_chunks, success_consecutive_chunks):
        return False
    recent = _recent_chunks(chunks, success_consecutive_chunks)
    return bool(recent) and all(bool(chunk.get("success_chunk")) for chunk in recent)


def _warning_state(chunks: Sequence[Dict[str, Any]]) -> str:
    recent = _recent_chunks(chunks, WARNING_WINDOW)
    if not recent:
        return "ok"
    if all(chunk.get("unique_smiles_ratio", 1.0) < 0.40 for chunk in recent):
        return "warning_diversity"
    if all(chunk.get("zero_score_fraction", 0.0) > 0.50 for chunk in recent):
        return "warning_diversity"
    commensal_values = [chunk.get("commensal_panel_mean_probability") for chunk in recent]
    if all(value is not None and value > 0.20 for value in commensal_values):
        return "warning_commensal_damage"
    return "ok"


def _monitor_commensal_panel(
    score_url: str,
    smiles_list: Sequence[str],
    objective: Dict[str, Any],
    monitor_strains: Sequence[str],
    timeout_seconds: int,
) -> Dict[str, Any]:
    if not smiles_list or not monitor_strains:
        return {
            "commensal_panel_mean_probability": None,
            "commensal_panel_max_probability": None,
            "commensal_monitor_count": 0,
        }
    monitor_objective = {
        "mode": "single_strain",
        "strains": list(monitor_strains),
        "tau": objective["tau"],
        "app_threshold": objective["app_threshold"],
        "min_nkill": objective["min_nkill"],
        "label": "commensal_monitor",
    }
    response = http_json_request(
        score_url,
        build_score_request(smiles_list, monitor_objective),
        timeout_seconds=timeout_seconds,
    )
    panel_means: List[float] = []
    panel_maxes: List[float] = []
    for item in response.get("items", []):
        selected = item.get("selected_probabilities", {})
        if not selected:
            continue
        values = [float(probability) for probability in selected.values()]
        panel_means.append(sum(values) / len(values))
        panel_maxes.append(max(values))
    if not panel_means:
        return {
            "commensal_panel_mean_probability": None,
            "commensal_panel_max_probability": None,
            "commensal_monitor_count": 0,
        }
    return {
        "commensal_panel_mean_probability": sum(panel_means) / len(panel_means),
        "commensal_panel_max_probability": max(panel_maxes),
        "commensal_monitor_count": len(panel_means),
    }


def _chunk_run_summary(
    *,
    chunk_index: int,
    total_steps_completed: int,
    chunk_dir: Path,
    settings: Dict[str, Any],
    command: Sequence[str],
    config_path: Path,
    objective_path: Path,
    scaffold_path: Path,
    api_health: Dict[str, Any],
    dry_run: bool,
) -> Dict[str, Any]:
    return {
        "chunk_index": chunk_index,
        "total_steps_completed": total_steps_completed,
        "run_dir": str(chunk_dir),
        "reinvent_root": str(settings["reinvent_root"]),
        "command": list(command),
        "rendered_config": str(config_path),
        "api_health": api_health,
        "scaffold_file": str(scaffold_path),
        "objective_file": str(objective_path),
        "dry_run": dry_run,
        "agent_file": str(settings["agent_file"]),
    }


def main() -> int:
    args = parse_args()
    settings = load_runtime_settings(args.env_file)
    spec = load_long_run_spec(args.experiment_spec)
    scaffold_path = resolve_path(args.scaffold_file)
    template_path = resolve_path(args.template)
    objective = spec["objective"]
    scaffold_records = validate_scaffold_file(scaffold_path)
    health = http_health_check(
        settings["api_health_url"],
        timeout_seconds=settings["request_timeout_seconds"],
    )

    run_dir = create_run_directory(
        args.output_root,
        experiment_name=spec["experiment_name"],
        run_id=args.run_id or None,
    )
    objective_path = run_dir / "objective.normalized.json"
    spec_path = run_dir / "long_run.spec.normalized.json"
    manifest_path = run_dir / "run_manifest.json"
    write_json_file(objective_path, objective)
    write_json_file(spec_path, spec)
    write_json_file(run_dir / "scaffold.validation.json", scaffold_records)

    current_agent_file = settings["agent_file"]
    chunk_entries: List[Dict[str, Any]] = []
    stop_reason = "max_total_steps_reached"

    manifest: Dict[str, Any] = {
        "experiment_name": spec["experiment_name"],
        "run_dir": str(run_dir),
        "template": str(template_path),
        "scaffold_file": str(scaffold_path),
        "objective_file": str(objective_path),
        "normalized_spec_file": str(spec_path),
        "api_health": health,
        "dry_run": args.dry_run,
        "chunks": chunk_entries,
        "stop_reason": None,
        "status": "ok",
    }
    write_json_file(manifest_path, manifest)

    for chunk_index in range(1, spec["max_chunks"] + 1):
        total_steps_completed = chunk_index * spec["chunk_steps"]
        chunk_dir = run_dir / f"chunk_{chunk_index:03d}"
        chunk_dir.mkdir(parents=True, exist_ok=True)

        chunk_settings = dict(settings)
        chunk_settings["agent_file"] = resolve_path(current_agent_file)
        chunk_settings["min_steps"] = spec["chunk_steps"]
        chunk_settings["max_steps"] = spec["chunk_steps"]
        chunk_settings["max_score"] = REINVENT_SIMPLE_MAX_SCORE

        config_path = chunk_dir / "reinvent.toml"
        context = build_reinvent_context(
            settings=chunk_settings,
            objective=objective,
            scaffold_path=scaffold_path,
            objective_path=objective_path,
            run_dir=chunk_dir,
        )
        config_path.write_text(render_template_file(template_path, context), encoding="utf-8")

        command = build_reinvent_command(chunk_settings, config_path)
        run_summary = _chunk_run_summary(
            chunk_index=chunk_index,
            total_steps_completed=total_steps_completed,
            chunk_dir=chunk_dir,
            settings=chunk_settings,
            command=command,
            config_path=config_path,
            objective_path=objective_path,
            scaffold_path=scaffold_path,
            api_health=health,
            dry_run=args.dry_run,
        )
        write_json_file(chunk_dir / "run.summary.json", run_summary)

        if args.dry_run:
            stop_reason = "dry_run_rendered_first_chunk"
            chunk_entries.append(
                {
                    **run_summary,
                    "success_chunk": False,
                    "status": "dry_run",
                }
            )
            break

        run_reinvent_command(chunk_settings, config_path, dry_run=False)
        csv_path = chunk_dir / "rl_run_1.csv"
        rows = read_reinvent_csv(csv_path)
        reward_summary = chunk_reward_summary(rows)
        selected = extract_top_unique_smiles(rows, top_n=spec["monitor_top_n"])
        selected_smiles = [item["smiles"] for item in selected]

        chunk_entry: Dict[str, Any] = {
            **run_summary,
            **reward_summary,
            "selected_top_n": spec["monitor_top_n"],
            "selected_unique_molecules": len(selected),
            "unique_murcko_count": unique_murcko_count(selected_smiles),
            "checkpoint_file": str(chunk_dir / "stage1.chkpt"),
            "success_chunk": (
                reward_summary["chunk_max_reward"] >= spec["success_high"]
                and reward_summary["chunk_top5_mean_reward"] >= spec["success_top5_mean"]
            ),
        }
        if spec["commensal_monitor_strains"]:
            chunk_entry.update(
                _monitor_commensal_panel(
                    score_url=settings["score_url"],
                    smiles_list=selected_smiles,
                    objective=objective,
                    monitor_strains=spec["commensal_monitor_strains"],
                    timeout_seconds=settings["request_timeout_seconds"],
                )
            )
        else:
            chunk_entry.update(
                {
                    "commensal_panel_mean_probability": None,
                    "commensal_panel_max_probability": None,
                    "commensal_monitor_count": 0,
                }
            )

        chunk_entries.append(chunk_entry)
        chunk_entry["status"] = _warning_state(chunk_entries)
        write_json_file(chunk_dir / "chunk.summary.json", chunk_entry)

        manifest["status"] = chunk_entry["status"]
        manifest["stop_reason"] = None
        write_json_file(manifest_path, manifest)

        checkpoint_file = chunk_dir / "stage1.chkpt"
        if not checkpoint_file.is_file():
            raise RuntimeError(f"Expected checkpoint file was not created: {checkpoint_file}")
        current_agent_file = checkpoint_file

        if _consecutive_success(
            chunk_entries,
            min_chunks=spec["min_chunks"],
            success_consecutive_chunks=spec["success_consecutive_chunks"],
        ):
            stop_reason = "success_threshold_reached"
            break

    manifest["stop_reason"] = stop_reason
    manifest["status"] = _warning_state(chunk_entries)
    write_json_file(manifest_path, manifest)
    sys.stdout.write(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
