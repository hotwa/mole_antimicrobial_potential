#!/usr/bin/env python3
"""Validate inputs, render a REINVENT4 config, and optionally launch RL."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.reinvent4_workflow import (  # noqa: E402
    build_reinvent_context,
    build_reinvent_command,
    create_run_directory,
    http_health_check,
    load_objective_spec,
    load_runtime_settings,
    render_template_file,
    run_reinvent_command,
    resolve_path,
    validate_scaffold_file,
    write_json_file,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True, help="Path to a REINVENT4 TOML template")
    parser.add_argument("--env-file", required=True, help="Path to local.env style settings")
    parser.add_argument("--objective-file", required=True, help="Objective JSON file")
    parser.add_argument("--scaffold-file", required=True, help="LibInvent scaffold .smi file")
    parser.add_argument("--experiment-name", required=True, help="Human-readable experiment name")
    parser.add_argument(
        "--output-root",
        default=str(REPO_ROOT / "workflows" / "reinvent4" / "results"),
        help="Root directory for rendered configs and outputs",
    )
    parser.add_argument("--run-id", default="", help="Optional run identifier")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and render everything, but do not launch REINVENT4",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    settings = load_runtime_settings(args.env_file)
    objective = load_objective_spec(args.objective_file)
    scaffold_records = validate_scaffold_file(args.scaffold_file)

    health = http_health_check(
        settings["api_health_url"],
        timeout_seconds=settings["request_timeout_seconds"],
    )

    run_dir = create_run_directory(
        args.output_root,
        experiment_name=args.experiment_name,
        run_id=args.run_id or None,
    )
    config_path = run_dir / "reinvent.toml"
    context = build_reinvent_context(
        settings=settings,
        objective=objective,
        scaffold_path=args.scaffold_file,
        objective_path=args.objective_file,
        run_dir=run_dir,
    )
    rendered = render_template_file(args.template, context)
    config_path.write_text(rendered, encoding="utf-8")
    write_json_file(run_dir / "objective.normalized.json", objective)
    write_json_file(run_dir / "scaffold.validation.json", scaffold_records)

    command = build_reinvent_command(settings, config_path)
    summary = {
        "run_dir": str(run_dir),
        "reinvent_root": str(settings["reinvent_root"]),
        "command": command,
        "rendered_config": str(config_path),
        "api_health": health,
        "scaffold_file": str(resolve_path(args.scaffold_file)),
        "objective_file": str(resolve_path(args.objective_file)),
        "dry_run": args.dry_run,
    }
    write_json_file(run_dir / "run.summary.json", summary)
    sys.stdout.write(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    run_reinvent_command(settings, config_path, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
