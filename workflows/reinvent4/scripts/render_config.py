#!/usr/bin/env python3
"""Render a REINVENT4 config from a local template and objective spec."""

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
    load_objective_spec,
    load_runtime_settings,
    render_template_file,
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
    parser.add_argument("--output-file", required=True, help="Rendered TOML output file")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    settings = load_runtime_settings(args.env_file)
    objective = load_objective_spec(args.objective_file)
    validate_scaffold_file(args.scaffold_file)

    output_file = Path(args.output_file).expanduser().resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)

    context = build_reinvent_context(
        settings=settings,
        objective=objective,
        scaffold_path=args.scaffold_file,
        objective_path=args.objective_file,
        run_dir=output_file.parent,
    )
    rendered = render_template_file(args.template, context)
    output_file.write_text(rendered, encoding="utf-8")

    normalized_objective_path = output_file.parent / "objective.normalized.json"
    write_json_file(normalized_objective_path, objective)

    summary = {
        "template": str(resolve_path(args.template)),
        "rendered_config": str(output_file),
        "objective_file": str(resolve_path(args.objective_file)),
        "normalized_objective_file": str(normalized_objective_path),
        "scaffold_file": str(resolve_path(args.scaffold_file)),
        "score_url": settings["score_url"],
    }
    sys.stdout.write(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
