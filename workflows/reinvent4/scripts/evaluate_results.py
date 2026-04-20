#!/usr/bin/env python3
"""Evaluate a REINVENT4 RL CSV using the local /score and /predict APIs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.reinvent4_workflow import (  # noqa: E402
    build_score_request,
    extract_top_unique_smiles,
    http_json_request,
    load_objective_spec,
    load_runtime_settings,
    read_reinvent_csv,
    summarize_predictions,
    write_json_file,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", required=True, help="Path to local.env style settings")
    parser.add_argument("--objective-file", required=True, help="Objective JSON file")
    parser.add_argument("--reinvent-csv", required=True, help="RL output CSV to analyze")
    parser.add_argument(
        "--output-dir",
        default="",
        help="Directory for evaluation artifacts (defaults next to the CSV)",
    )
    parser.add_argument("--top-n", type=int, default=50, help="Number of unique molecules to score")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    settings = load_runtime_settings(args.env_file)
    objective = load_objective_spec(args.objective_file)
    rows = read_reinvent_csv(args.reinvent_csv)
    selected = extract_top_unique_smiles(rows, top_n=args.top_n)

    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else Path(args.reinvent_csv).expanduser().resolve().parent / "evaluation"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    smiles = [item["smiles"] for item in selected]
    score_request = build_score_request(smiles, objective)
    score_response = http_json_request(
        settings["score_url"],
        score_request,
        timeout_seconds=settings["request_timeout_seconds"],
    )
    predict_base_url = settings["score_url"][: -len("/score")] if settings["score_url"].endswith("/score") else settings["score_url"].rstrip("/")
    predict_per_strain = http_json_request(
        f"{predict_base_url}/predict",
        {
            "smiles": smiles,
            "chem_id": score_request["chem_id"],
            "aggregate_scores": False,
            "app_threshold": objective["app_threshold"],
            "min_nkill": objective["min_nkill"],
        },
        timeout_seconds=settings["request_timeout_seconds"],
    )
    predict_aggregate = http_json_request(
        f"{predict_base_url}/predict",
        {
            "smiles": smiles,
            "chem_id": score_request["chem_id"],
            "aggregate_scores": True,
            "app_threshold": objective["app_threshold"],
            "min_nkill": objective["min_nkill"],
        },
        timeout_seconds=settings["request_timeout_seconds"],
    )

    summary_rows = summarize_predictions(selected, score_response, predict_aggregate)
    (output_dir / "top_smiles.smi").write_text(
        "\n".join(item["smiles"] for item in selected) + "\n",
        encoding="utf-8",
    )
    write_json_file(output_dir / "top_smiles.json", selected)
    write_json_file(output_dir / "score_response.json", score_response)
    write_json_file(output_dir / "predict_per_strain.json", predict_per_strain)
    write_json_file(output_dir / "predict_aggregate.json", predict_aggregate)
    write_tsv(output_dir / "summary.tsv", summary_rows)

    sys.stdout.write(
        json.dumps(
            {
                "output_dir": str(output_dir),
                "selected_count": len(selected),
                "summary_tsv": str(output_dir / "summary.tsv"),
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
