#!/usr/bin/env python3
"""Re-score per-strain predictions with a pathogen-selective panel.

Reads a CSV containing per-strain prediction probabilities (with a
``pred_id`` column formatted as ``chem_id:strain_name``) and computes
pathogen/commensal panel scores defined by a panel configuration file.

Usage::

    python scripts/pathogen_panel_rescore.py \\
        --panel-file workflows/reinvent4/inputs/objectives/pathogen_selective.panel.v1.json \\
        --input predictions_per_strain.csv \\
        --output predictions_scored.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402

from src.panel_scoring import load_panel_config, compute_panel_scores_from_dataframe  # noqa: E402


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Re-score per-strain predictions with a pathogen-selective panel.",
    )
    parser.add_argument(
        "--panel-file",
        required=True,
        help="Path to the panel configuration JSON file.",
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input CSV with pred_id and per-strain probability columns.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output CSV path. Defaults to stdout if omitted.",
    )
    parser.add_argument(
        "--probability-column",
        default="1",
        help="Name of the probability column in the input CSV (default: '1').",
    )
    parser.add_argument(
        "--lambda",
        dest="lambda_",
        type=float,
        default=None,
        help="Override lambda from panel config (commensal penalty weight).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Override app_threshold from panel config.",
    )
    parser.add_argument(
        "--tau",
        type=float,
        default=None,
        help="Override tau from panel config (sigmoid temperature).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    panel_path = Path(args.panel_file)
    if not panel_path.is_file():
        print(f"error: panel file not found: {panel_path}", file=sys.stderr)
        return 1

    panel = load_panel_config(panel_path)

    input_path = Path(args.input)
    if not input_path.is_file():
        print(f"error: input file not found: {input_path}", file=sys.stderr)
        return 1

    df = pd.read_csv(input_path)
    if "pred_id" not in df.columns:
        print(
            f"error: input CSV missing 'pred_id' column. Found: {list(df.columns)}",
            file=sys.stderr,
        )
        return 1

    scores_df = compute_panel_scores_from_dataframe(
        df,
        panel,
        pred_id_col="pred_id",
        probability_col=args.probability_column,
        threshold=args.threshold,
        tau=args.tau,
        lambda_=args.lambda_,
    )

    scores_df = scores_df.reset_index()
    output = df.merge(scores_df, on="chem_id", how="left")

    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output.to_csv(output_path, index=False)
    else:
        output.to_csv(sys.stdout, index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main())
