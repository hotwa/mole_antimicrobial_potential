#!/usr/bin/env python3
"""Analyze per-site decoration sizes from a REINVENT4 RL CSV."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import mean, median

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-file", required=True, help="REINVENT4 rl_run_*.csv file")
    parser.add_argument(
        "--output-dir",
        default="",
        help="Directory for analysis outputs. Defaults to the input file parent.",
    )
    parser.add_argument(
        "--rgroup-column",
        default="R-groups",
        help="Column containing LibInvent decorations separated by '|'.",
    )
    parser.add_argument("--score-column", default="Score", help="Column used for ranking molecules.")
    parser.add_argument(
        "--reward-column",
        default="Antimicrobial Reward",
        help="Column containing the raw antimicrobial reward.",
    )
    parser.add_argument(
        "--range-min",
        type=int,
        default=4,
        help="Reference lower bound for per-site heavy atom counts.",
    )
    parser.add_argument(
        "--range-max",
        type=int,
        default=27,
        help="Reference upper bound for per-site heavy atom counts.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="How many top-scoring molecules to summarize explicitly.",
    )
    return parser.parse_args()


def heavy_atom_count(fragment_smiles: str) -> int | None:
    mol = Chem.MolFromSmiles(fragment_smiles)
    if mol is None:
        return None
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


def parse_rgroups(raw_value: str) -> list[str]:
    return [part.strip() for part in raw_value.split("|") if part.strip()]


def analyze_row(
    row: dict[str, str],
    *,
    rgroup_column: str,
    score_column: str,
    reward_column: str,
    range_min: int,
    range_max: int,
) -> dict[str, object] | None:
    try:
        score = float(row[score_column])
    except (KeyError, TypeError, ValueError):
        return None

    reward_raw = row.get(reward_column, "")
    try:
        reward = float(reward_raw)
    except (TypeError, ValueError):
        reward = None

    fragments = parse_rgroups(row.get(rgroup_column, ""))
    site_counts = [heavy_atom_count(fragment) for fragment in fragments]
    while len(site_counts) < 4:
        site_counts.append(None)

    in_range = [
        count is not None and range_min <= count <= range_max for count in site_counts[:4]
    ]
    ge_min = [count is not None and count >= range_min for count in site_counts[:4]]

    result: dict[str, object] = {
        "step": row.get("step", ""),
        "score": score,
        "antimicrobial_reward": reward,
        "smiles": row.get("SMILES", ""),
        "input_scaffold": row.get("Input_Scaffold", ""),
        "scaffold": row.get("Scaffold", ""),
        "r_groups": row.get(rgroup_column, ""),
        "site_1_heavy_atoms": site_counts[0],
        "site_2_heavy_atoms": site_counts[1],
        "site_3_heavy_atoms": site_counts[2],
        "site_4_heavy_atoms": site_counts[3],
        "site_1_in_reference_range": int(in_range[0]),
        "site_2_in_reference_range": int(in_range[1]),
        "site_3_in_reference_range": int(in_range[2]),
        "site_4_in_reference_range": int(in_range[3]),
        "site_1_ge_min": int(ge_min[0]),
        "site_2_ge_min": int(ge_min[1]),
        "site_3_ge_min": int(ge_min[2]),
        "site_4_ge_min": int(ge_min[3]),
        "sites_in_reference_range_count": sum(in_range),
        "sites_ge_min_count": sum(ge_min),
        "all_four_sites_in_reference_range": int(all(in_range)),
        "all_four_sites_ge_min": int(all(ge_min)),
        "invalid_site_count": sum(count is None for count in site_counts[:4]),
    }
    return result


def summarize_site(values: list[int]) -> dict[str, float | int | None]:
    if not values:
        return {
            "count": 0,
            "min": None,
            "max": None,
            "mean": None,
            "median": None,
        }
    return {
        "count": len(values),
        "min": min(values),
        "max": max(values),
        "mean": mean(values),
        "median": median(values),
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input_file).expanduser().resolve()
    output_dir = (
        Path(args.output_dir).expanduser().resolve() if args.output_dir else input_path.parent.resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    analyzed_rows: list[dict[str, object]] = []
    with input_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            analyzed = analyze_row(
                row,
                rgroup_column=args.rgroup_column,
                score_column=args.score_column,
                reward_column=args.reward_column,
                range_min=args.range_min,
                range_max=args.range_max,
            )
            if analyzed is not None:
                analyzed_rows.append(analyzed)

    analyzed_rows.sort(key=lambda row: float(row["score"]), reverse=True)

    site_value_lists = {
        f"site_{index}_heavy_atoms": [
            int(row[f"site_{index}_heavy_atoms"])
            for row in analyzed_rows
            if row[f"site_{index}_heavy_atoms"] is not None
        ]
        for index in range(1, 5)
    }

    summary = {
        "input_file": str(input_path),
        "range_min": args.range_min,
        "range_max": args.range_max,
        "molecule_count": len(analyzed_rows),
        "all_four_sites_in_reference_range_fraction": (
            sum(int(row["all_four_sites_in_reference_range"]) for row in analyzed_rows) / len(analyzed_rows)
            if analyzed_rows
            else 0.0
        ),
        "all_four_sites_ge_min_fraction": (
            sum(int(row["all_four_sites_ge_min"]) for row in analyzed_rows) / len(analyzed_rows)
            if analyzed_rows
            else 0.0
        ),
        "invalid_site_count_histogram": {
            str(count): sum(1 for row in analyzed_rows if int(row["invalid_site_count"]) == count)
            for count in range(0, 5)
        },
        "expanded_site_count_histogram": {
            str(count): sum(1 for row in analyzed_rows if int(row["sites_ge_min_count"]) == count)
            for count in range(0, 5)
        },
        "in_range_site_count_histogram": {
            str(count): sum(
                1 for row in analyzed_rows if int(row["sites_in_reference_range_count"]) == count
            )
            for count in range(0, 5)
        },
        "site_statistics": {
            site_name: summarize_site(values) for site_name, values in site_value_lists.items()
        },
        "top_molecules": [
            {
                "step": row["step"],
                "score": row["score"],
                "antimicrobial_reward": row["antimicrobial_reward"],
                "smiles": row["smiles"],
                "site_heavy_atoms": [
                    row["site_1_heavy_atoms"],
                    row["site_2_heavy_atoms"],
                    row["site_3_heavy_atoms"],
                    row["site_4_heavy_atoms"],
                ],
                "sites_ge_min_count": row["sites_ge_min_count"],
                "sites_in_reference_range_count": row["sites_in_reference_range_count"],
            }
            for row in analyzed_rows[: args.top_n]
        ],
    }

    analysis_tsv = output_dir / "rgroup_site_analysis.tsv"
    summary_json = output_dir / "rgroup_site_summary.json"
    top_tsv = output_dir / "rgroup_site_top_molecules.tsv"

    write_tsv(analysis_tsv, analyzed_rows)
    write_tsv(top_tsv, analyzed_rows[: args.top_n])
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(analysis_tsv)
    print(summary_json)
    print(top_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
