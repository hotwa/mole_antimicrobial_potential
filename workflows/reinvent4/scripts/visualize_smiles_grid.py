#!/usr/bin/env python3
"""Render molecules from a CSV/TSV/SMI file into a grid image."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List

from rdkit import Chem
from rdkit.Chem import Draw


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-file", required=True, help="Source file: .csv, .tsv, or .smi")
    parser.add_argument("--output-file", required=True, help="Output image path, e.g. molecules.png")
    parser.add_argument("--smiles-column", default="SMILES", help="SMILES column for CSV/TSV inputs")
    parser.add_argument(
        "--legend-columns",
        nargs="*",
        default=[],
        help="Optional columns to include in each molecule legend",
    )
    parser.add_argument("--top-n", type=int, default=24, help="Maximum number of molecules to render")
    parser.add_argument("--mols-per-row", type=int, default=4, help="Molecules per row")
    parser.add_argument("--sub-img-size", type=int, nargs=2, default=[350, 250], help="Sub-image size")
    return parser.parse_args()


def load_rows(path: Path, smiles_column: str) -> List[dict]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        with path.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))
    if suffix == ".tsv":
        with path.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    if suffix == ".smi":
        rows = []
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                rows.append({smiles_column: parts[0], "line_number": str(line_number)})
        return rows
    raise ValueError(f"Unsupported input file type: {path}")


def main() -> int:
    args = parse_args()
    input_path = Path(args.input_file).expanduser().resolve()
    output_path = Path(args.output_file).expanduser().resolve()

    rows = load_rows(input_path, args.smiles_column)
    mols = []
    legends = []
    for row in rows:
        if len(mols) >= args.top_n:
            break
        smiles = row.get(args.smiles_column, "").strip()
        if not smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        mols.append(mol)
        legend_parts = [smiles]
        for column in args.legend_columns:
            if column in row and row[column] not in ("", None):
                legend_parts.append(f"{column}={row[column]}")
        legends.append("\n".join(legend_parts))

    if not mols:
        raise ValueError(f"No valid molecules found in {input_path}")

    image = Draw.MolsToGridImage(
        mols,
        molsPerRow=args.mols_per_row,
        subImgSize=tuple(args.sub_img_size),
        legends=legends,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    image.save(str(output_path))
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
