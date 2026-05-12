#!/usr/bin/env python3
"""
Run MolE antimicrobial predictions on reductive amination products.

Usage:
    pixi run python scripts/predict_reductive_amination_mole.py
"""

import asyncio
import sys
from pathlib import Path

import pandas as pd

# Add project root to path
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.batch_screening import screen_frame
from src.service import get_scheduler


INPUT_DIR = REPO_ROOT / "data" / "11.macrolide_reductive_amination"
OUTPUT_DIR = INPUT_DIR

FILES = [
    "tylosin_seed_amine_products.csv",
    "tylvalosin_seed_amine_products.csv",
]

SMILES_COL = "canonical_product_smiles"
CHEM_ID_COL = "product_id"

# MolE prediction parameters (same as esterification)
APP_THRESHOLD = 0.04374140128493309
MIN_NKILL = 10


def load_and_prepare(filepath: Path) -> pd.DataFrame:
    """Load CSV and prepare for MolE screening."""
    df = pd.read_csv(filepath, encoding="gb18030")
    # Rename columns to match MolE expected format
    df = df.rename(columns={
        SMILES_COL: "smiles",
        CHEM_ID_COL: "chem_id",
    })
    # Add input_order for stable sorting
    df["input_order"] = range(len(df))
    return df


async def predict_file(filepath: Path, output_path: Path) -> dict:
    """Run MolE prediction on a single file."""
    print(f"\nProcessing: {filepath.name}")

    # Load data
    df = load_and_prepare(filepath)
    print(f"  Loaded {len(df)} rows")

    # Run prediction
    scheduler = get_scheduler()
    result = await screen_frame(
        frame=df,
        aggregate_scores=True,
        app_threshold=APP_THRESHOLD,
        min_nkill=MIN_NKILL,
        scheduler=scheduler,
    )

    # Restore original column names
    result = result.rename(columns={
        "smiles": SMILES_COL,
        "chem_id": CHEM_ID_COL,
    })

    # Select and order output columns
    # Original columns from input (handle both origin_ and original_ variants)
    original_cols = [
        "product_id", "macrolide_name", "macrolide_input_smiles",
        "seed_metabolite_id", "seed_metabolite_name", "seed_canonical_neutral_smiles",
        "seed_n_primary_amine", "reaction_smarts", "aldehyde_match_atoms",
        "amine_match_atoms", "new_cn_bond_product_atoms", "product_smiles",
        "canonical_product_smiles", "product_mw", "product_heavy_atom_count",
        "reaction_status", "notes", "image_path",
        "pMIC_pred（μmol/ml）",
    ]
    # Add the origin/original pMIC column (name varies between files)
    for col in ["origin_pMIC_pred（μmol/ml）", "original_pMIC_pred（μmol/ml）"]:
        if col in result.columns:
            original_cols.append(col)
            break
    # MolE prediction columns
    mole_cols = [
        "apscore_total", "apscore_gnegative", "apscore_gpositive",
        "ginhib_total", "ginhib_gnegative", "ginhib_gpositive",
        "broad_spectrum",
    ]

    # Ensure all columns exist
    output_cols = []
    for col in original_cols + mole_cols:
        if col in result.columns:
            output_cols.append(col)
        else:
            print(f"  Warning: Column '{col}' not found in result")

    result = result[output_cols]

    # Save output
    result.to_csv(output_path, index=False, encoding="gb18030")
    print(f"  Saved to: {output_path.name}")
    print(f"  Shape: {result.shape}")

    # Summary stats
    stats = {
        "rows": len(result),
        "broad_spectrum_count": int(result["broad_spectrum"].sum()) if "broad_spectrum" in result.columns else 0,
        "apscore_total_min": float(result["apscore_total"].min()) if "apscore_total" in result.columns else None,
        "apscore_total_max": float(result["apscore_total"].max()) if "apscore_total" in result.columns else None,
        "ginhib_total_min": int(result["ginhib_total"].min()) if "ginhib_total" in result.columns else None,
        "ginhib_total_max": int(result["ginhib_total"].max()) if "ginhib_total" in result.columns else None,
    }
    print(f"  Stats: {stats}")
    return stats


async def main():
    """Main entry point."""
    print("=" * 60)
    print("MolE Antimicrobial Prediction for Reductive Amination Products")
    print("=" * 60)

    all_stats = {}
    for filename in FILES:
        input_path = INPUT_DIR / filename
        output_name = filename.replace(".csv", "_with_mole_predictions.csv")
        output_path = OUTPUT_DIR / output_name

        if not input_path.exists():
            print(f"Error: Input file not found: {input_path}")
            continue

        stats = await predict_file(input_path, output_path)
        all_stats[filename] = stats

    # Print summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for filename, stats in all_stats.items():
        print(f"\n{filename}:")
        for key, value in stats.items():
            print(f"  {key}: {value}")


if __name__ == "__main__":
    asyncio.run(main())
