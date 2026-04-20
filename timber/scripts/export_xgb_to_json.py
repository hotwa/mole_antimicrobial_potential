#!/usr/bin/env python3
"""
Export pickled XGBoost model to Timber-compatible JSON format.

Usage:
    cd timber
    pixi run export-json

Or directly:
    python scripts/export_xgb_to_json.py
"""

from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

try:
    import numpy as np
    import pandas as pd
    import xgboost as xgb
except ImportError as e:
    print(f"Error: Missing required package: {e}")
    print("Run: pixi install")
    sys.exit(1)


def load_pickle_model(pkl_path: Path) -> xgb.XGBClassifier:
    """Load pickled XGBoost sklearn wrapper."""
    print(f"Loading pickle model from: {pkl_path}")

    with open(pkl_path, "rb") as f:
        model = pickle.load(f)

    print(f"  Loaded object type: {type(model)}")

    # Verify it's an XGBClassifier
    if not isinstance(model, xgb.XGBClassifier):
        raise TypeError(
            f"Expected XGBClassifier, got {type(model)}. "
            "This model may not be compatible with Timber conversion."
        )

    return model


def export_to_json(model: xgb.XGBClassifier, output_path: Path) -> xgb.Booster:
    """Export XGBClassifier to JSON format for Timber."""
    print(f"\nExporting to JSON: {output_path}")

    # Get underlying booster
    booster = model.get_booster()
    print(f"  Booster type: {type(booster)}")
    print(f"  Number of features: {booster.num_features()}")
    print(f"  Number of trees: {booster.num_boosted_rounds()}")

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save as JSON (Timber compatible format)
    booster.save_model(str(output_path))

    print(f"  ✓ Successfully exported to: {output_path}")

    return booster


def save_model_info(model: xgb.XGBClassifier, booster: xgb.Booster, info_path: Path):
    """Save model metadata for reference."""
    info = {
        "xgboost_version": xgb.__version__,
        "model_type": str(type(model)),
        "booster_type": str(type(booster)),
        "num_features": int(booster.num_features()),
        "num_trees": int(booster.num_boosted_rounds()),
        "feature_names": booster.feature_names if booster.feature_names else None,
        "model_attributes": {
            "n_estimators": model.n_estimators if hasattr(model, "n_estimators") else None,
            "max_depth": model.max_depth if hasattr(model, "max_depth") else None,
            "learning_rate": model.learning_rate if hasattr(model, "learning_rate") else None,
        },
    }

    with open(info_path, "w") as f:
        json.dump(info, f, indent=2)

    print(f"  ✓ Model info saved to: {info_path}")


def main():
    # Paths (relative to timber/ directory)
    base_dir = Path(__file__).parent.parent
    repo_root = base_dir.parent

    pkl_path = repo_root / "data" / "03.model_evaluation" / "MolE-XGBoost-08.03.2024_14.20.pkl"
    output_dir = base_dir / "timber_output"
    json_path = output_dir / "mole_xgb.json"
    info_path = output_dir / "model_info.json"

    print("=" * 60)
    print("MolE XGBoost -> Timber JSON Export")
    print("=" * 60)

    # Verify source file exists
    if not pkl_path.exists():
        print(f"\nError: Model file not found: {pkl_path}")
        print("\nPlease ensure the model file exists:")
        print(f"  {pkl_path}")
        sys.exit(1)

    print(f"\nSource: {pkl_path}")
    print(f"Output: {output_dir}")

    # Load and export
    try:
        model = load_pickle_model(pkl_path)
        booster = export_to_json(model, json_path)
        save_model_info(model, booster, info_path)

        print("\n" + "=" * 60)
        print("Export completed successfully!")
        print("=" * 60)
        print(f"\nNext steps:")
        print(f"  1. Load into Timber:  timber load {json_path} --name mole-antimicrobial")
        print(f"  2. Validate:          timber validate mole-antimicrobial")
        print(f"  3. Test consistency:  pixi run test-consistency")
        print(f"  4. Serve:             timber serve mole-antimicrobial --port 8080")

    except Exception as e:
        print(f"\nError during export: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
