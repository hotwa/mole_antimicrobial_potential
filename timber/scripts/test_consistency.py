#!/usr/bin/env python3
"""
Test numerical consistency between original pickle model and Timber-compiled model.

This script:
1. Loads the original pickled XGBClassifier
2. Loads the Timber-exported JSON model via xgboost.Booster
3. Compares predictions on synthetic/random data
4. Reports max absolute difference (should be ~0 for identical models)

Usage:
    cd timber
    pixi run test-consistency
"""

from __future__ import annotations

import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
import xgboost as xgb


def generate_synthetic_features(n_samples: int = 100, n_features: int = 8100) -> np.ndarray:
    """Generate synthetic feature matrix for testing."""
    np.random.seed(42)
    # Use realistic range based on typical molecular embeddings
    return np.random.randn(n_samples, n_features).astype(np.float32)


def load_original_model(pkl_path: Path) -> xgb.Booster:
    """Load original pickled model and extract booster."""
    print(f"Loading original model: {pkl_path}")
    with open(pkl_path, "rb") as f:
        model = pickle.load(f)
    print(f"  Type: {type(model)}")

    # Extract booster from XGBClassifier to avoid sklearn version issues
    if hasattr(model, 'get_booster'):
        booster = model.get_booster()
        print(f"  Extracted booster type: {type(booster)}")
        return booster
    else:
        raise TypeError(f"Model does not have get_booster(): {type(model)}")


def load_timber_model(json_path: Path) -> xgb.Booster:
    """Load Timber-exported JSON model."""
    print(f"Loading Timber JSON model: {json_path}")
    booster = xgb.Booster()
    booster.load_model(str(json_path))
    print(f"  Type: {type(booster)}")
    print(f"  Features: {booster.num_features()}")
    return booster


def compare_predictions(
    original_model: xgb.Booster,
    timber_booster: xgb.Booster,
    X: np.ndarray,
) -> dict:
    """Compare predictions between original and Timber models."""
    print(f"\nComparing predictions on {X.shape[0]} samples...")

    # Original model prediction (Booster API to avoid sklearn version issues)
    # Get feature names from original model and apply to test data
    feature_names = original_model.feature_names
    dmat_orig = xgb.DMatrix(X, feature_names=feature_names)
    pred_orig = original_model.predict(dmat_orig, validate_features=False)
    print(f"  Original predictions shape: {pred_orig.shape}")

    # Timber model prediction (native Booster API)
    dmatrix = xgb.DMatrix(X, feature_names=feature_names)
    pred_timber = timber_booster.predict(dmatrix, validate_features=False)
    print(f"  Timber predictions shape: {pred_timber.shape}")

    # Calculate differences
    abs_diff = np.abs(pred_orig - pred_timber)
    max_diff = float(np.max(abs_diff))
    mean_diff = float(np.mean(abs_diff))

    results = {
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "max_absolute_difference": max_diff,
        "mean_absolute_difference": mean_diff,
        "predictions_match": max_diff < 1e-6,
    }

    return results, pred_orig, pred_timber


def main():
    base_dir = Path(__file__).parent.parent
    repo_root = base_dir.parent

    pkl_path = repo_root / "data" / "03.model_evaluation" / "MolE-XGBoost-08.03.2024_14.20.pkl"
    json_path = base_dir / "timber_output" / "mole_xgb.json"

    print("=" * 60)
    print("Timber Conversion Consistency Test")
    print("=" * 60)

    # Check files exist
    if not pkl_path.exists():
        print(f"\nError: Original model not found: {pkl_path}")
        sys.exit(1)

    if not json_path.exists():
        print(f"\nError: Timber JSON model not found: {json_path}")
        print("Run 'pixi run export-json' first!")
        sys.exit(1)

    # Load models
    try:
        original = load_original_model(pkl_path)
        timber = load_timber_model(json_path)

        # Get expected feature count
        n_features = timber.num_features()
        print(f"\nExpected input features: {n_features}")

        # Generate test data
        X_test = generate_synthetic_features(n_samples=50, n_features=n_features)
        print(f"Generated test data: {X_test.shape}")

        # Compare
        results, pred_orig, pred_timber = compare_predictions(original, timber, X_test)

        # Report
        print("\n" + "=" * 60)
        print("Results")
        print("=" * 60)
        print(f"Max absolute difference:  {results['max_absolute_difference']:.2e}")
        print(f"Mean absolute difference: {results['mean_absolute_difference']:.2e}")

        if results["predictions_match"]:
            print("\n✓ PASS: Predictions are numerically identical!")
            print("  The Timber conversion is consistent with the original model.")
        else:
            print("\n⚠ WARNING: Predictions differ!")
            print(f"  Max diff: {results['max_absolute_difference']}")
            print("  This may indicate an issue with the conversion.")

        # Show sample predictions
        print("\nSample predictions (first 5):")
        print(f"  Original: {pred_orig[:5]}")
        print(f"  Timber:   {pred_timber[:5]}")

        return 0 if results["predictions_match"] else 1

    except Exception as e:
        print(f"\nError during testing: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
