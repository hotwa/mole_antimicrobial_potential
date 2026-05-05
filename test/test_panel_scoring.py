"""Tests for the pathogen-selective panel scoring module."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from src.panel_scoring import (
    PanelConfig,
    compute_panel_scores,
    compute_panel_scores_from_dataframe,
    load_panel_config,
    sigmoid,
    validate_panel_config,
)


VALID_PANEL_RAW = {
    "mode": "pathogen_selective",
    "label": "test_panel",
    "app_threshold": 0.05,
    "tau": 0.02,
    "lambda": 0.5,
    "primary_pathogen_panel": ["StrainA (NT0001)", "StrainB (NT0002)"],
    "secondary_pathogen_panel": ["StrainC (NT0003)"],
    "commensal_sparing_panel": ["StrainD (NT0004)", "StrainE (NT0005)"],
}

ALL_STRAINS = [
    "StrainA (NT0001)",
    "StrainB (NT0002)",
    "StrainC (NT0003)",
    "StrainD (NT0004)",
    "StrainE (NT0005)",
]


class TestSigmoid(unittest.TestCase):
    """Test sigmoid numerical correctness."""

    def test_sigmoid_zero(self):
        result = sigmoid(np.array([0.0]))
        self.assertAlmostEqual(float(result[0]), 0.5, places=10)

    def test_sigmoid_large_positive(self):
        result = sigmoid(np.array([100.0]))
        self.assertAlmostEqual(float(result[0]), 1.0, places=10)

    def test_sigmoid_large_negative(self):
        result = sigmoid(np.array([-100.0]))
        self.assertAlmostEqual(float(result[0]), 0.0, places=10)

    def test_sigmoid_symmetry(self):
        """sigmoid(x) + sigmoid(-x) == 1 for all x."""
        for x in [0.1, 1.0, 5.0, 10.0]:
            pos = float(sigmoid(np.array([x]))[0])
            neg = float(sigmoid(np.array([-x]))[0])
            self.assertAlmostEqual(pos + neg, 1.0, places=10)

    def test_sigmoid_known_values(self):
        self.assertAlmostEqual(float(sigmoid(np.array([1.0]))[0]), 0.7310586, places=5)
        self.assertAlmostEqual(float(sigmoid(np.array([-1.0]))[0]), 0.2689414, places=5)

    def test_sigmoid_array(self):
        arr = np.array([0.0, 1.0, -1.0])
        result = sigmoid(arr)
        self.assertEqual(result.shape, (3,))
        self.assertAlmostEqual(float(result[0]), 0.5, places=10)

    def test_sigmoid_numerical_stability(self):
        """No overflow or underflow for extreme inputs."""
        result = sigmoid(np.array([1000.0, -1000.0]))
        self.assertAlmostEqual(float(result[0]), 1.0, places=10)
        self.assertAlmostEqual(float(result[1]), 0.0, places=10)


class TestValidatePanelConfig(unittest.TestCase):
    """Test panel config validation."""

    def test_valid_config(self):
        panel = validate_panel_config(VALID_PANEL_RAW)
        self.assertEqual(panel.label, "test_panel")
        self.assertEqual(panel.mode, "pathogen_selective")
        self.assertAlmostEqual(panel.app_threshold, 0.05)
        self.assertAlmostEqual(panel.tau, 0.02)
        self.assertAlmostEqual(panel.lambda_, 0.5)
        self.assertEqual(len(panel.primary_pathogen_panel), 2)
        self.assertEqual(len(panel.secondary_pathogen_panel), 1)
        self.assertEqual(len(panel.commensal_sparing_panel), 2)

    def test_pathogen_panel_property(self):
        panel = validate_panel_config(VALID_PANEL_RAW)
        self.assertEqual(
            panel.pathogen_panel,
            ["StrainA (NT0001)", "StrainB (NT0002)", "StrainC (NT0003)"],
        )

    def test_all_strains(self):
        panel = validate_panel_config(VALID_PANEL_RAW)
        self.assertEqual(panel.all_strains(), ALL_STRAINS)

    def test_missing_key_raises(self):
        raw = dict(VALID_PANEL_RAW)
        del raw["tau"]
        with self.assertRaises(ValueError) as ctx:
            validate_panel_config(raw)
        self.assertIn("tau", str(ctx.exception))

    def test_empty_primary_panel_raises(self):
        raw = dict(VALID_PANEL_RAW)
        raw["primary_pathogen_panel"] = []
        with self.assertRaises(ValueError) as ctx:
            validate_panel_config(raw)
        self.assertIn("primary_pathogen_panel", str(ctx.exception))

    def test_duplicate_strains_across_panels_raises(self):
        raw = dict(VALID_PANEL_RAW)
        raw["commensal_sparing_panel"] = ["StrainA (NT0001)"]
        with self.assertRaises(ValueError) as ctx:
            validate_panel_config(raw)
        self.assertIn("duplicate", str(ctx.exception))

    def test_duplicate_within_panel_raises(self):
        raw = dict(VALID_PANEL_RAW)
        raw["primary_pathogen_panel"] = ["StrainA (NT0001)", "StrainA (NT0001)"]
        with self.assertRaises(ValueError):
            validate_panel_config(raw)

    def test_empty_strain_name_raises(self):
        raw = dict(VALID_PANEL_RAW)
        raw["primary_pathogen_panel"] = ["StrainA (NT0001)", ""]
        with self.assertRaises(ValueError):
            validate_panel_config(raw)

    def test_non_list_panel_raises(self):
        raw = dict(VALID_PANEL_RAW)
        raw["commensal_sparing_panel"] = "not_a_list"
        with self.assertRaises(ValueError):
            validate_panel_config(raw)

    def test_negative_tau_raises(self):
        raw = dict(VALID_PANEL_RAW)
        raw["tau"] = -0.01
        with self.assertRaises(ValueError):
            validate_panel_config(raw)

    def test_empty_mode_raises(self):
        raw = dict(VALID_PANEL_RAW)
        raw["mode"] = ""
        with self.assertRaises(ValueError):
            validate_panel_config(raw)


class TestLoadPanelConfig(unittest.TestCase):
    """Test loading panel config from a JSON file."""

    def test_load_valid_file(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            json.dump(VALID_PANEL_RAW, f)
            f.flush()
            panel = load_panel_config(f.name)
        self.assertEqual(panel.label, "test_panel")

    def test_load_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            load_panel_config("/nonexistent/path.json")

    def test_load_invalid_json_raises(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            f.write("{invalid json")
            f.flush()
            with self.assertRaises(json.JSONDecodeError):
                load_panel_config(f.name)


class TestComputePanelScores(unittest.TestCase):
    """Test panel score computation for a single molecule."""

    def setUp(self):
        self.panel = validate_panel_config(VALID_PANEL_RAW)
        self.strains = ALL_STRAINS

    def test_basic_computation(self):
        probs = np.array([0.9, 0.8, 0.7, 0.1, 0.05])
        scores = compute_panel_scores(probs, self.strains, self.panel)
        self.assertIn("pathogen_soft", scores)
        self.assertIn("pathogen_hard", scores)
        self.assertIn("commensal_soft", scores)
        self.assertIn("selectivity_score", scores)

    def test_high_pathogen_low_commensal(self):
        """High pathogen probs, low commensal probs -> high selectivity."""
        probs = np.array([0.99, 0.99, 0.99, 0.01, 0.01])
        scores = compute_panel_scores(probs, self.strains, self.panel)
        self.assertGreater(scores["pathogen_soft"], 0.9)
        # sigmoid((0.01 - 0.05) / 0.02) = sigmoid(-2) ≈ 0.119
        self.assertLess(scores["commensal_soft"], 0.15)
        self.assertGreater(scores["selectivity_score"], 0.5)

    def test_low_pathogen_high_commensal(self):
        """Low pathogen probs, high commensal probs -> negative selectivity."""
        probs = np.array([0.01, 0.01, 0.01, 0.99, 0.99])
        scores = compute_panel_scores(probs, self.strains, self.panel)
        # sigmoid((0.01 - 0.05) / 0.02) = sigmoid(-2) ≈ 0.119
        self.assertLess(scores["pathogen_soft"], 0.15)
        self.assertGreater(scores["commensal_soft"], 0.9)
        self.assertLess(scores["selectivity_score"], 0.0)

    def test_all_zeros(self):
        probs = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        scores = compute_panel_scores(probs, self.strains, self.panel)
        # sigmoid((0 - 0.05) / 0.02) = sigmoid(-2.5) ≈ 0.076
        self.assertAlmostEqual(scores["pathogen_soft"], float(sigmoid(np.array([-2.5]))[0]), places=5)
        self.assertEqual(scores["pathogen_hard"], 0)

    def test_all_ones(self):
        probs = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        scores = compute_panel_scores(probs, self.strains, self.panel)
        self.assertGreater(scores["pathogen_soft"], 0.99)
        self.assertEqual(scores["pathogen_hard"], 3)
        self.assertGreater(scores["commensal_soft"], 0.99)

    def test_selectivity_formula(self):
        """Verify selectivity_score = pathogen_soft - lambda * commensal_soft."""
        probs = np.array([0.8, 0.7, 0.6, 0.3, 0.2])
        scores = compute_panel_scores(probs, self.strains, self.panel)
        expected = scores["pathogen_soft"] - 0.5 * scores["commensal_soft"]
        self.assertAlmostEqual(scores["selectivity_score"], expected, places=10)

    def test_hard_count_at_threshold(self):
        """Probs exactly at threshold should be counted as hard hits."""
        t = self.panel.app_threshold
        probs = np.array([t, t + 0.001, t - 0.001, 0.5, 0.5])
        scores = compute_panel_scores(probs, self.strains, self.panel)
        # t and t+0.001 should count, t-0.001 should not
        self.assertEqual(scores["pathogen_hard"], 2)

    def test_override_threshold(self):
        probs = np.array([0.5, 0.5, 0.5, 0.5, 0.5])
        scores_default = compute_panel_scores(probs, self.strains, self.panel)
        scores_high = compute_panel_scores(
            probs, self.strains, self.panel, threshold=0.6
        )
        # With higher threshold, fewer hard hits
        self.assertGreaterEqual(scores_default["pathogen_hard"], scores_high["pathogen_hard"])

    def test_override_lambda(self):
        probs = np.array([0.9, 0.8, 0.7, 0.3, 0.2])
        scores_l05 = compute_panel_scores(probs, self.strains, self.panel)
        scores_l0 = compute_panel_scores(probs, self.strains, self.panel, lambda_=0.0)
        # With lambda=0, commensal penalty is zero
        self.assertAlmostEqual(scores_l0["selectivity_score"], scores_l0["pathogen_soft"])
        # With lambda=0.5, selectivity should be lower
        self.assertGreater(scores_l0["selectivity_score"], scores_l05["selectivity_score"])

    def test_missing_pathogen_strain_raises(self):
        probs = np.array([0.5, 0.5])  # Only 2 values, missing strains
        with self.assertRaises(ValueError) as ctx:
            compute_panel_scores(probs, ["StrainA (NT0001)", "StrainB (NT0002)"], self.panel)
        self.assertIn("not found", str(ctx.exception))

    def test_missing_commensal_strain_raises(self):
        """Missing commensal strains should raise an error."""
        strains = ["StrainA (NT0001)", "StrainB (NT0002)", "StrainC (NT0003)", "StrainD (NT0004)"]
        probs = np.array([0.5, 0.5, 0.5, 0.5])
        with self.assertRaises(ValueError) as ctx:
            compute_panel_scores(probs, strains, self.panel)
        self.assertIn("commensal", str(ctx.exception))

    def test_2d_array_raises(self):
        with self.assertRaises(ValueError):
            compute_panel_scores(np.array([[0.5, 0.5]]), self.strains, self.panel)


class TestComputePanelScoresFromDataframe(unittest.TestCase):
    """Test panel scoring from a per-strain prediction DataFrame."""

    def setUp(self):
        self.panel = validate_panel_config(VALID_PANEL_RAW)

    def _make_df(self, probs_by_mol: dict[str, list[float]]) -> pd.DataFrame:
        """Build a pred_id DataFrame from {chem_id: [prob per strain]}."""
        rows = []
        for chem_id, probs in probs_by_mol.items():
            for strain, p in zip(ALL_STRAINS, probs):
                rows.append({"pred_id": f"{chem_id}:{strain}", "1": p})
        return pd.DataFrame(rows)

    def test_basic(self):
        df = self._make_df({"mol1": [0.9, 0.8, 0.7, 0.1, 0.05]})
        scores = compute_panel_scores_from_dataframe(df, self.panel)
        self.assertEqual(len(scores), 1)
        self.assertIn("pathogen_soft", scores.columns)
        self.assertIn("selectivity_score", scores.columns)

    def test_multiple_molecules(self):
        df = self._make_df({
            "mol1": [0.9, 0.8, 0.7, 0.1, 0.05],
            "mol2": [0.1, 0.1, 0.1, 0.9, 0.9],
        })
        scores = compute_panel_scores_from_dataframe(df, self.panel)
        self.assertEqual(len(scores), 2)
        self.assertIn("mol1", scores.index)
        self.assertIn("mol2", scores.index)

    def test_selectivity_ordering(self):
        """A pathogen-killing molecule should score higher than a commensal-killing one."""
        df = self._make_df({
            "pathogen_killer": [0.99, 0.99, 0.99, 0.01, 0.01],
            "commensal_killer": [0.01, 0.01, 0.01, 0.99, 0.99],
        })
        scores = compute_panel_scores_from_dataframe(df, self.panel)
        self.assertGreater(
            scores.loc["pathogen_killer", "selectivity_score"],
            scores.loc["commensal_killer", "selectivity_score"],
        )

    def test_missing_pred_id_column_raises(self):
        df = pd.DataFrame({"chem_id": ["mol1"], "1": [0.5]})
        with self.assertRaises(ValueError) as ctx:
            compute_panel_scores_from_dataframe(df, self.panel)
        self.assertIn("pred_id", str(ctx.exception))

    def test_missing_probability_column_raises(self):
        df = pd.DataFrame({"pred_id": ["mol1:StrainA (NT0001)"]})
        with self.assertRaises(ValueError) as ctx:
            compute_panel_scores_from_dataframe(df, self.panel, probability_col="prob")
        self.assertIn("prob", str(ctx.exception))

    def test_override_parameters(self):
        df = self._make_df({"mol1": [0.5, 0.5, 0.5, 0.5, 0.5]})
        scores_default = compute_panel_scores_from_dataframe(df, self.panel)
        scores_custom = compute_panel_scores_from_dataframe(
            df, self.panel, lambda_=0.0, threshold=0.6
        )
        # Different parameters should give different results
        self.assertNotAlmostEqual(
            float(scores_default["selectivity_score"].iloc[0]),
            float(scores_custom["selectivity_score"].iloc[0]),
        )


class TestPanelConfigWithRealStrains(unittest.TestCase):
    """Test with the actual panel config file from the project."""

    PANEL_PATH = Path(__file__).resolve().parents[1] / "workflows" / "reinvent4" / "inputs" / "objectives" / "pathogen_selective.panel.v1.json"

    def test_load_real_config(self):
        if not self.PANEL_PATH.exists():
            self.skipTest("Panel config file not found")
        panel = load_panel_config(self.PANEL_PATH)
        self.assertEqual(panel.label, "pathogen_selective_v1")
        self.assertEqual(len(panel.primary_pathogen_panel), 4)
        self.assertEqual(len(panel.secondary_pathogen_panel), 3)
        self.assertEqual(len(panel.commensal_sparing_panel), 8)
        self.assertEqual(len(panel.pathogen_panel), 7)
        self.assertEqual(len(panel.all_strains()), 15)

    def test_real_config_has_no_duplicate_strains(self):
        if not self.PANEL_PATH.exists():
            self.skipTest("Panel config file not found")
        panel = load_panel_config(self.PANEL_PATH)
        all_strains = panel.all_strains()
        self.assertEqual(len(all_strains), len(set(all_strains)))


if __name__ == "__main__":
    unittest.main()
