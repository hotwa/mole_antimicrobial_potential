#!/usr/bin/env python3
"""Unit tests for the dedicated REINVENT4 /score API."""

import math
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from fastapi.testclient import TestClient

import api_server


class FakePredictor:
    """Minimal predictor stub for API tests."""

    def get_status(self):
        return SimpleNamespace(
            service="mole-antimicrobial-prediction",
            version="1.0.0",
            loaded=True,
            device="cpu",
        )

    async def ensure_loaded(self) -> None:
        return None

    async def predict(self, input_data):
        chem_ids = [molecule.chem_id for molecule in input_data.molecules]
        probabilities = {
            "mol1": {
                "Strain A (NT1)": 0.8,
                "Strain B (NT2)": 0.2,
                "Strain C (NT3)": 0.04,
            },
            "mol2": {
                "Strain A (NT1)": 0.6,
                "Strain B (NT2)": 0.9,
                "Strain C (NT3)": 0.01,
            },
        }
        items = []
        for chem_id in chem_ids:
            for strain_name, probability in probabilities[chem_id].items():
                items.append(
                    {
                        "pred_id": f"{chem_id}:{strain_name}",
                        "antimicrobial_predictive_probability": probability,
                        "growth_inhibition": int(probability >= input_data.app_threshold),
                    }
                )
        return items


class ScoreApiTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.predictor = FakePredictor()
        self.patcher = patch("api_server.get_predictor", return_value=self.predictor)
        self.patcher.start()
        self.client = TestClient(api_server.app)

    def tearDown(self) -> None:
        self.patcher.stop()

    def test_single_target_score_uses_target_probability(self) -> None:
        response = self.client.post(
            "/score",
            json={
                "smiles": ["CCO", "CCN"],
                "chem_id": ["mol1", "mol2"],
                "objective": {
                    "mode": "single_strain",
                    "strain": "Strain A (NT1)",
                },
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["mode"], "reinvent_score")
        self.assertEqual(payload["objective"]["mode"], "single_strain")

        scores = {item["chem_id"]: item["score"] for item in payload["items"]}
        self.assertAlmostEqual(scores["mol1"], 0.8, places=6)
        self.assertAlmostEqual(scores["mol2"], 0.6, places=6)

    def test_multi_target_score_uses_normalized_weighted_mean(self) -> None:
        response = self.client.post(
            "/score",
            json={
                "smiles": ["CCO", "CCN"],
                "chem_id": ["mol1", "mol2"],
                "objective": {
                    "mode": "single_strain",
                    "strains": ["Strain A (NT1)", "Strain B (NT2)"],
                    "weights": [3, 1],
                },
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        items = {item["chem_id"]: item for item in payload["items"]}

        self.assertAlmostEqual(items["mol1"]["score"], 0.65, places=6)
        self.assertAlmostEqual(items["mol2"]["score"], 0.675, places=6)
        self.assertEqual(items["mol1"]["weights"], [0.75, 0.25])

    def test_broad_spectrum_soft_score_is_continuous_and_normalized(self) -> None:
        app_threshold = 0.1
        tau = 0.05
        response = self.client.post(
            "/score",
            json={
                "smiles": ["CCO", "CCN"],
                "chem_id": ["mol1", "mol2"],
                "app_threshold": app_threshold,
                "objective": {
                    "mode": "broad_spectrum_soft",
                    "tau": tau,
                },
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        items = {item["chem_id"]: item for item in payload["items"]}

        def sigmoid(value: float) -> float:
            return 1.0 / (1.0 + math.exp(-value))

        mol1_expected = (
            sigmoid((0.8 - app_threshold) / tau)
            + sigmoid((0.2 - app_threshold) / tau)
            + sigmoid((0.04 - app_threshold) / tau)
        ) / 3.0
        mol2_expected = (
            sigmoid((0.6 - app_threshold) / tau)
            + sigmoid((0.9 - app_threshold) / tau)
            + sigmoid((0.01 - app_threshold) / tau)
        ) / 3.0

        self.assertAlmostEqual(items["mol1"]["score"], mol1_expected, places=6)
        self.assertAlmostEqual(items["mol2"]["score"], mol2_expected, places=6)
        self.assertAlmostEqual(
            items["mol1"]["soft_inhibition_count"],
            items["mol1"]["score"] * 3,
            places=6,
        )

    def test_unknown_strain_returns_validation_error(self) -> None:
        response = self.client.post(
            "/score",
            json={
                "smiles": "CCO",
                "chem_id": "mol1",
                "objective": {
                    "mode": "single_strain",
                    "strain": "Missing strain",
                },
            },
        )

        self.assertEqual(response.status_code, 422, response.text)
        payload = response.json()
        self.assertIn("Unknown strains requested", payload["detail"])

    def test_site_reward_can_be_enabled_as_structural_auxiliary_term(self) -> None:
        response = self.client.post(
            "/score",
            json={
                "smiles": "NCCCO",
                "chem_id": "mol1",
                "objective": {
                    "mode": "single_strain",
                    "strain": "Strain A (NT1)",
                    "site_reward": {
                        "enabled": True,
                        "scaffold_smiles": "[*]CC[*]",
                        "range_min": 1,
                        "range_max": 3,
                        "alpha": 1.0,
                        "beta": 1.0,
                        "coverage_weight": 0.7,
                        "balance_weight": 0.3,
                        "lambda": 0.8,
                    },
                },
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        item = response.json()["items"][0]

        def sigmoid(value: float) -> float:
            return 1.0 / (1.0 + math.exp(-value))

        site_1 = sigmoid((1 - 1) / 1.0) * sigmoid((3 - 1) / 1.0)
        site_2 = sigmoid((2 - 1) / 1.0) * sigmoid((3 - 2) / 1.0)
        coverage = (site_1 + site_2) / 2.0
        balance = math.exp(-0.5 / (1.5 + 1e-6))
        site_reward = (0.7 * coverage) + (0.3 * balance)
        final_score = (0.8 * 0.8) + (0.2 * site_reward)

        self.assertAlmostEqual(item["mole_reward"], 0.8, places=6)
        self.assertAlmostEqual(item["site_reward"], site_reward, places=6)
        self.assertAlmostEqual(item["score"], final_score, places=6)
        self.assertEqual(item["site_heavy_atoms"], [1, 2])
        self.assertTrue(item["site_decomposition_success"])


if __name__ == "__main__":
    unittest.main()
