from __future__ import annotations

import asyncio
from pathlib import Path
from types import SimpleNamespace
from unittest import IsolatedAsyncioTestCase, TestCase
from unittest.mock import AsyncMock, Mock, patch

import numpy as np
import pandas as pd
from pydantic import ValidationError

from src.models import MoleculeInfo, MoleculeInput
from src.predictor import AntimicrobialPredictor


class MoleculeInputNormalizationTestCase(TestCase):
    def test_molecule_info_rejects_invalid_smiles_at_construction(self) -> None:
        with self.assertRaises(ValidationError):
            MoleculeInfo(smiles="not_a_smiles")

    def test_molecule_input_rejects_invalid_smiles_in_list(self) -> None:
        with self.assertRaises(ValidationError):
            MoleculeInput(smiles=["CCO", "not_a_smiles"])

    def test_molecule_input_normalize_assigns_default_chem_ids(self) -> None:
        normalized = MoleculeInput(
            smiles=["CCO", "CCN"],
            aggregate_scores=True,
        ).normalize()

        self.assertIsNone(normalized.smiles)
        self.assertIsNone(normalized.chem_id)
        self.assertEqual(
            [(molecule.smiles, molecule.chem_id) for molecule in normalized.molecules or []],
            [("CCO", "mol1"), ("CCN", "mol2")],
        )

    def test_molecule_input_normalize_rejects_duplicate_chem_ids(self) -> None:
        with self.assertRaisesRegex(ValueError, "chem_id values must be unique"):
            MoleculeInput(
                smiles=["CCO", "CCN"],
                chem_id=["dup", "dup"],
            ).normalize()


class PredictorNormalizationRegressionTestCase(IsolatedAsyncioTestCase):
    def _build_predictor(self) -> AntimicrobialPredictor:
        fake_probe = SimpleNamespace(
            preference="pickle",
            pickle_path=Path("/tmp/fake_model.pkl"),
            timber_model_dir=Path("/tmp/fake_timber"),
        )

        with patch("src.predictor.inspect_classifier_backends", return_value=fake_probe):
            predictor = AntimicrobialPredictor()

        predictor._model_loaded = True
        predictor.ensure_loaded = AsyncMock()
        predictor._get_mole_representation = Mock(
            return_value=pd.DataFrame({"feat": [1.0, 2.0]}, index=["x", "y"])
        )
        predictor._add_strains = Mock(
            side_effect=lambda _: pd.DataFrame(
                {"feature": [1.0, 2.0, 3.0, 4.0]},
                index=pd.Index(
                    [
                        "x:Strain A (NT1)",
                        "x:Strain B (NT2)",
                        "y:Strain A (NT1)",
                        "y:Strain B (NT2)",
                    ],
                    name="pred_id",
                ),
            )
        )
        predictor.model = Mock()
        predictor.model.predict_proba = Mock(
            return_value=np.array(
                [
                    [0.2, 0.8],
                    [0.7, 0.3],
                    [0.4, 0.6],
                    [0.1, 0.9],
                ]
            )
        )
        return predictor

    async def test_predict_accepts_prebuilt_molecule_info_list(self) -> None:
        predictor = self._build_predictor()
        input_data = MoleculeInput(
            molecules=[MoleculeInfo(smiles="CCO", chem_id="a")],
            aggregate_scores=False,
        )

        predictor._get_mole_representation = Mock(
            return_value=pd.DataFrame({"feat": [1.0]}, index=["a"])
        )
        predictor._add_strains = Mock(
            return_value=pd.DataFrame(
                {"feature": [1.0, 2.0]},
                index=pd.Index(
                    ["a:Strain A (NT1)", "a:Strain B (NT2)"],
                    name="pred_id",
                ),
            )
        )
        predictor.model.predict_proba = Mock(return_value=np.array([[0.2, 0.8], [0.7, 0.3]]))

        result = await predictor.predict(input_data)

        self.assertEqual(
            result,
            [
                {
                    "pred_id": "a:Strain A (NT1)",
                    "antimicrobial_predictive_probability": 0.8,
                    "growth_inhibition": 1,
                },
                {
                    "pred_id": "a:Strain B (NT2)",
                    "antimicrobial_predictive_probability": 0.3,
                    "growth_inhibition": 1,
                },
            ],
        )

    async def test_predict_defaults_single_smiles_chem_id_to_mol1(self) -> None:
        predictor = self._build_predictor()
        predictor._get_mole_representation = Mock(
            return_value=pd.DataFrame({"feat": [1.0]}, index=["mol1"])
        )
        predictor._add_strains = Mock(
            return_value=pd.DataFrame(
                {"feature": [1.0]},
                index=pd.Index(["mol1:Strain A (NT1)"], name="pred_id"),
            )
        )
        predictor.model.predict_proba = Mock(return_value=np.array([[0.2, 0.8]]))

        await predictor.predict(MoleculeInput(smiles="CCO", aggregate_scores=False))

        molecules = predictor._get_mole_representation.call_args.args[0]
        self.assertEqual([(molecule.smiles, molecule.chem_id) for molecule in molecules], [("CCO", "mol1")])

    async def test_predict_preserves_input_order_and_pred_id_order(self) -> None:
        predictor = self._build_predictor()

        result = await predictor.predict(
            MoleculeInput(
                smiles=["CCO", "CCN"],
                chem_id=["x", "y"],
                aggregate_scores=False,
            )
        )

        self.assertEqual(
            [item["pred_id"] for item in result],
            [
                "x:Strain A (NT1)",
                "x:Strain B (NT2)",
                "y:Strain A (NT1)",
                "y:Strain B (NT2)",
            ],
        )
