#!/usr/bin/env python3
"""Regression tests for MolE representation model reuse."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import gmean

from src.models import MoleculeInfo
from src.predictor import AntimicrobialPredictor


class PredictorMolEModelCacheTestCase(TestCase):
    def _build_predictor(self) -> AntimicrobialPredictor:
        fake_probe = SimpleNamespace(
            preference="pickle",
            pickle_path=Path("/tmp/fake_model.pkl"),
            timber_model_dir=Path("/tmp/fake_timber"),
        )

        with patch("src.predictor.inspect_classifier_backends", return_value=fake_probe):
            predictor = AntimicrobialPredictor()
        return predictor

    def _legacy_antimicrobial_potential(
        self,
        predictor: AntimicrobialPredictor,
        score_df: pd.DataFrame,
    ) -> pd.DataFrame:
        split = score_df["pred_id"].astype(str).str.rsplit(":", n=1, expand=True)
        score_df = score_df.copy()
        score_df["chem_id"] = split[0]
        score_df["strain_name"] = split[1]
        score_df["nt_number"] = score_df["strain_name"].str.extract(r".*?\((NT\d+)\)", expand=False)
        score_df["gram_stain"] = score_df["nt_number"].map(predictor._gram_dict)

        apscore_total = (
            score_df.groupby("chem_id")["1"].apply(gmean).to_frame().rename(columns={"1": "apscore_total"})
        )
        apscore_total["apscore_total"] = np.log(apscore_total["apscore_total"])

        apscore_gram = (
            score_df.groupby(["chem_id", "gram_stain"])["1"]
            .apply(gmean)
            .unstack()
            .rename(columns={"negative": "apscore_gnegative", "positive": "apscore_gpositive"})
        )
        apscore_gram["apscore_gnegative"] = np.log(apscore_gram["apscore_gnegative"])
        apscore_gram["apscore_gpositive"] = np.log(apscore_gram["apscore_gpositive"])

        inhibited_total = (
            score_df.groupby("chem_id")["growth_inhibition"]
            .sum()
            .to_frame()
            .rename(columns={"growth_inhibition": "ginhib_total"})
        )
        inhibited_gram = (
            score_df.groupby(["chem_id", "gram_stain"])["growth_inhibition"]
            .sum()
            .unstack()
            .rename(columns={"negative": "ginhib_gnegative", "positive": "ginhib_gpositive"})
        )

        return apscore_total.join(apscore_gram).join(inhibited_total).join(inhibited_gram)

    def test_mole_representation_model_is_loaded_once_and_reused(self) -> None:
        predictor = self._build_predictor()

        fake_model = object()
        fake_representation = object()

        with patch("src.predictor.load_pretrained_model", return_value=fake_model) as load_mock, patch(
            "src.predictor.process_representation", return_value=fake_representation
        ) as process_mock:
            predictor._get_mole_representation(
                [MoleculeInfo(smiles="CCO", chem_id="mol1")]
            )
            predictor._get_mole_representation(
                [MoleculeInfo(smiles="CCN", chem_id="mol2")]
            )

        self.assertEqual(load_mock.call_count, 1)
        self.assertEqual(process_mock.call_count, 2)
        for call in process_mock.call_args_list:
            self.assertIs(call.kwargs["model"], fake_model)

    def test_representation_settings_are_forwarded_to_process_representation(self) -> None:
        predictor = self._build_predictor()

        fake_model = object()
        fake_representation = object()

        with patch("src.predictor.load_pretrained_model", return_value=fake_model), patch(
            "src.predictor.process_representation", return_value=fake_representation
        ) as process_mock:
            predictor._get_mole_representation(
                [MoleculeInfo(smiles="CCO", chem_id="mol1")],
                num_graph_workers=4,
                graph_batch_size=256,
                prefetch_batches=6,
                deterministic_representation=True,
            )

        self.assertEqual(process_mock.call_count, 1)
        kwargs = process_mock.call_args.kwargs
        self.assertEqual(kwargs["num_graph_workers"], 4)
        self.assertEqual(kwargs["graph_batch_size"], 256)
        self.assertEqual(kwargs["prefetch_batches"], 6)
        self.assertTrue(kwargs["deterministic_representation"])

    def test_predictor_marks_dataframe_input_as_prevalidated_for_representation(self) -> None:
        predictor = self._build_predictor()

        fake_model = object()
        fake_representation = object()

        with patch("src.predictor.load_pretrained_model", return_value=fake_model), patch(
            "src.predictor.process_representation", return_value=fake_representation
        ) as process_mock:
            predictor._get_mole_representation(
                [MoleculeInfo(smiles="CCO", chem_id="mol1")]
            )

        self.assertEqual(process_mock.call_count, 1)
        self.assertFalse(
            process_mock.call_args.kwargs["validate_smiles"],
            "Predictor path should skip redundant read_smiles() RDKit validation for already-validated MoleculeInfo input",
        )

    def test_add_strains_matches_legacy_cross_join_layout(self) -> None:
        predictor = self._build_predictor()
        chemfeats = pd.DataFrame(
            [[1.0, 2.0], [3.0, 4.0]],
            index=pd.Index(["mol1", "mol2"], name="chem_id"),
            columns=["feat_a", "feat_b"],
        )
        predictor.strain_ohe = pd.DataFrame(
            [[1.0, 0.0], [0.0, 1.0]],
            index=["strainA", "strainB"],
            columns=["('strainA',)", "('strainB',)"],
        )

        chemfe = chemfeats.reset_index().rename(columns={"index": "chem_id"})
        chemfe["chem_id"] = chemfe["chem_id"].astype(str)
        sohe = predictor.strain_ohe.reset_index().rename(columns={"index": "strain_name"})
        legacy = chemfe.merge(sohe, how="cross")
        legacy["pred_id"] = legacy["chem_id"].str.cat(legacy["strain_name"], sep=":")
        legacy = legacy.set_index("pred_id").drop(columns=["chem_id", "strain_name"])
        legacy.columns = [str(column) for column in legacy.columns]

        current = predictor._add_strains(chemfeats)

        self.assertEqual(list(current.index), list(legacy.index))
        self.assertEqual(list(current.columns), list(legacy.columns))
        np.testing.assert_allclose(current.to_numpy(), legacy.to_numpy(), rtol=0.0, atol=0.0)

    def test_strain_feature_array_matches_add_strains_values(self) -> None:
        predictor = self._build_predictor()
        chemfeats = pd.DataFrame(
            [[1.0, 2.0], [3.0, 4.0]],
            index=pd.Index(["mol1", "mol2"], name="chem_id"),
            columns=["feat_a", "feat_b"],
        )
        predictor.strain_ohe = pd.DataFrame(
            [[1.0, 0.0], [0.0, 1.0]],
            index=["strainA", "strainB"],
            columns=["('strainA',)", "('strainB',)"],
        )

        current = predictor._strain_feature_array(chemfeats)
        expected = predictor._add_strains(chemfeats).to_numpy(dtype=np.float32, copy=False)

        self.assertIsInstance(current, np.ndarray)
        self.assertEqual(current.dtype, np.float32)
        np.testing.assert_allclose(current, expected, rtol=0.0, atol=0.0)

    def test_predict_passes_ndarray_into_classifier_backend(self) -> None:
        predictor = self._build_predictor()
        predictor._model_loaded = True
        predictor.strain_ohe = pd.DataFrame(
            [[1.0, 0.0], [0.0, 1.0]],
            index=["strainA", "strainB"],
            columns=["('strainA',)", "('strainB',)"],
        )
        predictor._gram_dict = {"NT1": "negative"}

        representation = pd.DataFrame(
            [[1.0, 2.0], [3.0, 4.0]],
            index=pd.Index(["mol1", "mol2"], name="chem_id"),
            columns=["feat_a", "feat_b"],
        )

        class _RecordingClassifier:
            def __init__(self):
                self.inputs = []

            def predict_proba(self, X):
                self.inputs.append(X)
                rows = len(X)
                return np.tile(np.array([[0.25, 0.75]], dtype=np.float32), (rows, 1))

        predictor.model = _RecordingClassifier()

        with patch.object(predictor, "_get_mole_representation", return_value=representation):
            import asyncio
            from src.models import MoleculeInput

            result = asyncio.run(
                predictor.predict(
                    MoleculeInput(smiles=["CCO", "CCN"], aggregate_scores=False).normalize()
                )
            )

        self.assertIsInstance(predictor.model.inputs[0], np.ndarray)
        self.assertEqual(len(result), 4)

    def test_predict_records_aggregate_profiling_fields(self) -> None:
        predictor = self._build_predictor()
        predictor._model_loaded = True
        predictor.strain_ohe = pd.DataFrame(
            [[1.0, 0.0], [0.0, 1.0]],
            index=["Strain A (NT1)", "Strain B (NT2)"],
            columns=["('strainA',)", "('strainB',)"],
        )
        predictor._gram_dict = {"NT1": "negative", "NT2": "positive"}

        representation = pd.DataFrame(
            [[1.0, 2.0], [3.0, 4.0]],
            index=pd.Index(["mol1", "mol2"], name="chem_id"),
            columns=["feat_a", "feat_b"],
        )
        representation.attrs["profiling"] = {"graph_items": 2, "graph_total_seconds": 0.25}

        class _RecordingClassifier:
            def predict_proba(self, X):
                rows = len(X)
                return np.tile(np.array([[0.25, 0.75]], dtype=np.float32), (rows, 1))

        predictor.model = _RecordingClassifier()

        batches = [
            {
                "batch_index": 0,
                "chem_ids": ["mol1"],
                "embedding_batch": representation.iloc[[0]],
                "profiling": {
                    "graph_items": 1,
                    "rdkit_parse_seconds": 0.10,
                    "add_hs_seconds": 0.20,
                    "atom_feature_seconds": 0.30,
                    "bond_feature_seconds": 0.40,
                    "graph_total_seconds": 1.00,
                    "dataloader_iter_seconds": 0.01,
                    "model_forward_seconds": 0.02,
                },
            },
            {
                "batch_index": 1,
                "chem_ids": ["mol2"],
                "embedding_batch": representation.iloc[[1]],
                "profiling": {
                    "graph_items": 1,
                    "rdkit_parse_seconds": 0.15,
                    "add_hs_seconds": 0.25,
                    "atom_feature_seconds": 0.35,
                    "bond_feature_seconds": 0.45,
                    "graph_total_seconds": 1.20,
                    "dataloader_iter_seconds": 0.03,
                    "model_forward_seconds": 0.04,
                },
            },
        ]

        with patch.object(predictor, "_iter_mole_representation_batches", return_value=iter(batches)):
            import asyncio
            from src.models import MoleculeInput

            result = asyncio.run(
                predictor.predict(
                    MoleculeInput(smiles=["CCO", "CCN"], aggregate_scores=True).normalize(),
                    enable_profiling=True,
                )
            )

        self.assertEqual(len(result), 2)
        profile = predictor.last_profile_snapshot()
        self.assertIsInstance(profile, dict)
        self.assertIn("prediction_frame_seconds", profile)
        self.assertIn("growth_inhibition_seconds", profile)
        self.assertIn("aggregate_scores_seconds", profile)
        self.assertIn("first_result_latency_seconds", profile)
        self.assertIn("streaming_batches", profile)
        self.assertEqual(len(profile["streaming_batches"]), 2)
        self.assertIn("representation_batch_production_seconds", profile["streaming_batches"][0])
        self.assertIn("classifier_batch_inference_seconds", profile["streaming_batches"][0])
        self.assertIn("aggregate_accumulate_seconds", profile["streaming_batches"][0])
        self.assertEqual(profile["graph_build"]["graph_items"], 2)
        self.assertAlmostEqual(profile["graph_build"]["graph_total_seconds"], 2.2)
        self.assertGreaterEqual(profile["prediction_frame_seconds"], 0.0)
        self.assertGreaterEqual(profile["growth_inhibition_seconds"], 0.0)
        self.assertGreaterEqual(profile["aggregate_scores_seconds"], 0.0)

    def test_aggregate_scores_from_matrix_matches_legacy_groupby_result(self) -> None:
        predictor = self._build_predictor()
        predictor._gram_dict = {
            "NT1": "negative",
            "NT2": "positive",
            "NT3": "negative",
            "NT4": "positive",
        }

        chem_ids = np.array(["mol2", "mol10"], dtype=object)
        strain_names = np.array(
            [
                "Strain A (NT1)",
                "Strain B (NT2)",
                "Strain C (NT3)",
                "Strain D (NT4)",
            ],
            dtype=object,
        )
        probability_matrix = np.array(
            [
                [0.20, 0.50, 0.80, 0.90],
                [0.70, 0.40, 0.60, 0.30],
            ],
            dtype=np.float32,
        )
        growth_matrix = (probability_matrix >= 0.50).astype(np.int64)

        score_df = pd.DataFrame(
            {
                "pred_id": [
                    f"{chem_id}:{strain_name}"
                    for chem_id in chem_ids
                    for strain_name in strain_names
                ],
                "1": probability_matrix.reshape(-1),
                "growth_inhibition": growth_matrix.reshape(-1),
            }
        )

        legacy = self._legacy_antimicrobial_potential(predictor, score_df)
        current = predictor._aggregate_scores_from_matrix(
            probability_matrix=probability_matrix,
            growth_inhibition_matrix=growth_matrix,
            chem_ids=chem_ids,
            strain_names=strain_names,
        )

        pd.testing.assert_frame_equal(current, legacy)

    def test_streaming_aggregate_matches_full_representation_baseline(self) -> None:
        predictor = self._build_predictor()
        predictor._model_loaded = True
        predictor.strain_ohe = pd.DataFrame(
            np.eye(4, dtype=np.float32),
            index=[
                "Strain A (NT1)",
                "Strain B (NT2)",
                "Strain C (NT3)",
                "Strain D (NT4)",
            ],
            columns=[
                "('Strain A (NT1)',)",
                "('Strain B (NT2)',)",
                "('Strain C (NT3)',)",
                "('Strain D (NT4)',)",
            ],
        )
        predictor._gram_dict = {
            "NT1": "negative",
            "NT2": "positive",
            "NT3": "negative",
            "NT4": "positive",
        }

        full_representation = pd.DataFrame(
            [[1.0, 5.0], [2.0, 1.0], [3.0, 4.0]],
            index=pd.Index(["mol2", "mol10", "mol1"], name="chem_id"),
            columns=["feat_a", "feat_b"],
        )

        class _FeatureAwareClassifier:
            def predict_proba(self, X):
                strain_index = np.argmax(X[:, 2:], axis=1)
                probability = np.clip(
                    0.05 + (X[:, 0] * 0.11) + (X[:, 1] * 0.03) + (strain_index * 0.07),
                    1e-6,
                    1.0 - 1e-6,
                )
                return np.column_stack([1.0 - probability, probability]).astype(np.float64, copy=False)

        predictor.model = _FeatureAwareClassifier()
        threshold = 0.55
        min_nkill = 2

        x_input = predictor._add_strains(full_representation)
        probabilities = predictor.model.predict_proba(x_input.to_numpy(dtype=np.float32, copy=False))[:, 1]
        strain_names = predictor.strain_ohe.index.astype(str).to_numpy(copy=False)
        probability_matrix = probabilities.reshape(len(full_representation), len(strain_names))
        growth_matrix = (probability_matrix >= threshold).astype(np.int64, copy=False)
        expected_df = predictor._aggregate_scores_from_matrix(
            probability_matrix=probability_matrix,
            growth_inhibition_matrix=growth_matrix,
            chem_ids=full_representation.index.astype(str).to_numpy(copy=False),
            strain_names=strain_names,
        )
        expected_df["broad_spectrum"] = (
            expected_df["ginhib_total"].to_numpy(copy=False) >= min_nkill
        ).astype(np.int64, copy=False)
        expected_records = expected_df.reset_index().to_dict(orient="records")

        staged_batches = [
            {
                "batch_index": 0,
                "chem_ids": ["mol2", "mol10"],
                "embedding_batch": full_representation.iloc[:2],
                "profiling": None,
            },
            {
                "batch_index": 1,
                "chem_ids": ["mol1"],
                "embedding_batch": full_representation.iloc[2:],
                "profiling": None,
            },
        ]

        predictor._get_mole_representation = lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("aggregate streaming path should not build a full representation frame first")
        )

        with patch.object(predictor, "_iter_mole_representation_batches", return_value=iter(staged_batches)):
            import asyncio
            from src.models import MoleculeInput

            result = asyncio.run(
                predictor.predict(
                    MoleculeInput(
                        molecules=[
                            MoleculeInfo(smiles="CCO", chem_id="mol2"),
                            MoleculeInfo(smiles="CCN", chem_id="mol10"),
                            MoleculeInfo(smiles="CCC", chem_id="mol1"),
                        ],
                        aggregate_scores=True,
                        app_threshold=threshold,
                        min_nkill=min_nkill,
                    ).normalize()
                )
            )

        self.assertEqual(result, expected_records)
        self.assertEqual(
            list(result[0].keys()),
            [
                "chem_id",
                "apscore_total",
                "apscore_gnegative",
                "apscore_gpositive",
                "ginhib_total",
                "ginhib_gnegative",
                "ginhib_gpositive",
                "broad_spectrum",
            ],
        )
        self.assertEqual(
            [row["chem_id"] for row in result],
            expected_df.index.tolist(),
        )

    def test_streaming_aggregate_path_does_not_wrap_entire_job_in_wait_for(self) -> None:
        predictor = self._build_predictor()
        predictor._model_loaded = True

        expected = [{"chem_id": "mol1", "broad_spectrum": 0}]
        predictor._predict_aggregate_streaming_sync = lambda *args, **kwargs: (expected, None)

        import asyncio
        from src.models import MoleculeInput

        with patch("src.predictor.asyncio.wait_for", side_effect=AssertionError("aggregate path should not use wait_for")):
            result = asyncio.run(
                predictor.predict(
                    MoleculeInput(smiles=["CCO"], aggregate_scores=True).normalize(),
                )
            )

        self.assertEqual(result, expected)
