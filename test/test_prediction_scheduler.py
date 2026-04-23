"""Tests for src/prediction_scheduler.py – Task 1 (failing first, then pass after impl)."""

from __future__ import annotations

import asyncio
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pandas as pd
import torch


# ---------------------------------------------------------------------------
# Helper: run a coroutine in an isolated event loop
# ---------------------------------------------------------------------------

def _run(coro):
    return asyncio.get_event_loop().run_until_complete(coro)


# ---------------------------------------------------------------------------
# Test: batch size scales with GPU memory budget
# ---------------------------------------------------------------------------

class TestBatchSizeSelection(unittest.TestCase):
    """select_batch_size must return larger batches on bigger GPUs."""

    # Using completely isolated static function.
    def test_22gb_gpu_smaller_than_32gb_gpu(self) -> None:
        """Use bytes_per_item large enough that neither budget saturates max_batch_size."""
        from src.prediction_scheduler import PredictionScheduler

        gpu22 = 22 * 1024 ** 3   # 22 GB
        gpu32 = 32 * 1024 ** 3   # 32 GB
        # 200 MB per item → 22 GB * 0.8 / 200 MB ≈ 88, 32 GB * 0.8 / 200 MB ≈ 128
        bytes_per_item = 200 * 1024 * 1024

        bs22 = PredictionScheduler._compute_batch_size(
            free_bytes=gpu22,
            bytes_per_item=bytes_per_item,
            target_fraction=0.80,
            max_batch_size=2048,
        )
        bs32 = PredictionScheduler._compute_batch_size(
            free_bytes=gpu32,
            bytes_per_item=bytes_per_item,
            target_fraction=0.80,
            max_batch_size=2048,
        )

        self.assertGreater(bs32, bs22, "32 GB GPU should produce a larger batch size than 22 GB GPU")
        self.assertGreater(bs22, 0)
        self.assertLessEqual(bs32, 2048)

    def test_batch_size_never_exceeds_max(self) -> None:
        from src.prediction_scheduler import PredictionScheduler

        bs = PredictionScheduler._compute_batch_size(
            free_bytes=128 * 1024 ** 3,   # absurdly large
            bytes_per_item=1,
            target_fraction=0.80,
            max_batch_size=512,
        )
        self.assertLessEqual(bs, 512)

    def test_large_gpu_not_capped_by_2048(self) -> None:
        from src.prediction_scheduler import DEFAULT_MAX_BATCH_SIZE
        self.assertGreater(DEFAULT_MAX_BATCH_SIZE, 2048, "DEFAULT_MAX_BATCH_SIZE must be raised to allow large GPUs to burst")

    def test_batch_size_at_least_one(self) -> None:
        from src.prediction_scheduler import PredictionScheduler

        bs = PredictionScheduler._compute_batch_size(
            free_bytes=1024,              # almost nothing
            bytes_per_item=500_000,
            target_fraction=0.80,
            max_batch_size=2048,
        )
        self.assertGreaterEqual(bs, 1)


# ---------------------------------------------------------------------------
# Test: OOM backoff halves batch size and retries
# ---------------------------------------------------------------------------

class TestOOMBackoff(unittest.TestCase):
    """Scheduler must halve batch size when the predictor raises OOM and retry."""

    def _oom_then_ok(self, call_count_box: list):
        """Return a mock predictor whose predict() raises OOM on first call."""

        async def _predict(input_data, **kwargs):
            call_count_box[0] += 1
            if call_count_box[0] == 1:
                raise torch.cuda.OutOfMemoryError("CUDA out of memory")  # type: ignore[attr-defined]
            # Second call succeeds with dummy data
            return [{"chem_id": "mol1", "apscore_total": -1.5}]

        predictor = mock.AsyncMock()
        predictor.predict.side_effect = _predict
        predictor.ensure_loaded = mock.AsyncMock()
        return predictor

    def test_oom_splits_batch_correctly(self) -> None:
        from src.prediction_scheduler import PredictionScheduler
        from src.models import MoleculeInfo

        call_args_box = []

        async def _predict(input_data, **kwargs):
            call_args_box.append(len(input_data.molecules))
            if len(input_data.molecules) > 2:
                raise torch.cuda.OutOfMemoryError("CUDA out of memory")  # type: ignore[attr-defined]
            return [{"chem_id": m.chem_id, "apscore_total": -1.5} for m in input_data.molecules]

        mock_predictor = mock.AsyncMock()
        mock_predictor.predict.side_effect = _predict
        mock_predictor.ensure_loaded = mock.AsyncMock()

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=4,
            max_batch_size=4,
            min_batch_size=1,
            target_memory_fraction=0.80,
            deterministic_representation=True,
        )

        molecules = [MoleculeInfo(smiles="CCO", chem_id=f"mol{i}") for i in range(4)]

        async def _run_predict():
            return await sched.predict_molecules(
                molecules=molecules,
                aggregate_scores=True,
                app_threshold=0.04374,
                min_nkill=10,
            )

        result = _run(_run_predict())
        # First call size 4 -> OOM. Batch size halves to 2.
        # Next call size 2 -> OK.
        # Next call size 2 -> OK.
        self.assertEqual(call_args_box, [4, 2, 2])
        self.assertEqual(sched.current_batch_size, 2)
        self.assertEqual(len(result), 4)

    def test_oom_halves_batch_and_retries(self) -> None:
        from src.prediction_scheduler import PredictionScheduler
        from src.models import MoleculeInfo

        call_count_box = [0]
        mock_predictor = self._oom_then_ok(call_count_box)

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=8,
            max_batch_size=8,
            min_batch_size=1,
            target_memory_fraction=0.80,
        )

        molecules = [MoleculeInfo(smiles="CCO", chem_id="mol1")]
        from src.models import MoleculeInput

        async def _run_predict():
            return await sched.predict_molecules(
                molecules=molecules,
                aggregate_scores=True,
                app_threshold=0.04374,
                min_nkill=10,
            )

        result = _run(_run_predict())
        self.assertEqual(call_count_box[0], 2, "Predictor must be called exactly twice (OOM + retry)")
        self.assertEqual(sched.current_batch_size, 4, "Batch size must be halved after OOM")
        self.assertIsInstance(result, list)

    def test_oom_min_batch_raises_after_exhaust(self) -> None:
        """If OOM persists at minimum batch size, the error must propagate."""
        from src.prediction_scheduler import PredictionScheduler
        from src.models import MoleculeInfo

        async def _always_oom(input_data, **kwargs):
            raise torch.cuda.OutOfMemoryError("CUDA out of memory")  # type: ignore[attr-defined]

        mock_predictor = mock.AsyncMock()
        mock_predictor.predict.side_effect = _always_oom
        mock_predictor.ensure_loaded = mock.AsyncMock()

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=2,
            max_batch_size=2,
            min_batch_size=2,   # no room to reduce
            target_memory_fraction=0.80,
        )

        molecules = [MoleculeInfo(smiles="CCO", chem_id="mol1")]

        async def _run_predict():
            return await sched.predict_molecules(
                molecules=molecules,
                aggregate_scores=True,
                app_threshold=0.04374,
                min_nkill=10,
            )

        with self.assertRaises(torch.cuda.OutOfMemoryError):  # type: ignore[attr-defined]
            _run(_run_predict())


# ---------------------------------------------------------------------------
# Test: MolE model loaded once and reused across batches
# ---------------------------------------------------------------------------

class TestModelReuse(unittest.TestCase):
    """ensure_loaded must be called once per scheduler, not once per batch."""

    def test_mole_model_loaded_once_across_batches(self) -> None:
        from src.prediction_scheduler import PredictionScheduler
        from src.models import MoleculeInfo

        async def _fake_predict(input_data, **kwargs):
            return [{"chem_id": m.chem_id, "apscore_total": -1.0} for m in input_data.molecules or []]

        mock_predictor = mock.AsyncMock()
        mock_predictor.predict.side_effect = _fake_predict
        mock_predictor.ensure_loaded = mock.AsyncMock()

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=4,
            max_batch_size=4,
            min_batch_size=1,
            target_memory_fraction=0.80,
            deterministic_representation=True,
        )

        molecules = [MoleculeInfo(smiles="CCO", chem_id=f"mol{i}") for i in range(10)]

        async def _run_two_batches():
            await sched.predict_molecules(
                molecules=molecules[:5],
                aggregate_scores=True,
                app_threshold=0.04374,
                min_nkill=10,
            )
            await sched.predict_molecules(
                molecules=molecules[5:],
                aggregate_scores=True,
                app_threshold=0.04374,
                min_nkill=10,
            )

        _run(_run_two_batches())

        mock_predictor.ensure_loaded.assert_called_once()

    def test_graph_construction_settings_are_forwarded_to_predictor(self) -> None:
        from src.prediction_scheduler import PredictionScheduler
        from src.models import MoleculeInfo

        mock_predictor = mock.AsyncMock()
        mock_predictor.ensure_loaded = mock.AsyncMock()
        mock_predictor.predict.return_value = [{"chem_id": "mol1", "apscore_total": -1.0}]

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=4,
            max_batch_size=4,
            min_batch_size=1,
            target_memory_fraction=0.80,
            deterministic_representation=True,
        )

        molecules = [MoleculeInfo(smiles="CCO", chem_id="mol1")]

        async def _run_predict():
            return await sched.predict_molecules(
                molecules=molecules,
                aggregate_scores=True,
                app_threshold=0.04374,
                min_nkill=10,
                num_graph_workers=3,
                graph_batch_size=64,
                prefetch_batches=5,
            )

        _run(_run_predict())

        _, kwargs = mock_predictor.predict.call_args
        self.assertEqual(kwargs["num_graph_workers"], 3)
        self.assertEqual(kwargs["graph_batch_size"], 64)
        self.assertEqual(kwargs["prefetch_batches"], 5)
        self.assertTrue(kwargs["deterministic_representation"])

    def test_auto_tune_warmup_uses_scheduler_graph_settings(self) -> None:
        from src.prediction_scheduler import PredictionScheduler

        mock_predictor = mock.AsyncMock()
        mock_predictor.ensure_loaded = mock.AsyncMock()
        mock_predictor.device = "cuda:1"
        mock_predictor.predict.return_value = [{"chem_id": "warmup1", "apscore_total": -1.0}]

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=None,
            max_batch_size=256,
            min_batch_size=1,
            target_memory_fraction=0.80,
            num_graph_workers=7,
            graph_batch_size=96,
            prefetch_batches=4,
            deterministic_representation=True,
        )

        with mock.patch("torch.cuda.is_available", return_value=True), \
             mock.patch("torch.cuda.reset_peak_memory_stats"), \
             mock.patch("torch.cuda.memory_allocated", return_value=100), \
             mock.patch("torch.cuda.max_memory_allocated", return_value=200), \
             mock.patch("torch.cuda.mem_get_info", return_value=(1024 * 1024 * 1024, 2 * 1024 * 1024 * 1024)):
            batch_size = _run(sched._auto_tune_batch_size(app_threshold=0.04374, min_nkill=10))

        self.assertGreaterEqual(batch_size, 1)
        _, kwargs = mock_predictor.predict.call_args
        self.assertEqual(kwargs["num_graph_workers"], 7)
        self.assertEqual(kwargs["graph_batch_size"], 96)
        self.assertEqual(kwargs["prefetch_batches"], 4)
        self.assertTrue(kwargs["deterministic_representation"])

    def test_predict_molecules_preserves_multi_chunk_result_order(self) -> None:
        from src.prediction_scheduler import PredictionScheduler
        from src.models import MoleculeInfo

        async def _fake_predict(input_data, **kwargs):
            return [
                {"chem_id": molecule.chem_id, "pred_id": f"{molecule.chem_id}:strain1"}
                for molecule in input_data.molecules or []
            ]

        mock_predictor = mock.AsyncMock()
        mock_predictor.predict.side_effect = _fake_predict
        mock_predictor.ensure_loaded = mock.AsyncMock()

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=2,
            max_batch_size=2,
            min_batch_size=1,
            target_memory_fraction=0.80,
        )

        molecules = [MoleculeInfo(smiles="CCO", chem_id=f"mol{i}") for i in range(5)]

        result = _run(
            sched.predict_molecules(
                molecules=molecules,
                aggregate_scores=False,
                app_threshold=0.04374,
                min_nkill=10,
            )
        )

        self.assertEqual(
            [item["pred_id"] for item in result],
            [f"mol{i}:strain1" for i in range(5)],
        )
        self.assertEqual(len(result), 5)
        self.assertEqual(len({item["chem_id"] for item in result}), 5)

    def test_predict_molecules_normalizes_request_once(self) -> None:
        from src.models import MoleculeInfo, MoleculeInput
        from src.prediction_scheduler import PredictionScheduler
        from src.predictor import AntimicrobialPredictor

        fake_probe = SimpleNamespace(
            preference="pickle",
            pickle_path=Path("/tmp/fake_model.pkl"),
            timber_model_dir=Path("/tmp/fake_timber"),
        )

        with mock.patch("src.predictor.inspect_classifier_backends", return_value=fake_probe):
            predictor = AntimicrobialPredictor()

        predictor._model_loaded = True
        predictor.ensure_loaded = mock.AsyncMock()
        predictor._get_mole_representation = mock.Mock(
            return_value=pd.DataFrame({"feat": [1.0]}, index=["mol1"])
        )
        predictor._add_strains = mock.Mock(
            return_value=pd.DataFrame(
                {"feature": [1.0]},
                index=["mol1:Strain A (NT1)"],
            )
        )
        predictor.model = mock.Mock()
        predictor.model.predict_proba = mock.Mock(return_value=np.array([[0.2, 0.8]]))

        sched = PredictionScheduler(
            predictor=predictor,
            initial_batch_size=4,
            max_batch_size=4,
            min_batch_size=1,
            target_memory_fraction=0.80,
        )

        original_normalize = MoleculeInput.normalize
        with mock.patch.object(
            MoleculeInput,
            "normalize",
            autospec=True,
            side_effect=original_normalize,
        ) as normalize_mock:
            _run(
                sched.predict_molecules(
                    molecules=[MoleculeInfo(smiles="CCO", chem_id="mol1")],
                    aggregate_scores=False,
                    app_threshold=0.04374,
                    min_nkill=10,
                )
            )

        self.assertEqual(normalize_mock.call_count, 1)

    def test_predict_molecules_records_last_profile_when_profiling_enabled(self) -> None:
        from src.prediction_scheduler import PredictionScheduler
        from src.models import MoleculeInfo

        mock_predictor = mock.AsyncMock()
        mock_predictor.ensure_loaded = mock.AsyncMock()
        mock_predictor.predict = mock.AsyncMock(return_value=[{"chem_id": "mol1", "apscore_total": -1.0}])
        mock_predictor.last_profile_snapshot = mock.Mock(
            return_value={
                "representation_seconds": 1.0,
                "prediction_frame_seconds": 0.2,
                "growth_inhibition_seconds": 0.1,
                "aggregate_scores_seconds": 0.3,
                "classifier_stage_seconds": 0.9,
                "classifier_workers": 3,
                "classifier_inflight_batches": 5,
                "xgboost_seconds": 0.5,
                "graph_build": {"graph_total_seconds": 0.25, "graph_items": 1},
            }
        )

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=4,
            max_batch_size=4,
            min_batch_size=1,
            target_memory_fraction=0.80,
        )

        _run(
            sched.predict_molecules(
                molecules=[MoleculeInfo(smiles="CCO", chem_id="mol1")],
                aggregate_scores=True,
                app_threshold=0.04374,
                min_nkill=10,
                enable_profiling=True,
                classifier_workers=3,
                classifier_inflight_batches=5,
            )
        )

        _, kwargs = mock_predictor.predict.call_args
        self.assertTrue(kwargs["enable_profiling"])
        self.assertEqual(kwargs["classifier_workers"], 3)
        self.assertEqual(kwargs["classifier_inflight_batches"], 5)
        self.assertEqual(
            sched.runtime_snapshot()["last_profile"]["graph_build"]["graph_items"],
            1,
        )
        self.assertAlmostEqual(
            sched.runtime_snapshot()["last_profile"]["aggregate_scores_seconds"],
            0.3,
        )
        self.assertAlmostEqual(
            sched.runtime_snapshot()["last_profile"]["classifier_stage_seconds"],
            0.9,
        )
        self.assertEqual(
            sched.runtime_snapshot()["last_profile"]["classifier_workers"],
            3,
        )
        self.assertEqual(
            sched.runtime_snapshot()["last_profile"]["classifier_inflight_batches"],
            5,
        )

    def test_predict_molecules_aggregates_profiles_across_internal_chunks(self) -> None:
        from src.prediction_scheduler import PredictionScheduler
        from src.models import MoleculeInfo

        async def _fake_predict(input_data, **kwargs):
            return [{"chem_id": m.chem_id, "apscore_total": -1.0} for m in input_data.molecules]

        mock_predictor = mock.AsyncMock()
        mock_predictor.ensure_loaded = mock.AsyncMock()
        mock_predictor.predict = mock.AsyncMock(side_effect=_fake_predict)
        mock_predictor.last_profile_snapshot = mock.Mock(
            side_effect=[
                {
                    "representation_seconds": 1.0,
                    "strain_expand_seconds": 0.1,
                    "prediction_frame_seconds": 0.2,
                    "growth_inhibition_seconds": 0.05,
                    "aggregate_scores_seconds": 0.3,
                    "classifier_stage_seconds": 0.95,
                    "classifier_workers": 2,
                    "classifier_inflight_batches": 4,
                    "xgboost_seconds": 0.5,
                    "graph_build": {
                        "graph_items": 2,
                        "rdkit_parse_seconds": 0.2,
                        "add_hs_seconds": 0.1,
                        "atom_feature_seconds": 0.3,
                        "bond_feature_seconds": 0.4,
                        "graph_total_seconds": 1.0,
                        "dataloader_iter_seconds": 1.1,
                        "model_forward_seconds": 0.2,
                        "graph_batch_size": 1024,
                        "graph_workers": 8,
                    },
                },
                {
                    "representation_seconds": 2.0,
                    "strain_expand_seconds": 0.2,
                    "prediction_frame_seconds": 0.3,
                    "growth_inhibition_seconds": 0.06,
                    "aggregate_scores_seconds": 0.4,
                    "classifier_stage_seconds": 1.26,
                    "classifier_workers": 2,
                    "classifier_inflight_batches": 4,
                    "xgboost_seconds": 0.6,
                    "graph_build": {
                        "graph_items": 2,
                        "rdkit_parse_seconds": 0.3,
                        "add_hs_seconds": 0.2,
                        "atom_feature_seconds": 0.4,
                        "bond_feature_seconds": 0.5,
                        "graph_total_seconds": 1.2,
                        "dataloader_iter_seconds": 1.3,
                        "model_forward_seconds": 0.3,
                        "graph_batch_size": 1024,
                        "graph_workers": 8,
                    },
                },
                {
                    "representation_seconds": 3.0,
                    "strain_expand_seconds": 0.3,
                    "prediction_frame_seconds": 0.4,
                    "growth_inhibition_seconds": 0.07,
                    "aggregate_scores_seconds": 0.5,
                    "classifier_stage_seconds": 1.57,
                    "classifier_workers": 2,
                    "classifier_inflight_batches": 4,
                    "xgboost_seconds": 0.7,
                    "graph_build": {
                        "graph_items": 1,
                        "rdkit_parse_seconds": 0.4,
                        "add_hs_seconds": 0.3,
                        "atom_feature_seconds": 0.5,
                        "bond_feature_seconds": 0.6,
                        "graph_total_seconds": 1.4,
                        "dataloader_iter_seconds": 1.5,
                        "model_forward_seconds": 0.4,
                        "graph_batch_size": 1024,
                        "graph_workers": 8,
                    },
                },
            ]
        )

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=2,
            max_batch_size=2,
            min_batch_size=1,
            target_memory_fraction=0.80,
        )

        _run(
            sched.predict_molecules(
                molecules=[MoleculeInfo(smiles="CCO", chem_id=f"mol{i}") for i in range(5)],
                aggregate_scores=True,
                app_threshold=0.04374,
                min_nkill=10,
                enable_profiling=True,
            )
        )

        profile = sched.runtime_snapshot()["last_profile"]
        self.assertEqual(profile["representation_seconds"], 6.0)
        self.assertAlmostEqual(profile["strain_expand_seconds"], 0.6)
        self.assertEqual(profile["xgboost_seconds"], 1.8)
        self.assertAlmostEqual(profile["prediction_frame_seconds"], 0.9)
        self.assertAlmostEqual(profile["growth_inhibition_seconds"], 0.18)
        self.assertAlmostEqual(profile["aggregate_scores_seconds"], 1.2)
        self.assertAlmostEqual(profile["classifier_stage_seconds"], 3.78)
        self.assertEqual(profile["classifier_workers"], 2)
        self.assertEqual(profile["classifier_inflight_batches"], 4)
        self.assertEqual(profile["graph_build"]["graph_items"], 5)
        self.assertAlmostEqual(profile["graph_build"]["graph_total_seconds"], 3.6)
        self.assertAlmostEqual(profile["graph_build"]["bond_feature_seconds"], 1.5)
        self.assertEqual(profile["graph_build"]["graph_batch_size"], 1024)
        self.assertEqual(profile["graph_build"]["graph_workers"], 8)


# ---------------------------------------------------------------------------
# Test: runtime metadata snapshot
# ---------------------------------------------------------------------------

class TestRuntimeSnapshot(unittest.TestCase):
    """Scheduler must expose a snapshot with device-level fields."""

    REQUIRED_FIELDS = {
        "gpu_name",
        "device",
        "torch_version",
        "total_memory_bytes",
        "free_memory_bytes",
        "selected_batch_size",
        "used_cuda",
        "classifier_workers",
        "classifier_inflight_batches",
    }

    def test_snapshot_contains_required_fields(self) -> None:
        from src.prediction_scheduler import PredictionScheduler

        mock_predictor = mock.AsyncMock()
        mock_predictor.ensure_loaded = mock.AsyncMock()

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=32,
            max_batch_size=256,
            min_batch_size=1,
            target_memory_fraction=0.80,
            classifier_workers=6,
            classifier_inflight_batches=9,
        )

        snap = sched.runtime_snapshot()
        for field in self.REQUIRED_FIELDS:
            self.assertIn(field, snap, f"Snapshot missing required field '{field}'")
        self.assertEqual(snap["classifier_workers"], 6)
        self.assertEqual(snap["classifier_inflight_batches"], 9)

    def test_snapshot_selected_batch_size_matches_current(self) -> None:
        from src.prediction_scheduler import PredictionScheduler

        mock_predictor = mock.AsyncMock()
        mock_predictor.ensure_loaded = mock.AsyncMock()

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=64,
            max_batch_size=256,
            min_batch_size=1,
            target_memory_fraction=0.80,
        )

        snap = sched.runtime_snapshot()
        self.assertEqual(snap["selected_batch_size"], sched.current_batch_size)

    def test_snapshot_used_cuda_reflects_torch(self) -> None:
        from src.prediction_scheduler import PredictionScheduler

        mock_predictor = mock.AsyncMock()
        mock_predictor.ensure_loaded = mock.AsyncMock()

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=32,
            max_batch_size=256,
            min_batch_size=1,
            target_memory_fraction=0.80,
        )

        snap = sched.runtime_snapshot()
        self.assertEqual(snap["used_cuda"], torch.cuda.is_available())

    def test_snapshot_reports_cpu_device_when_cuda_unavailable(self) -> None:
        from src.prediction_scheduler import PredictionScheduler

        mock_predictor = mock.AsyncMock()
        mock_predictor.ensure_loaded = mock.AsyncMock()

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=32,
            max_batch_size=256,
            min_batch_size=1,
            target_memory_fraction=0.80,
        )

        with mock.patch("torch.cuda.is_available", return_value=False):
            snap = sched.runtime_snapshot()

        self.assertEqual(snap["device"], "cpu")

    def test_snapshot_prefers_predictor_device_string(self) -> None:
        from src.prediction_scheduler import PredictionScheduler

        mock_predictor = mock.AsyncMock()
        mock_predictor.ensure_loaded = mock.AsyncMock()
        mock_predictor.device = "cuda:1"

        sched = PredictionScheduler(
            predictor=mock_predictor,
            initial_batch_size=32,
            max_batch_size=256,
            min_batch_size=1,
            target_memory_fraction=0.80,
        )

        with mock.patch("torch.cuda.is_available", return_value=True), \
             mock.patch("torch.cuda.mem_get_info", return_value=(0, 0)), \
             mock.patch("torch.cuda.get_device_name", return_value="fake-gpu"), \
             mock.patch("torch.cuda.current_device", return_value=0):
            snap = sched.runtime_snapshot()

        self.assertEqual(snap["device"], "cuda:1")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
