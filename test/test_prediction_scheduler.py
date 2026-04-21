"""Tests for src/prediction_scheduler.py – Task 1 (failing first, then pass after impl)."""

from __future__ import annotations

import asyncio
import unittest
from unittest import mock

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

    def _make_scheduler(self, free_bytes: int, total_bytes: int):
        from src.prediction_scheduler import PredictionScheduler

        sched = PredictionScheduler.__new__(PredictionScheduler)
        sched._warmup_bytes_per_item = 500_000          # 500 KB / item (mocked)
        sched._batch_size: int | None = None            # not yet cached
        sched.max_batch_size = 2048
        sched.target_memory_fraction = 0.80
        sched._free_bytes = free_bytes
        sched._total_bytes = total_bytes
        return sched

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

        async def _predict(input_data):
            call_count_box[0] += 1
            if call_count_box[0] == 1:
                raise torch.cuda.OutOfMemoryError("CUDA out of memory")  # type: ignore[attr-defined]
            # Second call succeeds with dummy data
            return [{"chem_id": "mol1", "apscore_total": -1.5}]

        predictor = mock.AsyncMock()
        predictor.predict.side_effect = _predict
        predictor.ensure_loaded = mock.AsyncMock()
        return predictor

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

        async def _always_oom(input_data):
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

        async def _fake_predict(input_data):
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
        )

        snap = sched.runtime_snapshot()
        for field in self.REQUIRED_FIELDS:
            self.assertIn(field, snap, f"Snapshot missing required field '{field}'")

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


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
