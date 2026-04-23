"""Single-GPU adaptive batch scheduler for MolE antimicrobial prediction.

Architecture
------------
- One :class:`PredictionScheduler` instance per process (one GPU per process).
- Auto-tunes ``current_batch_size`` from ``torch.cuda.mem_get_info()`` at
  construction time and updates it downward on :class:`torch.cuda.OutOfMemoryError`.
- Keeps the model resident between calls by delegating to a shared
  :class:`~src.predictor.AntimicrobialPredictor` whose ``ensure_loaded()`` is
  called exactly once per scheduler lifetime.
- Exposes :py:meth:`runtime_snapshot` so manifest writers can record which
  device and batch size were actually used.

Multi-GPU note
--------------
Multi-GPU scaling is intentionally out of scope for v1.  Start one process per
GPU and set ``CUDA_VISIBLE_DEVICES`` to pin each process to a specific device.
"""

from __future__ import annotations

import logging
import math
import copy
import inspect
from typing import Any, Dict, List, Optional, TYPE_CHECKING

import torch

if TYPE_CHECKING:
    from src.predictor import AntimicrobialPredictor
    from src.models import MoleculeInfo

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_TARGET_MEMORY_FRACTION: float = 0.80
DEFAULT_MAX_BATCH_SIZE: int = 16384
DEFAULT_MIN_BATCH_SIZE: int = 1
# Estimated bytes consumed by one molecule during MolE forward pass (8000-dim
# float32 embedding = 32 KB; the graph tensors and XGBoost feature matrix add
# more — 500 KB is a conservative default that avoids OOM on most GPUs).
DEFAULT_BYTES_PER_ITEM: int = 500_000


def _initial_profile_aggregate() -> Dict[str, Any]:
    return {
        "representation_seconds": 0.0,
        "strain_expand_seconds": 0.0,
        "xgboost_seconds": 0.0,
        "prediction_frame_seconds": 0.0,
        "growth_inhibition_seconds": 0.0,
        "aggregate_scores_seconds": 0.0,
        "result_records_seconds": 0.0,
        "graph_build": {
            "graph_items": 0,
            "rdkit_parse_seconds": 0.0,
            "add_hs_seconds": 0.0,
            "atom_feature_seconds": 0.0,
            "bond_feature_seconds": 0.0,
            "graph_total_seconds": 0.0,
            "dataloader_setup_seconds": 0.0,
            "dataloader_iter_seconds": 0.0,
            "model_forward_seconds": 0.0,
            "graph_batch_size": None,
            "graph_workers": None,
        },
    }


def _merge_profile_aggregate(aggregate: Dict[str, Any], profile: Dict[str, Any] | None) -> None:
    if not profile:
        return

    for key in (
        "representation_seconds",
        "strain_expand_seconds",
        "xgboost_seconds",
        "prediction_frame_seconds",
        "growth_inhibition_seconds",
        "aggregate_scores_seconds",
        "result_records_seconds",
    ):
        aggregate[key] += float(profile.get(key, 0.0))

    graph_profile = profile.get("graph_build")
    if not isinstance(graph_profile, dict):
        return

    aggregate_graph = aggregate["graph_build"]
    for key in (
        "graph_items",
        "rdkit_parse_seconds",
        "add_hs_seconds",
        "atom_feature_seconds",
        "bond_feature_seconds",
        "graph_total_seconds",
        "dataloader_setup_seconds",
        "dataloader_iter_seconds",
        "model_forward_seconds",
    ):
        if key == "graph_items":
            aggregate_graph[key] += int(graph_profile.get(key, 0))
        else:
            aggregate_graph[key] += float(graph_profile.get(key, 0.0))

    for key in ("graph_batch_size", "graph_workers"):
        if aggregate_graph.get(key) is None and key in graph_profile:
            aggregate_graph[key] = graph_profile.get(key)


class PredictionScheduler:
    """Adaptive single-GPU batch scheduler.

    Parameters
    ----------
    predictor:
        Shared :class:`~src.predictor.AntimicrobialPredictor` instance.
        ``ensure_loaded()`` is called once on the first
        :py:meth:`predict_molecules` call.
    initial_batch_size:
        Starting batch size.  If not given, it is auto-computed from
        ``torch.cuda.mem_get_info()`` at construction time.
    max_batch_size:
        Hard upper cap on batch size.
    min_batch_size:
        Smallest batch size to try before giving up on OOM.
    target_memory_fraction:
        Fraction of free GPU memory to aim for when auto-computing batch size.
    bytes_per_item:
        Expected GPU bytes consumed per molecule.  Used when
        ``initial_batch_size`` is *not* supplied (auto-tune mode).
    """

    def __init__(
        self,
        predictor: "AntimicrobialPredictor",
        initial_batch_size: Optional[int] = None,
        max_batch_size: int = DEFAULT_MAX_BATCH_SIZE,
        min_batch_size: int = DEFAULT_MIN_BATCH_SIZE,
        target_memory_fraction: float = DEFAULT_TARGET_MEMORY_FRACTION,
        num_graph_workers: int | str = "auto",
        graph_batch_size: int = 1024,
        prefetch_batches: int = 2,
        deterministic_representation: bool = False,
    ) -> None:
        self._predictor = predictor
        self.max_batch_size = max_batch_size
        self.min_batch_size = min_batch_size
        self.target_memory_fraction = target_memory_fraction
        self.num_graph_workers = num_graph_workers
        self.graph_batch_size = graph_batch_size
        self.prefetch_batches = prefetch_batches
        self.deterministic_representation = deterministic_representation
        self._last_profile: Optional[Dict[str, Any]] = None

        # Whether ensure_loaded() has already been awaited once
        self._model_ensured: bool = False

        if initial_batch_size is not None:
            self.current_batch_size: Optional[int] = max(min_batch_size, min(initial_batch_size, max_batch_size))
            logger.info(
                "PredictionScheduler initialised with fixed batch_size=%d",
                self.current_batch_size,
            )
        else:
            self.current_batch_size = None
            logger.info("PredictionScheduler initialised with dynamic auto-batching (will tune on first run)")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    async def predict_molecules(
        self,
        molecules: "List[MoleculeInfo] | None" = None,
        aggregate_scores: bool | None = None,
        app_threshold: float | None = None,
        min_nkill: int | None = None,
        input_data=None,
        already_normalized: bool = False,
        num_graph_workers: int | str | None = None,
        graph_batch_size: Optional[int] = None,
        prefetch_batches: Optional[int] = None,
        enable_profiling: bool = False,
        deterministic_representation: Optional[bool] = None,
    ) -> List[Dict[str, Any]]:
        """Run prediction over *molecules* with adaptive batching.

        The MolE model is loaded lazily the first time this method is called
        and then stays resident for every subsequent call on this scheduler.

        On :class:`torch.cuda.OutOfMemoryError` the current batch is split in
        half and retried until either the prediction succeeds or
        ``min_batch_size`` is reached, at which point the error propagates.
        """
        from src.models import MoleculeInput

        if not self._model_ensured:
            await self._predictor.ensure_loaded()
            self._model_ensured = True

        if self.current_batch_size is None:
            self.current_batch_size = await self._auto_tune_batch_size(app_threshold, min_nkill)
            logger.info(
                "Auto-tuned batch size to %d for device %s",
                self.current_batch_size,
                self._device_name(),
            )

        if input_data is not None:
            request = input_data if already_normalized else input_data.normalize()
        else:
            if molecules is None or aggregate_scores is None or app_threshold is None or min_nkill is None:
                raise ValueError(
                    "predict_molecules requires either input_data or all of molecules, aggregate_scores, app_threshold, min_nkill"
                )
            request = MoleculeInput(
                molecules=molecules,
                aggregate_scores=aggregate_scores,
                app_threshold=app_threshold,
                min_nkill=min_nkill,
            ).normalize()

        molecules = list(request.molecules or [])
        all_results: List[Dict[str, Any]] = []
        profile_aggregate = _initial_profile_aggregate() if enable_profiling else None

        _oom_type = getattr(torch.cuda, "OutOfMemoryError", None)
        unprocessed = molecules.copy()

        while unprocessed:
            batch_size = self.current_batch_size
            chunk = unprocessed[:batch_size]

            chunk_request = MoleculeInput.model_construct(
                molecules=chunk,
                smiles=None,
                chem_id=None,
                aggregate_scores=request.aggregate_scores,
                app_threshold=request.app_threshold,
                min_nkill=request.min_nkill,
            )

            try:
                result = await self._predictor.predict(
                    chunk_request,
                    num_graph_workers=(
                        self.num_graph_workers if num_graph_workers is None else num_graph_workers
                    ),
                    graph_batch_size=(
                        self.graph_batch_size if graph_batch_size is None else graph_batch_size
                    ),
                    prefetch_batches=(
                        self.prefetch_batches if prefetch_batches is None else prefetch_batches
                    ),
                    already_normalized=True,
                    enable_profiling=enable_profiling,
                    deterministic_representation=(
                        self.deterministic_representation
                        if deterministic_representation is None
                        else deterministic_representation
                    ),
                )
                if enable_profiling:
                    _merge_profile_aggregate(profile_aggregate, self._predictor.last_profile_snapshot())
                all_results.extend(result)
                unprocessed = unprocessed[len(chunk):]
            except Exception as exc:
                is_oom = _oom_type is not None and isinstance(exc, _oom_type)
                if not is_oom:
                    raise

                new_size = batch_size // 2
                if new_size < self.min_batch_size:
                    logger.error("OOM at minimum batch size %d; cannot reduce further.", batch_size)
                    raise

                logger.warning("OOM at batch size %d; retrying with batch size %d.", batch_size, new_size)
                self.current_batch_size = new_size

        if enable_profiling:
            self._last_profile = profile_aggregate

        return all_results

    def runtime_snapshot(self) -> Dict[str, Any]:
        """Return a serialisable dict describing the current GPU/batch state.

        The returned mapping always contains:

        - ``gpu_name``
        - ``device``
        - ``torch_version``
        - ``total_memory_bytes``
        - ``free_memory_bytes``
        - ``selected_batch_size``
        - ``used_cuda``
        """
        used_cuda = torch.cuda.is_available()
        device_str = self._device_name()
        predictor_status = None
        get_status = getattr(self._predictor, "get_status", None)
        if self._model_ensured and callable(get_status) and not inspect.iscoroutinefunction(get_status):
            predictor_status = get_status()

        if used_cuda:
            try:
                free_bytes, total_bytes = torch.cuda.mem_get_info()
            except Exception:
                free_bytes, total_bytes = 0, 0
            try:
                gpu_name = torch.cuda.get_device_name()
            except Exception:
                gpu_name = "unknown"
        else:
            free_bytes = 0
            total_bytes = 0
            gpu_name = "cpu"

        return {
            "gpu_name": gpu_name,
            "device": device_str,
            "torch_version": torch.__version__,
            "total_memory_bytes": total_bytes,
            "free_memory_bytes": free_bytes,
            "selected_batch_size": self.current_batch_size,
            "used_cuda": used_cuda,
            "deterministic_representation": self.deterministic_representation,
            "classifier_backend": getattr(predictor_status, "classifier_backend", None),
            "classifier_backend_preference": getattr(predictor_status, "classifier_backend_preference", None),
            "classifier_backend_path": getattr(predictor_status, "classifier_backend_path", None),
            "last_profile": copy.deepcopy(self._last_profile),
        }

    # ------------------------------------------------------------------
    # Static helper (pure function, easily unit-testable)
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_batch_size(
        free_bytes: int,
        bytes_per_item: int,
        target_fraction: float,
        max_batch_size: int,
    ) -> int:
        """Compute a batch size that uses at most *target_fraction* of *free_bytes*.

        Always returns at least 1.
        """
        usable = free_bytes * target_fraction
        computed = max(1, int(math.floor(usable / max(bytes_per_item, 1))))
        return min(computed, max_batch_size)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    async def _auto_tune_batch_size(self, app_threshold: float, min_nkill: int) -> int:
        """Run a warmup pass to estimate per-item memory and return a tuned batch size."""
        if not torch.cuda.is_available():
            return min(64, self.max_batch_size)

        from src.models import MoleculeInput, MoleculeInfo

        # Measure memory before
        torch.cuda.reset_peak_memory_stats()
        mem_before = torch.cuda.memory_allocated()

        # Warmup with 2 dummy items
        warmup_size = 2
        dummies = [
            MoleculeInfo(smiles="CCO", chem_id="warmup1"),
            MoleculeInfo(smiles="CCN", chem_id="warmup2"),
        ]
        request = MoleculeInput(
            molecules=dummies,
            aggregate_scores=True,
            app_threshold=app_threshold,
            min_nkill=min_nkill,
        )

        try:
            await self._predictor.predict(
                request,
                num_graph_workers=self.num_graph_workers,
                graph_batch_size=self.graph_batch_size,
                prefetch_batches=self.prefetch_batches,
                deterministic_representation=self.deterministic_representation,
            )
        except Exception as exc:
            logger.warning("Warmup prediction failed, falling back to static default: %s", exc)
            return min(64, self.max_batch_size)

        mem_after_peak = torch.cuda.max_memory_allocated()
        cost_for_batch = mem_after_peak - mem_before
        
        # Add conservative buffer and ensure it's at least DEFAULT_BYTES_PER_ITEM
        bytes_per_item = max(int(cost_for_batch / warmup_size), DEFAULT_BYTES_PER_ITEM)

        try:
            free_bytes, _total = torch.cuda.mem_get_info()
        except Exception:
            free_bytes = 0

        if free_bytes > 0:
            return self._compute_batch_size(
                free_bytes=free_bytes,
                bytes_per_item=bytes_per_item,
                target_fraction=self.target_memory_fraction,
                max_batch_size=self.max_batch_size,
            )
        
        return min(64, self.max_batch_size)

    def _device_name(self) -> str:
        predictor_device = getattr(self._predictor, "device", None)
        if isinstance(predictor_device, str) and predictor_device:
            return str(predictor_device)
        if torch.cuda.is_available():
            try:
                return f"cuda:{torch.cuda.current_device()}"
            except Exception:
                return "cuda"
        return "cpu"
