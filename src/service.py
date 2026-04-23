"""Predictor and scheduler singleton management."""

from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from src.predictor import AntimicrobialPredictor
    from src.prediction_scheduler import PredictionScheduler

_predictor: Optional["AntimicrobialPredictor"] = None
_scheduler: Optional["PredictionScheduler"] = None


def get_predictor() -> "AntimicrobialPredictor":
    """Return the singleton predictor instance."""
    global _predictor
    if _predictor is None:
        from src.predictor import AntimicrobialPredictor

        _predictor = AntimicrobialPredictor()
    return _predictor

def get_scheduler(
    initial_batch_size: Optional[int] = None,
    max_batch_size: Optional[int] = None,
    min_batch_size: Optional[int] = None,
    target_memory_fraction: Optional[float] = None,
    num_graph_workers: int | str = "auto",
    graph_batch_size: int = 1024,
    prefetch_batches: int = 2,
    deterministic_representation: bool = False,
    classifier_workers: int | str = "auto",
    classifier_inflight_batches: int | str = "auto",
) -> "PredictionScheduler":
    """Return the singleton prediction scheduler."""
    global _scheduler
    if _scheduler is None:
        from src.prediction_scheduler import (
            PredictionScheduler,
            DEFAULT_MAX_BATCH_SIZE,
            DEFAULT_MIN_BATCH_SIZE,
            DEFAULT_TARGET_MEMORY_FRACTION,
        )

        _scheduler = PredictionScheduler(
            predictor=get_predictor(),
            initial_batch_size=initial_batch_size,
            max_batch_size=max_batch_size if max_batch_size is not None else DEFAULT_MAX_BATCH_SIZE,
            min_batch_size=min_batch_size if min_batch_size is not None else DEFAULT_MIN_BATCH_SIZE,
            target_memory_fraction=target_memory_fraction if target_memory_fraction is not None else DEFAULT_TARGET_MEMORY_FRACTION,
            num_graph_workers=num_graph_workers,
            graph_batch_size=graph_batch_size,
            prefetch_batches=prefetch_batches,
            deterministic_representation=deterministic_representation,
            classifier_workers=classifier_workers,
            classifier_inflight_batches=classifier_inflight_batches,
        )
    return _scheduler


def create_scheduler(
    initial_batch_size: Optional[int] = None,
    max_batch_size: Optional[int] = None,
    min_batch_size: Optional[int] = None,
    target_memory_fraction: Optional[float] = None,
    num_graph_workers: int | str = "auto",
    graph_batch_size: int = 1024,
    prefetch_batches: int = 2,
    deterministic_representation: bool = False,
    classifier_workers: int | str = "auto",
    classifier_inflight_batches: int | str = "auto",
) -> "PredictionScheduler":
    """Create a fresh scheduler bound to the shared predictor singleton."""
    from src.prediction_scheduler import (
        PredictionScheduler,
        DEFAULT_MAX_BATCH_SIZE,
        DEFAULT_MIN_BATCH_SIZE,
        DEFAULT_TARGET_MEMORY_FRACTION,
    )

    return PredictionScheduler(
        predictor=get_predictor(),
        initial_batch_size=initial_batch_size,
        max_batch_size=max_batch_size if max_batch_size is not None else DEFAULT_MAX_BATCH_SIZE,
        min_batch_size=min_batch_size if min_batch_size is not None else DEFAULT_MIN_BATCH_SIZE,
        target_memory_fraction=target_memory_fraction if target_memory_fraction is not None else DEFAULT_TARGET_MEMORY_FRACTION,
        num_graph_workers=num_graph_workers,
        graph_batch_size=graph_batch_size,
        prefetch_batches=prefetch_batches,
        deterministic_representation=deterministic_representation,
        classifier_workers=classifier_workers,
        classifier_inflight_batches=classifier_inflight_batches,
    )
