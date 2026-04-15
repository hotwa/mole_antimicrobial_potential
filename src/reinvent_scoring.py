"""Reward shaping helpers for REINVENT4 integration."""

from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Sequence

from src.models import ReinventScoreRequest, ScoreObjective


def _sigmoid(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-value))


def _geometric_mean(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("values must be non-empty")
    if any(value <= 0 for value in values):
        return 0.0
    return math.exp(sum(math.log(value) for value in values) / len(values))


def _parse_predictions(
    raw_items: Sequence[Mapping[str, Any]],
) -> Dict[str, Dict[str, float]]:
    """Group predictor output as chem_id -> strain_name -> probability."""
    grouped: Dict[str, Dict[str, float]] = {}
    for item in raw_items:
        pred_id = item.get("pred_id")
        probability = item.get("antimicrobial_predictive_probability")
        if not isinstance(pred_id, str) or ":" not in pred_id:
            raise ValueError("Each prediction item must contain pred_id='chem_id:strain'")
        if probability is None:
            raise ValueError("Each prediction item must contain antimicrobial_predictive_probability")
        chem_id, strain_name = pred_id.split(":", 1)
        grouped.setdefault(chem_id, {})[strain_name] = float(probability)
    return grouped


def _resolve_panel(
    objective: ScoreObjective,
    available_strains: Sequence[str],
) -> List[str]:
    selected = objective.resolved_strains()
    if objective.mode == "single_strain":
        if not selected:
            raise ValueError("single_strain objective requires strain or strains")
    elif objective.mode == "broad_spectrum_soft":
        if not selected:
            return list(available_strains)
    else:  # pragma: no cover - defensive guard
        raise ValueError(f"Unsupported objective mode: {objective.mode}")

    assert selected is not None
    missing = [strain for strain in selected if strain not in available_strains]
    if missing:
        raise ValueError(f"Unknown strains requested: {missing}")
    return selected


def _score_single_strain(
    chem_id: str,
    probabilities: Mapping[str, float],
    panel: Sequence[str],
    objective: ScoreObjective,
) -> Dict[str, Any]:
    weights = objective.normalized_weights(len(panel))
    selected_probabilities = [float(probabilities[strain]) for strain in panel]
    weighted_probability = sum(
        weight * probability
        for weight, probability in zip(weights, selected_probabilities)
    )
    return {
        "chem_id": chem_id,
        "score": weighted_probability,
        "objective_mode": objective.mode,
        "selected_strains": list(panel),
        "weights": weights,
        "selected_probabilities": dict(zip(panel, selected_probabilities)),
        "weighted_probability": weighted_probability,
        "panel_mean_probability": sum(selected_probabilities) / len(selected_probabilities),
        "panel_gmean_probability": _geometric_mean(selected_probabilities),
    }


def _score_broad_spectrum_soft(
    chem_id: str,
    probabilities: Mapping[str, float],
    panel: Sequence[str],
    request: ReinventScoreRequest,
) -> Dict[str, Any]:
    selected_probabilities = [float(probabilities[strain]) for strain in panel]
    soft_hits = [
        _sigmoid((probability - request.app_threshold) / request.objective.tau)
        for probability in selected_probabilities
    ]
    soft_count = sum(soft_hits)
    soft_ratio = soft_count / len(panel)
    return {
        "chem_id": chem_id,
        "score": soft_ratio,
        "objective_mode": request.objective.mode,
        "selected_strains": list(panel),
        "selected_probabilities": dict(zip(panel, selected_probabilities)),
        "soft_inhibition_values": dict(zip(panel, soft_hits)),
        "soft_inhibition_count": soft_count,
        "soft_inhibition_ratio": soft_ratio,
        "panel_mean_probability": sum(selected_probabilities) / len(selected_probabilities),
        "panel_gmean_probability": _geometric_mean(selected_probabilities),
        "app_threshold": request.app_threshold,
        "tau": request.objective.tau,
        "min_nkill_reference": request.min_nkill,
    }


def score_reinvent_predictions(
    raw_items: Sequence[Mapping[str, Any]],
    request: ReinventScoreRequest,
) -> List[Dict[str, Any]]:
    """
    Convert per-strain predictions into continuous REINVENT4 rewards.

    All returned `score` values are normalized to the interval [0, 1].
    """

    grouped = _parse_predictions(raw_items)
    scored_items: List[Dict[str, Any]] = []

    for chem_id, probabilities in grouped.items():
        panel = _resolve_panel(request.objective, list(probabilities))
        if request.objective.mode == "single_strain":
            scored = _score_single_strain(
                chem_id=chem_id,
                probabilities=probabilities,
                panel=panel,
                objective=request.objective,
            )
        else:
            scored = _score_broad_spectrum_soft(
                chem_id=chem_id,
                probabilities=probabilities,
                panel=panel,
                request=request,
            )
        scored_items.append(scored)

    return scored_items
