"""Selective pathogen panel scoring logic.

Computes pathogen/commensal soft and hard scores from per-strain
antimicrobial probabilities and a panel configuration file.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PanelConfig:
    """Parsed pathogen-selective panel configuration."""

    label: str
    mode: str
    app_threshold: float
    tau: float
    lambda_: float
    primary_pathogen_panel: List[str]
    secondary_pathogen_panel: List[str]
    commensal_sparing_panel: List[str]

    @property
    def pathogen_panel(self) -> List[str]:
        """Union of primary and secondary pathogen panels."""
        return list(dict.fromkeys(
            self.primary_pathogen_panel + self.secondary_pathogen_panel
        ))

    def all_strains(self) -> List[str]:
        """All strains across all panels, in definition order."""
        return list(dict.fromkeys(
            self.primary_pathogen_panel
            + self.secondary_pathogen_panel
            + self.commensal_sparing_panel
        ))


def load_panel_config(path: str | Path) -> PanelConfig:
    """Load and validate a panel configuration JSON file."""
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return validate_panel_config(raw)


def validate_panel_config(raw: Dict) -> PanelConfig:
    """Validate a raw panel config dict and return a PanelConfig."""
    required = (
        "mode",
        "label",
        "app_threshold",
        "tau",
        "lambda",
        "primary_pathogen_panel",
        "secondary_pathogen_panel",
        "commensal_sparing_panel",
    )
    missing = [key for key in required if key not in raw]
    if missing:
        raise ValueError(f"panel config missing required keys: {missing}")

    mode = str(raw["mode"]).strip()
    if not mode:
        raise ValueError("panel config 'mode' must be a non-empty string")

    label = str(raw["label"]).strip()
    if not label:
        raise ValueError("panel config 'label' must be a non-empty string")

    app_threshold = float(raw["app_threshold"])
    tau = float(raw["tau"])
    lambda_ = float(raw["lambda"])

    if tau <= 0:
        raise ValueError(f"tau must be positive, got {tau}")

    primary = _validate_strain_list(raw["primary_pathogen_panel"], "primary_pathogen_panel")
    secondary = _validate_strain_list(raw["secondary_pathogen_panel"], "secondary_pathogen_panel")
    commensal = _validate_strain_list(raw["commensal_sparing_panel"], "commensal_sparing_panel")

    all_strains = primary + secondary + commensal
    if len(set(all_strains)) != len(all_strains):
        seen: set[str] = set()
        dupes = [s for s in all_strains if s in seen or seen.add(s)]  # type: ignore[func-returns-value]
        raise ValueError(f"panel config contains duplicate strains: {dupes}")

    return PanelConfig(
        label=label,
        mode=mode,
        app_threshold=app_threshold,
        tau=tau,
        lambda_=lambda_,
        primary_pathogen_panel=primary,
        secondary_pathogen_panel=secondary,
        commensal_sparing_panel=commensal,
    )


def _validate_strain_list(value: object, field_name: str) -> List[str]:
    if not isinstance(value, list) or len(value) == 0:
        raise ValueError(f"'{field_name}' must be a non-empty list of strain names")
    cleaned = [str(item).strip() for item in value]
    if any(not item for item in cleaned):
        raise ValueError(f"'{field_name}' entries must be non-empty strings")
    if len(set(cleaned)) != len(cleaned):
        raise ValueError(f"'{field_name}' must contain unique strain names")
    return cleaned


def sigmoid(x: np.ndarray) -> np.ndarray:
    """Numerically stable sigmoid function."""
    x = np.asarray(x, dtype=np.float64)
    positive = x >= 0
    result = np.empty_like(x)
    result[positive] = 1.0 / (1.0 + np.exp(-x[positive]))
    exp_x = np.exp(x[~positive])
    result[~positive] = exp_x / (1.0 + exp_x)
    return result


def compute_panel_scores(
    probabilities: np.ndarray,
    strain_names: Sequence[str],
    panel: PanelConfig,
    *,
    threshold: float | None = None,
    tau: float | None = None,
    lambda_: float | None = None,
) -> Dict[str, float]:
    """Compute panel scores for a single molecule from per-strain probabilities.

    Parameters
    ----------
    probabilities : 1-D array of per-strain probabilities, aligned with *strain_names*.
    strain_names : strain name for each column in *probabilities*.
    panel : parsed panel configuration.
    threshold : override ``app_threshold`` from *panel*.
    tau : override ``tau`` from *panel*.
    lambda_ : override ``lambda_`` from *panel*.

    Returns
    -------
    dict with keys ``pathogen_soft``, ``pathogen_hard``, ``commensal_soft``,
    ``selectivity_score``.
    """
    probabilities = np.asarray(probabilities, dtype=np.float64)
    if probabilities.ndim != 1:
        raise ValueError("probabilities must be a 1-D array")

    t = threshold if threshold is not None else panel.app_threshold
    tau_val = tau if tau is not None else panel.tau
    lambda_val = lambda_ if lambda_ is not None else panel.lambda_

    strain_to_idx: Dict[str, int] = {}
    for i, name in enumerate(strain_names):
        strain_to_idx[str(name)] = i

    pathogen_strains = panel.pathogen_panel
    commensal_strains = panel.commensal_sparing_panel

    missing_pathogens = [s for s in pathogen_strains if s not in strain_to_idx]
    if missing_pathogens:
        raise ValueError(
            f"pathogen strains not found in data: {missing_pathogens}"
        )

    pathogen_indices = [strain_to_idx[s] for s in pathogen_strains]
    p_pathogen = probabilities[pathogen_indices]

    pathogen_soft = float(np.mean(sigmoid((p_pathogen - t) / tau_val)))
    pathogen_hard = int(np.sum(p_pathogen >= t))

    if commensal_strains:
        missing_commensals = [s for s in commensal_strains if s not in strain_to_idx]
        if missing_commensals:
            raise ValueError(
                f"commensal strains not found in data: {missing_commensals}"
            )
        commensal_indices = [strain_to_idx[s] for s in commensal_strains]
        p_commensal = probabilities[commensal_indices]
        commensal_soft = float(np.mean(sigmoid((p_commensal - t) / tau_val)))
    else:
        commensal_soft = 0.0

    selectivity_score = pathogen_soft - lambda_val * commensal_soft

    return {
        "pathogen_soft": pathogen_soft,
        "pathogen_hard": pathogen_hard,
        "commensal_soft": commensal_soft,
        "selectivity_score": selectivity_score,
    }


def compute_panel_scores_from_dataframe(
    df: pd.DataFrame,
    panel: PanelConfig,
    *,
    pred_id_col: str = "pred_id",
    probability_col: str = "1",
    threshold: float | None = None,
    tau: float | None = None,
    lambda_: float | None = None,
) -> pd.DataFrame:
    """Compute panel scores for all molecules in a per-strain prediction dataframe.

    The input *df* must contain a ``pred_id`` column formatted as
    ``chem_id:strain_name`` and a probability column.

    Returns a DataFrame indexed by ``chem_id`` with columns
    ``pathogen_soft``, ``pathogen_hard``, ``commensal_soft``, ``selectivity_score``.
    """
    if pred_id_col not in df.columns:
        raise ValueError(f"input dataframe missing '{pred_id_col}' column")
    if probability_col not in df.columns:
        raise ValueError(f"input dataframe missing '{probability_col}' column")

    work = df[[pred_id_col, probability_col]].copy()
    split = work[pred_id_col].astype(str).str.rsplit(":", n=1, expand=True)
    if split.shape[1] != 2:
        raise ValueError(
            f"'{pred_id_col}' column must be in 'chem_id:strain_name' format"
        )
    work["chem_id"] = split[0]
    work["strain_name"] = split[1]

    pivot = work.pivot_table(
        index="chem_id",
        columns="strain_name",
        values=probability_col,
        aggfunc="first",
    )
    strain_names = list(pivot.columns)
    matrix = pivot.to_numpy(dtype=np.float64)

    results: list[dict[str, float]] = []
    for row_idx in range(matrix.shape[0]):
        scores = compute_panel_scores(
            matrix[row_idx],
            strain_names,
            panel,
            threshold=threshold,
            tau=tau,
            lambda_=lambda_,
        )
        results.append(scores)

    return pd.DataFrame(results, index=pivot.index)
