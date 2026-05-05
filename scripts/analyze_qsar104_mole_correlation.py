#!/usr/bin/env python
"""Run QSAR-104 MolE/pMIC correlation analysis.

This script intentionally reuses the production predictor singleton and keeps
all experiment-specific reshaping, statistics, plotting, and reporting here.
"""

from __future__ import annotations

import argparse
import asyncio
import json
import math
import platform
import re
import sys
from pathlib import Path
from typing import Any
import warnings

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.models import MoleculeInfo, MoleculeInput  # noqa: E402
from src.service import get_predictor  # noqa: E402


DEFAULT_INPUT = REPO_ROOT / "data" / "09.qsar_smiles" / "QSAR-104-SMILES.xlsx"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "data" / "09.qsar_smiles" / "analysis"
NT_NUMBER_PATTERN = re.compile(r".*?\((NT\d+)\)")
ECOLI_PATTERN = re.compile(r"escherichia\s+coli|e\.\s*coli|ecoli", re.IGNORECASE)
AGGREGATE_COLUMNS = ["apscore_total", "apscore_gnegative", "apscore_gpositive"]
DEFAULT_APP_THRESHOLD = 0.04374140128493309
OUTPUT_FILENAMES = {
    "matrix": "qsar104_mole_prediction_matrix.csv",
    "matrix_metadata": "qsar104_mole_prediction_matrix_with_metadata.csv",
    "correlations": "qsar104_pmic_correlations.csv",
    "ecoli_correlations": "qsar104_ecoli_correlations.csv",
    "feature_summary": "qsar104_feature_summary.csv",
    "manifest": "qsar104_analysis_manifest.json",
    "report": "EXPERIMENT_REPORT.md",
    "trend_original": "qsar104_trend_original_order_aligned_zscore.png",
    "trend_sorted": "qsar104_trend_pmic_sorted_aligned_zscore.png",
    "spearman_barplot": "qsar104_spearman_correlation_barplot.png",
    "key_scatterplots": "qsar104_key_scatterplots.png",
    "feature_heatmap": "qsar104_mole_feature_spearman_heatmap.png",
}


def parse_sheet(value: str | int) -> str | int:
    """Allow sheet names while treating integer-looking values as indices."""
    if isinstance(value, int):
        return value
    text = str(value)
    if text.lstrip("-").isdigit():
        return int(text)
    return text


def parse_auto_or_int(value: str) -> str | int:
    if value == "auto":
        return value
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be 'auto' or an integer >= 1")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze correlations between QSAR-104 pMIC and MolE predictions."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--sheet", type=parse_sheet, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--app-threshold", type=float, default=DEFAULT_APP_THRESHOLD)
    parser.add_argument("--min-nkill", type=int, default=10)
    parser.add_argument("--prediction-batch-size", type=int, default=16)
    parser.add_argument("--graph-batch-size", type=int, default=256)
    parser.add_argument("--num-graph-workers", type=parse_auto_or_int, default="auto")
    parser.add_argument("--prefetch-batches", type=int, default=2)
    parser.add_argument("--classifier-workers", type=parse_auto_or_int, default="auto")
    parser.add_argument("--classifier-inflight-batches", type=parse_auto_or_int, default="auto")
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260424)
    parser.add_argument(
        "--device",
        choices=["auto", "cpu", "cuda", "cuda:0"],
        default="auto",
    )
    parser.add_argument("--deterministic-representation", action="store_true")
    return parser.parse_args()


def read_qsar_excel(path: Path, sheet: str | int) -> pd.DataFrame:
    """Read the QSAR-104 workbook and enforce the analysis input contract."""
    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {path}")

    df = pd.read_excel(path, sheet_name=sheet, usecols=["SMILES", "name", "pMIC"])
    required = ["SMILES", "name", "pMIC"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Input file is missing required columns: {missing}")

    df = df.loc[:, required].copy()
    df["SMILES"] = df["SMILES"].astype("string").str.strip()
    df["name"] = df["name"].astype("string").str.strip()
    df["pMIC"] = pd.to_numeric(df["pMIC"], errors="raise")

    for column in required:
        if df[column].isna().any():
            bad_rows = (df.index[df[column].isna()] + 2).tolist()
            raise ValueError(f"Column {column!r} contains missing values at Excel rows {bad_rows}")
    for column in ["SMILES", "name"]:
        blank_mask = df[column].astype(str).str.len() == 0
        if blank_mask.any():
            bad_rows = (df.index[blank_mask] + 2).tolist()
            raise ValueError(f"Column {column!r} contains blank values at Excel rows {bad_rows}")

    if df["name"].duplicated().any():
        duplicates = sorted(df.loc[df["name"].duplicated(keep=False), "name"].astype(str).unique())
        raise ValueError(f"Duplicate name values are not allowed: {duplicates}")

    df.insert(0, "input_order", np.arange(1, len(df) + 1, dtype=int))
    df.insert(1, "chem_id", [f"QSAR{index:03d}" for index in range(1, len(df) + 1)])
    return df


async def predict_per_strain(
    qsar_df: pd.DataFrame, args: argparse.Namespace
) -> tuple[pd.DataFrame, list[str], dict[str, str]]:
    """Run MolE/XGBoost per-strain predictions for the QSAR molecules."""
    predictor = get_predictor()
    if args.device != "auto":
        predictor.device = args.device
    await predictor.ensure_loaded()

    if predictor.strain_ohe is None:
        raise RuntimeError("Predictor loaded without strain OHE data")
    if predictor._gram_dict is None:
        raise RuntimeError("Predictor loaded without Gram stain metadata")

    strain_names = [str(value) for value in predictor.strain_ohe.index.tolist()]
    if len(strain_names) != 40:
        raise ValueError(f"Expected 40 strains from predictor, got {len(strain_names)}")

    records: list[dict[str, Any]] = []
    total = len(qsar_df)
    for start in range(0, total, args.prediction_batch_size):
        stop = min(start + args.prediction_batch_size, total)
        batch = qsar_df.iloc[start:stop]
        molecules = [
            MoleculeInfo(smiles=str(row.SMILES), chem_id=str(row.chem_id))
            for row in batch.itertuples(index=False)
        ]
        input_data = MoleculeInput(
            molecules=molecules,
            aggregate_scores=False,
            app_threshold=args.app_threshold,
            min_nkill=args.min_nkill,
        )
        batch_records = await predictor.predict(
            input_data,
            num_graph_workers=args.num_graph_workers,
            graph_batch_size=args.graph_batch_size,
            prefetch_batches=args.prefetch_batches,
            classifier_workers=args.classifier_workers,
            classifier_inflight_batches=args.classifier_inflight_batches,
            deterministic_representation=args.deterministic_representation,
        )
        records.extend(batch_records)
        print(f"Predicted molecules {start + 1}-{stop} of {total}", flush=True)

    pred_df = pd.DataFrame.from_records(records)
    expected_rows = len(qsar_df) * len(strain_names)
    if len(pred_df) != expected_rows:
        raise AssertionError(f"Expected {expected_rows} per-strain rows, got {len(pred_df)}")
    if "pred_id" not in pred_df.columns:
        raise ValueError("Predictor output is missing pred_id")

    split = pred_df["pred_id"].astype(str).str.rsplit(":", n=1, expand=True)
    if split.shape[1] != 2:
        raise ValueError("Unable to split pred_id values into chem_id and strain_name")
    pred_df["chem_id"] = split[0]
    pred_df["strain_name"] = split[1]
    return pred_df, strain_names, dict(predictor._gram_dict)


def nt_number_for_strain(strain_name: str) -> str:
    match = NT_NUMBER_PATTERN.search(str(strain_name))
    return match.group(1) if match else ""


def gram_for_strain(strain_name: str, gram_dict: dict[str, str]) -> str:
    nt_number = nt_number_for_strain(strain_name)
    gram = gram_dict.get(nt_number, "")
    return str(gram).strip().lower()


def build_prediction_matrix(
    qsar_df: pd.DataFrame,
    pred_df: pd.DataFrame,
    strain_names: list[str],
    gram_dict: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build the 44-column matrix and feature metadata."""
    value_column = "antimicrobial_predictive_probability"
    required = {"chem_id", "strain_name", value_column}
    missing = required.difference(pred_df.columns)
    if missing:
        raise ValueError(f"Prediction output missing columns: {sorted(missing)}")

    pivot = pred_df.pivot(index="chem_id", columns="strain_name", values=value_column)
    chem_order = qsar_df["chem_id"].astype(str).tolist()
    pivot = pivot.reindex(index=chem_order, columns=strain_names)
    if pivot.isna().any().any():
        missing_locations = int(pivot.isna().sum().sum())
        raise ValueError(f"Prediction matrix contains {missing_locations} missing values")

    probabilities = pivot.to_numpy(dtype=float, copy=True)
    if not np.isfinite(probabilities).all():
        raise ValueError("Prediction probabilities contain non-finite values")
    if (probabilities <= 0).any():
        raise ValueError("Prediction probabilities must be > 0 to compute log apscores")

    gram_groups = [gram_for_strain(strain_name, gram_dict) for strain_name in strain_names]
    negative_mask = np.array([gram == "negative" for gram in gram_groups], dtype=bool)
    positive_mask = np.array([gram == "positive" for gram in gram_groups], dtype=bool)
    if not negative_mask.any():
        raise ValueError("No Gram-negative strains available for apscore_gnegative")
    if not positive_mask.any():
        raise ValueError("No Gram-positive strains available for apscore_gpositive")

    log_probabilities = np.log(probabilities)
    aggregate_df = pd.DataFrame(
        {
            "apscore_total": log_probabilities.mean(axis=1),
            "apscore_gnegative": log_probabilities[:, negative_mask].mean(axis=1),
            "apscore_gpositive": log_probabilities[:, positive_mask].mean(axis=1),
        },
        index=pivot.index,
    )

    matrix = pd.concat(
        [
            qsar_df.set_index("chem_id").loc[chem_order, ["pMIC"]],
            pivot,
            aggregate_df,
        ],
        axis=1,
    ).reset_index(drop=True)
    expected_columns = ["pMIC", *strain_names, *AGGREGATE_COLUMNS]
    matrix = matrix.loc[:, expected_columns]

    matrix_with_metadata = pd.concat(
        [
            qsar_df.loc[:, ["input_order", "chem_id", "name", "SMILES"]].reset_index(drop=True),
            matrix,
        ],
        axis=1,
    )

    feature_rows: list[dict[str, Any]] = []
    for strain_name, gram in zip(strain_names, gram_groups):
        feature_rows.append(
            {
                "feature": strain_name,
                "feature_type": "strain_probability",
                "gram_group": gram,
                "nt_number": nt_number_for_strain(strain_name),
                "is_ecoli": bool(ECOLI_PATTERN.search(strain_name)),
                "expected_correlation_with_pMIC": "positive",
            }
        )
    feature_rows.extend(
        [
            {
                "feature": "apscore_total",
                "feature_type": "aggregate_apscore",
                "gram_group": "all",
                "nt_number": "",
                "is_ecoli": False,
                "expected_correlation_with_pMIC": "positive",
            },
            {
                "feature": "apscore_gnegative",
                "feature_type": "aggregate_apscore",
                "gram_group": "negative",
                "nt_number": "",
                "is_ecoli": False,
                "expected_correlation_with_pMIC": "positive",
            },
            {
                "feature": "apscore_gpositive",
                "feature_type": "aggregate_apscore",
                "gram_group": "positive",
                "nt_number": "",
                "is_ecoli": False,
                "expected_correlation_with_pMIC": "positive",
            },
        ]
    )
    feature_meta = pd.DataFrame(feature_rows)

    if matrix.shape != (len(qsar_df), 44):
        raise AssertionError(f"Expected matrix shape {(len(qsar_df), 44)}, got {matrix.shape}")
    if matrix.columns[0] != "pMIC":
        raise AssertionError("First matrix column must be pMIC")
    if matrix.columns[-3:].tolist() != AGGREGATE_COLUMNS:
        raise AssertionError(f"Last matrix columns must be {AGGREGATE_COLUMNS}")
    if len(feature_meta) != 43:
        raise AssertionError(f"Expected 43 MolE feature metadata rows, got {len(feature_meta)}")

    return matrix, matrix_with_metadata, feature_meta


def safe_corr(
    x: np.ndarray, y: np.ndarray, method: str
) -> tuple[float, float]:
    """Calculate a correlation while returning NaN for undefined cases."""
    finite_mask = np.isfinite(x) & np.isfinite(y)
    x_values = x[finite_mask]
    y_values = y[finite_mask]
    if len(x_values) < 3 or np.nanstd(x_values) == 0 or np.nanstd(y_values) == 0:
        return math.nan, math.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if method == "pearson":
            result = stats.pearsonr(x_values, y_values)
            return float(result.statistic), float(result.pvalue)
        if method == "spearman":
            result = stats.spearmanr(x_values, y_values)
            return float(result.statistic), float(result.pvalue)
        if method == "kendall":
            result = stats.kendalltau(x_values, y_values)
            return float(result.statistic), float(result.pvalue)
    raise ValueError(f"Unsupported correlation method: {method}")


def bootstrap_corr_ci(
    x: np.ndarray,
    y: np.ndarray,
    method: str,
    iterations: int,
    seed: int,
) -> tuple[float, float]:
    """Bootstrap a percentile 95% confidence interval for paired correlations."""
    finite_mask = np.isfinite(x) & np.isfinite(y)
    x_values = x[finite_mask]
    y_values = y[finite_mask]
    if (
        iterations <= 0
        or len(x_values) < 3
        or np.nanstd(x_values) == 0
        or np.nanstd(y_values) == 0
    ):
        return math.nan, math.nan

    rng = np.random.default_rng(seed)
    n = len(x_values)
    samples: list[float] = []
    for _ in range(iterations):
        indices = rng.integers(0, n, size=n)
        statistic, _ = safe_corr(x_values[indices], y_values[indices], method)
        if np.isfinite(statistic):
            samples.append(statistic)

    if not samples:
        return math.nan, math.nan
    low, high = np.percentile(np.asarray(samples, dtype=float), [2.5, 97.5])
    return float(low), float(high)


def add_fdr_qvalues(df: pd.DataFrame, p_col: str, q_col: str) -> None:
    q_values = np.full(len(df), np.nan, dtype=float)
    p_values = df[p_col].to_numpy(dtype=float)
    valid = np.isfinite(p_values)
    if valid.any():
        q_values[valid] = multipletests(p_values[valid], method="fdr_bh")[1]
    df[q_col] = q_values


def calculate_correlations(
    matrix: pd.DataFrame,
    feature_meta: pd.DataFrame,
    bootstrap_iterations: int,
    bootstrap_seed: int,
) -> pd.DataFrame:
    """Calculate pMIC correlations and direction-normalized ranks."""
    p_mic = matrix["pMIC"].to_numpy(dtype=float)
    rows: list[dict[str, Any]] = []

    for feature_index, meta_row in feature_meta.reset_index(drop=True).iterrows():
        feature = str(meta_row["feature"])
        values = matrix[feature].to_numpy(dtype=float)
        finite_mask = np.isfinite(p_mic) & np.isfinite(values)
        n = int(finite_mask.sum())

        pearson_r, pearson_p = safe_corr(p_mic, values, "pearson")
        spearman_rho, spearman_p = safe_corr(p_mic, values, "spearman")
        kendall_tau, kendall_p = safe_corr(p_mic, values, "kendall")
        pearson_low, pearson_high = bootstrap_corr_ci(
            p_mic,
            values,
            "pearson",
            bootstrap_iterations,
            bootstrap_seed + feature_index * 2,
        )
        spearman_low, spearman_high = bootstrap_corr_ci(
            p_mic,
            values,
            "spearman",
            bootstrap_iterations,
            bootstrap_seed + feature_index * 2 + 1,
        )

        expected = str(meta_row["expected_correlation_with_pMIC"])
        direction_multiplier = -1.0 if expected == "negative" else 1.0
        directional_pearson = (
            float(pearson_r * direction_multiplier) if np.isfinite(pearson_r) else math.nan
        )
        directional_spearman = (
            float(spearman_rho * direction_multiplier) if np.isfinite(spearman_rho) else math.nan
        )
        if expected == "negative":
            direction_matches = bool(np.isfinite(spearman_rho) and spearman_rho <= 0)
        else:
            direction_matches = bool(np.isfinite(spearman_rho) and spearman_rho >= 0)

        rows.append(
            {
                "feature": feature,
                "n": n,
                "pearson_r": pearson_r,
                "pearson_p": pearson_p,
                "pearson_ci95_low": pearson_low,
                "pearson_ci95_high": pearson_high,
                "spearman_rho": spearman_rho,
                "spearman_p": spearman_p,
                "spearman_ci95_low": spearman_low,
                "spearman_ci95_high": spearman_high,
                "kendall_tau": kendall_tau,
                "kendall_p": kendall_p,
                "feature_type": meta_row["feature_type"],
                "gram_group": meta_row["gram_group"],
                "nt_number": meta_row["nt_number"],
                "is_ecoli": bool(meta_row["is_ecoli"]),
                "expected_correlation_with_pMIC": expected,
                "directional_pearson_r": directional_pearson,
                "directional_spearman_rho": directional_spearman,
                "direction_matches_expected": direction_matches,
            }
        )

    corr = pd.DataFrame(rows)
    add_fdr_qvalues(corr, "pearson_p", "pearson_q")
    add_fdr_qvalues(corr, "spearman_p", "spearman_q")
    add_fdr_qvalues(corr, "kendall_p", "kendall_q")

    q_column_order = [
        "feature",
        "n",
        "pearson_r",
        "pearson_p",
        "pearson_q",
        "pearson_ci95_low",
        "pearson_ci95_high",
        "spearman_rho",
        "spearman_p",
        "spearman_q",
        "spearman_ci95_low",
        "spearman_ci95_high",
        "kendall_tau",
        "kendall_p",
        "kendall_q",
        "feature_type",
        "gram_group",
        "nt_number",
        "is_ecoli",
        "expected_correlation_with_pMIC",
        "directional_pearson_r",
        "directional_spearman_rho",
        "direction_matches_expected",
    ]
    corr = corr.loc[:, q_column_order]

    corr["rank_by_directional_spearman"] = (
        corr["directional_spearman_rho"]
        .rank(method="min", ascending=False, na_option="bottom")
        .astype(int)
    )
    corr["rank_by_abs_spearman"] = (
        corr["spearman_rho"]
        .abs()
        .rank(method="min", ascending=False, na_option="bottom")
        .astype(int)
    )
    corr = corr.sort_values(
        ["rank_by_directional_spearman", "feature"], kind="stable"
    ).reset_index(drop=True)

    if len(corr) != 43:
        raise AssertionError(f"Expected 43 correlation rows, got {len(corr)}")
    ecoli_count = int(corr["is_ecoli"].sum())
    if ecoli_count != 2:
        raise AssertionError(f"Expected exactly 2 E. coli rows, got {ecoli_count}")
    return corr


def build_feature_summary(matrix: pd.DataFrame, feature_meta: pd.DataFrame) -> pd.DataFrame:
    """Summarize pMIC and all MolE feature distributions."""
    metadata_by_feature = feature_meta.set_index("feature").to_dict(orient="index")
    rows: list[dict[str, Any]] = []
    for feature in ["pMIC", *feature_meta["feature"].astype(str).tolist()]:
        values = matrix[feature].to_numpy(dtype=float)
        series = pd.Series(values[np.isfinite(values)])
        row: dict[str, Any] = {
            "feature": feature,
            "n": int(series.count()),
            "mean": float(series.mean()),
            "std": float(series.std()),
            "min": float(series.min()),
            "q25": float(series.quantile(0.25)),
            "median": float(series.median()),
            "q75": float(series.quantile(0.75)),
            "max": float(series.max()),
        }
        row.update(
            metadata_by_feature.get(
                feature,
                {
                    "feature_type": "observed_response",
                    "gram_group": "",
                    "nt_number": "",
                    "is_ecoli": False,
                    "expected_correlation_with_pMIC": "",
                },
            )
        )
        rows.append(row)
    return pd.DataFrame(rows)


def directional_feature_values(matrix: pd.DataFrame, feature: str) -> pd.Series:
    return matrix[feature]


def plot_label(feature: str) -> str:
    return feature


def zscore(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    std = float(values.std(ddof=0))
    if not np.isfinite(std) or std == 0:
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - float(values.mean())) / std


def plot_trend(
    matrix_with_metadata: pd.DataFrame,
    feature_meta: pd.DataFrame,
    output_path: Path,
    sort_by_pmic: bool,
) -> None:
    plot_df = matrix_with_metadata.copy()
    if sort_by_pmic:
        plot_df = plot_df.sort_values("pMIC", ascending=True, kind="stable").reset_index(drop=True)
        x_values = np.arange(1, len(plot_df) + 1)
        x_label = "Molecules sorted by ascending pMIC"
    else:
        x_values = plot_df["input_order"].to_numpy(dtype=int)
        x_label = "Input order"

    strain_features = feature_meta.loc[
        feature_meta["feature_type"] == "strain_probability", "feature"
    ].astype(str).tolist()
    ecoli_features = feature_meta.loc[feature_meta["is_ecoli"], "feature"].astype(str).tolist()

    fig, ax = plt.subplots(figsize=(16, 8))
    for feature in strain_features:
        if feature in ecoli_features:
            continue
        ax.plot(
            x_values,
            zscore(plot_df[feature]),
            color="#9ca3af",
            alpha=0.22,
            linewidth=0.8,
            zorder=1,
        )

    ecoli_colors = ["#f97316", "#ea580c"]
    for index, feature in enumerate(ecoli_features):
        ax.plot(
            x_values,
            zscore(plot_df[feature]),
            color=ecoli_colors[index % len(ecoli_colors)],
            linewidth=2.0,
            label=feature,
            zorder=3,
        )

    aggregate_styles = {
        "apscore_total": ("#16a34a", "--"),
        "apscore_gnegative": ("#15803d", "--"),
        "apscore_gpositive": ("#65a30d", "--"),
    }
    for feature, (color, linestyle) in aggregate_styles.items():
        ax.plot(
            x_values,
            zscore(plot_df[feature]),
            color=color,
            linestyle=linestyle,
            linewidth=2.0,
            label=feature,
            zorder=3,
        )

    ax.plot(
        x_values,
        zscore(plot_df["pMIC"]),
        color="black",
        linewidth=3.0,
        label="pMIC",
        zorder=4,
    )
    ax.set_title(
        "Direction-aligned z-score trend by pMIC order"
        if sort_by_pmic
        else "Direction-aligned z-score trend by original order"
    )
    ax.set_xlabel(x_label)
    ax.set_ylabel("Z-score (higher = stronger activity after direction alignment)")
    ax.axhline(0, color="#111827", linewidth=0.8, alpha=0.35)
    ax.grid(axis="y", color="#e5e7eb", linewidth=0.8)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_spearman_barplot(corr: pd.DataFrame, output_path: Path) -> None:
    plot_df = corr.sort_values("directional_spearman_rho", ascending=True).copy()
    colors = np.where(
        plot_df["is_ecoli"],
        "#f97316",
        np.where(plot_df["feature_type"] == "aggregate_apscore", "#16a34a", "#9ca3af"),
    )
    labels = [plot_label(feature) for feature in plot_df["feature"].astype(str)]
    y = np.arange(len(plot_df))

    fig, ax = plt.subplots(figsize=(12, 14))
    ax.barh(y, plot_df["directional_spearman_rho"], color=colors, alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color="#111827", linewidth=0.9)
    ax.set_xlabel("Directional Spearman rho (higher = stronger aligned association)")
    ax.set_title("pMIC vs MolE feature directional Spearman correlations")
    for idx, (_, row) in enumerate(plot_df.iterrows()):
        if np.isfinite(row["spearman_q"]) and row["spearman_q"] < 0.05:
            offset = 0.012 if row["directional_spearman_rho"] >= 0 else -0.025
            ax.text(
                row["directional_spearman_rho"] + offset,
                idx,
                "*",
                va="center",
                ha="left" if offset > 0 else "right",
                fontsize=12,
                color="#111827",
            )
    ax.grid(axis="x", color="#e5e7eb", linewidth=0.8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_key_scatterplots(
    matrix: pd.DataFrame,
    corr: pd.DataFrame,
    output_path: Path,
) -> None:
    ecoli_features = corr.loc[corr["is_ecoli"], "feature"].astype(str).tolist()
    requested_features: list[str] = [*ecoli_features]
    top_feature = str(corr.sort_values("rank_by_directional_spearman").iloc[0]["feature"])
    if top_feature not in requested_features:
        requested_features.append(top_feature)
    for feature in ["apscore_gnegative", "apscore_total"]:
        if feature not in requested_features:
            requested_features.append(feature)

    n_panels = len(requested_features)
    n_cols = 2
    n_rows = math.ceil(n_panels / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(13, 5 * n_rows), squeeze=False)
    for ax, feature in zip(axes.ravel(), requested_features):
        row = corr.loc[corr["feature"] == feature].iloc[0]
        y_values = directional_feature_values(matrix, feature)
        sns.regplot(
            x=matrix["pMIC"],
            y=y_values,
            ax=ax,
            scatter_kws={"s": 26, "alpha": 0.76, "color": "#2563eb"},
            line_kws={"color": "#dc2626", "linewidth": 1.5},
            ci=95,
        )
        ax.set_xlabel("pMIC")
        ax.set_ylabel(plot_label(feature))
        ax.set_title(
            f"{plot_label(feature)}\n"
            f"Spearman rho={row['spearman_rho']:.3f}, "
            f"q={row['spearman_q']:.3g}, "
            f"directional={row['directional_spearman_rho']:.3f}"
        )
        ax.grid(color="#e5e7eb", linewidth=0.8)

    for ax in axes.ravel()[n_panels:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_feature_heatmap(
    matrix: pd.DataFrame,
    feature_meta: pd.DataFrame,
    output_path: Path,
) -> None:
    aligned = pd.DataFrame(
        {
            plot_label(feature): directional_feature_values(matrix, feature)
            for feature in feature_meta["feature"].astype(str)
        }
    )
    spearman = aligned.corr(method="spearman")
    fig, ax = plt.subplots(figsize=(18, 16))
    sns.heatmap(
        spearman,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.05,
        cbar_kws={"label": "Spearman correlation"},
    )
    ax.set_title("Spearman correlation heatmap across direction-aligned MolE features")
    ax.tick_params(axis="x", labelrotation=90, labelsize=6)
    ax.tick_params(axis="y", labelsize=6)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def generate_plots(
    matrix: pd.DataFrame,
    matrix_with_metadata: pd.DataFrame,
    feature_meta: pd.DataFrame,
    corr: pd.DataFrame,
    output_paths: dict[str, Path],
) -> None:
    sns.set_theme(style="whitegrid")
    plot_trend(matrix_with_metadata, feature_meta, output_paths["trend_original"], sort_by_pmic=False)
    plot_trend(matrix_with_metadata, feature_meta, output_paths["trend_sorted"], sort_by_pmic=True)
    plot_spearman_barplot(corr, output_paths["spearman_barplot"])
    plot_key_scatterplots(matrix, corr, output_paths["key_scatterplots"])
    plot_feature_heatmap(matrix, feature_meta, output_paths["feature_heatmap"])


def finite_range(values: pd.Series) -> tuple[float, float]:
    finite = values.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.min(finite)), float(np.max(finite))


def markdown_table(df: pd.DataFrame, columns: list[str], max_rows: int | None = None) -> str:
    table_df = df.loc[:, columns].copy()
    if max_rows is not None:
        table_df = table_df.head(max_rows)
    for column in table_df.columns:
        if pd.api.types.is_float_dtype(table_df[column]):
            table_df[column] = table_df[column].map(
                lambda value: "" if pd.isna(value) else f"{value:.4g}"
            )
        else:
            table_df[column] = table_df[column].map(
                lambda value: "" if pd.isna(value) else str(value)
            )

    headers = [str(column) for column in table_df.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in table_df.itertuples(index=False):
        values = [str(value).replace("\n", " ") for value in row]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def make_manifest(
    args: argparse.Namespace,
    qsar_df: pd.DataFrame,
    strain_names: list[str],
    corr: pd.DataFrame,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    top10_cols = [
        "feature",
        "feature_type",
        "gram_group",
        "is_ecoli",
        "spearman_rho",
        "spearman_q",
        "directional_spearman_rho",
        "rank_by_directional_spearman",
    ]
    return {
        "input_path": str(args.input.resolve()),
        "sheet": args.sheet,
        "output_directory": str(args.output_dir.resolve()),
        "n_molecules": int(len(qsar_df)),
        "n_strains": int(len(strain_names)),
        "n_mole_features": 43,
        "app_threshold": args.app_threshold,
        "min_nkill": args.min_nkill,
        "prediction_batch_size": args.prediction_batch_size,
        "graph_batch_size": args.graph_batch_size,
        "num_graph_workers": args.num_graph_workers,
        "prefetch_batches": args.prefetch_batches,
        "classifier_workers": args.classifier_workers,
        "classifier_inflight_batches": args.classifier_inflight_batches,
        "bootstrap_iterations": args.bootstrap_iterations,
        "bootstrap_seed": args.bootstrap_seed,
        "device": args.device,
        "deterministic_representation": args.deterministic_representation,
        "python_version": sys.version,
        "platform": platform.platform(),
        "top_10_directional_spearman": corr.sort_values("rank_by_directional_spearman")
        .head(10)
        .loc[:, top10_cols]
        .to_dict(orient="records"),
        "output_paths": {key: str(path.resolve()) for key, path in output_paths.items()},
    }


def generate_report(
    args: argparse.Namespace,
    matrix: pd.DataFrame,
    matrix_with_metadata: pd.DataFrame,
    feature_summary: pd.DataFrame,
    corr: pd.DataFrame,
    output_paths: dict[str, Path],
) -> str:
    p_mic_min, p_mic_max = finite_range(matrix["pMIC"])
    probability_columns = [
        column for column in matrix.columns if column not in ["pMIC", *AGGREGATE_COLUMNS]
    ]
    probability_min = float(matrix[probability_columns].min().min())
    probability_max = float(matrix[probability_columns].max().max())
    probability_summary = feature_summary.loc[
        feature_summary["feature_type"] == "strain_probability", ["mean", "std", "min", "max"]
    ]
    significant_spearman = int((corr["spearman_q"] < 0.05).sum())
    significant_pearson = int((corr["pearson_q"] < 0.05).sum())
    significant_kendall = int((corr["kendall_q"] < 0.05).sum())
    ecoli = corr.loc[corr["is_ecoli"]].sort_values("feature")
    output_file_list = "\n".join(
        f"- `{path.relative_to(REPO_ROOT)}`" for path in output_paths.values()
    )
    top_table = markdown_table(
        corr.sort_values("rank_by_directional_spearman"),
        [
            "rank_by_directional_spearman",
            "feature",
            "feature_type",
            "gram_group",
            "spearman_rho",
            "spearman_q",
            "directional_spearman_rho",
            "direction_matches_expected",
        ],
        max_rows=10,
    )
    ecoli_table = markdown_table(
        ecoli,
        [
            "feature",
            "pearson_r",
            "pearson_q",
            "pearson_ci95_low",
            "pearson_ci95_high",
            "spearman_rho",
            "spearman_q",
            "spearman_ci95_low",
            "spearman_ci95_high",
            "directional_spearman_rho",
            "rank_by_directional_spearman",
        ],
    )
    aggregate_table = markdown_table(
        corr.loc[corr["feature"].isin(AGGREGATE_COLUMNS)].sort_values(
            "rank_by_directional_spearman"
        ),
        [
            "feature",
            "spearman_rho",
            "spearman_q",
            "directional_spearman_rho",
            "rank_by_directional_spearman",
        ],
    )

    command = "pixi run python scripts/analyze_qsar104_mole_correlation.py"
    if args.device != "auto":
        command += f" --device {args.device}"
    if args.prediction_batch_size != 16:
        command += f" --prediction-batch-size {args.prediction_batch_size}"
    if args.graph_batch_size != 256:
        command += f" --graph-batch-size {args.graph_batch_size}"

    return f"""# QSAR-104 MolE/pMIC 相关性分析报告

## 实验目的

本实验评估 QSAR-104 数据集中 104 个分子的真实 E. coli pMIC 与 MolE/XGBoost 模型输出之间是否存在探索性相关趋势。比较对象包括 40 个菌株的 `antimicrobial_predictive_probability` 和 3 个聚合 antimicrobial potential score：`apscore_total`、`apscore_gnegative`、`apscore_gpositive`。

## 实验数据

- 输入文件：`{args.input.relative_to(REPO_ROOT)}`
- 工作表：`{args.sheet}`
- 分子数量：{len(matrix_with_metadata)}
- pMIC 范围：{p_mic_min:.4g} 到 {p_mic_max:.4g}
- MolE/XGBoost 菌株数量：{len(probability_columns)}
- E. coli 菌株：Escherichia coli ED1a (NT5078) 与 Escherichia coli IAI1 (NT5077)

## 实验方案

1. 使用 `src.service.get_predictor()` 获取生产预测器单例，不修改 `src/predictor.py`、`src/models.py` 或 `src/service.py`。
2. 读取 `SMILES`、`name`、`pMIC` 三列，保留原始顺序，并生成 `QSAR001` 到 `QSAR104` 的 `chem_id`。
3. 对每个分子运行 MolE 表征和 XGBoost per-strain 预测，得到 40 个菌株概率。
4. 构建 44 列主矩阵：`pMIC`、40 个菌株概率、3 个原始 `apscore_*`。
5. 在脚本内按 `mean(log(probability))` 计算 `apscore_total`、`apscore_gnegative` 和 `apscore_gpositive`。
6. 对 43 个 MolE 特征分别计算 Pearson、Spearman、Kendall 相关性，使用 Benjamini-Hochberg 方法计算 FDR q 值，并用 bootstrap 估计 Pearson/Spearman 95% CI。
7. 对方向进行统一：菌株概率和 `apscore_*` 都保持原方向。`apscore_* = mean(log(probability))`，由于 probability < 1 通常为负值；数值越接近 0 表示几何平均预测概率越高。

## 实验过程

- 运行命令：`{command}`
- 输出目录：`{args.output_dir.relative_to(REPO_ROOT)}`
- bootstrap 次数：{args.bootstrap_iterations}
- bootstrap seed：{args.bootstrap_seed}
- APP 阈值：{args.app_threshold}
- min_nkill：{args.min_nkill}

输出文件：

{output_file_list}

## 结果解读

按方向统一后的 Spearman 相关性排序，前 10 个特征如下：

{top_table}

FDR 显著性数量：

- Pearson q < 0.05：{significant_pearson}
- Spearman q < 0.05：{significant_spearman}
- Kendall q < 0.05：{significant_kendall}

`apscore_*` 汇总特征结果如下，CSV 和图中均保留原始 `apscore_*` 方向：

{aggregate_table}

### E. coli 专项解读

两个 E. coli 菌株的相关性结果如下：

{ecoli_table}

E. coli 行用于检验本数据集真实 pMIC 与模型中最接近目标病原的 per-strain 概率是否一致。若其方向统一 Spearman 排名靠前，说明模型输出与本 QSAR 系列的 E. coli 活性有较强单调趋势；若排名靠后或 FDR 不显著，则说明该系列中的真实 pMIC 变化不能被这两个训练面板 E. coli 菌株概率稳定解释。

### 分数方向

菌株概率的方向是“越高表示预测抑制越强”，因此预期与 pMIC 正相关。`apscore_total`、`apscore_gnegative`、`apscore_gpositive` 是 `mean(log(probability))`；因为 probability 位于 0 到 1，log 后通常是负值。按当前代码公式解释，`apscore_*` 越接近 0 表示几何平均预测概率越高，越负表示几何平均预测概率越低。因此图和方向统一排名不再对 `apscore_*` 取负，CSV 也保留原始 `apscore_*` 数值。

### 预测值分布

40 个菌株预测概率的总体范围为 {probability_min:.4g} 到 {probability_max:.4g}。菌株概率分布摘要：平均值范围 {probability_summary["mean"].min():.4g} 到 {probability_summary["mean"].max():.4g}，标准差范围 {probability_summary["std"].min():.4g} 到 {probability_summary["std"].max():.4g}。如果某些菌株概率分布过窄，相关系数会受到限制，即使分子真实 pMIC 存在变化，模型输出也可能无法呈现高相关。

## 结论边界

本实验是探索性相关性分析，不能证明 MolE/XGBoost 能够准确回归 QSAR-104 的 pMIC，也不能把相关性解释为因果关系。pMIC 来自特定 E. coli 实验条件，而 MolE/XGBoost 输出来自 Maier 训练面板和预训练分子表征；菌株差异、实验体系差异、化合物系列偏置、预测概率校准不足和样本量有限都会影响结果解释。建议把这些输出作为后续模型诊断和特征比较依据，而不是作为单独的抗菌活性判定标准。
"""


def write_outputs(
    args: argparse.Namespace,
    matrix: pd.DataFrame,
    matrix_with_metadata: pd.DataFrame,
    corr: pd.DataFrame,
    feature_summary: pd.DataFrame,
    manifest: dict[str, Any],
    report: str,
    output_paths: dict[str, Path],
) -> None:
    matrix.to_csv(output_paths["matrix"], index=False)
    matrix_with_metadata.to_csv(output_paths["matrix_metadata"], index=False)
    corr.to_csv(output_paths["correlations"], index=False)
    ecoli = corr.loc[corr["is_ecoli"]].sort_values("feature")
    if len(ecoli) != 2:
        raise AssertionError(f"Expected exactly 2 E. coli rows, got {len(ecoli)}")
    ecoli.to_csv(output_paths["ecoli_correlations"], index=False)
    feature_summary.to_csv(output_paths["feature_summary"], index=False)
    output_paths["manifest"].write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    output_paths["report"].write_text(report, encoding="utf-8")


async def async_main() -> None:
    args = parse_args()
    args.input = args.input.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    output_paths = {
        key: args.output_dir / filename for key, filename in OUTPUT_FILENAMES.items()
    }

    qsar_df = read_qsar_excel(args.input, args.sheet)
    pred_df, strain_names, gram_dict = await predict_per_strain(qsar_df, args)
    matrix, matrix_with_metadata, feature_meta = build_prediction_matrix(
        qsar_df,
        pred_df,
        strain_names,
        gram_dict,
    )
    corr = calculate_correlations(
        matrix,
        feature_meta,
        args.bootstrap_iterations,
        args.bootstrap_seed,
    )
    feature_summary = build_feature_summary(matrix, feature_meta)
    generate_plots(matrix, matrix_with_metadata, feature_meta, corr, output_paths)
    manifest = make_manifest(args, qsar_df, strain_names, corr, output_paths)
    report = generate_report(args, matrix, matrix_with_metadata, feature_summary, corr, output_paths)
    write_outputs(args, matrix, matrix_with_metadata, corr, feature_summary, manifest, report, output_paths)

    print(f"Wrote QSAR-104 MolE correlation analysis to: {args.output_dir}")
    print(f"Main matrix: {output_paths['matrix']}")
    print(f"Report: {output_paths['report']}")


def main() -> None:
    asyncio.run(async_main())


if __name__ == "__main__":
    main()
