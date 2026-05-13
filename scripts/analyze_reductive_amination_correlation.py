#!/usr/bin/env python3
"""
Analyze correlations between MolE predictions and pMIC_pred for reductive amination products.

Usage:
    pixi run python scripts/analyze_reductive_amination_correlation.py
"""

import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = REPO_ROOT / "data" / "11.macrolide_reductive_amination"
OUTPUT_DIR = INPUT_DIR / "analysis"

FILES = [
    "tylosin_seed_amine_products_with_mole_predictions.csv",
    "tylvalosin_seed_amine_products_with_mole_predictions.csv",
]

PMIC_COL = "pMIC_pred（μmol/ml）"
MOLE_COLS = [
    "apscore_total",
    "apscore_gnegative",
    "apscore_gpositive",
    "ginhib_total",
    "ginhib_gnegative",
    "ginhib_gpositive",
    "broad_spectrum",
]


def safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> tuple[float, float]:
    """Calculate correlation, returning NaN for undefined cases."""
    finite_mask = np.isfinite(x) & np.isfinite(y)
    x_vals = x[finite_mask]
    y_vals = y[finite_mask]
    if len(x_vals) < 3 or np.std(x_vals) == 0 or np.std(y_vals) == 0:
        return np.nan, np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if method == "pearson":
            r, p = stats.pearsonr(x_vals, y_vals)
        elif method == "spearman":
            r, p = stats.spearmanr(x_vals, y_vals)
        elif method == "kendall":
            r, p = stats.kendalltau(x_vals, y_vals)
        else:
            raise ValueError(f"Unknown method: {method}")
    return float(r), float(p)


def analyze_file(filepath: Path) -> pd.DataFrame:
    """Calculate correlations for a single file."""
    df = pd.read_csv(filepath, encoding="gb18030")

    if PMIC_COL not in df.columns:
        # Try alternative column name
        alt_col = "origin_pMIC_pred（μmol/ml）"
        if alt_col in df.columns:
            pmic_col = alt_col
        else:
            raise ValueError(f"pMIC column not found in {filepath.name}")
    else:
        pmic_col = PMIC_COL

    pmic = df[pmic_col].to_numpy(dtype=float)

    rows = []
    for col in MOLE_COLS:
        if col not in df.columns:
            continue
        values = df[col].to_numpy(dtype=float)

        pearson_r, pearson_p = safe_corr(pmic, values, "pearson")
        spearman_rho, spearman_p = safe_corr(pmic, values, "spearman")
        kendall_tau, kendall_p = safe_corr(pmic, values, "kendall")

        rows.append({
            "feature": col,
            "n": int(np.sum(np.isfinite(pmic) & np.isfinite(values))),
            "pearson_r": pearson_r,
            "pearson_p": pearson_p,
            "spearman_rho": spearman_rho,
            "spearman_p": spearman_p,
            "kendall_tau": kendall_tau,
            "kendall_p": kendall_p,
        })

    return pd.DataFrame(rows)


def plot_correlation_matrix(df: pd.DataFrame, output_path: Path, title: str) -> None:
    """Plot correlation heatmap between pMIC and MolE features."""
    # Select relevant columns
    cols_to_plot = [PMIC_COL] + MOLE_COLS
    cols_available = [c for c in cols_to_plot if c in df.columns]

    if len(cols_available) < 2:
        return

    corr_matrix = df[cols_available].corr(method="spearman")

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(
        corr_matrix,
        ax=ax,
        annot=True,
        fmt=".3f",
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
    )
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_scatter_grid(df: pd.DataFrame, output_path: Path, title: str) -> None:
    """Plot scatter plots of pMIC vs each MolE feature."""
    pmic_col = PMIC_COL if PMIC_COL in df.columns else "origin_pMIC_pred（μmol/ml）"

    mole_cols_available = [c for c in MOLE_COLS if c in df.columns]
    n_cols = 3
    n_rows = (len(mole_cols_available) + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    axes = axes.flatten() if n_rows > 1 else [axes] if n_rows == 1 else axes.flatten()

    for idx, col in enumerate(mole_cols_available):
        ax = axes[idx]
        x = df[pmic_col].to_numpy(dtype=float)
        y = df[col].to_numpy(dtype=float)

        finite_mask = np.isfinite(x) & np.isfinite(y)
        x_plot = x[finite_mask]
        y_plot = y[finite_mask]

        ax.scatter(x_plot, y_plot, alpha=0.6, s=30, color="#2563eb")

        # Add regression line
        if len(x_plot) > 2:
            z = np.polyfit(x_plot, y_plot, 1)
            p = np.poly1d(z)
            x_line = np.linspace(x_plot.min(), x_plot.max(), 100)
            ax.plot(x_line, p(x_line), color="#dc2626", linewidth=2)

        # Calculate correlation
        r, p_val = safe_corr(x_plot, y_plot, "spearman")
        ax.set_xlabel("pMIC_pred")
        ax.set_ylabel(col)
        ax.set_title(f"{col}\nSpearman rho={r:.3f}, p={p_val:.3g}")
        ax.grid(alpha=0.3)

    # Hide unused axes
    for idx in range(len(mole_cols_available), len(axes)):
        axes[idx].axis("off")

    fig.suptitle(title, fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    """Main entry point."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Reductive Amination Products: MolE vs pMIC Correlation Analysis")
    print("=" * 60)

    all_results = []

    for filename in FILES:
        filepath = INPUT_DIR / filename
        if not filepath.exists():
            print(f"Warning: File not found: {filepath}")
            continue

        print(f"\nAnalyzing: {filename}")

        # Calculate correlations
        corr_df = analyze_file(filepath)
        corr_df["file"] = filename
        all_results.append(corr_df)

        # Print results
        print("\nCorrelation Results:")
        print("-" * 80)
        for _, row in corr_df.iterrows():
            sig = ""
            if row["spearman_p"] < 0.001:
                sig = "***"
            elif row["spearman_p"] < 0.01:
                sig = "**"
            elif row["spearman_p"] < 0.05:
                sig = "*"

            print(f"  {row['feature']:25s} | "
                  f"Pearson r={row['pearson_r']:+.4f} (p={row['pearson_p']:.4g}) | "
                  f"Spearman rho={row['spearman_rho']:+.4f} (p={row['spearman_p']:.4g}) {sig}")

        # Generate plots
        df = pd.read_csv(filepath, encoding="gb18030")

        prefix = filename.replace("_with_mole_predictions.csv", "")
        plot_correlation_matrix(
            df,
            OUTPUT_DIR / f"{prefix}_correlation_heatmap.png",
            f"{prefix}: Spearman Correlation Matrix"
        )
        plot_scatter_grid(
            df,
            OUTPUT_DIR / f"{prefix}_scatter_plots.png",
            f"{prefix}: pMIC vs MolE Features"
        )

    # Combine results
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined_path = OUTPUT_DIR / "combined_correlations.csv"
        combined.to_csv(combined_path, index=False)
        print(f"\n\nCombined results saved to: {combined_path}")

    print("\n" + "=" * 60)
    print("Analysis Complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
