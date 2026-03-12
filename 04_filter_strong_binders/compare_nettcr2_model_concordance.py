#!/usr/bin/env python3
"""
Compare TCR-peptide model binding probabilities with NetTCR-2.0 predictions.

Since no ground truth exists for these pairs, we assess model concordance:
- Spearman ρ and Pearson r (rank and linear correlation)
- Scatter plot and Bland-Altman plot

Inputs:
  - strong_binders_best_cdr3b.csv (binding_probability)
  - nettcr2_test_predictions.csv (prediction)

Outputs:
  - visualizations/model_vs_nettcr2_concordance.png (and .pdf)
  - model_vs_nettcr2_merged.csv (merged data for reference)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

SCRIPT_DIR = Path(__file__).resolve().parent
OUR_CSV = SCRIPT_DIR / "strong_binders_best_cdr3b.csv"
NETTCR2_CSV = SCRIPT_DIR / "nettcr2_test_predictions.csv"
OUT_DIR = SCRIPT_DIR / "visualizations"
MERGED_CSV = SCRIPT_DIR / "model_vs_nettcr2_merged.csv"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def set_pub_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.alpha": 0.25,
            "grid.color": "#9CA3AF",
            "grid.linestyle": "--",
        }
    )


def main() -> None:
    our = pd.read_csv(OUR_CSV)
    nettcr2 = pd.read_csv(NETTCR2_CSV)

    # Align columns: our has best_cdr3b, nettcr2 has CDR3b
    our_aligned = our[["peptide", "best_cdr3b", "binding_probability"]].rename(
        columns={"best_cdr3b": "CDR3b"}
    )
    nettcr2_aligned = nettcr2[["peptide", "CDR3b", "prediction"]].rename(
        columns={"prediction": "nettcr2_prediction"}
    )

    df = our_aligned.merge(
        nettcr2_aligned, on=["peptide", "CDR3b"], how="inner"
    )
    if len(df) == 0:
        raise ValueError("No overlapping (peptide, CDR3b) pairs between files.")
    print(f"Merged {len(df)} pairs on (peptide, CDR3b)")

    prob_ours = df["binding_probability"].values
    prob_nettcr2 = df["nettcr2_prediction"].values

    spearman_r, spearman_p = stats.spearmanr(prob_ours, prob_nettcr2)
    pearson_r, pearson_p = stats.pearsonr(prob_ours, prob_nettcr2)

    print(f"Spearman ρ = {spearman_r:.4f}  (p = {spearman_p:.2e})")
    print(f"Pearson r  = {pearson_r:.4f}  (p = {pearson_p:.2e})")

    df.to_csv(MERGED_CSV, index=False)
    print(f"Merged data saved to {MERGED_CSV}")

    set_pub_style()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    # --- Panel A: Scatter ---
    ax = axes[0]
    ax.scatter(prob_nettcr2, prob_ours, alpha=0.6, s=25, c="#2563EB", edgecolors="white", linewidths=0.5)
    lims = [0, 1]
    ax.plot(lims, lims, "k--", alpha=0.5, label="y=x")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal")
    ax.set_xlabel("NetTCR-2.0 prediction")
    ax.set_ylabel("TCR-peptide model probability")
    ax.set_title(f"A. Scatter (ρ = {spearman_r:.3f})")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)

    # --- Panel B: Bland-Altman ---
    ax = axes[1]
    mean_vals = (prob_ours + prob_nettcr2) / 2
    diff = prob_ours - prob_nettcr2
    mean_diff = np.mean(diff)
    std_diff = np.std(diff)
    ax.scatter(mean_vals, diff, alpha=0.6, s=25, c="#7C3AED", edgecolors="white", linewidths=0.5)
    ax.axhline(mean_diff, color="#DC2626", linestyle="-", label=f"Mean Δ = {mean_diff:.3f}")
    ax.axhline(mean_diff + 1.96 * std_diff, color="#DC2626", linestyle="--", alpha=0.7)
    ax.axhline(mean_diff - 1.96 * std_diff, color="#DC2626", linestyle="--", alpha=0.7)
    ax.axhline(0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Mean (ours + NetTCR-2.0) / 2")
    ax.set_ylabel("Difference (ours − NetTCR-2.0)")
    ax.set_title("B. Bland–Altman")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        out_path = OUT_DIR / f"model_vs_nettcr2_concordance.{ext}"
        plt.savefig(out_path, bbox_inches="tight")
        print(f"Saved {out_path}")
    plt.close()

    print("\nDone. Use Spearman ρ to report model concordance (not validation).")


if __name__ == "__main__":
    main()
