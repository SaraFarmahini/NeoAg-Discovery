#!/usr/bin/env python3
"""Plot distribution of best CDR3β binding scores from strong_binders_best_cdr3b.csv."""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
# Organized layout: strong_binders_best_cdr3b in 04_filter_strong_binders
PROJECT_ROOT = SCRIPT_DIR.parent if (SCRIPT_DIR.parent / "01_nsclc_mutation_data").exists() else SCRIPT_DIR
DIR_04 = PROJECT_ROOT / "04_filter_strong_binders"
CSV = (DIR_04 / "strong_binders_best_cdr3b.csv") if (DIR_04 / "strong_binders_best_cdr3b.csv").exists() else (SCRIPT_DIR / "strong_binders_best_cdr3b.csv")
OUT_DIR = SCRIPT_DIR / "visualizations"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    df = pd.read_csv(CSV)
    if "binding_score" not in df.columns:
        print("No binding_score column in", CSV)
        return
    scores = df["binding_score"].dropna()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Histogram
    ax = axes[0]
    ax.hist(scores, bins=25, color="steelblue", edgecolor="white", alpha=0.85)
    ax.axvline(scores.mean(), color="darkred", linestyle="--", linewidth=2, label=f"Mean = {scores.mean():.2f}")
    ax.axvline(0, color="gray", linestyle="-", linewidth=1)
    ax.set_xlabel("Binding score (logit)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Distribution of best CDR3β binding scores", fontsize=13)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    # KDE-style density
    ax = axes[1]
    ax.hist(scores, bins=25, density=True, color="steelblue", edgecolor="white", alpha=0.6, label="Histogram")
    try:
        scores.plot.kde(ax=ax, color="darkblue", linewidth=2, label="KDE")
    except Exception:
        pass
    ax.axvline(0, color="gray", linestyle="-", linewidth=1)
    ax.set_xlabel("Binding score (logit)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title("Score density", fontsize=13)
    ax.legend()
    ax.grid(alpha=0.3)

    fig.suptitle(f"Strong binders: best TCR per peptide (n = {len(scores)})", fontsize=12, y=1.02)
    fig.tight_layout()
    out_path = OUT_DIR / "cdr3b_binding_score_distribution.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved", out_path)


if __name__ == "__main__":
    main()
