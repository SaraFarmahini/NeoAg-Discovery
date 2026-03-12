#!/usr/bin/env python3
"""
Publication-quality visualization of binding score distributions:
- Panel A: Logit score distribution (histogram + KDE)
- Panel B: Binding probability distribution (histogram + KDE)

Input: strong_binders_best_cdr3b.csv (requires binding_score and binding_probability columns)
Output: step6_NSCLC/visualizations/strong_binders_logit_probability_distributions.(png|pdf)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

SCRIPT_DIR = Path(__file__).resolve().parent
# Organized layout: strong_binders_best_cdr3b in 04_filter_strong_binders
PROJECT_ROOT = SCRIPT_DIR.parent if (SCRIPT_DIR.parent / "01_nsclc_mutation_data").exists() else SCRIPT_DIR
DIR_04 = PROJECT_ROOT / "04_filter_strong_binders"
CSV_PATH = (DIR_04 / "strong_binders_best_cdr3b.csv") if (DIR_04 / "strong_binders_best_cdr3b.csv").exists() else (SCRIPT_DIR / "strong_binders_best_cdr3b.csv")
OUT_DIR = SCRIPT_DIR / "visualizations"
OUT_DIR.mkdir(parents=True, exist_ok=True)
THESIS_IMG_DIR = (
    SCRIPT_DIR.parents[1] / "Thesis_Template_DE_EN" / "BA_MA_English" / "images" / "c0"
)
THESIS_IMG_DIR.mkdir(parents=True, exist_ok=True)


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


def plot_histogram_kde(
    ax: plt.Axes,
    values: np.ndarray,
    xlabel: str,
    title: str,
    color: str = "#2563EB",
    mean_color: str = "#7C3AED",
) -> None:
    """Histogram + KDE with mean and median lines."""
    values = np.array(values, dtype=float)
    values = values[~np.isnan(values)]

    ax.hist(
        values,
        bins=30,
        color=color,
        alpha=0.75,
        edgecolor="white",
        linewidth=0.5,
        density=True,
        label="Histogram",
    )
    try:
        kde = gaussian_kde(values, bw_method=0.15)
        x_min, x_max = values.min(), values.max()
        x_pad = max(0.1 * (x_max - x_min), 0.1)
        x = np.linspace(x_min - x_pad, x_max + x_pad, 300)
        ax.plot(x, kde(x), color=mean_color, linewidth=2.5, label="KDE")
    except Exception:
        pass

    mean_v = float(np.mean(values))
    med_v = float(np.median(values))
    ax.axvline(mean_v, color=mean_color, linestyle="--", linewidth=1.6, label=f"Mean={mean_v:.2f}")
    ax.axvline(med_v, color="#1F2937", linestyle=":", linewidth=1.6, label=f"Median={med_v:.2f}")

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.legend(frameon=False, loc="upper right", fontsize=8)
    ax.set_ylim(0, None)
    ax.grid(True, alpha=0.25, linestyle="--")


def main() -> None:
    set_pub_style()
    df = pd.read_csv(CSV_PATH)

    if "binding_score" not in df.columns or "binding_probability" not in df.columns:
        raise ValueError(
            f"CSV must have 'binding_score' and 'binding_probability' columns. Found: {list(df.columns)}"
        )

    logits = df["binding_score"].dropna().values
    probs = df["binding_probability"].dropna().values
    n = len(logits)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.5))

    plot_histogram_kde(
        axes[0],
        logits,
        xlabel="Binding score (logit)",
        title="A) Logit score distribution",
        color="#2563EB",
    )
    plot_histogram_kde(
        axes[1],
        probs,
        xlabel="Binding probability",
        title="B) Probability distribution",
        color="#0F766E",
    )

    fig.suptitle(
        f"NSCLC strong-binder TCR–peptide scores (n = {n:,})",
        fontsize=12,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout()

    stem = "strong_binders_logit_probability_distributions"
    for out_dir in (OUT_DIR, THESIS_IMG_DIR):
        fig.savefig(out_dir / f"{stem}.png", bbox_inches="tight", facecolor="white")
        fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved {OUT_DIR / stem}.png and .pdf")

    # Second visualization: logit vs probability scatter with sigmoid curve
    fig2, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(logits, probs, alpha=0.5, s=20, color="#2563EB", edgecolors="white", linewidth=0.3)
    x_sigmoid = np.linspace(logits.min() - 0.5, logits.max() + 0.5, 200)
    y_sigmoid = 1 / (1 + np.exp(-x_sigmoid))
    ax.plot(x_sigmoid, y_sigmoid, color="#7C3AED", linewidth=2.5, label="Sigmoid")
    ax.axhline(0.5, color="#9CA3AF", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline(0, color="#9CA3AF", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.set_xlabel("Binding score (logit)")
    ax.set_ylabel("Binding probability")
    ax.set_title("Logit–probability transformation")
    ax.legend(frameon=False, loc="lower right")
    ax.set_xlim(x_sigmoid.min(), x_sigmoid.max())
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25, linestyle="--")
    fig2.tight_layout()

    stem2 = "strong_binders_logit_vs_probability"
    for out_dir in (OUT_DIR, THESIS_IMG_DIR):
        fig2.savefig(out_dir / f"{stem2}.png", bbox_inches="tight", facecolor="white")
        fig2.savefig(out_dir / f"{stem2}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig2)
    print(f"Saved {OUT_DIR / stem2}.png and .pdf")
    print(f"Thesis images: {THESIS_IMG_DIR} ({stem}.png, {stem2}.png)")


if __name__ == "__main__":
    main()
