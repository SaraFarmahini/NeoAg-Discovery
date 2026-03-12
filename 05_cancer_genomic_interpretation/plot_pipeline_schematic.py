#!/usr/bin/env python3
"""
Generate a clean schematic of the NSCLC pipeline data flow.
Shows only counts at each stage, no tool versions.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / "visualizations"
OUT_DIR.mkdir(parents=True, exist_ok=True)
THESIS_IMG_DIR = SCRIPT_DIR.parent.parent / "Thesis_Template_DE_EN" / "BA_MA_English" / "images" / "c0"
# Organized layout: data in 01, 02, 04
PROJECT_ROOT = SCRIPT_DIR.parent if (SCRIPT_DIR.parent / "01_nsclc_mutation_data").exists() else SCRIPT_DIR
DIR_01 = PROJECT_ROOT / "01_nsclc_mutation_data"
DIR_02 = PROJECT_ROOT / "02_generate_mutant_peptides"
DIR_04 = PROJECT_ROOT / "04_filter_strong_binders"

# Colour palette for pipeline overview (teal/blue gradient)
PIPELINE_COLORS = [
    "#0b1d78",
    "#0045a5",
    "#0069c0",
    "#008ac5",
    "#00a9b5",
    "#00c698",
    "#1fe074",  # green
]


def load_counts():
    """Load counts from data files."""
    cosmic = len(pd.read_csv(DIR_01 / "Cosmic_ResistanceMutations_v101_GRCh38.tsv", sep="\t")) if (DIR_01 / "Cosmic_ResistanceMutations_v101_GRCh38.tsv").exists() else 5199
    n_fasta = len(list((DIR_01 / "protein_sequences").glob("*.fasta"))) if (DIR_01 / "protein_sequences").exists() else 70
    primary = len(pd.read_csv(DIR_02 / "neoantigens_primary.csv")) if (DIR_02 / "neoantigens_primary.csv").exists() else 382
    df_9 = pd.read_csv(DIR_02 / "neoantigen_9mers.csv") if (DIR_02 / "neoantigen_9mers.csv").exists() else None
    n_9mer_pairs = len(df_9) if df_9 is not None else 3418
    n_unique_9mers = df_9["neoantigen_peptide"].nunique() if df_9 is not None else 1739
    sb = pd.read_csv(DIR_04 / "strong_binders_peptide_wt_metrics_unique.csv") if (DIR_04 / "strong_binders_peptide_wt_metrics_unique.csv").exists() else None
    n_strong = len(sb) if sb is not None else 167
    return [
        ("COSMIC\nmutations", cosmic),
        ("Protein\nsequences", n_fasta),
        ("Primary\nneoantigens", primary),
        ("9-mer\npairs", n_9mer_pairs),
        ("Unique\n9-mers", n_unique_9mers),
        ("Strong\nbinders", n_strong),
    ]


def load_counts_overview():
    """Load counts for full pipeline overview (5 stages; first block COSMIC mutations omitted)."""
    n_fasta = len(list((DIR_01 / "protein_sequences").glob("*.fasta"))) if (DIR_01 / "protein_sequences").exists() else 70
    primary = len(pd.read_csv(DIR_02 / "neoantigens_primary.csv")) if (DIR_02 / "neoantigens_primary.csv").exists() else 382
    df_9 = pd.read_csv(DIR_02 / "neoantigen_9mers.csv") if (DIR_02 / "neoantigen_9mers.csv").exists() else None
    n_9mer_pairs = len(df_9) if df_9 is not None else 3418
    n_unique_9mers = df_9["neoantigen_peptide"].nunique() if df_9 is not None else 1739
    sb = pd.read_csv(DIR_04 / "strong_binders_peptide_wt_metrics_unique.csv") if (DIR_04 / "strong_binders_peptide_wt_metrics_unique.csv").exists() else None
    n_strong = len(sb) if sb is not None else 167
    return [
        ("Protein\nsequences", n_fasta),
        ("Primary\nneoantigens", primary),
        ("9-mer\npairs", n_9mer_pairs),
        ("Unique\n9-mers", n_unique_9mers),
        ("Strong\nbinders", n_strong),
    ]


def plot_schematic(stages, out_path, copy_to_thesis=True, colors=None):
    """Draw a clean horizontal flow schematic."""
    n = len(stages)
    labels = [s[0] for s in stages]
    values = [s[1] for s in stages]

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 4)
    ax.set_aspect("equal")
    ax.axis("off")

    box_w, box_h = 1.4, 1.8
    if colors is None:
        colors = plt.cm.Blues(np.linspace(0.4, 0.9, n))
    else:
        colors = [colors[i % len(colors)] for i in range(n)]
    x_centers = np.linspace(1.2, 12.8, n)

    for i, (label, val) in enumerate(stages):
        x = x_centers[i] - box_w / 2
        y = 1.1
        rect = mpatches.FancyBboxPatch(
            (x, y), box_w, box_h,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            facecolor=colors[i],
            edgecolor="white",
            linewidth=1.5,
        )
        ax.add_patch(rect)
        ax.text(x_centers[i], y + box_h * 0.62, f"{val:,}", ha="center", va="center", fontsize=14, fontweight="bold", color="white")
        ax.text(x_centers[i], y + box_h * 0.25, label, ha="center", va="center", fontsize=8, color="white", linespacing=1.2)

        # Arrow to next
        if i < n - 1:
            ax.annotate(
                "",
                xy=(x_centers[i + 1] - box_w / 2 - 0.05, y + box_h / 2),
                xytext=(x + box_w + 0.05, y + box_h / 2),
                arrowprops=dict(arrowstyle="->", color="#555", lw=2),
            )

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved {out_path}")
    if copy_to_thesis:
        import shutil
        THESIS_IMG_DIR.mkdir(parents=True, exist_ok=True)
        dest = THESIS_IMG_DIR / out_path.name
        shutil.copy(out_path, dest)
        print(f"Copied to {dest}")


def main():
    stages = load_counts()
    plot_schematic(stages, OUT_DIR / "nsclc_pipeline_schematic.png", copy_to_thesis=False)

    # Pipeline overview with new colour palette (for thesis)
    overview_stages = load_counts_overview()
    plot_schematic(
        overview_stages,
        OUT_DIR / "nsclc_pipeline_overview.png",
        copy_to_thesis=True,
        colors=PIPELINE_COLORS,
    )


if __name__ == "__main__":
    main()
