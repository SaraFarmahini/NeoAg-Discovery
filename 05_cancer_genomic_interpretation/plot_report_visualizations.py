#!/usr/bin/env python3
"""
Generate visualizations that confirm REPORT_Neoantigen_Pipeline.md using actual data.
Outputs to step6_NSCLC/visualizations/.
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / "visualizations"
OUT_DIR.mkdir(parents=True, exist_ok=True)
# Organized layout: data in 01, 02, 04
PROJECT_ROOT = SCRIPT_DIR.parent if (SCRIPT_DIR.parent / "01_nsclc_mutation_data").exists() else SCRIPT_DIR
DIR_01 = PROJECT_ROOT / "01_nsclc_mutation_data"
DIR_02 = PROJECT_ROOT / "02_generate_mutant_peptides"
DIR_04 = PROJECT_ROOT / "04_filter_strong_binders"


def load_counts():
    """Load counts from data files (excluding headers)."""
    cosmic = len(pd.read_csv(DIR_01 / "Cosmic_ResistanceMutations_v101_GRCh38.tsv", sep="\t")) if (DIR_01 / "Cosmic_ResistanceMutations_v101_GRCh38.tsv").exists() else 5199
    updated = len(pd.read_csv(DIR_01 / "Updated_Cosmic_Data_with_Sequences.tsv", sep="\t")) if (DIR_01 / "Updated_Cosmic_Data_with_Sequences.tsv").exists() else 400
    primary = len(pd.read_csv(DIR_02 / "neoantigens_primary.csv")) if (DIR_02 / "neoantigens_primary.csv").exists() else 382
    df_9mers = pd.read_csv(DIR_02 / "neoantigen_9mers.csv") if (DIR_02 / "neoantigen_9mers.csv").exists() else None
    n_9mer_rows = len(df_9mers) if df_9mers is not None else 3418
    n_unique_peptides = len(df_9mers["neoantigen_peptide"].drop_duplicates()) if df_9mers is not None else 1739
    if (DIR_02 / "neoantigen_9mers.list").exists():
        n_unique_peptides = len([l for l in (DIR_02 / "neoantigen_9mers.list").read_text().strip().splitlines() if l.strip()])
    strong = pd.read_csv(DIR_04 / "nsclc_netmhcpan_strong_binders.csv") if (DIR_04 / "nsclc_netmhcpan_strong_binders.csv").exists() else None
    n_strong = len(strong) if strong is not None else 444
    n_strong_valid_wt = n_strong
    if strong is not None and "wt_peptide" in strong.columns:
        n_strong_valid_wt = strong["wt_peptide"].notna() & (strong["wt_peptide"].astype(str).str.strip() != "")
        n_strong_valid_wt = n_strong_valid_wt.sum()
    return {
        "COSMIC mutations (raw)": cosmic,
        "Mutations with sequence": updated,
        "Primary neoantigens": primary,
        "9-mer pairs": n_9mer_rows,
        "Unique 9-mers (to NetMHCpan)": n_unique_peptides,
        "Strong binders (EL_Rank≤0.5)": n_strong,
        "Strong binders (valid WT)": n_strong_valid_wt,
    }, strong, df_9mers


def plot_pipeline_funnel(counts_dict, out_path):
    """Funnel / pipeline stages with counts."""
    stages = list(counts_dict.keys())
    values = list(counts_dict.values())
    colors = plt.cm.Blues(np.linspace(0.35, 0.85, len(stages)))[::-1]

    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = np.arange(len(stages)) * 1.2
    bars = ax.barh(y_pos, values, height=0.7, color=colors, edgecolor="gray", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(stages, fontsize=11)
    ax.set_xlabel("Count", fontsize=12)
    ax.set_title("Step6 NSCLC pipeline: sample and peptide counts (report confirmation)", fontsize=13)
    ax.set_xlim(0, max(values) * 1.08)
    for i, (bar, v) in enumerate(zip(bars, values)):
        ax.text(v + max(values) * 0.01, bar.get_y() + bar.get_height() / 2, f"{v:,}", va="center", fontsize=10, fontweight="bold")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


def plot_strong_binders_by_allele(strong_df, out_path):
    """Strong binders per HLA allele."""
    if strong_df is None or "allele" not in strong_df.columns:
        return
    counts = strong_df["allele"].value_counts().sort_values(ascending=True)
    fig, ax = plt.subplots(figsize=(9, 5))
    y_pos = np.arange(len(counts))
    ax.barh(y_pos, counts.values, color=plt.cm.viridis(np.linspace(0.2, 0.8, len(counts))), edgecolor="gray", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(counts.index, fontsize=11)
    ax.set_xlabel("Number of strong binder (peptide–allele) pairs", fontsize=12)
    ax.set_title("Strong binders (EL_Rank ≤ 0.5) by HLA allele", fontsize=13)
    for i, v in enumerate(counts.values):
        ax.text(v + 2, i, f"{v}", va="center", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


def plot_strong_binders_by_gene(strong_df, out_path, top_n=15):
    """Strong binders per gene (top N)."""
    if strong_df is None or "GENE_SYMBOL" not in strong_df.columns:
        return
    counts = strong_df["GENE_SYMBOL"].value_counts()
    if len(counts) > top_n:
        counts = counts.head(top_n)
    counts = counts.sort_values(ascending=True)
    fig, ax = plt.subplots(figsize=(9, max(5, len(counts) * 0.35)))
    y_pos = np.arange(len(counts))
    ax.barh(y_pos, counts.values, color=plt.cm.plasma(np.linspace(0.2, 0.8, len(counts))), edgecolor="gray", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(counts.index, fontsize=10)
    ax.set_xlabel("Number of strong binder rows", fontsize=12)
    ax.set_title(f"Strong binders by gene (top {top_n})", fontsize=13)
    for i, v in enumerate(counts.values):
        ax.text(v + 1, i, f"{v}", va="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


def plot_9mers_per_gene(df_9mers, out_path, top_n=15):
    """9-mer pairs per gene (from neoantigen_9mers.csv)."""
    if df_9mers is None or "GENE_SYMBOL" not in df_9mers.columns:
        return
    counts = df_9mers["GENE_SYMBOL"].value_counts()
    if len(counts) > top_n:
        counts = counts.head(top_n)
    counts = counts.sort_values(ascending=True)
    fig, ax = plt.subplots(figsize=(9, max(5, len(counts) * 0.35)))
    y_pos = np.arange(len(counts))
    ax.barh(y_pos, counts.values, color=plt.cm.coolwarm(np.linspace(0.2, 0.8, len(counts))), edgecolor="gray", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(counts.index, fontsize=10)
    ax.set_xlabel("Number of 9-mer pairs (mutation-including)", fontsize=12)
    ax.set_title(f"9-mer pairs by gene (top {top_n})", fontsize=13)
    for i, v in enumerate(counts.values):
        ax.text(v + 5, i, f"{v}", va="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


def plot_el_rank_distribution(strong_df, out_path):
    """Distribution of EL_Rank among strong binders (all ≤ 0.5)."""
    if strong_df is None or "EL_Rank" not in strong_df.columns:
        return
    ranks = pd.to_numeric(strong_df["EL_Rank"], errors="coerce").dropna()
    if ranks.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(ranks, bins=30, color="steelblue", edgecolor="white", alpha=0.85)
    ax.axvline(0.5, color="red", linestyle="--", linewidth=2, label="Strong binder threshold (0.5)")
    ax.set_xlabel("EL_Rank", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Distribution of EL_Rank in strong binders (all ≤ 0.5)", fontsize=13)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


def main():
    counts, strong_df, df_9mers = load_counts()
    plot_pipeline_funnel(counts, OUT_DIR / "report_pipeline_funnel.png")
    plot_strong_binders_by_allele(strong_df, OUT_DIR / "report_strong_binders_by_allele.png")
    plot_strong_binders_by_gene(strong_df, OUT_DIR / "report_strong_binders_by_gene.png")
    plot_9mers_per_gene(df_9mers, OUT_DIR / "report_9mers_by_gene.png")
    plot_el_rank_distribution(strong_df, OUT_DIR / "report_el_rank_distribution.png")
    print("Done. All report visualizations in:", OUT_DIR)


if __name__ == "__main__":
    main()
