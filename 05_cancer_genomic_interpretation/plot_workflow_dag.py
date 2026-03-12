#!/usr/bin/env python3
"""
Generate a visualization of the Snakemake workflow (rule DAG).

Outputs:
  - visualizations/workflow_dag.png  (matplotlib + networkx)
  - workflow_rulegraph.dot           (Graphviz format; use: dot -Tpng workflow_rulegraph.dot -o workflow_dag.png)

Run from project root:
  python 05_cancer_genomic_interpretation/plot_workflow_dag.py
"""

from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / "visualizations"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DOT_FILE = SCRIPT_DIR.parent / "workflow_rulegraph.dot"

# Rule graph: edges (source -> target) in dependency order
EDGES = [
    ("fetch_primary", "generate_9mers"),
    ("generate_9mers", "predict_mhc"),
    ("predict_mhc", "peptide_wt"),
    ("predict_mhc", "copy_strong_binders"),
    ("peptide_wt", "blosum_metrics"),
    ("blosum_metrics", "metrics_unique"),
    ("blosum_metrics", "score_tcr"),
    ("copy_strong_binders", "plot_report"),
    ("generate_9mers", "plot_report"),
    ("fetch_primary", "plot_report"),
    ("metrics_unique", "plot_schematic"),
    ("generate_9mers", "plot_schematic"),
    ("fetch_primary", "plot_schematic"),
    ("score_tcr", "plot_cdr3b"),
    ("score_tcr", "plot_logit"),
    ("plot_report", "all"),
    ("plot_schematic", "all"),
    ("plot_cdr3b", "all"),
    ("plot_logit", "all"),
    ("score_tcr", "all"),
]


# Shorter display labels so text fits inside nodes
LABELS = {
    "fetch_primary": "fetch_primary",
    "generate_9mers": "generate_9mers",
    "predict_mhc": "predict_mhc",
    "copy_strong_binders": "copy_strong",
    "peptide_wt": "peptide_wt",
    "blosum_metrics": "blosum",
    "metrics_unique": "metrics_uni",
    "score_tcr": "score_tcr",
    "plot_report": "plot_report",
    "plot_schematic": "schematic",
    "plot_cdr3b": "plot_cdr3b",
    "plot_logit": "plot_logit",
    "all": "all",
}

# Pipeline stages for colouring and legend (rule -> stage index 0..4)
STAGE_NAMES = [
    "1. Mutation data",
    "2. Generate peptides",
    "3. MHC prediction",
    "4. Filter strong binders",
    "5. Interpretation",
]
RULE_TO_STAGE = {
    "fetch_primary": 0,
    "generate_9mers": 1,
    "predict_mhc": 2,
    "copy_strong_binders": 2,
    "peptide_wt": 3,
    "blosum_metrics": 3,
    "metrics_unique": 3,
    "score_tcr": 3,
    "plot_report": 4,
    "plot_schematic": 4,
    "plot_cdr3b": 4,
    "plot_logit": 4,
    "all": 4,
}
# Colours per stage (distinct, print-friendly)
STAGE_COLORS = [
    "#2166ac",   # 1. Mutation data – dark blue
    "#1b9e77",   # 2. Generate peptides – teal
    "#d95f02",   # 3. MHC prediction – orange
    "#7570b3",   # 4. Filter strong binders – purple
    "#e7298a",   # 5. Interpretation – magenta
]


def main():
    try:
        import matplotlib.pyplot as plt
        import networkx as nx
    except ImportError:
        print("Install matplotlib and networkx: pip install matplotlib networkx")
        return 1

    G = nx.DiGraph()
    G.add_edges_from(EDGES)

    pos = nx.spring_layout(G, seed=42, k=2.5, iterations=50)

    labels = {n: LABELS.get(n, n) for n in G.nodes()}
    node_colors = [STAGE_COLORS[RULE_TO_STAGE.get(n, 4)] for n in G.nodes()]
    node_edge_colors = [STAGE_COLORS[RULE_TO_STAGE.get(n, 4)] for n in G.nodes()]
    # Darken edge for visibility
    import matplotlib.colors as mcolors
    node_edges = [mcolors.to_rgba(c, alpha=0.9) for c in node_edge_colors]
    for i, c in enumerate(node_edges):
        r, g, b, _ = c
        node_edges[i] = (min(1, r * 0.6), min(1, g * 0.6), min(1, b * 0.6), 1)

    fig, ax = plt.subplots(figsize=(15, 10))
    nx.draw_networkx_nodes(
        G, pos, node_color=node_colors, edgecolors=node_edges, linewidths=2,
        node_size=5200, ax=ax,
    )
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight="bold", ax=ax, font_color="white")
    nx.draw_networkx_edges(
        G, pos, edge_color="gray", arrows=True, arrowsize=22,
        connectionstyle="arc3,rad=0.1", ax=ax, width=1.5,
    )
    ax.set_title("NSCLC Neoantigen Pipeline — Snakemake workflow by stage", fontsize=14)

    # Legend: one patch per stage
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=STAGE_COLORS[i], edgecolor="black", linewidth=1, label=STAGE_NAMES[i])
        for i in range(len(STAGE_NAMES))
    ]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=9, framealpha=0.95)
    ax.axis("off")
    plt.tight_layout()
    out_png = OUT_DIR / "workflow_dag.png"
    plt.savefig(out_png, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved {out_png}")

    lines = ["digraph workflow {", "    rankdir=TB;", "    node [shape=box, style=rounded];"]
    for u, v in EDGES:
        lines.append(f'    "{u}" -> "{v}";')
    lines.append("}")
    DOT_FILE.write_text("\n".join(lines))
    print(f"Saved {DOT_FILE} (use: dot -Tpng {DOT_FILE.name} -o workflow_dag.png)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
