# NSCLC Neoantigen Discovery Pipeline

A five-stage computational workflow for translating lung cancer mutation data into genomic insights.

[**View Repository**](https://github.com/SaraFarmahini/NeoAg-Discovery)

### Pipeline Diagram

```
     COSMIC mutations
            ↓
   Neoantigen generation
            ↓
 NetMHCpan binding prediction
            ↓
  Strong binder filtering
            ↓
  TCR interaction scoring
            ↓
  Candidate neoantigens
```

### Biological Motivation

Tumor-specific mutations can generate **neoantigens**—peptides recognized by T cells as non-self. Identifying these neoantigens is critical for personalized cancer immunotherapy and vaccine design. This pipeline addresses that need by turning NSCLC mutation data into ranked, MHC-presented peptide candidates with associated TCR-binding scores.

---

This repository hosts a modular pipeline designed to identify high-affinity neoantigens from Non-Small Cell Lung Cancer (NSCLC) variants. By integrating COSMIC mutation data, MHC binding predictions (NetMHCpan 4.2), and TCR–peptide scoring, the pipeline provides a streamlined path from raw genomic data to interpreted visualizations.

---

## Pipeline Overview

The project is structured as a sequential flow, where each stage builds upon the output of the last. While scripts can be run individually, the entire workflow is orchestrated via Snakemake for reproducibility and efficiency.

### The Logic

1. **Mutation Data:** Standardizing COSMIC variants and fetching wildtype protein sequences.
2. **Peptide Generation:** Translating mutations into full-length primary neoantigens and sliding-window 9-mers.
3. **MHC Binding:** Predicting affinity using NetMHCpan 4.2 to identify "strong binders" (EL_Rank ≤ 0.5).
4. **Refinement & Scoring:** Applying BLOSUM distance metrics and TCR-binding probability models (NetTCR-2).
5. **Interpretation:** Generating diagnostic plots and funnel reports to visualize the selection process.

---

## Key Results

The pipeline produces interpretable outputs at each stage:

| Output | Description |
|--------|-------------|
| **Neoantigen / EL_Rank distribution** | Funnel and binding-affinity distribution across alleles and genes (`report_pipeline_funnel.png`, `report_el_rank_distribution.png`) |
| **TCR interaction score distribution** | CDR3β–peptide binding score and logit/probability distributions (`cdr3b_binding_score_distribution.png`, `strong_binders_logit_vs_probability.png`) |
| **Top candidate neoantigens** | Ranked table of strong binders with best TCR scores (`strong_binders_best_cdr3b.csv`) |

All figures are written to `05_cancer_genomic_interpretation/visualizations/` after running the workflow.

---

## Workflow Execution

One command runs the full pipeline from mutation data to visualizations:

```bash
snakemake --cores 4
```

Snakemake runs stages in dependency order and skips steps whose outputs already exist. For a dry-run to preview the plan: `snakemake -n`.

---

## Skills Demonstrated

- **Cancer genomics analysis** — COSMIC mutation data, transcript-level variants, protein sequence handling
- **Neoantigen discovery pipelines** — End-to-end workflow from mutations to candidate peptides
- **Immunoinformatics tools** — NetMHCpan (MHC binding), BLOSUM metrics, NetTCR-2 comparison
- **PyTorch deep learning models** — TCR–peptide binding scoring with custom checkpoints
- **Reproducible workflows using Snakemake** — Dependency-aware automation and incremental runs

---

## Getting Started

### Prerequisites

- **Python 3.x:** pandas, requests, torch, matplotlib, networkx
- **NetMHCpan 4.2:** Must be installed locally ([NetMHCpan license](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.2)).
- **Data:** Due to licensing, COSMIC mutation files are not included. Please place your `Cosmic_ResistanceMutations_v101_GRCh38.tsv` in `01_nsclc_mutation_data/`.

### Execution

The easiest way to run the entire pipeline is through Snakemake. This ensures dependencies are met and prevents re-running completed steps.

```bash
# Run the full pipeline using 4 cores
snakemake --cores 4

# Perform a dry-run to see the execution plan
snakemake -n
```

---

## Project Architecture

```
NeoAg-Discovery/
├── Snakefile              # The "brain" of the project (workflow automation)
├── 01_nsclc_mutation_data/      # Raw COSMIC inputs and cached protein FASTA files
├── 02_generate_mutant_peptides/ # Scripts for 9-mer generation & Ensembl API calls
├── 03_predict_mhc_binding/      # NetMHCpan wrappers and EL_Rank filtering
├── 04_filter_strong_binders/    # BLOSUM metrics, TCR model checkpoints (PyTorch), and NetTCR-2
└── 05_cancer_genomic_interpretation/  # The "Results" hub: visualizers and the final report
```

### Key Workflow Rules

| Rule | Function | Output |
|------|----------|--------|
| `fetch_primary` | Maps COSMIC IDs to protein sequences | `neoantigens_primary.csv` |
| `predict_mhc` | Runs NetMHCpan on 9-mer candidates | `nsclc_netmhcpan_out.xls` |
| `score_tcr` | Ranks peptides by TCR binding probability | `strong_binders_best_cdr3b.csv` |
| `plot_report` | Generates the pipeline "funnel" visualization | `report_pipeline_funnel.png` |

---

## Visualizing the Workflow

We've included a utility to visualize the logic of the pipeline. To see how the data flows between rules, run:

```bash
python 05_cancer_genomic_interpretation/plot_workflow_dag.py
```

Alternatively, if you have Graphviz installed:

```bash
snakemake --rulegraph | dot -Tpng -o workflow_dag.png
```

---

## Portability & Research Use

This pipeline was built with modularity in mind. Each folder contains its own README and requirements where specific environment needs differ (e.g., Stage 4's deep learning dependencies).

**Notes for Researchers:**

- **COSMIC Licensing:** Ensure you have a valid academic or commercial license before using the data scripts in Stage 1.
- **Performance:** For large datasets, use the `step3_run_netmhcpan_chunked.py` script provided in Stage 3 to handle memory overhead.
