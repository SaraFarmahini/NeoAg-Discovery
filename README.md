# NSCLC Neoantigen Pipeline

**Repository:** [github.com/SaraFarmahini/NeoAg-Discovery](https://github.com/SaraFarmahini/NeoAg-Discovery)

This project is organized as a **five-stage pipeline** from NSCLC mutation data to cancer genomic interpretation. Run stages in order; each stage uses outputs from the previous one.

---

## How the project is organized

### Pipeline flow

```
01  NSCLC mutation data
    ↓
02  Generate mutant peptides
    ↓
03  Predict MHC binding affinity
    ↓
04  Filter strong binders
    ↓
05  Cancer genomic interpretation
```

### Folder structure (tree)

```
NeoAg-Discovery/
├── README.md
├── Snakefile                        ← Snakemake workflow (run: snakemake --cores 4)
│
├── 01_nsclc_mutation_data/          ← NSCLC mutation data
│   ├── Cosmic_ResistanceMutations_v101_GRCh38.tsv
│   ├── Updated_Cosmic_Data_with_Sequences.tsv
│   └── protein_sequences/           (cached .fasta per transcript)
│
├── 02_generate_mutant_peptides/     ← Generate mutant peptides
│   ├── step1_fetch_proteins_and_primary_neoantigens.py
│   ├── step2_generate_neoantigen_9mers.py
│   ├── requirements_step1.txt
│   ├── README_step1.md
│   ├── neoantigens_primary.csv
│   ├── neoantigen_9mers.csv
│   └── neoantigen_9mers.list
│
├── 03_predict_mhc_binding/          ← Predict MHC binding affinity
│   ├── step3_run_netmhcpan.py
│   ├── step3_run_netmhcpan_chunked.py
│   ├── run_netmhcpan_cluster.slurm
│   ├── check_netmhcpan_progress.sh
│   ├── create_threshold_files.py
│   ├── README_step3_netmhcpan.md
│   ├── README_cluster.md
│   ├── nsclc_netmhcpan_out.xls
│   └── nsclc_netmhcpan_strong_binders.csv
│
├── 04_filter_strong_binders/        ← Filter strong binders
│   ├── nsclc_netmhcpan_strong_binders.csv   (from step 3)
│   ├── strong_binders_peptide_wt.csv
│   ├── strong_binders_peptide_wt_metrics.csv
│   ├── strong_binders_peptide_wt_metrics_unique.csv
│   ├── strong_binders_best_cdr3b.csv
│   ├── compute_blosum_distance.py
│   ├── score_peptide_cdr3.py
│   ├── run_score_peptide_cdr3.slurm
│   ├── check_promiscuous_cdr3.py
│   ├── ds_cdr3b_unique.csv
│   ├── tuned_model/                 (TCR scoring: final_model_new_hpo_bce_logits.py, finetune/best_model.pt)
│   ├── result1/                     (filter variants, extra outputs)
│   ├── README_strong_binders_duplicates.md
│   └── NetTCR-2: run_nettcr2_predict.py, score_nettcr2_test_pairs.py, compare_nettcr2_model_concordance.py, ...
│
└── 05_cancer_genomic_interpretation/  ← Cancer genomic interpretation
    ├── REPORT_Neoantigen_Pipeline.md
    ├── plot_report_visualizations.py
    ├── plot_pipeline_schematic.py
    ├── plot_workflow_dag.py            ← workflow DAG (Snakemake rules)
    ├── plot_cdr3b_score_distribution.py
    ├── plot_logit_probability_distributions.py
    ├── promiscuity_check_report1.txt
    └── visualizations/              (generated .png: report_*.png, nsclc_pipeline_*.png, cdr3b_binding_score_distribution.png, strong_binders_logit_*.png, workflow_dag.png)
```

At project root you also get **`workflow_rulegraph.dot`** after running `plot_workflow_dag.py` (for Graphviz).

### Quick reference

| Stage | Folder | What it does |
|-------|--------|----------------|
| **1** | `01_nsclc_mutation_data/` | COSMIC mutations + protein sequences |
| **2** | `02_generate_mutant_peptides/` | Primary neoantigens → 9-mer peptides |
| **3** | `03_predict_mhc_binding/` | NetMHCpan → strong binders (EL_Rank ≤ 0.5) |
| **4** | `04_filter_strong_binders/` | Peptide–WT metrics, TCR scoring, NetTCR-2 |
| **5** | `05_cancer_genomic_interpretation/` | Report + plots and visualizations |

Scripts detect this layout and read/write paths relative to the project root.

---

## Running with Snakemake

The pipeline can be run as a single workflow so Snakemake runs steps in the correct order, skips already-completed steps, and tracks dependencies.

**Run the full pipeline (from project root):**

```bash
snakemake --cores 4
```

**What Snakemake does:**

- Runs steps in dependency order: mutation data → peptides → MHC prediction → strong binders → visualizations.
- Avoids recomputing: if an output file is already present and inputs unchanged, that step is skipped.
- Tracks dependencies: changing an input (e.g. COSMIC TSV) triggers only the downstream rules that need it.

**Rule summary:**

| Rule | Input(s) | Output(s) |
|------|----------|-----------|
| `fetch_primary` | COSMIC TSV | `neoantigens_primary.csv` |
| `generate_9mers` | primary CSV | `neoantigen_9mers.csv`, `neoantigen_9mers.list` |
| `predict_mhc` | 9-mer list + CSV | NetMHCpan `.xls`, `nsclc_netmhcpan_strong_binders.csv` |
| `copy_strong_binders` | strong binders (03) | strong binders copy in 04 |
| `peptide_wt` | strong binders | `strong_binders_peptide_wt.csv` |
| `blosum_metrics` | peptide_wt | `strong_binders_peptide_wt_metrics.csv` |
| `metrics_unique` | metrics | `strong_binders_peptide_wt_metrics_unique.csv` |
| `score_tcr` | metrics + CDR3 list | `strong_binders_best_cdr3b.csv` |
| `plot_report` | COSMIC, primary, 9mers, strong | `report_pipeline_funnel.png`, `report_strong_binders_by_allele.png`, `report_strong_binders_by_gene.png`, `report_9mers_by_gene.png`, `report_el_rank_distribution.png` |
| `plot_schematic` | primary, 9mers, metrics_unique | `nsclc_pipeline_schematic.png`, `nsclc_pipeline_overview.png` |
| `plot_cdr3b` | best_cdr3b | `cdr3b_binding_score_distribution.png` |
| `plot_logit` | best_cdr3b | `strong_binders_logit_probability_distributions.png`, `strong_binders_logit_vs_probability.png` |

**Prerequisites:** Python (pandas, requests, torch), NetMHCpan 4.2 for `predict_mhc`, and optionally the `sequence_encodings` module for `blosum_metrics` (see `04_filter_strong_binders/compute_blosum_distance.py`). TCR scoring uses `04_filter_strong_binders/tuned_model/final_model_new_hpo_bce_logits.py` and `tuned_model/finetune/best_model.pt`.

**Run only up to a given stage (example: stop after strong binders CSV):**

```bash
snakemake --cores 4 04_filter_strong_binders/nsclc_netmhcpan_strong_binders.csv
```

**Dry run (show what would be run):**

```bash
snakemake -n
```

---

## Visualize the workflow

You can generate a **rule DAG** of the Snakemake pipeline in two ways.

### 1. Python script (matplotlib + networkx)

From the project root:

```bash
python 05_cancer_genomic_interpretation/plot_workflow_dag.py
```

This creates:

- **`05_cancer_genomic_interpretation/visualizations/workflow_dag.png`** — diagram of all rules and dependencies.
- **`workflow_rulegraph.dot`** — Graphviz source (see below).

Requires: `matplotlib`, `networkx` (`pip install matplotlib networkx`).

### 2. Snakemake + Graphviz (optional)

If Graphviz is installed (`brew install graphviz` or `conda install graphviz`):

```bash
snakemake --rulegraph | dot -Tpng -o workflow_dag.png
```

Or use the `.dot` file produced by the script above:

```bash
dot -Tpng workflow_rulegraph.dot -o workflow_dag.png
```

The diagram shows rules (e.g. `fetch_primary`, `generate_9mers`, `predict_mhc`) and arrows for dependencies (e.g. `fetch_primary` → `generate_9mers` → `predict_mhc` → … → `all`).

---

## Stage 1: NSCLC mutation data  
**Folder:** `01_nsclc_mutation_data/`

- **Input:** COSMIC resistance mutations (`Cosmic_ResistanceMutations_v101_GRCh38.tsv`)
- **Outputs:** `Updated_Cosmic_Data_with_Sequences.tsv`, `protein_sequences/*.fasta`
- No script to run here; place or generate the COSMIC TSV and (optionally) pre-fetched protein sequences in this folder.

---

## Stage 2: Generate mutant peptides  
**Folder:** `02_generate_mutant_peptides/`

- **Step 1:** Fetch wildtype proteins and build primary neoantigens (full-length).
- **Step 2:** Generate mutation-including 9-mers for MHC prediction.

**Run (from project root or from inside `02_generate_mutant_peptides/`):**

```bash
cd 02_generate_mutant_peptides
pip install pandas requests   # if needed
python step1_fetch_proteins_and_primary_neoantigens.py
python step2_generate_neoantigen_9mers.py
```

**Outputs:** `neoantigens_primary.csv`, `neoantigen_9mers.csv`, `neoantigen_9mers.list`

---

## Stage 3: Predict MHC binding affinity  
**Folder:** `03_predict_mhc_binding/`

- Runs **NetMHCpan** on the 9-mer list and filters to strong binders (EL_Rank ≤ 0.5).

**Run:**

```bash
cd 03_predict_mhc_binding
python step3_run_netmhcpan.py
# Optional: --limit 100 for testing; --netmhcpan_dir /path/to/netMHCpan-4.2
```

**Outputs:** `nsclc_netmhcpan_out.xls`, `nsclc_netmhcpan_strong_binders.csv` (and a copy in `04_filter_strong_binders/` for the next stage)

---

## Stage 4: Filter strong binders  
**Folder:** `04_filter_strong_binders/`

- Build peptide–WT tables, distance metrics (BLOSUM, Boman), and TCR–peptide scoring (best CDR3β per peptide).
- Optional: NetTCR-2 comparison scripts.

**Typical order:**

1. Build peptide–WT and metrics from strong binders (e.g. derive `strong_binders_peptide_wt.csv` from `nsclc_netmhcpan_strong_binders.csv` if needed, then run `compute_blosum_distance.py`).
2. Run TCR scoring: `score_peptide_cdr3.py` (uses `tuned_model/` and `ds_cdr3b_unique.csv`).
3. Optional: `check_promiscuous_cdr3.py`, NetTCR-2 comparison.

**Outputs:** `strong_binders_peptide_wt_metrics.csv`, `strong_binders_best_cdr3b.csv`, etc.

---

## Stage 5: Cancer genomic interpretation  
**Folder:** `05_cancer_genomic_interpretation/`

- Report, funnel plots, and visualizations that confirm pipeline counts and distributions.

**Run:**

```bash
cd 05_cancer_genomic_interpretation
python plot_report_visualizations.py
python plot_pipeline_schematic.py
python plot_cdr3b_score_distribution.py
python plot_logit_probability_distributions.py
```

**Outputs:** `visualizations/*.png` (and optionally thesis image folder if configured)

---

## Using this project on GitHub (e.g. as NGS / bioinformatics experience)

This repo is suitable to publish on GitHub as a portfolio project. It demonstrates:

- **NGS / cancer genomics:** COSMIC mutation data, transcript-level variants, protein sequences
- **Computational pipelines:** Snakemake, dependency-aware workflows, reproducibility
- **Immunoinformatics:** Neoantigen peptide generation, MHC binding (NetMHCpan), TCR–peptide scoring
- **Python:** Data handling (pandas), REST APIs (Ensembl), optional ML (PyTorch for TCR model)

**Before pushing:**

1. **Data and licensing**
   - **COSMIC:** Mutation data is from [COSMIC](https://cancer.sanger.ac.uk/cosmic); their [license](https://cancer.sanger.ac.uk/cosmic/license) may restrict redistribution. To keep the repo code-only, add the COSMIC TSV to `.gitignore` (see commented lines) and document in the README: *“COSMIC data not included; download from COSMIC after registration and place in `01_nsclc_mutation_data/`.”*
   - **NetMHCpan:** Requires a separate [license](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.2); do not commit the NetMHCpan binary or data. Document the install path in the README or config.
   - **TCR model:** If `tuned_model/` checkpoints (e.g. `best_model.pt`) are large or private, add `*.pt` and `tuned_model/finetune/` to `.gitignore` and add a short note on how to obtain or retrain the model.

2. **What to commit**
   - All Python and Snakemake code, README, and small config files.
   - Optionally: small example outputs (e.g. a few rows of CSVs or one schematic PNG) so viewers see the pipeline runs; keep the repo under GitHub’s recommended size by ignoring large CSVs, `.xls`, and `protein_sequences/` if needed (see `.gitignore`).

3. **Repository description**
   - Example: *“NSCLC neoantigen discovery pipeline: COSMIC mutations → peptide generation → NetMHCpan → TCR scoring (Snakemake).”*
   - Add topics such as: `bioinformatics`, `ngs`, `neoantigen`, `snakemake`, `immunoinformatics`, `cancer-genomics`.

**Clone / push:**  
`git clone https://github.com/SaraFarmahini/NeoAg-Discovery.git`  
To push from an existing local copy: `git remote add origin https://github.com/SaraFarmahini/NeoAg-Discovery.git`, then `git push -u origin main`.

---

## Folder layout (summary)

| Folder | Purpose |
|--------|--------|
| `01_nsclc_mutation_data/` | COSMIC TSV, Updated Cosmic table, protein_sequences |
| `02_generate_mutant_peptides/` | Step 1 & 2 scripts, primary neoantigens, 9-mer list/CSV |
| `03_predict_mhc_binding/` | NetMHCpan scripts, raw output, strong binders CSV |
| `04_filter_strong_binders/` | Strong-binder tables, metrics, TCR scoring, NetTCR-2 |
| `05_cancer_genomic_interpretation/` | Report, plot scripts, visualizations |
