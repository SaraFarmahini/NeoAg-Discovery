# Step6 NSCLC: Complete Pipeline Report — From COSMIC Mutations to Strong-Binder Neoantigens and TCR Matching

This report documents the full workflow used in `step6_NSCLC` to derive **neoantigen peptides** and their **wildtype counterparts** from COSMIC mutation data, extract **9-mers** spanning the mutation, predict **MHC strong binders** with NetMHCpan, and prepare them for **TCR–peptide scoring**.

---

## 1. Overview

| Stage | Input | Output | Purpose |
|--------|--------|--------|---------|
| **1** | COSMIC TSV | Primary neoantigens (full-length) | Get wildtype and mutant protein sequences per mutation |
| **2** | Primary neoantigens | 9-mer peptides (mutation-including) | Sliding 9-mers that contain the mutation and differ from WT |
| **3** | 9-mer list | NetMHCpan predictions → strong binders | Filter to MHC-I strong binders (EL_Rank ≤ 0.5) |
| **4** | Strong binders | Peptide–WT tables, distance metrics, visualizations | BLOSUM/Boman/aliphatic distances; QC and figures |
| **5** | Strong binders + CDR3β list | Best-matching TCR per peptide | Finetuned model scores peptide–CDR3 pairs |

Counts achieved in this run are summarised in **Section 7**.

---

## 2. Data Source: COSMIC

**Input file:** `Cosmic_ResistanceMutations_v101_GRCh38.tsv`

- **Source:** COSMIC (Catalogue Of Somatic Mutations In Cancer), resistance mutations dataset (GRCh38), version 101.
- **Content:** Somatic mutations associated with drug resistance (e.g. NSCLC, EGFR, KIT, etc.), with sample and gene identifiers, transcript accessions, and protein-level mutation annotations.
- **Key columns used:**
  - `TRANSCRIPT_ACCESSION` — Ensembl transcript ID (e.g. ENST00000275493.6)
  - `HGVSP` — protein-level mutation (e.g. ENSP00000275493.2:p.Thr790Met)
  - `MUTATION_ID`, `GENE_SYMBOL`, `MUTATION_AA`, `MUTATION_CDS`, and optionally `MUTATION_SOMATIC_STATUS`
- **Raw row count (this run):** 5,199 mutation rows (5,200 lines including header).

These mutations are the starting point for defining **wildtype** and **neoantigen** protein sequences.

---

## 3. Step 1: From COSMIC Mutations to Primary Neoantigens (Full-Length)

**Script:** `step1_fetch_proteins_and_primary_neoantigens.py`

### 3.1 Procedure

1. **Load and clean the COSMIC TSV**
   - Drop duplicates on `(GENE_SYMBOL, MUTATION_AA, MUTATION_CDS)` (keep first).
   - Optionally remove rows with `?` in `MUTATION_AA` or `MUTATION_CDS`.
   - Optionally keep only rows with `MUTATION_SOMATIC_STATUS == 'Confirmed somatic variant'`.
2. **Parse HGVSP** (e.g. `p.Thr790Met`) to obtain:
   - Wildtype amino acid (1-letter), mutant amino acid (1-letter), and 1-based position.
3. **Fetch wildtype protein sequences** for each unique `TRANSCRIPT_ACCESSION`:
   - From Ensembl REST API: `GET /sequence/id/{transcript_id}?type=protein&format=fasta`.
   - Cached under `protein_sequences/<TRANSCRIPT_BASE>.fasta`.
4. **Build mutant sequence** by applying the single amino-acid substitution at the parsed position to the wildtype sequence (validation: wildtype residue must match at that position).
5. **Outputs:**
   - **Updated_Cosmic_Data_with_Sequences.tsv** — mutation table plus a `protein_sequence` column (wildtype) for each row with a successfully fetched sequence.
   - **neoantigens_primary.csv** — one row per mutation with:
     - Full-length **wildtype_protein_sequence** and **neoantigen_protein_sequence** (mutant),
     - `mutation_pos`, `wt_aa`, `mut_aa`, `GENE_SYMBOL`, `TRANSCRIPT_ACCESSION`, `HGVSP`, `MUTATION_ID`, etc.

No 9-mer slicing is done here; only full-length wildtype and neoantigen proteins are defined.

### 3.2 Counts (this run)

- **Updated_Cosmic_Data_with_Sequences.tsv:** 400 rows (mutations with a valid protein sequence).
- **neoantigens_primary.csv:** 382 primary neoantigen rows (full-length wildtype + mutant pairs).

---

## 4. Step 2: Primary Neoantigens → 9-mers That Include the Mutation

**Script:** `step2_generate_neoantigen_9mers.py`

### 4.1 Rationale

- MHC class I typically presents **8–11 mer** peptides; we use **9-mers**.
- Only 9-mers that **contain the mutated position** and **differ from the wildtype** at that window are true neoantigens; others would be identical to wildtype.

### 4.2 Procedure

1. **Input:** `neoantigens_primary.csv` (columns: `wildtype_protein_sequence`, `neoantigen_protein_sequence`, `mutation_pos`, plus gene/transcript/mutation metadata).
2. **For each primary row:**  
   Sliding 9-mer windows are generated such that:
   - The window contains the mutation position (1-based `mutation_pos`).
   - The 9-mer from the **mutant** sequence is different from the 9-mer from the **wildtype** sequence at the same start position (so the mutation is reflected in the peptide).
3. **Outputs:**
   - **neoantigen_9mers.csv** — one row per (neoantigen_peptide, wt_peptide) 9-mer pair, with:
     - `neoantigen_peptide`, `wt_peptide`, `peptide_length` (9), `mutation_pos`, `peptide_start_pos`, `wt_aa`, `mut_aa`, and mutation identifiers.
   - **neoantigen_9mers.list** — a plain list of **unique** neoantigen 9-mer sequences (one per line), used as the peptide input for NetMHCpan.

### 4.3 Counts (this run)

- **neoantigen_9mers.csv:** 3,418 rows (each row is a distinct 9-mer pair; the same peptide can appear in multiple overlapping windows for a given mutation).
- **neoantigen_9mers.list:** 1,739 **unique** neoantigen 9-mer peptides (deduplicated) sent to NetMHCpan.

---

## 5. Step 3: NetMHCpan and Strong Binder Selection

**Script:** `step3_run_netmhcpan.py`  
**Cluster job (optional):** `run_netmhcpan_cluster.slurm`

### 5.1 NetMHCpan Run

- **Input peptides:** `neoantigen_9mers.list` (1,739 unique 9-mers).
- **Tool:** NetMHCpan 4.2.
- **MHC alleles (default):** 9 alleles — HLA-A01:01, HLA-A02:01, HLA-A03:01, HLA-A24:02, HLA-B07:02, HLA-B08:01, HLA-B15:01, HLA-C07:01, HLA-C07:02.
- **Output:** **nsclc_netmhcpan_out.xls** — predictions for each (peptide × allele) pair (Score, Rank, BA_score, BA_Rank, etc.). The XLS has 2 header lines and 1,739 data rows (one per peptide), with blocks of columns per allele.

### 5.2 Strong Binder Definition and Parsing

- **Criterion:** **EL_Rank ≤ 0.5** (percentile rank of the predicted binding, strong binder threshold).
- The script parses the XLS (handling column names **Rank** / **Score** as well as **EL_Rank** / **EL-score**), builds a long-format table of (peptide, allele, EL_Rank, …), and keeps rows with `EL_Rank ≤ 0.5`.
- These strong-binder rows are **merged** with `neoantigen_9mers.csv` on peptide (as `neoantigen_peptide`) to attach wildtype peptide and mutation metadata. Missing `wt_peptide` is filled from `neoantigen_9mers.csv` where available.
- **Output:** **nsclc_netmhcpan_strong_binders.csv** — strong binders with columns: peptide, allele, EL_Rank, EL_score, BA_Rank, BA_score, neoantigen_peptide, wt_peptide, peptide_length, mutation_pos, peptide_start_pos, wt_aa, mut_aa, GENE_SYMBOL, TRANSCRIPT_ACCESSION, MUTATION_ID, HGVSP, MUTATION_AA.

### 5.3 Counts (this run)

- **Peptides sent to NetMHCpan:** 1,739 unique 9-mers.
- **NetMHCpan output:** 1,739 peptides × 9 alleles = 15,651 (peptide, allele) predictions in the XLS.
- **Strong binders (EL_Rank ≤ 0.5):** 449 rows (some peptides are strong binders for more than one allele).
- After **removing rows with missing wildtype** (e.g. peptides not in `neoantigen_9mers.csv` or with invalid residues): **444 strong-binder rows** retained for downstream analysis.

So: **1,739 peptides** were tested; **444 peptide–allele strong-binder rows** (with valid peptide and wildtype) were kept for the rest of the pipeline.

---

## 6. Downstream Steps (After Strong Binders)

### 6.1 Peptide–Wildtype Tables

- **strong_binders_peptide_wt.csv** — two columns: `peptide`, `wt` (wildtype 9-mer), one row per strong-binder row (444 rows after removing missing WT).
- **strong_binders_peptide_wt_metrics.csv** — same rows plus **BLOSUM62 distance**, **aliphatic index distance**, and **Boman distance** between peptide and WT (computed by `compute_blosum_distance.py`).

### 6.2 Visualizations

- **Paired distance plots** (same style as step5): ECDF and scatter of the three distance metrics across strong-binder peptide–WT pairs.
- Generated with the step5 script:  
  `plot_distance_metrics.py --input strong_binders_peptide_wt_metrics.csv --out-dir step6_NSCLC/visualizations`  
  Outputs: `paired_distance_ecdf.png`, `paired_distance_plot.png`, `paired_distance_ranked.png`.

### 6.3 TCR–Peptide Scoring (Best-Matching CDR3β)

- **Model:** Architecture from `final_model_new_hpo_bce_logits.py`; weights from the **blosum_boman** finetuned checkpoint (`finetune_outputs_ablation/ablation_blosum_boman/best_model.pt`).
- **Inputs:**  
  - Peptides and deltas: `strong_binders_peptide_wt_metrics.csv` (peptide, wt, blosum_distance, boman_distance).  
  - TCR repertoire: `ds_cdr3b_unique.csv` (197,851 unique CDR3β sequences).
- **Script:** `score_peptide_cdr3.py` — for each peptide row, scores all (peptide, CDR3β) pairs in batches, and records the **best_cdr3b** (highest binding score) and **binding_score**.
- **Output:** e.g. **strong_binders_best_cdr3b.csv** (all strong-binder columns plus `best_cdr3b`, `binding_score`).
- **Quality check:** `check_promiscuous_cdr3.py` — flags if the same CDR3β is best for every (or too many) peptides, as a simple check for promiscuous binders.
- **Cluster:** `run_score_peptide_cdr3.slurm` runs the full scoring on GPU and then the promiscuity check.

---

## 7. Summary: Sample and Peptide Counts

| Step | Description | Count |
|------|-------------|--------|
| **COSMIC input** | Mutation rows in `Cosmic_ResistanceMutations_v101_GRCh38.tsv` | 5,199 |
| **After cleaning & sequence fetch** | Rows in `Updated_Cosmic_Data_with_Sequences.tsv` | 400 |
| **Primary neoantigens** | Rows in `neoantigens_primary.csv` (full-length WT + mutant) | 382 |
| **9-mer pairs** | Rows in `neoantigen_9mers.csv` (mutation-including 9-mers) | 3,418 |
| **Unique 9-mers to NetMHCpan** | Lines in `neoantigen_9mers.list` | 1,739 |
| **NetMHCpan predictions** | Peptides × alleles in `nsclc_netmhcpan_out.xls` | 1,739 × 9 |
| **Strong binders (EL_Rank ≤ 0.5)** | Rows before removing missing WT | 449 |
| **Strong binders (with valid WT)** | Rows used in peptide–WT and TCR steps | **444** |
| **Unique CDR3β for TCR scoring** | Rows in `ds_cdr3b_unique.csv` | 197,851 |

**Summary sentence:** From **5,199 COSMIC resistance mutations**, we obtained **400 mutations** with valid protein sequences and **382 primary neoantigens**. These yielded **3,418 mutation-including 9-mer pairs** and **1,739 unique 9-mer peptides** run through NetMHCpan across **9 HLA alleles**. **444 strong-binder (peptide–allele) rows** with valid wildtype were retained and used for distance metrics, visualizations, and for finding the best-matching CDR3β per peptide from **197,851 unique TCRs**.

---

## 8. Report visualizations (data-confirmed figures)

The following figures are generated from the actual data in this directory by `plot_report_visualizations.py` and confirm the counts and distributions described in this report.

| Figure | Description |
|--------|--------------|
| **visualizations/report_pipeline_funnel.png** | Pipeline funnel: counts at each stage (COSMIC raw → mutations with sequence → primary neoantigens → 9-mer pairs → unique 9-mers → strong binders → strong binders with valid WT). |
| **visualizations/report_strong_binders_by_allele.png** | Number of strong-binder (peptide–allele) rows per HLA allele. |
| **visualizations/report_strong_binders_by_gene.png** | Number of strong-binder rows per gene (top 15 genes). |
| **visualizations/report_9mers_by_gene.png** | Number of 9-mer pairs per gene (top 15), from `neoantigen_9mers.csv`. |
| **visualizations/report_el_rank_distribution.png** | Histogram of EL_Rank in strong binders (all ≤ 0.5), with threshold line. |

To regenerate: `python3 plot_report_visualizations.py` (from `step6_NSCLC`).

---

## 9. File Reference (Key Outputs)

| File | Description |
|------|-------------|
| `Cosmic_ResistanceMutations_v101_GRCh38.tsv` | Input COSMIC mutations |
| `Updated_Cosmic_Data_with_Sequences.tsv` | Mutations + wildtype protein sequence |
| `protein_sequences/*.fasta` | Cached wildtype protein per transcript |
| `neoantigens_primary.csv` | Full-length wildtype and neoantigen per mutation |
| `neoantigen_9mers.csv` | 9-mer (neoantigen_peptide, wt_peptide) + metadata |
| `neoantigen_9mers.list` | Unique 9-mer list for NetMHCpan |
| `nsclc_netmhcpan_out.xls` | NetMHCpan raw predictions |
| `nsclc_netmhcpan_strong_binders.csv` | Strong binders (EL_Rank ≤ 0.5) + metadata |
| `strong_binders_peptide_wt.csv` | peptide, wt only |
| `strong_binders_peptide_wt_metrics.csv` | + blosum/aliphatic/boman distances |
| `strong_binders_best_cdr3b.csv` | + best_cdr3b and binding_score per peptide row |
| `visualizations/paired_distance_*.png` | Distance metric figures (BLOSUM/Boman/aliphatic) |
| `visualizations/report_*.png` | Report-confirmation figures (funnel, by allele, by gene, EL_Rank) |
| `plot_report_visualizations.py` | Script to regenerate report figures from data |

This report reflects the pipeline as implemented in **step6_NSCLC** and the counts from the current data files in that directory.
