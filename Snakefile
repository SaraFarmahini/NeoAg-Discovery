"""
NSCLC Neoantigen Pipeline — Snakemake workflow

Pipeline order:
  mutation data → mutant peptide generation → MHC prediction → strong binder filtering → genomic interpretation

Run from project root:
  snakemake --cores 4

Prerequisites:
  - Python 3 with pandas, requests (step 1–2), torch (step 4 TCR scoring)
  - NetMHCpan 4.2 for rule predict_mhc (set netmhcpan_dir in config or use default)
  - sequence_encodings module for rule blosum_metrics (BLOSUM62); see compute_blosum_distance.py
"""

# ---------------------------------------------------------------------------
# Paths (config or defaults)
# ---------------------------------------------------------------------------
DIR_01 = "01_nsclc_mutation_data"
DIR_02 = "02_generate_mutant_peptides"
DIR_03 = "03_predict_mhc_binding"
DIR_04 = "04_filter_strong_binders"
DIR_05 = "05_cancer_genomic_interpretation"

COSMIC_TSV = f"{DIR_01}/Cosmic_ResistanceMutations_v101_GRCh38.tsv"
PRIMARY_CSV = f"{DIR_02}/neoantigens_primary.csv"
NEO_9MER_CSV = f"{DIR_02}/neoantigen_9mers.csv"
NEO_9MER_LIST = f"{DIR_02}/neoantigen_9mers.list"
NETMHCPAN_XLS = f"{DIR_03}/nsclc_netmhcpan_out.xls"
STRONG_BINDERS_CSV = f"{DIR_03}/nsclc_netmhcpan_strong_binders.csv"
PEPTIDE_WT_CSV = f"{DIR_04}/strong_binders_peptide_wt.csv"
METRICS_CSV = f"{DIR_04}/strong_binders_peptide_wt_metrics.csv"
BEST_CDR3_CSV = f"{DIR_04}/strong_binders_best_cdr3b.csv"
VIS_DIR = f"{DIR_05}/visualizations"

# ---------------------------------------------------------------------------
# Final outputs (rule all)
# ---------------------------------------------------------------------------
REPORT_FIGS = [
    f"{VIS_DIR}/report_pipeline_funnel.png",
    f"{VIS_DIR}/report_strong_binders_by_allele.png",
    f"{VIS_DIR}/report_strong_binders_by_gene.png",
    f"{VIS_DIR}/report_9mers_by_gene.png",
    f"{VIS_DIR}/report_el_rank_distribution.png",
]
SCHEMATIC_FIGS = [
    f"{VIS_DIR}/nsclc_pipeline_schematic.png",
    f"{VIS_DIR}/nsclc_pipeline_overview.png",
]
CDR3_FIG = f"{VIS_DIR}/cdr3b_binding_score_distribution.png"
LOGIT_FIGS = [
    f"{VIS_DIR}/strong_binders_logit_probability_distributions.png",
    f"{VIS_DIR}/strong_binders_logit_vs_probability.png",
]

rule all:
    input:
        BEST_CDR3_CSV,
        REPORT_FIGS,
        SCHEMATIC_FIGS,
        CDR3_FIG,
        LOGIT_FIGS,


# ---------------------------------------------------------------------------
# Stage 1–2: Mutation data → primary neoantigens → 9-mers
# ---------------------------------------------------------------------------
rule fetch_primary:
    """Step 1: COSMIC TSV → protein sequences + primary neoantigens (full-length)."""
    input:
        COSMIC_TSV,
    output:
        PRIMARY_CSV,
    shell:
        "python {DIR_02}/step1_fetch_proteins_and_primary_neoantigens.py"


rule generate_9mers:
    """Step 2: Primary neoantigens → 9-mer peptides (mutation-including) + list for NetMHCpan."""
    input:
        PRIMARY_CSV,
    output:
        NEO_9MER_CSV,
        NEO_9MER_LIST,
    shell:
        "python {DIR_02}/step2_generate_neoantigen_9mers.py"


# ---------------------------------------------------------------------------
# Stage 3: MHC binding prediction (NetMHCpan)
# ---------------------------------------------------------------------------
rule predict_mhc:
    """Step 3: Run NetMHCpan on 9-mers; filter to strong binders (EL_Rank ≤ 0.5)."""
    input:
        peptide_list=NEO_9MER_LIST,
        neo_csv=NEO_9MER_CSV,
    output:
        xls=NETMHCPAN_XLS,
        strong=STRONG_BINDERS_CSV,
    shell:
        "python {DIR_03}/step3_run_netmhcpan.py"


rule copy_strong_binders:
    """Copy strong binders CSV to 04 for scripts that read from 04."""
    input:
        STRONG_BINDERS_CSV,
    output:
        f"{DIR_04}/nsclc_netmhcpan_strong_binders.csv",
    shell:
        "cp {input} {output}"


# ---------------------------------------------------------------------------
# Stage 4: Strong binder filtering (peptide–WT → metrics → TCR scoring)
# ---------------------------------------------------------------------------
rule peptide_wt:
    """Build peptide/wt table from strong binders for distance metrics."""
    input:
        STRONG_BINDERS_CSV,
    output:
        PEPTIDE_WT_CSV,
    run:
        import pandas as pd
        df = pd.read_csv(input[0])
        if "wt_peptide" not in df.columns:
            raise ValueError("Strong binders CSV must have wt_peptide column")
        out = df[["peptide", "wt_peptide"]].drop_duplicates()
        out = out.rename(columns={"wt_peptide": "wt"})
        out.to_csv(output[0], index=False)


rule blosum_metrics:
    """Compute BLOSUM/Boman/aliphatic distances for peptide–WT pairs."""
    input:
        PEPTIDE_WT_CSV,
    output:
        METRICS_CSV,
    shell:
        "python {DIR_04}/compute_blosum_distance.py --input {input} --output {output}"


rule metrics_unique:
    """One row per unique peptide (dedup of metrics for schematic)."""
    input:
        METRICS_CSV,
    output:
        f"{DIR_04}/strong_binders_peptide_wt_metrics_unique.csv",
    run:
        import pandas as pd
        df = pd.read_csv(input[0])
        df.drop_duplicates(subset=["peptide"], keep="first").to_csv(output[0], index=False)


rule score_tcr:
    """TCR–peptide scoring: best CDR3β per peptide (uses tuned_model + ds_cdr3b_unique.csv)."""
    input:
        peptides=METRICS_CSV,
        cdr3=f"{DIR_04}/ds_cdr3b_unique.csv",
    output:
        BEST_CDR3_CSV,
    params:
        model_py=f"{DIR_04}/tuned_model/final_model_new_hpo_bce_logits.py",
        checkpoint=f"{DIR_04}/tuned_model/finetune/best_model.pt",
    shell:
        "python {DIR_04}/score_peptide_cdr3.py --peptides {input.peptides} --cdr3_list {input.cdr3} --output {output}"


# ---------------------------------------------------------------------------
# Stage 5: Cancer genomic interpretation (visualizations)
# ---------------------------------------------------------------------------
rule plot_report:
    """Report visualizations: funnel, by allele, by gene, EL_Rank distribution."""
    input:
        cosmic=COSMIC_TSV,
        primary=PRIMARY_CSV,
        neo_csv=NEO_9MER_CSV,
        neo_list=NEO_9MER_LIST,
        strong=f"{DIR_04}/nsclc_netmhcpan_strong_binders.csv",
    output:
        REPORT_FIGS,
    shell:
        "python {DIR_05}/plot_report_visualizations.py"


rule plot_schematic:
    """Pipeline schematic and overview figures."""
    input:
        cosmic=COSMIC_TSV,
        primary=PRIMARY_CSV,
        neo_csv=NEO_9MER_CSV,
        metrics_unique=f"{DIR_04}/strong_binders_peptide_wt_metrics_unique.csv",
    output:
        SCHEMATIC_FIGS,
    shell:
        "python {DIR_05}/plot_pipeline_schematic.py"


rule plot_cdr3b:
    """CDR3β binding score distribution."""
    input:
        BEST_CDR3_CSV,
    output:
        CDR3_FIG,
    shell:
        "python {DIR_05}/plot_cdr3b_score_distribution.py"


rule plot_logit:
    """Logit and probability distribution figures."""
    input:
        BEST_CDR3_CSV,
    output:
        LOGIT_FIGS,
    shell:
        "python {DIR_05}/plot_logit_probability_distributions.py"
