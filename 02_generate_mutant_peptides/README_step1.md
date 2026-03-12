# Step 1: Protein sequences and primary neoantigens (full-length)

Input: `Cosmic_ResistanceMutations_v101_GRCh38.tsv`  
Output:
- **Updated_Cosmic_Data_with_Sequences.tsv** — mutations + `protein_sequence` column (wildtype full length).
- **protein_sequences/*.fasta** — cached wildtype protein per transcript.
- **neoantigens_primary.csv** — primary form: for each mutation, **wildtype_protein_sequence** (full) and **neoantigen_protein_sequence** (full, with mutation applied). No 9-mer splitting.

## Run

```bash
cd /path/to/step6_NSCLC
pip install pandas requests   # if needed
python step1_fetch_proteins_and_primary_neoantigens.py
```

Optional:
- `--no_clean` — skip dedup and drop of `?` in MUTATION_AA/MUTATION_CDS.
- `--somatic_only` — keep only `MUTATION_SOMATIC_STATUS == 'Confirmed somatic variant'`.
- `--delay 0.3` — pause between Ensembl requests (default 0.2 s).
- `--limit 100` — fetch only first 100 unique transcripts (for testing).

## Requirements

- Python 3.8+
- pandas, requests
