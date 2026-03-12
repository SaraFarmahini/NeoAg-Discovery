#!/usr/bin/env python3
"""
Step 2: From primary neoantigens (full-length) → only 9-mers that INCLUDE the mutation.

Reads neoantigens_primary.csv. For each row, generates 9-mer windows that span the
mutation position and that differ from wildtype (mut_pep != wt_pep). Only these
are true neoantigen 9-mers; 9-mers that don't contain the mutation are wildtype.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def generate_mutant_kmers(wt_seq: str, mut_seq: str, mut_pos: int, k: int = 9):
    """Yield (start_1based, mut_9mer, wt_9mer) for 9-mers that contain the mutation and differ from WT."""
    if not wt_seq or not mut_seq or len(wt_seq) != len(mut_seq):
        return
    idx = mut_pos - 1
    if idx < 0 or idx >= len(mut_seq):
        return
    start_min = max(0, idx - k + 1)
    start_max = min(idx, len(mut_seq) - k)
    for start in range(start_min, start_max + 1):
        mut_pep = mut_seq[start : start + k]
        wt_pep = wt_seq[start : start + k]
        if len(mut_pep) == k and mut_pep != wt_pep:
            yield start + 1, mut_pep, wt_pep


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate only 9-mers that include the mutation (neoantigen 9-mers)."
    )
    parser.add_argument(
        "--input",
        default="neoantigens_primary.csv",
        help="Primary neoantigens CSV from step1",
    )
    parser.add_argument(
        "--out_csv",
        default="neoantigen_9mers.csv",
        help="Output CSV: one row per neoantigen 9-mer (mutation included)",
    )
    parser.add_argument(
        "--out_list",
        default="neoantigen_9mers.list",
        help="Output peptide list (one 9-mer per line) for NetMHCpan",
    )
    parser.add_argument("--k", type=int, default=9, help="Peptide length (default 9)")
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    # Support organized layout with numbered stage folders
    project_root = base_dir.parent if (base_dir.parent / "01_nsclc_mutation_data").exists() else base_dir
    if (project_root / "02_generate_mutant_peptides").exists():
        input_path = project_root / "02_generate_mutant_peptides" / "neoantigens_primary.csv"
        out_csv = project_root / "02_generate_mutant_peptides" / args.out_csv
        out_list = project_root / "02_generate_mutant_peptides" / args.out_list
    else:
        input_path = base_dir / args.input
        out_csv = base_dir / args.out_csv
        out_list = base_dir / args.out_list
    if not input_path.exists():
        print(f"Error: not found: {input_path}. Run step1 first.", file=sys.stderr)
        return 1

    df = pd.read_csv(input_path)
    for col in ["wildtype_protein_sequence", "neoantigen_protein_sequence", "mutation_pos"]:
        if col not in df.columns:
            print(f"Error: missing column {col}", file=sys.stderr)
            return 1

    rows = []
    for _, row in df.iterrows():
        wt = str(row["wildtype_protein_sequence"]).strip()
        mut = str(row["neoantigen_protein_sequence"]).strip()
        pos = int(row["mutation_pos"])
        if not wt or not mut:
            continue
        for start_1based, mut_pep, wt_pep in generate_mutant_kmers(wt, mut, pos, k=args.k):
            rows.append({
                "neoantigen_peptide": mut_pep,
                "wt_peptide": wt_pep,
                "peptide_length": args.k,
                "mutation_pos": pos,
                "peptide_start_pos": start_1based,
                "wt_aa": row["wt_aa"],
                "mut_aa": row["mut_aa"],
                "MUTATION_ID": row["MUTATION_ID"],
                "TRANSCRIPT_ACCESSION": row["TRANSCRIPT_ACCESSION"],
                "GENE_SYMBOL": row["GENE_SYMBOL"],
                "HGVSP": row["HGVSP"],
                "MUTATION_AA": row.get("MUTATION_AA"),
            })

    out_df = pd.DataFrame(rows)
    if out_df.empty:
        print("No neoantigen 9-mers generated. Check primary CSV.", file=sys.stderr)
        return 1

    out_df.to_csv(out_csv, index=False)
    print(f"Wrote {len(out_df)} neoantigen 9-mers (mutation included) -> {out_csv.name}")

    peptides = out_df["neoantigen_peptide"].drop_duplicates().tolist()
    out_list.write_text("\n".join(peptides) + "\n")
    print(f"Wrote {len(peptides)} unique peptides -> {out_list.name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
