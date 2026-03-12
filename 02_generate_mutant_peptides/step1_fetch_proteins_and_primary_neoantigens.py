#!/usr/bin/env python3
"""
Step 1: COSMIC mutations → protein sequences and primary neoantigen sequences.

Reads Cosmic_ResistanceMutations_v101_GRCh38.tsv, optionally cleans it, fetches
wildtype protein sequences from Ensembl for each transcript, parses HGVSP to get
mutation position and amino acids, and produces:

1. Updated_Cosmic_Data_with_Sequences.tsv  — mutations + protein_sequence column
2. protein_sequences/<TRANSCRIPT>.fasta    — cached wildtype protein FASTA per transcript
3. neoantigens_primary.csv                 — primary form: full-length wildtype and
   neoantigen (mutant) protein sequences per mutation, NO 9-mer splitting.

Requires: pandas, requests
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import pandas as pd

try:
    import requests
except ImportError:
    print("Install requests: pip install requests", file=sys.stderr)
    sys.exit(1)


# -----------------------------------------------------------------------------
# Amino acid 3-letter → 1-letter
# -----------------------------------------------------------------------------
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*", "Stop": "*",
}


def parse_hgvsp(hgvsp: str) -> dict | None:
    """
    Parse HGVSP string (e.g. ENSP00000275493.2:p.Thr790Met or p.Thr790Met).
    Returns dict with keys: wt (1-letter), mut (1-letter), pos (1-based int).
    """
    if pd.isna(hgvsp) or not str(hgvsp).strip():
        return None
    text = str(hgvsp).strip()
    if ":" in text:
        _, text = text.split(":", 1)
    match = re.search(r"p\.([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})", text)
    if not match:
        return None
    wt_raw, pos_raw, mut_raw = match.groups()
    pos = int(pos_raw)
    wt = AA3_TO_1.get(wt_raw.capitalize(), wt_raw.upper() if len(wt_raw) == 1 else None)
    mut = AA3_TO_1.get(mut_raw.capitalize(), mut_raw.upper() if len(mut_raw) == 1 else None)
    if wt is None or mut is None or len(wt) != 1 or len(mut) != 1:
        return None
    return {"wt": wt, "mut": mut, "pos": pos}


def normalize_protein_sequence(raw: str) -> str:
    """Extract a single continuous amino acid sequence from FASTA or multiline text."""
    if pd.isna(raw):
        return ""
    lines = [line.strip() for line in str(raw).splitlines() if line.strip()]
    lines = [line for line in lines if not line.startswith(">")]
    return "".join(lines)


def fetch_protein_from_ensembl(transcript_id: str, timeout: int = 30) -> str | None:
    """
    Fetch protein sequence for an Ensembl transcript ID.
    Uses REST: GET /sequence/id/{id}?type=protein&format=fasta
    Returns the sequence string (no header) or None on failure.
    """
    base_id = transcript_id.split(".")[0]
    url = f"https://rest.ensembl.org/sequence/id/{base_id}"
    params = {"type": "protein", "format": "fasta"}
    headers = {"Content-Type": "text/x-fasta"}
    try:
        r = requests.get(url, params=params, headers=headers, timeout=timeout)
        if r.status_code != 200:
            return None
        return normalize_protein_sequence(r.text)
    except requests.RequestException:
        return None


def apply_mutation(seq: str, pos_1based: int, wt_aa: str, mut_aa: str) -> str | None:
    """
    Apply a single amino acid mutation at 1-based position.
    Returns mutant sequence or None if validation fails (wrong wt at position).
    """
    if not seq or pos_1based < 1 or pos_1based > len(seq):
        return None
    idx = pos_1based - 1
    if seq[idx] != wt_aa:
        return None
    return seq[:idx] + mut_aa + seq[idx + 1 :]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Step 1: Fetch protein sequences and build primary neoantigen sequences (full-length, no 9-mer split)."
    )
    parser.add_argument(
        "--input",
        default="Cosmic_ResistanceMutations_v101_GRCh38.tsv",
        help="Input COSMIC TSV path",
    )
    parser.add_argument(
        "--out_tsv",
        default="Updated_Cosmic_Data_with_Sequences.tsv",
        help="Output TSV with protein_sequence column",
    )
    parser.add_argument(
        "--out_primary",
        default="neoantigens_primary.csv",
        help="Output CSV: primary wildtype and neoantigen full sequences (no 9-mers)",
    )
    parser.add_argument(
        "--fasta_dir",
        default="protein_sequences",
        help="Directory to cache FASTA files (one per transcript)",
    )
    parser.add_argument(
        "--no_clean",
        action="store_true",
        help="Skip cleaning (no dedup, no drop of ? or somatic filter)",
    )
    parser.add_argument(
        "--somatic_only",
        action="store_true",
        help="Keep only rows with MUTATION_SOMATIC_STATUS == 'Confirmed somatic variant'",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.2,
        help="Seconds between Ensembl API requests (default 0.2)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Max number of unique transcripts to fetch (0 = no limit)",
    )
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    # Support organized layout: 01_nsclc_mutation_data, 02_generate_mutant_peptides, ...
    project_root = base_dir.parent if (base_dir.parent / "01_nsclc_mutation_data").exists() else base_dir
    if (project_root / "01_nsclc_mutation_data").exists():
        dir_01 = project_root / "01_nsclc_mutation_data"
        dir_02 = project_root / "02_generate_mutant_peptides"
        input_path = dir_01 / "Cosmic_ResistanceMutations_v101_GRCh38.tsv"
        out_tsv_path = dir_01 / "Updated_Cosmic_Data_with_Sequences.tsv"
        fasta_dir = dir_01 / "protein_sequences"
        out_primary_path = dir_02 / "neoantigens_primary.csv"
    else:
        input_path = base_dir / args.input
        out_tsv_path = base_dir / args.out_tsv
        fasta_dir = base_dir / args.fasta_dir
        out_primary_path = base_dir / args.out_primary
    if not input_path.exists():
        print(f"Error: input not found: {input_path}", file=sys.stderr)
        return 1

    fasta_dir = Path(fasta_dir)
    fasta_dir.mkdir(parents=True, exist_ok=True)

    # ----- Load and optionally clean -----
    df = pd.read_csv(input_path, sep="\t", low_memory=False)
    df.columns = df.columns.str.strip()
    n_raw = len(df)

    if not args.no_clean:
        # Drop duplicates by mutation identity (keep first)
        df = df.drop_duplicates(
            subset=["GENE_SYMBOL", "MUTATION_AA", "MUTATION_CDS"],
            keep="first",
        )
        # Drop rows with ? in MUTATION_AA or MUTATION_CDS
        for col in ["MUTATION_AA", "MUTATION_CDS"]:
            if col in df.columns:
                df = df[~df[col].astype(str).str.contains(r"\?", na=False, regex=True)]
        if args.somatic_only and "MUTATION_SOMATIC_STATUS" in df.columns:
            df = df[df["MUTATION_SOMATIC_STATUS"].astype(str).str.strip() == "Confirmed somatic variant"]
        print(f"After cleaning: {len(df)} rows (from {n_raw})")

    # Require columns
    for col in ["TRANSCRIPT_ACCESSION", "HGVSP", "MUTATION_ID", "GENE_SYMBOL"]:
        if col not in df.columns:
            print(f"Error: missing column {col}", file=sys.stderr)
            return 1

    # Drop rows with missing or unparseable HGVSP
    df = df[df["HGVSP"].notna() & (df["HGVSP"].astype(str).str.strip() != "")]
    parsed = df["HGVSP"].map(parse_hgvsp)
    df = df[parsed.notna()].copy()
    parsed = parsed[parsed.notna()]
    print(f"Rows with valid HGVSP: {len(df)}")

    # ----- Unique transcripts and fetch sequences -----
    transcripts = df["TRANSCRIPT_ACCESSION"].astype(str).dropna().unique().tolist()
    if args.limit > 0:
        transcripts = transcripts[: args.limit]
        print(f"Limited to {len(transcripts)} transcripts")
    sequences: dict[str, str] = {}
    for i, tid in enumerate(transcripts):
        fasta_path = fasta_dir / f"{tid.split('.')[0]}.fasta"
        if fasta_path.exists():
            seq = normalize_protein_sequence(fasta_path.read_text())
            if seq:
                sequences[tid] = seq
        else:
            seq = fetch_protein_from_ensembl(tid)
            if seq:
                sequences[tid] = seq
                fasta_path.write_text(f">{tid}\n{seq}\n")
            if args.delay > 0:
                time.sleep(args.delay)
        if (i + 1) % 50 == 0:
            print(f"  Fetched/cached {i + 1}/{len(transcripts)} transcripts")
    print(f"Protein sequences obtained: {len(sequences)} unique transcripts")

    # ----- Add protein_sequence to dataframe -----
    df = df.copy()
    df["protein_sequence"] = df["TRANSCRIPT_ACCESSION"].astype(str).map(sequences)
    df_with_seq = df[df["protein_sequence"].notna() & (df["protein_sequence"].astype(str).str.len() > 0)]
    df_with_seq.to_csv(out_tsv_path, sep="\t", index=False)
    print(f"Wrote {len(df_with_seq)} rows -> {out_tsv_path.name}")

    # ----- Build primary neoantigens: full-length wt and mutant -----
    primary_rows = []
    for _, row in df_with_seq.iterrows():
        seq = row["protein_sequence"]
        if not seq or not isinstance(seq, str):
            continue
        p = parse_hgvsp(row["HGVSP"])
        if not p:
            continue
        mut_seq = apply_mutation(seq, p["pos"], p["wt"], p["mut"])
        if mut_seq is None:
            continue
        primary_rows.append({
            "MUTATION_ID": row["MUTATION_ID"],
            "TRANSCRIPT_ACCESSION": row["TRANSCRIPT_ACCESSION"],
            "GENE_SYMBOL": row["GENE_SYMBOL"],
            "HGVSP": row["HGVSP"],
            "MUTATION_AA": row.get("MUTATION_AA"),
            "wt_aa": p["wt"],
            "mut_aa": p["mut"],
            "mutation_pos": p["pos"],
            "wildtype_protein_sequence": seq,
            "neoantigen_protein_sequence": mut_seq,
            "sequence_length": len(seq),
        })
    primary_df = pd.DataFrame(primary_rows)
    primary_df.to_csv(out_primary_path, index=False)
    print(f"Wrote {len(primary_df)} primary neoantigens (full-length wt + mutant) -> {out_primary_path.name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
