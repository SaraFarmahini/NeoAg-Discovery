#!/usr/bin/env python3
"""
Step 3: Run NetMHCpan on neoantigen 9-mers (same setup as step6_NSCLC_test/Thesis/test).

Uses peptide list from step2 output, runs NetMHCpan, filters to strong binders (EL_Rank <= 0.5),
and merges with neoantigen_9mers.csv to produce strong binders with full metadata.
"""

import argparse
import subprocess
import sys
import shutil
from pathlib import Path

import pandas as pd


def _run_netmhcpan(
    netmhcpan_bin: Path,
    netmhc_home: Path,
    peptide_list: Path,
    alleles: str,
    out_xls: Path,
    mode: int,
    synlist: Path,
) -> None:
    tmp_template = netmhc_home / "tmp" / "netMHCpan_XXXXXX"
    tmp_template.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(netmhcpan_bin),
        "-rdir", str(netmhc_home),
        "-tdir", str(tmp_template),
        "-syn", str(synlist),
        "-hlapseudo", str(netmhc_home / "data" / "MHC_pseudo.dat"),
        "-version", str(netmhc_home / "data" / "version"),
        "-allname", str(netmhc_home / "data" / "allelenames"),
        "-mode", str(mode),
        "-BA", "-p", "-f", str(peptide_list),
        "-a", alleles,
        "-xls", "-xlsfile", str(out_xls),
    ]
    subprocess.run(cmd, check=True)


def _parse_netmhcpan_xls(raw_path: Path, out_path: Path, neo_csv: Path) -> int:
    with raw_path.open("r") as f:
        line1 = f.readline().rstrip("\n")
        line2 = f.readline().rstrip("\n")
    cols1 = line1.split("\t")
    cols2 = line2.split("\t")
    alleles = []
    current = ""
    for c in cols1:
        if c.strip():
            current = c.strip()
        alleles.append(current)
    raw = pd.read_csv(raw_path, sep="\t", skiprows=2, header=None)
    raw.columns = cols2
    col_allele = {i: alleles[i] for i in range(len(alleles))}
    records = []
    for _, row in raw.iterrows():
        pep = row.get("Peptide")
        if pd.isna(pep):
            continue
        for allele in sorted(set(a for a in alleles if a)):
            idxs = [i for i, a in col_allele.items() if a == allele]
            values = {cols2[i]: row.iloc[i] for i in idxs if i < len(cols2)}
            if not values:
                continue
            # NetMHCpan XLS may use "Rank"/"Score" or "EL_Rank"/"EL-score" per allele
            el_rank = values.get("EL_Rank") or values.get("Rank")
            el_score = values.get("EL-score") or values.get("Score")
            ba_rank = values.get("BA_Rank")
            ba_score = values.get("BA-score")
            records.append({
                "peptide": pep,
                "allele": allele,
                "EL_Rank": el_rank,
                "EL_score": el_score,
                "BA_Rank": ba_rank,
                "BA_score": ba_score,
            })
    long_df = pd.DataFrame(records)
    rank_vals = pd.to_numeric(long_df["EL_Rank"], errors="coerce")
    strong = long_df[rank_vals <= 0.5].copy()
    neo_df = pd.read_csv(neo_csv)
    merged = strong.merge(
        neo_df, how="left", left_on="peptide", right_on="neoantigen_peptide"
    )
    # Fill missing wt_peptide from neoantigen_9mers (one wt per neoantigen_peptide)
    if "wt_peptide" in neo_df.columns and "wt_peptide" in merged.columns:
        pep_to_wt = neo_df.set_index("neoantigen_peptide")["wt_peptide"].to_dict()
        missing = merged["wt_peptide"].isna() | (merged["wt_peptide"].astype(str).str.strip() == "")
        merged.loc[missing, "wt_peptide"] = merged.loc[missing, "peptide"].map(pep_to_wt)
    # Normalise column names for output (keep both if present)
    out_cols = [
        "peptide", "allele", "EL_Rank", "EL_score", "BA_Rank", "BA_score",
        "neoantigen_peptide", "wt_peptide", "peptide_length", "mutation_pos",
        "peptide_start_pos", "wt_aa", "mut_aa",
        "GENE_SYMBOL", "TRANSCRIPT_ACCESSION", "MUTATION_ID", "HGVSP", "MUTATION_AA",
    ]
    for c in out_cols:
        if c not in merged.columns:
            merged[c] = ""
    merged = merged[[c for c in out_cols if c in merged.columns]]
    merged.to_csv(out_path, index=False)
    return len(merged)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run NetMHCpan on neoantigen 9-mers.")
    parser.add_argument(
        "--peptide_list",
        default="neoantigen_9mers.list",
        help="Peptide list (one per line) from step2",
    )
    parser.add_argument(
        "--neo_csv",
        default="neoantigen_9mers.csv",
        help="Neoantigen 9-mers CSV for merge",
    )
    parser.add_argument(
        "--netmhcpan_dir",
        default="../step6_NSCLC_test/netMHCpan-4.2",
        help="Path to NetMHCpan-4.2 directory",
    )
    parser.add_argument(
        "--platform",
        default="Darwin_arm64",
        help="Platform subfolder (Darwin_arm64 or Linux_x86_64)",
    )
    parser.add_argument(
        "--netmhcpan_bin",
        default="",
        help="Override: full path to netMHCpan binary",
    )
    parser.add_argument(
        "--alleles",
        default=(
            "HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,"
            "HLA-B07:02,HLA-B08:01,HLA-B15:01,HLA-C07:01,HLA-C07:02"
        ),
    )
    parser.add_argument("--mode", type=int, default=2, choices=[0, 1, 2])
    parser.add_argument("--synlist", default="synlist_cedar.bin")
    parser.add_argument("--out_xls", default="nsclc_netmhcpan_out.xls")
    parser.add_argument("--out_csv", default="nsclc_netmhcpan_strong_binders.csv")
    parser.add_argument("--limit", type=int, default=0, help="Use only first N peptides (0 = all, for testing)")
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    # Support organized layout: read from 02_generate_mutant_peptides, write to 03_predict_mhc_binding
    project_root = base_dir.parent if (base_dir.parent / "01_nsclc_mutation_data").exists() else base_dir
    if (project_root / "02_generate_mutant_peptides").exists() and (project_root / "03_predict_mhc_binding").exists():
        dir_02 = project_root / "02_generate_mutant_peptides"
        dir_03 = project_root / "03_predict_mhc_binding"
        peptide_list = dir_02 / "neoantigen_9mers.list"
        neo_csv = dir_02 / "neoantigen_9mers.csv"
        out_xls = dir_03 / args.out_xls
        out_csv = dir_03 / args.out_csv
    else:
        peptide_list = base_dir / args.peptide_list
        neo_csv = base_dir / args.neo_csv
        out_xls = base_dir / args.out_xls
        out_csv = base_dir / args.out_csv

    if not peptide_list.exists():
        print(f"Error: peptide list not found: {peptide_list}", file=sys.stderr)
        return 1
    if not neo_csv.exists():
        print(f"Error: neoantigen CSV not found: {neo_csv}", file=sys.stderr)
        return 1

    netmhc_home = (base_dir / args.netmhcpan_dir).resolve()
    if not netmhc_home.exists():
        print(f"Error: NetMHCpan dir not found: {netmhc_home}", file=sys.stderr)
        return 1
    if args.netmhcpan_bin:
        netmhcpan_bin = Path(args.netmhcpan_bin).expanduser().resolve()
    else:
        netmhcpan_bin = netmhc_home / args.platform / "bin" / "netMHCpan-4.2"
    if not netmhcpan_bin.exists():
        print(f"Error: binary not found: {netmhcpan_bin}", file=sys.stderr)
        return 1
    synlist = netmhc_home / "data" / args.synlist
    if not synlist.exists():
        print(f"Error: synlist not found: {synlist}", file=sys.stderr)
        return 1

    peptides = [p for p in peptide_list.read_text().strip().splitlines() if p.strip()]
    if args.limit > 0:
        peptides = peptides[: args.limit]
        print(f"Limited to first {len(peptides)} peptides (--limit {args.limit})")
    n_pep = len(peptides)
    if n_pep == 0:
        print("Error: no peptides in list", file=sys.stderr)
        return 1
    # Write possibly limited list to a temp file for NetMHCpan
    run_list = base_dir / "netmhcpan_peptide_list.txt"
    run_list.write_text("\n".join(peptides) + "\n")
    print(f"Running NetMHCpan on {n_pep} peptides -> {out_xls.name} (may take several minutes)")
    _run_netmhcpan(
        netmhcpan_bin, netmhc_home, run_list,
        args.alleles, out_xls, args.mode, synlist,
    )
    strong_count = _parse_netmhcpan_xls(out_xls, out_csv, neo_csv)
    print(f"Wrote strong binders (EL_Rank <= 0.5) -> {out_csv.name} (rows={strong_count})")
    # In organized layout, copy to 04 for downstream steps
    dir_04 = project_root / "04_filter_strong_binders"
    if dir_04.exists() and out_csv.parent == project_root / "03_predict_mhc_binding":
        dest = dir_04 / out_csv.name
        shutil.copy(out_csv, dest)
        print(f"Copied to {dest} for step 4")
    return 0


if __name__ == "__main__":
    sys.exit(main())
