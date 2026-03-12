#!/usr/bin/env python3
"""
Run NetMHCpan in chunks (e.g. 200 peptides per chunk) then merge results.
Use this on your Mac if the full run gets stuck, or on the cluster if you only
have the Mac binary (not recommended on cluster; Linux build is better).

Produces the same outputs as step3_run_netmhcpan.py: nsclc_netmhcpan_out.xls
and nsclc_netmhcpan_strong_binders.csv.
"""

import argparse
import sys
from pathlib import Path

# Reuse step3 logic
from step3_run_netmhcpan import _run_netmhcpan, _parse_netmhcpan_xls


def main() -> int:
    parser = argparse.ArgumentParser(description="Run NetMHCpan in chunks and merge.")
    parser.add_argument("--peptide_list", default="neoantigen_9mers.list")
    parser.add_argument("--neo_csv", default="neoantigen_9mers.csv")
    parser.add_argument("--netmhcpan_dir", default="../step6_NSCLC_test/netMHCpan-4.2")
    parser.add_argument("--platform", default="Darwin_arm64")
    parser.add_argument("--netmhcpan_bin", default="")
    parser.add_argument(
        "--alleles",
        default=(
            "HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,"
            "HLA-B07:02,HLA-B08:01,HLA-B15:01,HLA-C07:01,HLA-C07:02"
        ),
    )
    parser.add_argument("--mode", type=int, default=2)
    parser.add_argument("--synlist", default="synlist_cedar.bin")
    parser.add_argument("--out_xls", default="nsclc_netmhcpan_out.xls")
    parser.add_argument("--out_csv", default="nsclc_netmhcpan_strong_binders.csv")
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=200,
        help="Peptides per chunk (smaller = more stable, more chunks)",
    )
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    peptide_list = base_dir / args.peptide_list
    neo_csv = base_dir / args.neo_csv
    out_xls = base_dir / args.out_xls
    out_csv = base_dir / args.out_csv

    if not peptide_list.exists():
        print(f"Error: peptide list not found: {peptide_list}", file=sys.stderr)
        return 1
    if not neo_csv.exists():
        print(f"Error: neo CSV not found: {neo_csv}", file=sys.stderr)
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
    if not peptides:
        print("Error: no peptides in list", file=sys.stderr)
        return 1

    chunk_size = max(1, args.chunk_size)
    chunks = [
        peptides[i : i + chunk_size]
        for i in range(0, len(peptides), chunk_size)
    ]
    n_chunks = len(chunks)
    print(f"Running NetMHCpan in {n_chunks} chunks of up to {chunk_size} peptides each ({len(peptides)} total)")

    chunk_xls_files = []
    for i, chunk in enumerate(chunks):
        run_list = base_dir / f"netmhcpan_chunk_{i}.list"
        chunk_xls = base_dir / f"nsclc_netmhcpan_chunk_{i}.xls"
        run_list.write_text("\n".join(chunk) + "\n")
        print(f"  Chunk {i + 1}/{n_chunks}: {len(chunk)} peptides -> {chunk_xls.name}")
        _run_netmhcpan(
            netmhcpan_bin, netmhc_home, run_list,
            args.alleles, chunk_xls, args.mode, synlist,
        )
        chunk_xls_files.append(chunk_xls)
        run_list.unlink(missing_ok=True)

    # Merge XLS: first 2 lines (header) from first chunk, then data rows from all chunks
    with out_xls.open("w") as out:
        with chunk_xls_files[0].open() as f:
            out.write(f.readline())
            out.write(f.readline())
        for xls_path in chunk_xls_files:
            with xls_path.open() as f:
                f.readline()
                f.readline()
                for line in f:
                    out.write(line)
    for p in chunk_xls_files:
        p.unlink(missing_ok=True)

    strong_count = _parse_netmhcpan_xls(out_xls, out_csv, neo_csv)
    print(f"Wrote {out_xls.name} and {out_csv.name} (strong binders: {strong_count})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
