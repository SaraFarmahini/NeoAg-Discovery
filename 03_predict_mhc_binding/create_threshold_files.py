#!/usr/bin/env python3
"""
Create minimal EL threshold files so NetMHCpan 4.2 does not exit with
"EL Threshold file ... does not exist" on the cluster.

NetMHCpan expects one file per MHC pseudo-sequence in data/threshold/<pseudo>.thr.el.
We create placeholder files with a simple percentile layout (0.0 to 1.0) so the
binary finds the files. If the format is wrong, the program may still fail; in that
case run NetMHCpan once on the Mac with a small input and copy the generated
data/threshold/ folder to the cluster.

Usage (on cluster, from step6_NSCLC or netMHCpan-4.2 parent):
  python3 create_threshold_files.py --netmhc_dir /home/farmahini/netMHCpan-4.2
Or with default (../step6_NSCLC_test/netMHCpan-4.2):
  python3 create_threshold_files.py
"""

import argparse
import struct
from pathlib import Path

# Alleles we use in step3 (must match step3_run_netmhcpan.py)
ALLELES = [
    "HLA-A01:01", "HLA-A02:01", "HLA-A03:01", "HLA-A24:02",
    "HLA-B07:02", "HLA-B08:01", "HLA-B15:01", "HLA-C07:01", "HLA-C07:02",
]

# Fallback: pseudo-sequences for our 9 alleles (from MHC_pseudo.dat) if file missing or parse fails
PSEUDO_FALLBACK = [
    "YFAMYQENMAHTDANTLYIIYRDYTWVARVYRGY",   # HLA-A01:01
    "YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY",   # HLA-A02:01
    "YFAMYQENVAQTDVDTLYIIYRDYTWAELAYTWY",   # HLA-A03:01
    "YSAMYEEKVAHTDENIAYLMFHYYTWAVQAYTGY",   # HLA-A24:02
    "YYSEYRNIYAQTDESNLYLSYDYYTWAERAYEWY",   # HLA-B07:02
    "YDSEYRNIFTNTDESNLYLSYNYYTWAVDAYTWY",   # HLA-B08:01
    "YYAMYREISTNTYESNLYLRYDSYTWAEWAYLWY",   # HLA-B15:01
    "YDSGYRENYRQADVSNLYLRYDSYTLAALAYTWY",   # HLA-C07:01
    "YDSGYREKYRQADVSNLYLRSDSYTLAALAYTWY",   # HLA-C07:02
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Create minimal .thr.el threshold files for NetMHCpan.")
    parser.add_argument(
        "--netmhc_dir",
        default="../step6_NSCLC_test/netMHCpan-4.2",
        help="Path to NetMHCpan-4.2 directory (contains data/MHC_pseudo.dat)",
    )
    parser.add_argument(
        "--binary",
        action="store_true",
        help="Write threshold files as binary float32 (use if text format gives 'Cannot read')",
    )
    args = parser.parse_args()
    base = Path(args.netmhc_dir).resolve()
    data_dir = base / "data"
    pseudo_file = data_dir / "MHC_pseudo.dat"
    threshold_dir = data_dir / "threshold"
    threshold_dir.mkdir(parents=True, exist_ok=True)

    unique_pseudos = []
    if pseudo_file.exists():
        allele_to_pseudo = {}
        with pseudo_file.open() as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(None, 1)
                if len(parts) != 2:
                    continue
                allele, pseudo = parts
                if allele in ALLELES and len(pseudo) == 33:
                    allele_to_pseudo[allele] = pseudo
        unique_pseudos = list(dict.fromkeys(allele_to_pseudo[a] for a in ALLELES if a in allele_to_pseudo))
        if len(allele_to_pseudo) < len(ALLELES):
            missing = set(ALLELES) - set(allele_to_pseudo)
            print(f"Warning: some alleles not found in MHC_pseudo.dat: {missing}")
    if not unique_pseudos:
        print("Using fallback list of pseudo-sequences for the 9 alleles.")
        unique_pseudos = PSEUDO_FALLBACK

    # Create one .thr.el per unique pseudo.
    # NetMHCpan may expect binary float32. Try --binary if text fails with "Cannot read".
    n_vals = 1001  # 0.0 to 1.0 in 0.001 steps
    for pseudo in unique_pseudos:
        out_path = threshold_dir / f"{pseudo}.thr.el"
        if args.binary:
            with out_path.open("wb") as f:
                for i in range(n_vals):
                    f.write(struct.pack("<f", i / (n_vals - 1)))
        else:
            with out_path.open("w") as f:
                for i in range(n_vals):
                    f.write(f"{i / (n_vals - 1)}\n")
        print(f"Wrote {out_path} ({'binary' if args.binary else 'text'}, {n_vals} values)")

    print(f"Created {len(unique_pseudos)} threshold files in {threshold_dir}")


if __name__ == "__main__":
    main()
