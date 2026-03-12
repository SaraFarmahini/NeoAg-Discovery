#!/usr/bin/env python3
"""
Quality check: detect promiscuous CDR3β that appear as "best" for many peptides.
If the same CDR3β is best for every (or nearly every) peptide, it may be a
promiscuous binder rather than peptide-specific.

Usage:
  python check_promiscuous_cdr3.py strong_binders_best_cdr3b.csv
  python check_promiscuous_cdr3.py strong_binders_best_cdr3b.csv --top 20 --output report.txt
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Check for promiscuous best CDR3β in scoring results.")
    parser.add_argument("csv", type=str, help="Output CSV from score_peptide_cdr3.py (must have best_cdr3b)")
    parser.add_argument("--top", type=int, default=15, help="Show top N CDR3s by number of peptides they are best for")
    parser.add_argument("--warn-if-any-above", type=int, default=50, help="Warn if any CDR3 is best for more than this many peptides")
    parser.add_argument("--output", "-o", type=str, default="promiscuity_check_report.txt", help="Save report to this file (default: promiscuity_check_report.txt)")
    args = parser.parse_args()

    path = Path(args.csv)
    if not path.exists():
        print(f"Error: file not found: {path}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(path)
    if "best_cdr3b" not in df.columns:
        print("Error: CSV must have column 'best_cdr3b'. Columns:", list(df.columns), file=sys.stderr)
        sys.exit(1)

    n_peptides = len(df)
    counts = df["best_cdr3b"].value_counts()

    # Promiscuous: same CDR3 best for all (or nearly all) peptides
    max_count = counts.iloc[0] if len(counts) else 0
    n_unique_best = len(counts)

    lines = []
    lines.append(f"Total peptide rows: {n_peptides}")
    lines.append(f"Unique CDR3β that are 'best' for at least one peptide: {n_unique_best}")
    lines.append("")

    if max_count == n_peptides and n_unique_best == 1:
        lines.append("*** PROMISCUOUS BINDER WARNING ***")
        lines.append("The same CDR3β is the best match for every single peptide.")
        lines.append("This may indicate a promiscuous TCR; consider checking the model or data.")
        lines.append(f"CDR3β: {counts.index[0]}")
        report = "\n".join(lines)
        print(report)
        out_path = Path(args.output)
        out_path.write_text(report, encoding="utf-8")
        print(f"Saved report to {out_path}", file=sys.stderr)
        sys.exit(2)

    if max_count >= args.warn_if_any_above:
        lines.append("*** WARNING ***")
        lines.append(f"At least one CDR3β is best for {max_count} peptides (>{args.warn_if_any_above}).")
        lines.append("Top promiscuous candidates:")
        for cdr3, c in counts.head(5).items():
            lines.append(f"  {c:5d} peptides  {cdr3}")
        lines.append("")

    lines.append(f"Top {args.top} CDR3β by number of peptides they are 'best' for:")
    lines.append("-" * 60)
    for cdr3, c in counts.head(args.top).items():
        pct = 100.0 * c / n_peptides
        lines.append(f"  {c:5d} ({pct:5.1f}%)  {cdr3}")
    lines.append("")
    lines.append("Check complete. Many peptides sharing the same best CDR3 may warrant further inspection.")

    report = "\n".join(lines)
    print(report)

    out_path = Path(args.output)
    out_path.write_text(report, encoding="utf-8")
    print(f"Saved report to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
