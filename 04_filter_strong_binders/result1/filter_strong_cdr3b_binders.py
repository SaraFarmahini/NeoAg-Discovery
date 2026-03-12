#!/usr/bin/env python3
"""
Filter strong CDR3β binders from best CDR3β scoring results.

Similar to NetMHCpan's EL rank ≤ 0.5 threshold for strong MHC binders,
this script filters CDR3β-peptide pairs based on binding score thresholds.

Options:
- Probability threshold: filter by binding probability (sigmoid of logit)
- Logit threshold: filter by raw logit score
- Percentile threshold: filter by top X% of scores
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def sigmoid(x: np.ndarray) -> np.ndarray:
    """Convert logits to probabilities."""
    return 1 / (1 + np.exp(-np.clip(x, -500, 500)))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter strong CDR3β binders based on binding score threshold"
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input CSV with binding_score column (e.g., strong_binders_best_cdr3b1.csv)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output CSV path (default: input with '_strong_cdr3b' suffix)",
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=["probability", "logit", "percentile"],
        default="probability",
        help="Filtering method: probability (prob > threshold), logit (logit > threshold), or percentile (top X%%)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Threshold value: for probability/logit methods, minimum value; for percentile, top X%% (e.g., 0.25 = top 25%%)",
    )
    parser.add_argument(
        "--min_logit",
        type=float,
        default=None,
        help="Optional: also apply minimum logit threshold (e.g., 0.0 to require positive logits)",
    )
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    if args.output:
        output_path = Path(args.output).resolve()
    else:
        output_path = input_path.parent / f"{input_path.stem}_strong_cdr3b{input_path.suffix}"

    # Load data
    df = pd.read_csv(input_path)
    if "binding_score" not in df.columns:
        raise ValueError(f"Input CSV must have 'binding_score' column. Found columns: {list(df.columns)}")

    scores = df["binding_score"].values
    n_total = len(df)

    # Apply filtering
    if args.method == "probability":
        probs = sigmoid(scores)
        mask = probs > args.threshold
        threshold_desc = f"probability > {args.threshold:.3f}"
    elif args.method == "logit":
        mask = scores > args.threshold
        threshold_desc = f"logit > {args.threshold:.2f}"
    elif args.method == "percentile":
        # Top X% means scores above (100-X)th percentile
        percentile_threshold = np.percentile(scores, (1 - args.threshold) * 100)
        mask = scores > percentile_threshold
        threshold_desc = f"top {args.threshold*100:.1f}%% (logit > {percentile_threshold:.2f})"
    else:
        raise ValueError(f"Unknown method: {args.method}")

    # Optional minimum logit filter
    if args.min_logit is not None:
        mask = mask & (scores > args.min_logit)
        threshold_desc += f" and logit > {args.min_logit:.2f}"

    filtered_df = df[mask].copy()

    # Add probability column for reference
    if "binding_probability" not in filtered_df.columns:
        filtered_df["binding_probability"] = sigmoid(filtered_df["binding_score"].values)

    # Save results
    filtered_df.to_csv(output_path, index=False)

    # Print summary
    print(f"Input: {n_total} CDR3β-peptide pairs")
    print(f"Filter: {threshold_desc}")
    print(f"Output: {len(filtered_df)} strong CDR3β binders ({len(filtered_df)/n_total*100:.1f}%%)")
    print(f"\nScore statistics (filtered):")
    print(f"  Logit: mean={filtered_df['binding_score'].mean():.2f}, median={filtered_df['binding_score'].median():.2f}")
    print(f"  Probability: mean={filtered_df['binding_probability'].mean():.3f}, median={filtered_df['binding_probability'].median():.3f}")
    print(f"\nSaved to: {output_path}")


if __name__ == "__main__":
    main()
