#!/usr/bin/env python3
"""
Run NetTCR-2.0 (chain b: CDR3β only) to predict binding scores for nettcr2_test.csv.

NetTCR-2.0 does not provide pretrained weights for inference-only; it trains on a
training file then predicts on a test file. This script:

1. Prepares nettcr2_test.csv in NetTCR-2.0 format (CDR3b, peptide)
2. Uses a training file: --train_file, or all.csv (step5 neoantigen fine-tuning data), or ds.csv, or NetTCR-2.0 sample
3. Runs the original nettcr.py from NetTCR-2.0 repo

REQUIREMENTS:
- Clone NetTCR-2.0: git clone https://github.com/mnielLab/NetTCR-2.0.git
- Install: tensorflow, keras, numpy, pandas, scikit-learn

Usage:
  python run_nettcr2_predict.py [--nettcr2_repo /path/to/NetTCR-2.0] [--train_file /path/to/train.csv]
"""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
STEP6_DIR = SCRIPT_DIR
TEST_CSV = STEP6_DIR / "nettcr2_test.csv"
OUT_PREDICTIONS = STEP6_DIR / "nettcr2_test_predictions.csv"
# Project root (parent of step6_NSCLC)
PROJECT_ROOT = STEP6_DIR.parent
DEFAULT_NETTCR2_REPO = PROJECT_ROOT / "NetTCR-2.0"
# Default training: neoantigen fine-tuning dataset (peptide, CDR3B, binder)
ALL_CSV = PROJECT_ROOT / "step5_Finetuning" / "data" / "all.csv"
DS_CSV = PROJECT_ROOT / "step1_random_negatives_generation" / "data" / "ds.csv"


def prepare_test_file() -> Path:
    """Ensure test file has CDR3b and peptide columns (NetTCR-2.0 format)."""
    df = pd.read_csv(TEST_CSV)
    if "best_cdr3b" in df.columns and "CDR3b" not in df.columns:
        df = df.rename(columns={"best_cdr3b": "CDR3b"})
    # NetTCR-2.0 chain b only needs CDR3b and peptide
    required = ["CDR3b", "peptide"]
    if not all(c in df.columns for c in required):
        raise ValueError(f"Test file needs columns {required}. Found: {list(df.columns)}")
    out = STEP6_DIR / "nettcr2_test_formatted.csv"
    df[["CDR3b", "peptide"]].to_csv(out, index=False)
    return out


def prepare_all_csv_for_nettcr() -> Path:
    """Prepare all.csv (peptide, CDR3B, binder) for NetTCR-2.0 (expects CDR3b)."""
    df = pd.read_csv(ALL_CSV)
    df = df.rename(columns={"CDR3B": "CDR3b"})
    df = df[["peptide", "CDR3b", "binder"]].dropna()
    df["binder"] = df["binder"].astype(int)
    out = STEP6_DIR / "nettcr2_train_from_all.csv"
    df.to_csv(out, index=False)
    return out


def create_train_from_ds(n_rows: int = 10000) -> Path:
    """Create training file from ds.csv for NetTCR-2.0 chain b."""
    df = pd.read_csv(DS_CSV, low_memory=False, nrows=n_rows * 2)
    df = df[["cdr3.beta", "antigen.epitope", "label"]].dropna()
    df = df.rename(columns={"cdr3.beta": "CDR3b", "antigen.epitope": "peptide", "label": "binder"})
    df["binder"] = df["binder"].astype(int)
    df = df.drop_duplicates(subset=["CDR3b", "peptide"]).head(n_rows)
    out = STEP6_DIR / "nettcr2_train_from_ds.csv"
    df.to_csv(out, index=False)
    return out


def main():
    parser = argparse.ArgumentParser(description="Run NetTCR-2.0 predictions on nettcr2_test.csv")
    parser.add_argument(
        "--nettcr2_repo",
        type=Path,
        default=DEFAULT_NETTCR2_REPO,
        help=f"Path to NetTCR-2.0 repo (default: {DEFAULT_NETTCR2_REPO})",
    )
    parser.add_argument(
        "--train_file",
        type=Path,
        default=None,
        help="Training CSV (CDR3b, peptide, binder). If not set, uses all.csv (step5 fine-tuning data).",
    )
    parser.add_argument(
        "--train_rows",
        type=int,
        default=15000,
        help="Rows to use when creating train from ds.csv (default: 15000)",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=50,
        help="Training epochs (default: 50)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUT_PREDICTIONS,
        help="Output predictions file",
    )
    args = parser.parse_args()

    nettcr_repo = Path(args.nettcr2_repo)
    nettcr_py = nettcr_repo / "nettcr.py"
    if not nettcr_py.exists():
        print(
            f"ERROR: NetTCR-2.0 not found at {nettcr_repo}.\n"
            "Clone it with: git clone https://github.com/mnielLab/NetTCR-2.0.git",
            file=sys.stderr,
        )
        sys.exit(1)

    print("Preparing test file...")
    test_file = prepare_test_file()
    print(f"  Test file: {test_file} ({len(pd.read_csv(test_file))} rows)")

    if args.train_file and args.train_file.exists():
        train_file = args.train_file
        print(f"  Using provided train file: {train_file}")
    elif ALL_CSV.exists():
        print(f"Preparing train file from {ALL_CSV} (neoantigen fine-tuning data)...")
        train_file = prepare_all_csv_for_nettcr()
        print(f"  Train file: {train_file} ({len(pd.read_csv(train_file))} rows)")
    elif DS_CSV.exists():
        print(f"Creating train file from {DS_CSV}...")
        train_file = create_train_from_ds(n_rows=args.train_rows)
        print(f"  Train file: {train_file} ({len(pd.read_csv(train_file))} rows)")
    else:
        # Fallback: use NetTCR-2.0's sample_train.csv (bundled with the repo)
        sample_train = nettcr_repo / "test" / "sample_train.csv"
        if not sample_train.exists():
            print(
                f"ERROR: No training data found. Tried:\n"
                f"  - all.csv at {ALL_CSV}\n"
                f"  - ds.csv at {DS_CSV}\n"
                f"  - NetTCR-2.0 sample at {sample_train}\n"
                f"Provide --train_file with CDR3b, peptide, binder columns.",
                file=sys.stderr,
            )
            sys.exit(1)
        train_file = sample_train
        print(f"  Using NetTCR-2.0 sample_train.csv ({len(pd.read_csv(train_file))} rows)")

    cmd = [
        sys.executable,
        str(nettcr_py),
        "-tr", str(train_file),
        "-te", str(test_file),
        "-c", "b",
        "-o", str(args.output),
        "-e", str(args.epochs),
    ]
    print(f"\nRunning: {' '.join(cmd)}")
    print("(NetTCR-2.0 will train on the training file, then predict on the test file.)\n")

    result = subprocess.run(cmd, cwd=str(nettcr_repo))
    if result.returncode != 0:
        sys.exit(result.returncode)

    print(f"\nPredictions saved to: {args.output}")


if __name__ == "__main__":
    main()
