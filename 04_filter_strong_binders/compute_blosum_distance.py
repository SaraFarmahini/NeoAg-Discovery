#!/usr/bin/env python3
"""
Compute BLOSUM62 distance for neoantigen/wildtype peptide pairs.

Distance definition:
- For aligned positions, distance = (self_match_score - blosum_score)
- For length differences, apply a fixed gap penalty per extra residue
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent

# Try external sequence_encodings; otherwise use built-in BLOSUM62
try:
    FALLBACK_SAMPLE_DIRS = [
        SCRIPT_DIR / "../step2_model_selection/sample",
        SCRIPT_DIR / "../step3_HPO/HPO_scripts/sample",
        SCRIPT_DIR / "../step0_dataCollection/clustering/cluster_120/sample",
        SCRIPT_DIR / "../step0_dataCollection/clustering/cluster/sample",
        SCRIPT_DIR / "../../step2_model_selection/sample",
        SCRIPT_DIR / "../../step3_HPO/HPO_scripts/sample",
    ]
    for candidate in FALLBACK_SAMPLE_DIRS:
        if (candidate / "sequence_encodings.py").exists():
            sys.path.insert(0, str(candidate.resolve()))
            break
    from sequence_encodings import BLOSUM62, AA_TO_IDX
except ImportError:
    # Built-in: standard 20 amino acids and BLOSUM62 (order A R N D C Q E G H I L K M F P S T W Y V)
    AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
    AA_TO_IDX = {a: i for i, a in enumerate(AA_ORDER)}
    # BLOSUM62 matrix (20x20, row/col = AA_ORDER)
    _blosum62 = [
        [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
        [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
        [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
        [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
        [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
        [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
        [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
        [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
        [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
        [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
        [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
        [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
        [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
        [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
        [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
        [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
        [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
        [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
        [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
        [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
    ]
    BLOSUM62 = np.array(_blosum62, dtype=np.float64)


GAP_PENALTY = -4  # common default for BLOSUM62
BO_MAN_SI = {
    "A": 1.81,
    "C": 1.28,
    "D": 8.72,
    "E": 6.63,
    "F": -2.27,
    "G": 0.00,
    "H": 4.66,
    "I": -3.13,
    "K": 5.55,
    "L": -2.81,
    "M": -1.48,
    "N": 6.64,
    "P": 0.00,
    "Q": 5.54,
    "R": 14.92,
    "S": 3.40,
    "T": 2.57,
    "V": -1.69,
    "W": -2.09,
    "Y": 1.39,
}


def blosum_distance(seq_a: str, seq_b: str) -> float:
    seq_a = str(seq_a).upper().strip()
    seq_b = str(seq_b).upper().strip()

    min_len = min(len(seq_a), len(seq_b))
    distance = 0.0

    for i in range(min_len):
        aa = seq_a[i]
        bb = seq_b[i]
        if aa not in AA_TO_IDX or bb not in AA_TO_IDX:
            # Unknown residue: treat as gap
            distance += -GAP_PENALTY
            continue
        score = BLOSUM62[AA_TO_IDX[aa], AA_TO_IDX[bb]]
        self_score = BLOSUM62[AA_TO_IDX[aa], AA_TO_IDX[aa]]
        distance += float(self_score - score)

    # Gap penalty for extra residues
    if len(seq_a) != len(seq_b):
        gap_len = abs(len(seq_a) - len(seq_b))
        distance += gap_len * (-GAP_PENALTY)

    return float(distance)


def main():
    parser = argparse.ArgumentParser(description="Compute BLOSUM62 distance for neo/wt pairs")
    parser.add_argument("--input", type=str, required=True, help="Input CSV with neoantigen,wildtype columns")
    parser.add_argument("--output", type=str, required=True, help="Output CSV with blosum_distance column")
    args = parser.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)

    df = pd.read_csv(in_path)
    if {"neoantigen", "wildtype"}.issubset(df.columns):
        neo_col, wt_col = "neoantigen", "wildtype"
    elif {"peptide", "wt"}.issubset(df.columns):
        neo_col, wt_col = "peptide", "wt"
    else:
        raise ValueError("Input CSV must have columns: neoantigen,wildtype or peptide,wt")

    df["blosum_distance"] = [
        blosum_distance(n, w) for n, w in zip(df[neo_col], df[wt_col])
    ]

    # Aliphatic index distance (Kyte/Doolittle style coefficients)
    # AI = 100 * (X(Ala) + 2.9*X(Val) + 3.9*(X(Ile)+X(Leu)))
    def aliphatic_index(seq: str) -> float:
        seq = str(seq).upper().strip()
        if not seq:
            return 0.0
        counts = {aa: seq.count(aa) for aa in ["A", "V", "I", "L"]}
        length = len(seq)
        x_a = counts["A"] / length
        x_v = counts["V"] / length
        x_i = counts["I"] / length
        x_l = counts["L"] / length
        return 100.0 * (x_a + 2.9 * x_v + 3.9 * (x_i + x_l))

    df["aliphatic_index_distance"] = [
        abs(aliphatic_index(n) - aliphatic_index(w))
        for n, w in zip(df[neo_col], df[wt_col])
    ]

    # Boman potential distance (average solubility value per residue)
    def boman_index(seq: str) -> float:
        seq = str(seq).upper().strip()
        if not seq:
            return 0.0
        total = 0.0
        count = 0
        for aa in seq:
            if aa in BO_MAN_SI:
                total += BO_MAN_SI[aa]
                count += 1
        return total / count if count > 0 else 0.0

    df["boman_distance"] = [
        abs(boman_index(n) - boman_index(w))
        for n, w in zip(df[neo_col], df[wt_col])
    ]

    df.to_csv(out_path, index=False)
    print(f"Saved {len(df)} rows to: {out_path}")


if __name__ == "__main__":
    main()
