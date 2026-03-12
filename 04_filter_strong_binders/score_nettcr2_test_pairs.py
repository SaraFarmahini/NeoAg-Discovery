#!/usr/bin/env python3
"""
Score (peptide, CDR3b) pairs from nettcr2_test.csv using the finetuned model.

Model: tuned_model/final_model_new_hpo_bce_logits.py
Checkpoint: tuned_model/finetune/best_model.pt
Trained on: step5_Finetuning/data/all.csv

The model requires blosum_distance and boman_distance (neoantigen–wildtype dissimilarity).
These are obtained by merging nettcr2_test peptides with strong_binders_peptide_wt_metrics_unique.csv.

Usage:
  python score_nettcr2_test_pairs.py [--input nettcr2_test.csv] [--output nettcr2_test_model_scores.csv]
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MODEL_PY = SCRIPT_DIR / "tuned_model" / "final_model_new_hpo_bce_logits.py"
DEFAULT_CHECKPOINT = SCRIPT_DIR / "tuned_model" / "finetune" / "best_model.pt"
DEFAULT_INPUT = SCRIPT_DIR / "nettcr2_test.csv"
DEFAULT_OUTPUT = SCRIPT_DIR / "nettcr2_test_model_scores.csv"
PEPTIDE_METRICS_CSV = SCRIPT_DIR / "strong_binders_peptide_wt_metrics_unique.csv"

MAX_TCR_LEN = 50
MAX_PEPTIDE_LEN = 30


def sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -500, 500)))


def load_model_and_encoder(model_py: Path, checkpoint_path: Path, device: torch.device):
    """Load TransformerModel (delta_feature_dim=2) and AAIndexEncoder; load weights from checkpoint."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("model_module", str(model_py))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)

    encoder = m.AAIndexEncoder()
    model = m.TransformerModel(delta_feature_dim=2).to(device)
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    if isinstance(ckpt, dict):
        state = ckpt.get("model_state_dict", ckpt.get("state_dict", ckpt))
    else:
        state = ckpt
    model_state = model.state_dict()
    loaded = {k: v for k, v in state.items() if k in model_state and v.shape == model_state[k].shape}
    model.load_state_dict(loaded, strict=False)
    model.eval()
    return encoder, model


def run(
    model_py: Path,
    checkpoint: Path,
    input_csv: Path,
    output_csv: Path,
    peptide_metrics_csv: Path,
    batch_size: int = 64,
    device_name: str = "cuda" if torch.cuda.is_available() else "cpu",
):
    device = torch.device(device_name)

    # Load pairs
    pairs = pd.read_csv(input_csv)
    pairs.columns = [c.strip() for c in pairs.columns]
    cdr3_col = "best_cdr3b" if "best_cdr3b" in pairs.columns else "CDR3b"
    if "peptide" not in pairs.columns or cdr3_col not in pairs.columns:
        raise ValueError(
            f"Input must have 'peptide' and '{cdr3_col}'. Found: {list(pairs.columns)}"
        )

    # Merge with peptide metrics for blosum_distance, boman_distance
    if not peptide_metrics_csv.exists():
        raise FileNotFoundError(
            f"Peptide metrics file not found: {peptide_metrics_csv}\n"
            "Required for delta features (blosum_distance, boman_distance)."
        )
    metrics = pd.read_csv(peptide_metrics_csv)
    metrics.columns = [c.strip() for c in metrics.columns]
    for col in ["peptide", "blosum_distance", "boman_distance"]:
        if col not in metrics.columns:
            raise ValueError(
                f"Peptide metrics must have '{col}'. Found: {list(metrics.columns)}"
            )

    df = pairs.merge(
        metrics[["peptide", "blosum_distance", "boman_distance"]],
        on="peptide",
        how="left",
    )
    missing = df["blosum_distance"].isna()
    if missing.any():
        n_miss = missing.sum()
        print(f"Warning: {n_miss} peptides missing delta features; using zeros.")
        df.loc[missing, "blosum_distance"] = 0.0
        df.loc[missing, "boman_distance"] = 0.0

    df["blosum_distance"] = df["blosum_distance"].astype(np.float32)
    df["boman_distance"] = df["boman_distance"].astype(np.float32)

    # Delta normalization (from peptide metrics)
    delta_mean = metrics[["blosum_distance", "boman_distance"]].mean().values.astype(np.float32)
    delta_std = metrics[["blosum_distance", "boman_distance"]].std().values.astype(np.float32)
    delta_std = np.where(delta_std < 1e-6, 1.0, delta_std)

    print(f"Loading model from {checkpoint}...")
    encoder, model = load_model_and_encoder(model_py, checkpoint, device)
    print(f"Scoring {len(df)} (peptide, CDR3b) pairs...")

    peptides = df["peptide"].astype(str).str.strip().tolist()
    cdr3s = df[cdr3_col].astype(str).str.strip().tolist()
    deltas = (
        (df[["blosum_distance", "boman_distance"]].values.astype(np.float32) - delta_mean)
        / delta_std
    )

    # Encode all sequences
    pep_enc = np.zeros((len(peptides), MAX_PEPTIDE_LEN, 5), dtype=np.float32)
    pep_mask = np.ones((len(peptides), MAX_PEPTIDE_LEN), dtype=np.float32)
    tcr_enc = np.zeros((len(cdr3s), MAX_TCR_LEN, 5), dtype=np.float32)
    tcr_mask = np.ones((len(cdr3s), MAX_TCR_LEN), dtype=np.float32)

    for i, (p, c) in enumerate(zip(peptides, cdr3s)):
        pep_enc[i] = encoder.encode(p, MAX_PEPTIDE_LEN)
        pep_mask[i] = encoder.create_mask(p, MAX_PEPTIDE_LEN)
        tcr_enc[i] = encoder.encode(c, MAX_TCR_LEN)
        tcr_mask[i] = encoder.create_mask(c, MAX_TCR_LEN)

    # Score in batches
    logits = []
    with torch.no_grad():
        for start in range(0, len(df), batch_size):
            end = min(start + batch_size, len(df))
            B = end - start
            tcr_t = torch.from_numpy(tcr_enc[start:end]).to(device)
            tcr_m = torch.from_numpy(tcr_mask[start:end]).to(device)
            pep_t = torch.from_numpy(pep_enc[start:end]).to(device)
            pep_m = torch.from_numpy(pep_mask[start:end]).to(device)
            delta_t = torch.from_numpy(deltas[start:end]).to(device)

            out = model(tcr_t, pep_t, tcr_m, pep_m, delta_t)
            logits.append(out.cpu().numpy().ravel())

    logits = np.concatenate(logits)
    probs = sigmoid(logits)

    out_df = df.copy()
    out_df["binding_score"] = logits
    out_df["binding_probability"] = probs

    out_df.to_csv(output_csv, index=False)
    print(f"Saved {len(out_df)} rows to {output_csv}")
    return output_csv


def main():
    parser = argparse.ArgumentParser(
        description="Score nettcr2_test.csv pairs with the TCR-peptide model"
    )
    parser.add_argument(
        "--model_py",
        type=Path,
        default=DEFAULT_MODEL_PY,
        help="Path to model definition",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=DEFAULT_CHECKPOINT,
        help="Path to checkpoint",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Input CSV with peptide and best_cdr3b",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output CSV with binding_score and binding_probability",
    )
    parser.add_argument(
        "--peptide_metrics",
        type=Path,
        default=PEPTIDE_METRICS_CSV,
        help="CSV with peptide, blosum_distance, boman_distance",
    )
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--device", type=str, default="")
    args = parser.parse_args()

    base = SCRIPT_DIR
    for p in ("model_py", "checkpoint", "input", "output", "peptide_metrics"):
        v = getattr(args, p)
        if v and not v.is_absolute():
            setattr(args, p, base / v)

    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    run(
        model_py=args.model_py,
        checkpoint=args.checkpoint,
        input_csv=args.input,
        output_csv=args.output,
        peptide_metrics_csv=args.peptide_metrics,
        batch_size=args.batch_size,
        device_name=device,
    )


if __name__ == "__main__":
    main()
