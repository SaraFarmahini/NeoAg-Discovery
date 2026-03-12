#!/usr/bin/env python3
"""
For each peptide in strong_binders_peptide_wt_metrics.csv, find the CDR3β from
ds_cdr3b_unique.csv with the highest binding score using the finetuned
blosum_boman model in step6 tuned_model.

Default model/checkpoint:
  tuned_model/final_model_new_hpo_bce_logits.py
  tuned_model/finetune/best_model.pt

Usage:
  python score_peptide_cdr3.py \\
    --peptides strong_binders_peptide_wt_metrics.csv \\
    --cdr3_list ds_cdr3b_unique.csv \\
    --output strong_binders_best_cdr3b.csv \\
    [--batch_size 512] [--device cuda]

Performance note:
  CDR3 encodings are precomputed once and reused for all peptides.
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MODEL_PY = SCRIPT_DIR / "tuned_model" / "final_model_new_hpo_bce_logits.py"
DEFAULT_CHECKPOINT = SCRIPT_DIR / "tuned_model" / "finetune" / "best_model.pt"


def load_model_and_encoder(model_py: Path, checkpoint_path: Path, device: torch.device):
    """Load TransformerModel (delta_feature_dim=2) and AAIndexEncoder from model_py; load weights from checkpoint."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("model_module", str(model_py))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)

    encoder = m.AAIndexEncoder()
    # Ablation blosum_boman uses only blosum + boman -> delta_feature_dim=2
    model = m.TransformerModel(delta_feature_dim=2).to(device)
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    if isinstance(ckpt, dict):
        if "model_state_dict" in ckpt:
            state = ckpt["model_state_dict"]
        elif "state_dict" in ckpt:
            state = ckpt["state_dict"]
        else:
            state = ckpt
    else:
        state = ckpt
    # Filter to matching keys (in case checkpoint has extra)
    model_state = model.state_dict()
    loaded = {k: v for k, v in state.items() if k in model_state and v.shape == model_state[k].shape}
    model.load_state_dict(loaded, strict=False)
    model.eval()
    return encoder, model


def run(
    model_py: Path,
    checkpoint: Path,
    peptides_csv: Path,
    cdr3_csv: Path,
    output_csv: Path,
    batch_size: int = 512,
    device_name: str = "cuda" if torch.cuda.is_available() else "cpu",
    cdr3_column: str = None,
    limit_peptides: int = 0,
    limit_cdr3: int = 0,
):
    device = torch.device(device_name)
    encoder, model = load_model_and_encoder(model_py, checkpoint, device)

    max_tcr_len = 50
    max_peptide_len = 30

    # Peptides: must have peptide, blosum_distance, boman_distance
    pep_df = pd.read_csv(peptides_csv)
    pep_df.columns = [c.strip() for c in pep_df.columns]
    for col in ["peptide", "blosum_distance", "boman_distance"]:
        if col not in pep_df.columns:
            raise ValueError(f"Peptides CSV must have column '{col}'. Got: {list(pep_df.columns)}")

    # Delta normalization (blosum, boman only) from peptide data
    delta_cols = ["blosum_distance", "boman_distance"]
    delta_mean = pep_df[delta_cols].mean().values.astype(np.float32)
    delta_std = pep_df[delta_cols].std().values.astype(np.float32)
    delta_std = np.where(delta_std < 1e-6, 1.0, delta_std)

    # CDR3 list
    cdr3_df = pd.read_csv(cdr3_csv)
    cdr3_df.columns = [c.strip() for c in cdr3_df.columns]
    if cdr3_column is None:
        cand = [c for c in cdr3_df.columns if "cdr3" in c.lower() and "beta" in c.lower()]
        cdr3_column = cand[0] if cand else cdr3_df.columns[0]
    if cdr3_column not in cdr3_df.columns:
        raise ValueError(f"CDR3 CSV must have column '{cdr3_column}'. Got: {list(cdr3_df.columns)}")
    cdr3_list = cdr3_df[cdr3_column].astype(str).str.strip().tolist()
    cdr3_list = [s for s in cdr3_list if s and s.lower() != "nan"]
    if limit_cdr3 > 0:
        cdr3_list = cdr3_list[: limit_cdr3]
    n_cdr3 = len(cdr3_list)
    if limit_peptides > 0:
        pep_df = pep_df.head(limit_peptides).copy()
    print(f"Loaded {len(pep_df)} peptide rows, {n_cdr3} CDR3β sequences. Batch size={batch_size}.")

    # Pre-encode CDR3 pool once (major speedup vs re-encoding for every peptide).
    print("Pre-encoding CDR3β pool...")
    all_tcr_enc = np.zeros((n_cdr3, max_tcr_len, 5), dtype=np.float32)
    all_tcr_mask = np.ones((n_cdr3, max_tcr_len), dtype=np.float32)
    for i, cdr3 in enumerate(cdr3_list):
        all_tcr_enc[i] = encoder.encode(cdr3, max_tcr_len)
        all_tcr_mask[i] = encoder.create_mask(cdr3, max_tcr_len)
        if (i + 1) % 50000 == 0 or (i + 1) == n_cdr3:
            print(f"  Encoded {i + 1} / {n_cdr3} CDR3β sequences.")

    results = []
    with torch.no_grad():
        for idx, row in pep_df.iterrows():
            peptide = str(row["peptide"]).strip()
            if not peptide:
                results.append({**row.to_dict(), "best_cdr3b": "", "binding_score": np.nan})
                continue

            pep_enc = encoder.encode(peptide, max_peptide_len)
            pep_mask = encoder.create_mask(peptide, max_peptide_len)
            delta = np.array([
                float(row["blosum_distance"]),
                float(row["boman_distance"]),
            ], dtype=np.float32)
            delta = (delta - delta_mean) / delta_std
            delta_t = torch.tensor(delta, device=device, dtype=torch.float32).unsqueeze(0)

            best_score = -float("inf")
            best_cdr3 = ""

            for start in range(0, n_cdr3, batch_size):
                batch_cdr3 = cdr3_list[start : start + batch_size]
                B = len(batch_cdr3)
                tcr_t = torch.from_numpy(all_tcr_enc[start : start + B]).to(device)
                tcr_mask_t = torch.from_numpy(all_tcr_mask[start : start + B]).to(device)
                pep_t = torch.from_numpy(pep_enc).unsqueeze(0).expand(B, -1, -1).to(device)
                pep_mask_t = torch.from_numpy(pep_mask).unsqueeze(0).expand(B, -1).to(device)
                delta_batch = delta_t.expand(B, -1)

                logits = model(tcr_t, pep_t, tcr_mask_t, pep_mask_t, delta_batch)
                scores = logits.cpu().numpy()

                for i, cdr3 in enumerate(batch_cdr3):
                    if scores[i] > best_score:
                        best_score = float(scores[i])
                        best_cdr3 = cdr3

            out_row = {**row.to_dict(), "best_cdr3b": best_cdr3, "binding_score": best_score}
            results.append(out_row)

            if (len(results)) % 50 == 0 or len(results) == len(pep_df):
                print(f"  Processed {len(results)} / {len(pep_df)} peptides.")

    out_df = pd.DataFrame(results)
    out_df.to_csv(output_csv, index=False)
    print(f"Saved {len(out_df)} rows to {output_csv}")
    return output_csv


def main():
    parser = argparse.ArgumentParser(description="Score peptide–CDR3 pairs and write best CDR3 per peptide.")
    parser.add_argument("--model_py", type=str, default=str(DEFAULT_MODEL_PY), help="Path to final_model_new_hpo_bce_logits.py")
    parser.add_argument("--checkpoint", type=str, default=str(DEFAULT_CHECKPOINT), help="Path to finetuned blosum_boman checkpoint")
    parser.add_argument("--peptides", type=str, default="strong_binders_peptide_wt_metrics.csv", help="CSV with peptide, blosum_distance, boman_distance")
    parser.add_argument("--cdr3_list", type=str, default="ds_cdr3b_unique.csv", help="CSV with CDR3β column")
    parser.add_argument("--output", type=str, default="strong_binders_best_cdr3b.csv", help="Output CSV")
    parser.add_argument("--batch_size", type=int, default=512, help="Batch size over CDR3s")
    parser.add_argument("--device", type=str, default="", help="Device: cuda or cpu (default: cuda if available)")
    parser.add_argument("--cdr3_column", type=str, default=None, help="CDR3 column name (default: auto-detect)")
    parser.add_argument("--limit_peptides", type=int, default=0, help="Use only first N peptide rows (0 = all, for testing)")
    parser.add_argument("--limit_cdr3", type=int, default=0, help="Use only first N CDR3s (0 = all, for testing)")
    args = parser.parse_args()

    base = SCRIPT_DIR
    model_py = Path(args.model_py)
    if not model_py.is_absolute():
        model_py = base / args.model_py
    checkpoint = Path(args.checkpoint)
    if not checkpoint.is_absolute():
        checkpoint = base / args.checkpoint
    # Paths from Snakemake are relative to cwd (project root); resolve so we don't double 04_filter_strong_binders
    peptides_csv = Path(args.peptides).resolve() if (os.path.sep in args.peptides or Path(args.peptides).is_absolute()) else base / args.peptides
    cdr3_csv = Path(args.cdr3_list).resolve() if (os.path.sep in args.cdr3_list or Path(args.cdr3_list).is_absolute()) else base / args.cdr3_list
    output_csv = Path(args.output).resolve() if (os.path.sep in args.output or Path(args.output).is_absolute()) else base / args.output

    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using model_py={model_py}")
    print(f"Using checkpoint={checkpoint}")
    run(
        model_py=model_py,
        checkpoint=checkpoint,
        peptides_csv=peptides_csv,
        cdr3_csv=cdr3_csv,
        output_csv=output_csv,
        batch_size=args.batch_size,
        device_name=device,
        cdr3_column=args.cdr3_column,
        limit_peptides=args.limit_peptides,
        limit_cdr3=args.limit_cdr3,
    )


if __name__ == "__main__":
    main()
