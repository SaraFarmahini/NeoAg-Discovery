#!/bin/bash
#SBATCH --job-name=score_nettcr2_pairs
#SBATCH --output=score_nettcr2_pairs_%j.out
#SBATCH --error=score_nettcr2_pairs_%j.err
#SBATCH --partition=V100-32GB,A100-PCI,A100-40GB,A100-80GB,A100-RP,H100,H100-RP,H100-PCI,H200,H200-PCI,B200,RTX3090,RTXA6000,L40S,L40S-DSA
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:1

set -euo pipefail

WORKDIR="${WORKDIR:-/home/farmahini}"
STEP6="${WORKDIR}/step6_NSCLC"

cd "${STEP6}"

python3 score_nettcr2_test_pairs.py \
  --model_py "${STEP6}/tuned_model/final_model_new_hpo_bce_logits.py" \
  --checkpoint "${STEP6}/tuned_model/finetune/best_model.pt" \
  --input nettcr2_test.csv \
  --output nettcr2_test_model_scores.csv \
  --peptide_metrics strong_binders_peptide_wt_metrics_unique.csv \
  --batch_size 64 \
  --device cuda

echo ""
echo "=== Done. Scores saved to nettcr2_test_model_scores.csv ==="
