#!/bin/bash
#SBATCH --job-name=score_peptide_cdr3
#SBATCH --output=score_peptide_cdr3_%j.out
#SBATCH --error=score_peptide_cdr3_%j.err
#SBATCH --partition=V100-32GB,A100-PCI,A100-40GB,A100-80GB,A100-RP,H100,H100-RP,H100-PCI,H200,H200-PCI,B200,RTX3090,RTXA6000,L40S,L40S-DSA
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1

set -euo pipefail

# Full run can exceed 1h depending on node/GPU load; keep a safe time budget.
WORKDIR="${WORKDIR:-/home/farmahini}"
STEP6="${WORKDIR}/step6_NSCLC"

cd "${STEP6}"

python3 score_peptide_cdr3.py \
  --model_py "${STEP6}/tuned_model/final_model_new_hpo_bce_logits.py" \
  --checkpoint "${STEP6}/tuned_model/finetune/best_model.pt" \
  --peptides strong_binders_peptide_wt_metrics.csv \
  --cdr3_list ds_cdr3b_unique.csv \
  --output strong_binders_best_cdr3b.csv \
  --batch_size 512 \
  --device cuda

echo ""
echo "=== Promiscuous CDR3β check ==="
python3 check_promiscuous_cdr3.py strong_binders_best_cdr3b.csv --top 15 || true
