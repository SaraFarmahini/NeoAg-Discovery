#!/bin/bash
#SBATCH --job-name=nettcr2_predict
#SBATCH --output=nettcr2_predict_%j.out
#SBATCH --error=nettcr2_predict_%j.err
#SBATCH --partition=V100-32GB,A100-PCI,A100-40GB,A100-80GB,A100-RP,H100,H100-RP,H100-PCI,H200,H200-PCI,B200,RTX3090,RTXA6000,L40S,L40S-DSA
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --mail-user=sara.farmahini_farahani@dfki.de
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# Paths – adjust to your cluster layout
WORKDIR="${WORKDIR:-/home/farmahini}"
# Project root containing step6_NSCLC, step1_random_negatives_generation, etc.
PROJECT_ROOT="${PROJECT_ROOT:-${WORKDIR}}"
STEP6_DIR="${PROJECT_ROOT}/step6_NSCLC"
# NetTCR-2.0 repo (cloned from https://github.com/mnielLab/NetTCR-2.0)
NETTCR2_REPO="${NETTCR2_REPO:-${PROJECT_ROOT}/NetTCR-2.0}"

cd "${STEP6_DIR}"
echo "Working directory: $(pwd)"
echo "Project root: ${PROJECT_ROOT}"
echo "NetTCR-2.0 repo: ${NETTCR2_REPO}"

# Clone NetTCR-2.0 if not present
if [[ ! -f "${NETTCR2_REPO}/nettcr.py" ]]; then
  echo "Cloning NetTCR-2.0 into ${PROJECT_ROOT}..."
  mkdir -p "${PROJECT_ROOT}"
  (cd "${PROJECT_ROOT}" && git clone https://github.com/mnielLab/NetTCR-2.0.git) || true
  if [[ ! -f "${NETTCR2_REPO}/nettcr.py" ]]; then
    echo "ERROR: NetTCR-2.0 not found at ${NETTCR2_REPO}. Clone manually:"
    echo "  cd ${PROJECT_ROOT} && git clone https://github.com/mnielLab/NetTCR-2.0.git"
    exit 1
  fi
fi

# Activate venv if present, otherwise use system Python
if [[ -d "${PROJECT_ROOT}/venv" ]] && [[ -f "${PROJECT_ROOT}/venv/bin/activate" ]]; then
  echo "Activating venv..."
  source "${PROJECT_ROOT}/venv/bin/activate"
elif [[ -d "${STEP6_DIR}/venv" ]] && [[ -f "${STEP6_DIR}/venv/bin/activate" ]]; then
  source "${STEP6_DIR}/venv/bin/activate"
fi

# Install TensorFlow + deps if needed (NetTCR-2.0 uses Keras/TF)
pip install tensorflow "numpy<2.4" pandas scikit-learn --quiet 2>/dev/null || true

# Run prediction (trains on ds.csv subset, then predicts on nettcr2_test)
EPOCHS="${EPOCHS:-50}"
TRAIN_ROWS="${TRAIN_ROWS:-15000}"

echo ""
echo "=== Running NetTCR-2.0 prediction on nettcr2_test.csv ==="
python3 run_nettcr2_predict.py \
  --nettcr2_repo "${NETTCR2_REPO}" \
  --train_rows "${TRAIN_ROWS}" \
  --epochs "${EPOCHS}" \
  --output nettcr2_test_predictions.csv

echo ""
echo "=== Done. Predictions saved to nettcr2_test_predictions.csv ==="
