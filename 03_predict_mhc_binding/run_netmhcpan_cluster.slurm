#!/bin/bash
#SBATCH --job-name=nsclc_netmhcpan
#SBATCH --output=nsclc_netmhcpan_%j.out
#SBATCH --error=nsclc_netmhcpan_%j.err
#SBATCH --partition=V100-32GB,A100-PCI,A100-40GB,A100-80GB,A100-RP,H100,H100-RP,H100-PCI,H200,H200-PCI,B200,RTX3090,RTXA6000,L40S,L40S-DSA
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# No GPU needed for NetMHCpan

set -euo pipefail

# Paths – adjust to your cluster layout
WORKDIR="${WORKDIR:-/home/farmahini}"
STEP6_DIR="${WORKDIR}/step6_NSCLC"
# NetMHCpan 4.2 directory on cluster (must contain data/ and Linux_x86_64/bin/netMHCpan-4.2)
NETMHC_DIR="${NETMHC_DIR:-${WORKDIR}/netMHCpan-4.2}"

cd "${STEP6_DIR}"

# Find Linux NetMHCpan binary (DTU may use bin/ or put binary directly in platform folder)
NETMHC_BIN=""
NETMHC_PLATFORM=""
for platform in Linux_x86_64 Linux; do
  for subpath in "bin/netMHCpan-4.2" "netMHCpan-4.2"; do
    candidate="${NETMHC_DIR}/${platform}/${subpath}"
    if [[ -f "${candidate}" ]]; then
      chmod +x "${candidate}" 2>/dev/null || true
      if [[ -x "${candidate}" ]]; then
        NETMHC_BIN="${candidate}"
        NETMHC_PLATFORM="${platform}"
        break 2
      fi
    fi
  done
done
if [[ -z "${NETMHC_BIN}" ]]; then
  echo "Error: No Linux NetMHCpan binary found under ${NETMHC_DIR}"
  echo "  You have: $(ls -d ${NETMHC_DIR}/*/ 2>/dev/null | xargs -I{} basename {} | tr '\n' ' ')"
  echo "  Check: ls -la ${NETMHC_DIR}/Linux_x86_64/"
  echo "  Required: Linux_x86_64/bin/netMHCpan-4.2 or Linux_x86_64/netMHCpan-4.2 (and chmod +x)"
  exit 1
fi
echo "Using NetMHCpan binary: ${NETMHC_BIN}"

# NetMHCpan looks for helper binaries (nnalign_gaps_...) in rdir/bin/, but Linux tarball puts them in Linux_x86_64/bin/. Symlink so rdir/bin exists.
PLATFORM_BIN="${NETMHC_DIR}/${NETMHC_PLATFORM}/bin"
if [[ ! -d "${NETMHC_DIR}/bin" ]]; then
  if [[ -d "${PLATFORM_BIN}" ]]; then
    ln -s "${PLATFORM_BIN}" "${NETMHC_DIR}/bin"
    echo "Created ${NETMHC_DIR}/bin -> ${PLATFORM_BIN}"
  else
    echo "Error: ${PLATFORM_BIN} not found; NetMHCpan needs rdir/bin/ for helper binaries."
    exit 1
  fi
fi
# Ensure all binaries in bin/ are executable (fix "Permission denied" on helper programs)
chmod +x "${NETMHC_DIR}/bin/"* 2>/dev/null || true

# Binary looks for data/threshold/<pseudo>.thr.el. Do NOT overwrite real files (~2300 bytes).
export NETMHCpan="${NETMHC_DIR}"
export NMHOME="${NETMHC_DIR}"
mkdir -p "${NETMHC_DIR}/data/threshold"
REAL_THR="${NETMHC_DIR}/data/threshold/YFAMYQENMAHTDANTLYIIYRDYTWVARVYRGY.thr.el"
REAL_SIZE=0
[[ -f "${REAL_THR}" ]] && REAL_SIZE=$(stat -c%s "${REAL_THR}" 2>/dev/null || echo 0)
if [[ "${REAL_SIZE}" -gt 2000 && "${REAL_SIZE}" -lt 3000 ]]; then
  echo "Real threshold files present (~${REAL_SIZE} bytes); skipping placeholder creation."
else
  echo "Creating placeholder threshold files (do not overwrite if you have real DTU data/threshold)..."
  python3 "${STEP6_DIR}/create_threshold_files.py" --netmhc_dir "${NETMHC_DIR}" --binary
fi
if [[ ! -f "${NETMHC_DIR}/data/threshold/YYSEYRNIYAQTDESNLYLSYDYYTWAERAYEWY.thr.el" ]]; then
  echo "Error: threshold file missing. Restore data/threshold/ from DTU Linux package or run create_threshold_files.py."
  exit 1
fi
chmod 755 "${NETMHC_DIR}/data/threshold" 2>/dev/null || true

# Synlist: try synlist_cedar.bin, then synlist.bin, then synlist (no extension; from DTU data upload)
SYNLIST=""
for f in synlist_cedar.bin synlist.bin synlist; do
  if [[ -f "${NETMHC_DIR}/data/${f}" ]]; then
    SYNLIST="${f}"
    echo "Using synlist file: ${f}"
    break
  fi
done
if [[ -z "${SYNLIST}" ]]; then
  echo "Error: no synlist file found in ${NETMHC_DIR}/data/ (tried synlist_cedar.bin, synlist.bin, synlist)"
  exit 1
fi

# Python with pandas (use module load if your cluster has one)
# module load python/3.10  # uncomment and adjust if needed

# Mode 2 (Neoepitopes) requires CEDAR threshold files in data/threshold/ that may not be in the Linux data. Use mode 0 (Antigen presentation) to get EL/BA predictions without them.
python3 step3_run_netmhcpan.py \
  --netmhcpan_dir "${NETMHC_DIR}" \
  --platform "${NETMHC_PLATFORM}" \
  --netmhcpan_bin "${NETMHC_BIN}" \
  --synlist "${SYNLIST}" \
  --mode 0 \
  --peptide_list neoantigen_9mers.list \
  --neo_csv neoantigen_9mers.csv \
  --out_xls nsclc_netmhcpan_out.xls \
  --out_csv nsclc_netmhcpan_strong_binders.csv

echo "Done. Check nsclc_netmhcpan_out.xls and nsclc_netmhcpan_strong_binders.csv"
