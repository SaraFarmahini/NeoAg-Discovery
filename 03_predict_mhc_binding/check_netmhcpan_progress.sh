#!/bin/bash
# Quick check: is the NetMHCpan (step3) run still going?
# Usage: ./check_netmhcpan_progress.sh   (run every few minutes to see progress)

STEP3_PATTERN="step3_run_netmhcpan"
TMP_DIR="/Users/sarafarmahinifarahani/Downloads/2026/TCR_peptide_transformer/step6_NSCLC_test/netMHCpan-4.2/tmp"
OUT_DIR="/Users/sarafarmahinifarahani/Downloads/2026/TCR_peptide_transformer/step6_NSCLC"

echo "=== NetMHCpan run status ==="
if pgrep -f "$STEP3_PATTERN" >/dev/null; then
  echo "Status:  RUNNING (step3 process found)"
  ps -o pid,etime,state -p $(pgrep -f "$STEP3_PATTERN" | head -1) 2>/dev/null || true
else
  echo "Status:  step3 process not running (finished or never started)"
fi

if ls "$OUT_DIR"/nsclc_netmhcpan_out.xls "$OUT_DIR"/nsclc_netmhcpan_strong_binders.csv 2>/dev/null; then
  echo "Output:  DONE - output files present"
  wc -l "$OUT_DIR"/nsclc_netmhcpan_strong_binders.csv 2>/dev/null
else
  echo "Output:  not yet (no nsclc_netmhcpan_*.xls / *.csv)"
fi

if [ -d "$TMP_DIR" ]; then
  echo "Tmp dir: $(du -sh "$TMP_DIR" 2>/dev/null | cut -f1) total"
  echo "0.dat sizes (bytes): $(find "$TMP_DIR" -name "0.dat" -exec wc -c {} \; 2>/dev/null | awk '{s+=$1} END {print s+0}')"
fi
echo "============================"
