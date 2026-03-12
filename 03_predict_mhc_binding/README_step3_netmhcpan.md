# Step 3: Run NetMHCpan on neoantigen 9-mers

Same setup as `step6_NSCLC_test/Thesis/test/run_netmhcpan_cosmic.py`.

## Prerequisites

- Step 2 done: `neoantigen_9mers.list` and `neoantigen_9mers.csv` in this folder.
- NetMHCpan 4.2: by default uses `../step6_NSCLC_test/netMHCpan-4.2` (Darwin_arm64).

## Run (full, ~1739 peptides × 9 alleles — can take 10–30+ minutes)

```bash
cd /path/to/step6_NSCLC
python3 step3_run_netmhcpan.py
```

Optional:

- `--limit 100` — run on first 100 peptides only (for testing).
- `--netmhcpan_dir /path/to/netMHCpan-4.2` — if NetMHCpan is elsewhere.
- `--platform Linux_x86_64` — if on Linux.

## Outputs

- **nsclc_netmhcpan_out.xls** — raw NetMHCpan output.
- **nsclc_netmhcpan_strong_binders.csv** — peptides with EL_Rank ≤ 0.5, merged with `neoantigen_9mers.csv` (metadata: gene, transcript, mutation, wt_peptide, etc.).

## Run in background (recommended for full list)

```bash
cd /path/to/step6_NSCLC
nohup python3 step3_run_netmhcpan.py > netmhcpan.log 2>&1 &
tail -f netmhcpan.log
```
