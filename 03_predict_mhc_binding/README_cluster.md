# Running the NSCLC pipeline on a cluster

You can run the full pipeline (or just NetMHCpan step 3) on a cluster so the long NetMHCpan job doesn’t block your laptop.

## What to upload

### 1. Step 6 NSCLC directory

Upload the whole `step6_NSCLC` folder to the cluster (e.g. under `$HOME` or your project dir). It should contain at least:

- `step1_fetch_proteins_and_primary_neoantigens.py`
- `step2_generate_neoantigen_9mers.py`
- `step3_run_netmhcpan.py`
- `run_netmhcpan_cluster.slurm`
- Input data you already have:
  - `Cosmic_ResistanceMutations_v101_GRCh38.tsv` (for a full run from step 1), **or**
  - For **only step 3 (NetMHCpan)** you only need:
    - `neoantigen_9mers.list`
    - `neoantigen_9mers.csv`
    - (optional) `protein_sequences/` if you care about keeping paths consistent; step 3 only needs the list and CSV above.)

### 2. NetMHCpan 4.2 for Linux

On your **Mac** you have NetMHCpan under `step6_NSCLC_test/netMHCpan-4.2` with **Darwin_arm64** only.  
Clusters are usually **Linux x86_64**, so you need a **Linux build** of NetMHCpan 4.2.

**Important:** The Mac and Linux builds are **separate downloads** from DTU. The Mac tarball does **not** contain any Linux binary. You must request the Linux package in addition.

- Go to **https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/** and sign in (academic account).
- Request the **Linux** (or **Linux x86_64**) package. You will get a **different** tarball that contains a `Linux` or `Linux_x86_64` folder with the Linux binary.
- On the cluster you already have `~/netMHCpan-4.2/` with `data/` and `Darwin_arm64/`. After downloading the Linux tarball:
  1. Unpack it (on your machine or on the cluster).
  2. Copy **only** the `Linux_x86_64` or `Linux` folder from the unpacked archive into `~/netMHCpan-4.2/`. Do **not** overwrite `data/` (you keep the existing one).
  3. You should end up with: `~/netMHCpan-4.2/data/`, `~/netMHCpan-4.2/Darwin_arm64/`, and `~/netMHCpan-4.2/Linux_x86_64/` (or `Linux/`).
- Direct link to request **Linux** build:  
  **https://services.healthtech.dtu.dk/cgi-bin/sw_request?packageversion=4.2b&platform=Linux&software=netMHCpan&version=4.2**

So on the cluster you should have something like:

```text
$HOME/
  step6_NSCLC/
    step3_run_netmhcpan.py
    neoantigen_9mers.list
    neoantigen_9mers.csv
    run_netmhcpan_cluster.slurm
    ...
  netMHCpan-4.2/
    data/
      synlist_cedar.bin
      MHC_pseudo.dat
      ...
    Linux_x86_64/
      bin/
        netMHCpan-4.2
```

(You can also put `netMHCpan-4.2` somewhere else and set `NETMHC_DIR` in the SLURM script; see below.)

### 3. Python

Step 3 needs **Python 3** and **pandas**. Use the cluster’s Python (e.g. `module load python/3.10` or the system `python3`) and install pandas if needed (`pip install pandas` in your env).

---

## Running on the cluster

### Option A: Only NetMHCpan (step 3) – you already have 9-mers

If you already have `neoantigen_9mers.list` and `neoantigen_9mers.csv` (e.g. from your Mac):

1. Upload `step6_NSCLC` (with those two files + `step3_run_netmhcpan.py`) and `netMHCpan-4.2` (Linux version) as above.
2. In the SLURM script, set:
   - `WORKDIR` to the directory that contains `step6_NSCLC` (e.g. `/home/farmahini`).
   - `NETMHC_DIR` if NetMHCpan is not in `$WORKDIR/netMHCpan-4.2`.
3. Submit:

```bash
cd $HOME/step6_NSCLC
sbatch run_netmhcpan_cluster.slurm
```

Outputs will be in the same directory: `nsclc_netmhcpan_out.xls`, `nsclc_netmhcpan_strong_binders.csv`.

### Option B: Full pipeline (steps 1 → 2 → 3) on the cluster

1. Upload `step6_NSCLC` with:
   - `Cosmic_ResistanceMutations_v101_GRCh38.tsv`
   - `step1_*.py`, `step2_*.py`, `step3_*.py`
   - `run_netmhcpan_cluster.slurm`
2. Upload NetMHCpan 4.2 (Linux) as above.
3. Run step 1 (fetch proteins, primary neoantigens):

   ```bash
   cd $HOME/step6_NSCLC
   python3 step1_fetch_proteins_and_primary_neoantigens.py
   ```

4. Run step 2 (generate 9-mers):

   ```bash
   python3 step2_generate_neoantigen_9mers.py
   ```

5. Run step 3 via SLURM:

   ```bash
   sbatch run_netmhcpan_cluster.slurm
   ```

---

## SLURM script options

In `run_netmhcpan_cluster.slurm` you can set:

- **WORKDIR** – directory that contains `step6_NSCLC` (default: `/home/farmahini`).
- **NETMHC_DIR** – path to the NetMHCpan 4.2 directory (default: `$WORKDIR/netMHCpan-4.2`).

The script uses `--platform Linux_x86_64` so it will use the Linux binary. Time and memory in the script are set to 8 hours and 16 GB; adjust `#SBATCH --time` and `#SBATCH --mem` if your cluster limits differ.

---

## Summary

| Item | Where |
|------|--------|
| NetMHCpan in use | **NetMHCpan 4.2** (same version as on Mac) |
| Mac (current) | `step6_NSCLC_test/netMHCpan-4.2` → `Darwin_arm64/bin/netMHCpan-4.2` |
| Cluster | Upload `netMHCpan-4.2` with **Linux_x86_64** build; path set by `NETMHC_DIR` (e.g. `$HOME/netMHCpan-4.2`). |

Yes, you can upload the code and data and run it on the cluster; use the Linux NetMHCpan 4.2 binary and the provided SLURM script for step 3.

---

## "I only have the same files as on Mac (Darwin_arm64)"

The DTU NetMHCpan 4.2 download gives **one platform per request**. The Mac tarball contains only `Darwin_arm64`; it does **not** include any Linux binary. To run on the cluster you must:

1. Go to **https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/** and sign in with the same (academic) account you used for the Mac version.
2. Request the **Linux** package (it is a separate download). Direct link:  
   **https://services.healthtech.dtu.dk/cgi-bin/sw_request?packageversion=4.2b&platform=Linux&software=netMHCpan&version=4.2**
3. Download the Linux tarball, unpack it, and copy the **Linux** or **Linux_x86_64** folder into your existing `~/netMHCpan-4.2/` on the cluster (keep your current `data/` and `Darwin_arm64/`; you are only adding the new folder).
4. Re-run `sbatch run_netmhcpan_cluster.slurm`.

---

## "EL Threshold file ... thr.el does not exist" or "Cannot read EL Threshold file"

NetMHCpan 4.2 expects **real** `.thr.el` files in `data/threshold/` (~2300 bytes each, specific binary format). The cluster script **no longer overwrites** existing threshold files: if a real file is present (size 2000–3000 bytes), placeholder creation is skipped.

- **If you get "Cannot read"**: Our placeholder files (4004 bytes) have the wrong format and were overwriting real ones. **Restore** `data/threshold/` from your DTU Linux data package. If you have the data in `step6_NSCLC/data/` locally, copy **only** the threshold folder to the cluster so you don’t overwrite `synlist_cedar.bin` (the cluster Linux package usually has it):  
  `scp -r step6_NSCLC/data/threshold farmahini@<CLUSTER>:~/netMHCpan-4.2/data/`  
  Do not run `create_threshold_files.py` if you have the official threshold data.
- **If the directory is empty**: Run `create_threshold_files.py --binary` only as a last resort; the binary may still reject the format. Prefer restoring from the DTU package or copying from a working install.

---

## "nnalign_gaps_... not found" or "rdir/bin/ ... not found"

NetMHCpan expects helper binaries in **rdir/bin/** (e.g. `~/netMHCpan-4.2/bin/`), but the Linux tarball puts them in `Linux_x86_64/bin/`. The SLURM script now creates a symlink automatically: `netMHCpan-4.2/bin` → `netMHCpan-4.2/Linux_x86_64/bin`. If you run without the script, create it manually:

```bash
ln -s ~/netMHCpan-4.2/Linux_x86_64/bin ~/netMHCpan-4.2/bin
```

---

## If you still cannot get the Linux build

- **Option A – Ask cluster admins:** Some HPC centres install NetMHCpan as a module (e.g. `module load netmhcpan`). If yours does, set `NETMHC_DIR` to that install path and ensure it has a `Linux` or `Linux_x86_64` folder with `bin/netMHCpan-4.2`.
- **Option B – Run on your Mac in chunks:** Use the chunked script so each run is smaller and less likely to hang:
  ```bash
  cd /path/to/step6_NSCLC
  python3 step3_run_netmhcpan_chunked.py --chunk_size 200
  ```
  This runs NetMHCpan in batches of 200 peptides, then merges results into `nsclc_netmhcpan_out.xls` and `nsclc_netmhcpan_strong_binders.csv`. You can run it in the background or in a terminal you leave open. First kill any stuck run: `pkill -f step3_run_netmhcpan`.
