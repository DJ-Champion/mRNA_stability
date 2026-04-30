# mRNA_stability
Code for experiments and interpretation on mRNA stability

# RNAfold Pipeline

A staged pipeline for extracting transcript regions, computing RNAfold MFE / shuffle-based z-scores, and preparing results for downstream half-life analysis.

## Pipeline shape

```
01_extract.py        → extracted_regions/{extracted_<region>.fa, manifest.tsv, run_manifest.yaml}
02_stratify.sh       → lists/{tier_<n>/, tier_<n>.txt, .lengths.tsv}
03_calibrate.sh      → calibration/<timestamp>/recommendations.tsv
04_submit.sh         → SLURM array per tier with calibrated resources
05_collate.py        → results/combined_<region>.tsv (one table per metric source)
```

Each step's output is the next step's input. Any step is regenerable from the previous step's output.

## Quick start

```bash
# 1. Extract regions from a GFF/FASTA pair (creates 7 multifastas + manifest.tsv)
./01_extract.py config.yaml

# 2. Stratify into 10 length tiers per region (reads manifest.tsv, no re-scanning)
./bin/02_stratify.sh

# 3. Calibrate resources (run on a compute node)
srun --pty -p aoraki --time=01:00:00 -c 1 --mem=4G bash -c './bin/03_calibrate.sh'

# 4. Submit each tier (one SLURM array per tier)
./bin/04_submit.sh           # submits all tiers using calibration recommendations
# or manually:
N_SHUFFLES=1000 FASTA_LIST=lists/tier_4.txt SLURM_TIME=00:05:00 SLURM_MEM_PER_CPU=300MB ./bin/04_submit.sh

# 5. Find any sequences that didn't complete and resubmit
./bin/find_missing.sh > lists/redo.txt
FASTA_LIST=lists/redo.txt ./bin/04_submit.sh
```

## Configuration layers

- **`config.yaml`** — biological inputs only (genome, annotation, gene list, regions to extract)
- **`config.main.sh`** — cluster/scheduling concerns only (paths to binaries, SLURM resources, tier policy)

Either can be edited without touching the other.

## Files

| Path | Purpose |
|---|---|
| `config/config.yaml` | Extraction config (genome, GFF, gene list, regions) |
| `config/config.main.sh` | Cluster config (SLURM, tool paths, tier bounds) |
| `bin/01_extract.py` | Extract regions from GFF + genome → multifastas + manifest |
| `bin/02_stratify.sh` | Bin sequences into 10 length tiers (uses manifest) |
| `bin/03_calibrate.sh` | Time representative sequences, recommend SLURM resources |
| `bin/04_submit.sh` | Submit SLURM array job(s) for one or all tiers |
| `bin/rnafold_worker.sh` | Per-sequence worker (called by SLURM tasks) |
| `bin/find_missing.sh` | Identify sequences whose output is missing/incomplete |
| `slurm/rnafold.sbatch` | SLURM array job template |

## Notes

- **Idempotency**: every step skips work that's already complete. Re-running any step is safe.
- **Resume**: `find_missing.sh` produces a list of unfinished sequences; feed it back through `submit.sh`.
- **Tier policy**: `N_SHUFFLES` per tier is set in `config.sh`. Long sequences get fewer shuffles or MFE-only.
- **Re-tiering**: change `TIER_BOUNDS` in `config.sh` and re-run `02_stratify.sh` (seconds — uses manifest).
