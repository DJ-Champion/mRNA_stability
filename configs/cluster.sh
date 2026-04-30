# configs/cluster.sh
# Cluster / scheduling concerns only. Tool-specific settings live in
# configs/tools/<tool>.sh; biological inputs live in configs/datasets/<name>.yaml.
#
# Every value is overridable via env on the command line:
#   SLURM_PARTITION=other ./bin/04_submit.sh -d human_liver -t rnafold

# --- Output root ---
# Override RUNS_ROOT to keep run outputs off the repo (e.g. on scratch FS):
#   RUNS_ROOT=/scratch/me/mrna_runs ./bin/...
# Default: $PROJECT_ROOT/runs (set by lib/paths.sh).

# --- Length tier definitions ---
# 9 boundaries -> 10 tiers. Geometric (~1.8x ratio).
# Shared across all tools — re-tiering means re-running 02_stratify.sh.
: "${TIER_BOUNDS:=80 150 280 500 900 1600 2800 5000 10000}"

# --- SLURM defaults ---
# Used as fallbacks when no calibration recommendation is available
# (e.g. --list mode in 04_submit.sh). Calibrated submissions ignore these.
: "${SLURM_JOB_NAME:=mrna}"
: "${SLURM_PARTITION:=aoraki}"
: "${SLURM_TIME:=03:00:00}"
: "${SLURM_MEM_PER_CPU:=1000MB}"
: "${SLURM_CPUS_PER_TASK:=1}"
: "${CONCURRENT_LIMIT:=10}"

# --- Calibration parameters (tool-agnostic) ---
: "${SAMPLES_PER_TIER:=3}"
: "${SAFETY_FACTOR_TIME:=2}"
: "${SAFETY_FACTOR_MEM:=2}"

export TIER_BOUNDS \
       SLURM_JOB_NAME SLURM_PARTITION SLURM_TIME SLURM_MEM_PER_CPU \
       SLURM_CPUS_PER_TASK CONCURRENT_LIMIT \
       SAMPLES_PER_TIER SAFETY_FACTOR_TIME SAFETY_FACTOR_MEM