#!/bin/bash
# config.sh
# Cluster / scheduling concerns. Source this from any pipeline step.
# Every value is overridable via env:
#   N_SHUFFLES=10 FASTA_LIST=lists/tier_9.txt ./bin/04_submit.sh

# --- Tools (absolute paths — adjust once, here) ---
: "${RNAFOLD_BIN:=/home/chado47p/Software/ViennaRNA/bin/RNAfold}"
: "${ESL_SHUFFLE_BIN:=/home/chado47p/Software/hmmer3.4/bin/esl-shuffle}"

# Make tools visible to sub-shells without per-task conda activation:
rnafold_dir="$(dirname "$RNAFOLD_BIN")"
esl_shuffle_dir="$(dirname "$ESL_SHUFFLE_BIN")"
export PATH="$rnafold_dir:$esl_shuffle_dir:$PATH"

# --- Inputs / outputs ---
: "${MULTIFASTA_DIR:=extracted_regions}"
: "${MANIFEST_TSV:=${MULTIFASTA_DIR}/manifest.tsv}"
: "${TIER_ROOT:=lists}"
: "${RESULTS_DIR:=results/rnafold}"
: "${RAW_SHUFFLES_DIR:=results/rnafold_raw_shuffles}"
: "${ERRORS_DIR:=errors/rnafold}"
: "${TMP_DIR:=tmp/rnafold}"
: "${SLURM_LOG_DIR:=slurm_logs}"
: "${CALIBRATION_ROOT:=calibration}"

# --- Tier definitions ---
# 9 boundaries → 10 tiers. Geometric (~1.8x ratio).
: "${TIER_BOUNDS:=80 150 280 500 900 1600 2800 5000 10000}"

# Per-tier N_SHUFFLES policy (tier 1..10). 0 = MFE-only, no shuffles.
# Long sequences get progressively fewer shuffles to keep wallclock tractable.
: "${TIER_SHUFFLES:=1000 1000 1000 1000 1000 1000 1000 100 10 0}"

# --- Shuffle parameters ---
: "${N_SHUFFLES:=1000}"             # default; per-tier values come from TIER_SHUFFLES
: "${MIN_VALID_PERCENT:=90}"

# --- SLURM resources ---
: "${SLURM_JOB_NAME:=rnafold}"
: "${SLURM_PARTITION:=aoraki}"
: "${SLURM_TIME:=03:00:00}"
: "${SLURM_MEM_PER_CPU:=1000MB}"
: "${SLURM_CPUS_PER_TASK:=1}"
: "${CONCURRENT_LIMIT:=10}"

# --- Worker ---
: "${WORKER_SCRIPT:=bin/rnafold_worker.sh}"

# --- Calibration ---
: "${SAMPLES_PER_TIER:=3}"
: "${SAFETY_FACTOR_TIME:=2}"
: "${SAFETY_FACTOR_MEM:=2}"

export RNAFOLD_BIN ESL_SHUFFLE_BIN \
       MULTIFASTA_DIR MANIFEST_TSV TIER_ROOT \
       RESULTS_DIR RAW_SHUFFLES_DIR ERRORS_DIR TMP_DIR SLURM_LOG_DIR CALIBRATION_ROOT \
       TIER_BOUNDS TIER_SHUFFLES \
       N_SHUFFLES MIN_VALID_PERCENT \
       SLURM_JOB_NAME SLURM_PARTITION SLURM_TIME SLURM_MEM_PER_CPU SLURM_CPUS_PER_TASK CONCURRENT_LIMIT \
       WORKER_SCRIPT \
       SAMPLES_PER_TIER SAFETY_FACTOR_TIME SAFETY_FACTOR_MEM
