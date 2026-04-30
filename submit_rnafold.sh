#!/bin/bash
set -euo pipefail
source ./config/config.main.sh
mkdir -p "$SLURM_LOG_DIR"

last_idx=$(( TOTAL_TASKS - 1 ))
echo "Submitting array 0-${last_idx}%${CONCURRENT_LIMIT} (N_SHUFFLES=$N_SHUFFLES, list=$FASTA_LIST)"

sbatch \
    --job-name="$SLURM_JOB_NAME" \
    --partition="$SLURM_PARTITION" \
    --time="$SLURM_TIME" \
    --mem-per-cpu="$SLURM_MEM_PER_CPU" \
    --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    --array="0-${last_idx}%${CONCURRENT_LIMIT}" \
    slurm/rnafold.sbatch