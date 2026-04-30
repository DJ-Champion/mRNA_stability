#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=slurm_logs/rnafold_%A_%a.out
#SBATCH --error=slurm_logs/rnafold_%A_%a.err

# Resource directives (--time, --mem-per-cpu, --cpus-per-task, --partition,
# --job-name, --array) are passed on the sbatch CLI from 04_submit.sh and
# override any #SBATCH lines here.

set -euo pipefail
source ./config.sh

mkdir -p "$SLURM_LOG_DIR"

# FASTA_LIST is exported by 04_submit.sh; pick the line for this array index.
[[ -s "${FASTA_LIST:-}" ]] || { echo "ERROR: FASTA_LIST not set or empty" >&2; exit 1; }

line=$((SLURM_ARRAY_TASK_ID + 1))
fasta=$(sed -n "${line}{p;q;}" "$FASTA_LIST")
[[ -n "$fasta" ]] || { echo "ERROR: no path on line $line of $FASTA_LIST" >&2; exit 1; }

echo "JobID=$SLURM_JOB_ID ArrayTaskID=$SLURM_ARRAY_TASK_ID Fasta=$fasta N_SHUFFLES=$N_SHUFFLES"

srun "$WORKER_SCRIPT" "$fasta"
