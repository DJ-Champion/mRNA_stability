#!/bin/bash
# slurm/array.sbatch
# Generic SLURM array launcher for the pipeline.
# All #SBATCH directives are passed dynamically via the sbatch command in
# 04_submit.sh — this file deliberately has none.
#
# Required env (set by --export in 04_submit.sh):
#   DATASET, TOOL, FASTA_LIST, PROJECT_ROOT
#
# Optional per-tier overrides (also via --export, set by tier_params hook
# in configs/tools/<tool>.sh):
#   e.g. N_SHUFFLES for rnafold
#
# Everything else (WORKER_SCRIPT, tool binary paths, etc.) is derived
# locally by sourcing lib/paths.sh. This makes the launcher robust against
# partial env propagation: as long as DATASET, TOOL, and PROJECT_ROOT
# arrive, we can rebuild the full environment.

set -euo pipefail

# --- 1. Sanity-check the bootstrap env --------------------------------------

: "${PROJECT_ROOT:?PROJECT_ROOT not in env (check --export in 04_submit.sh)}"
: "${DATASET:?DATASET not in env}"
: "${TOOL:?TOOL not in env}"
: "${FASTA_LIST:?FASTA_LIST not in env}"
: "${SLURM_ARRAY_TASK_ID:?Not running inside a SLURM array job}"

# --- 2. Self-sufficient path resolution -------------------------------------
# Source the same machinery the workers use. This re-derives WORKER_SCRIPT,
# RESULTS_DIR, ERRORS_DIR, etc., regardless of what propagated through
# --export=ALL.

source "$PROJECT_ROOT/lib/paths.sh"
parse_pipeline_args     # no CLI args; just initialises globals
resolve_paths

# --- 3. Map array task -> FASTA path ----------------------------------------

LINE_NUM=$(( SLURM_ARRAY_TASK_ID + 1 ))
FASTA_FILE=$(sed -n "${LINE_NUM}p" "$FASTA_LIST")

if [[ -z "$FASTA_FILE" ]]; then
    echo "ERROR: line $LINE_NUM is empty in $FASTA_LIST" >&2
    exit 1
fi
if [[ ! -s "$FASTA_FILE" ]]; then
    echo "ERROR: FASTA file missing or empty: $FASTA_FILE" >&2
    exit 1
fi

# --- 4. Diagnostic banner (lands in <tool>_<jobid>_<task>.out) --------------

echo "=========================================================="
echo "Job ID      : ${SLURM_ARRAY_JOB_ID:-?}_${SLURM_ARRAY_TASK_ID}"
echo "Host        : $(hostname)"
echo "Started     : $(date -Iseconds)"
echo "CPUs alloc  : ${SLURM_CPUS_PER_TASK:-1}"
echo "Mem alloc   : ${SLURM_MEM_PER_CPU:-?}MB/cpu"
echo "Dataset     : $DATASET"
echo "Tool        : $TOOL"
echo "Worker      : $WORKER_SCRIPT"
echo "FASTA list  : $FASTA_LIST  (line $LINE_NUM)"
echo "Input       : $FASTA_FILE"
[[ -n "${N_SHUFFLES:-}" ]] && echo "N_SHUFFLES  : $N_SHUFFLES"
echo "=========================================================="

# --- 5. Hand off to the tool's worker ---------------------------------------
# exec replaces this shell so SLURM signals (e.g. timeout, cancel) reach
# the worker directly. Worker re-sources paths.sh itself, so it remains
# usable as a standalone script too.

exec "$WORKER_SCRIPT" "$FASTA_FILE"