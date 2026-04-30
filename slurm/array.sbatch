#!/bin/bash
# array.sbatch
# Generic SLURM array launcher for the pipeline.
# Note: #SBATCH directives are passed dynamically via the sbatch command in 04_submit.sh.

# Exit immediately on error, undefined variables, or pipeline failures.
set -euo pipefail

# The following variables are exported via --export=ALL from 04_submit.sh:
# DATASET, TOOL, FASTA_LIST, PROJECT_ROOT, WORKER_SCRIPT, etc.

# 1. Map the Slurm Array ID to a file in the list.
# SLURM_ARRAY_TASK_ID goes from 0 to N-1 (as defined in 04_submit.sh).
# sed needs a 1-based line number.
LINE_NUM=$(( SLURM_ARRAY_TASK_ID + 1 ))

# 2. Extract the specific FASTA file path for this array task.
FASTA_FILE=$(sed -n "${LINE_NUM}p" "$FASTA_LIST")

if [[ ! -s "$FASTA_FILE" ]]; then
    echo "ERROR: Could not read FASTA file at line $LINE_NUM in $FASTA_LIST" >&2
    exit 1
fi

echo "=========================================================="
echo "Host        : $(hostname)"
echo "Task ID     : $SLURM_ARRAY_TASK_ID"
echo "Dataset     : $DATASET"
echo "Tool        : $TOOL"
echo "Input       : $FASTA_FILE"
echo "Worker      : $WORKER_SCRIPT"
echo "=========================================================="

# 3. Execute the worker script.
# We use 'exec' so the worker process replaces this bash script,
# cleanly receiving all signals from SLURM. 
# lib/paths.sh will be sourced inside the worker to set up specific dirs.
exec "$WORKER_SCRIPT" "$FASTA_FILE"