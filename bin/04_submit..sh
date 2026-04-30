#!/bin/bash
# 04_submit.sh
# Submit SLURM array jobs. Two modes:
#
# (1) Submit ALL tiers using calibration recommendations:
#       ./bin/04_submit.sh
#
# (2) Submit a specific list (overrides), e.g. for resubmissions:
#       FASTA_LIST=lists/redo.txt N_SHUFFLES=1000 \
#         SLURM_TIME=00:30:00 SLURM_MEM_PER_CPU=500MB ./bin/04_submit.sh

set -euo pipefail
source "$(dirname "$0")/../config/config.main.sh"

mkdir -p "$SLURM_LOG_DIR"

submit_array() {
    local fasta_list="$1"
    local n_shuffles="$2"
    local slurm_time="$3"
    local slurm_mem="$4"
    local tag="${5:-}"

    [[ -s "$fasta_list" ]] || { echo "Empty/missing list: $fasta_list — skipping"; return; }

    local total
    total=$(wc -l < "$fasta_list")
    local last_idx=$((total - 1))
    local array_spec="0-${last_idx}%${CONCURRENT_LIMIT}"

    echo "Submitting: list=$fasta_list  n=$total  shuffles=$n_shuffles  time=$slurm_time  mem=$slurm_mem  ${tag:+[tag=$tag]}"

    sbatch \
        --job-name="${SLURM_JOB_NAME}${tag:+_$tag}" \
        --partition="$SLURM_PARTITION" \
        --time="$slurm_time" \
        --mem-per-cpu="$slurm_mem" \
        --cpus-per-task="$SLURM_CPUS_PER_TASK" \
        --array="$array_spec" \
        --export=ALL,FASTA_LIST="$fasta_list",N_SHUFFLES="$n_shuffles" \
        slurm/rnafold.sbatch
}

# --- Mode 1: explicit list via env ---
if [[ "${FASTA_LIST:-}" != "" && -s "${FASTA_LIST}" ]] \
   && [[ "${FASTA_LIST}" != "lists/tier_"*.txt || -n "${FORCE_LIST:-}" ]]; then
    # User specified an explicit list (resubmission, custom run)
    submit_array "$FASTA_LIST" "$N_SHUFFLES" "$SLURM_TIME" "$SLURM_MEM_PER_CPU" "custom"
    exit 0
fi

# --- Mode 2: submit all tiers from calibration ---
recs="${CALIBRATION_ROOT}/latest/recommendations.tsv"
if [[ ! -s "$recs" ]]; then
    echo "ERROR: no calibration recommendations at $recs" >&2
    echo "Run 03_calibrate.sh first, or submit explicitly with FASTA_LIST=...\n" >&2
    exit 1
fi

echo "Reading calibration recommendations from: $recs"
echo

# recommendations.tsv columns: tier  n_shuffles  n_samples  max_len  max_wall_s  max_rss_mb  rec_time  rec_mem
while IFS=$'\t' read -r tier n_shuf n_samples max_len max_wall max_rss rec_time rec_mem; do
    [[ "$tier" == "tier" ]] && continue   # header
    list="$TIER_ROOT/tier_${tier}.txt"
    submit_array "$list" "$n_shuf" "$rec_time" "$rec_mem" "t${tier}"
done < "$recs"

echo
echo "All tier submissions queued. Monitor: squeue -u $USER --name=${SLURM_JOB_NAME}*"
