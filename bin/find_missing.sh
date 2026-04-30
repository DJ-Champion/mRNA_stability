#!/bin/bash
# find_missing.sh
# Print paths to FASTAs whose RNAfold output is missing or incomplete.
# Pipe to a file and feed back through 04_submit.sh for resubmission.
#
# Usage:
#   ./bin/find_missing.sh                       # check all tier lists
#   ./bin/find_missing.sh lists/tier_5.txt      # check one specific list
#
# Example resume workflow:
#   ./bin/find_missing.sh > lists/redo.txt
#   FASTA_LIST=lists/redo.txt N_SHUFFLES=1000 \
#     SLURM_TIME=00:30:00 SLURM_MEM_PER_CPU=500MB ./bin/04_submit.sh

set -euo pipefail
source "$(dirname "$0")/../config/config.main.sh"

check_one_list() {
    local list="$1"
    [[ -s "$list" ]] || return 0

    while read -r fasta; do
        local seq_name
        seq_name=$(basename "$fasta" .fa)
        local out="$RESULTS_DIR/results_${seq_name}.csv"
        if [[ ! -s "$out" ]] || (( $(wc -l < "$out") < 2 )); then
            echo "$fasta"
        fi
    done < "$list"
}

if (( $# > 0 )); then
    for list in "$@"; do
        check_one_list "$list"
    done
else
    for i in $(seq 1 10); do
        check_one_list "$TIER_ROOT/tier_${i}.txt"
    done
fi
