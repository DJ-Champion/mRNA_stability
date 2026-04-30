#!/bin/bash
# find_missing.sh
# Print paths to FASTAs whose tool output is missing or incomplete.
# Pipe to a file and feed back through 04_submit.sh for resubmission.
#
# Usage:
#   ./bin/find_missing.sh -d human_liver -t rnafold                       # check all tiers
#   ./bin/find_missing.sh -d human_liver -t rnafold --list lists/tier_5.txt  # one list
#
# Resume workflow:
#   ./bin/find_missing.sh -d human_liver -t rnafold > redo.txt
#   ./bin/04_submit.sh -d human_liver -t rnafold --list redo.txt

set -euo pipefail
source "$(dirname "$0")/../lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

[[ -n "$TOOL" ]] || { echo "ERROR: --tool required" >&2; exit 1; }

check_one_list() {
    local list="$1"
    [[ -s "$list" ]] || return 0
    while IFS= read -r fasta; do
        local seq_name out
        seq_name=$(basename "$fasta" .fa)
        out="$RESULTS_DIR/results_${seq_name}.csv"
        if [[ ! -s "$out" ]] || (( $(awk 'END{print NR}' "$out") < 2 )); then
            echo "$fasta"
        fi
    done < "$list"
}

if [[ -n "$EXPLICIT_LIST" ]]; then
    check_one_list "$EXPLICIT_LIST"
else
    for i in $(seq 1 10); do
        check_one_list "$LISTS_DIR/tier_${i}.txt"
    done
fi