#!/bin/bash
# tools/rnaup_avoidance/worker.sh
# One tranche FASTA (M windows) -> $RESULTS_DIR/results_<tranche>.csv.
# Heavy lifting and column definitions: tools/rnaup_avoidance/avoidance.py.
#
# Usage: tools/rnaup_avoidance/worker.sh <tranche_fasta>

set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$(cd "$SRC_DIR/../.." && pwd)/lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

fasta_file="${1:?Usage: $0 <tranche_fasta>}"
[[ -r "$fasta_file" ]] || { echo "ERROR: cannot read $fasta_file" >&2; exit 1; }
command -v "$RNAUP_BIN" >/dev/null 2>&1 || { echo "ERROR: RNAup not found: $RNAUP_BIN" >&2; exit 1; }
command -v "$ESL_SHUFFLE_BIN" >/dev/null 2>&1 || { echo "ERROR: esl-shuffle not found: $ESL_SHUFFLE_BIN" >&2; exit 1; }
[[ -r "$NCRNA_FASTA" ]] || { echo "ERROR: ncRNA FASTA not found: $NCRNA_FASTA" >&2; exit 1; }

mkdir -p "$RESULTS_DIR" "$ERRORS_DIR" "$TMP_DIR"

tranche=$(basename "$fasta_file" .fa)
if [[ -s "$RESULTS_DIR/results_${tranche}.csv" ]]; then
    echo "Skip $tranche (done)"
    exit 0
fi
error_file="$ERRORS_DIR/errors_${tranche}.log"

# Job-local scratch: RNAup may drop files in its cwd despite --no_output_file.
work=$(mktemp -d "${TMPDIR:-$TMP_DIR}/${tranche}.XXXXXX")
trap 'rm -rf "$work"' EXIT
export TMPDIR="$work"

# Progress goes to the SLURM log; the errors file is only written on failure.
if ! python3 "$SRC_DIR/avoidance.py" "$fasta_file" 2>"$work/log.txt"; then
    cat "$work/log.txt" >> "$error_file"
    cat "$work/log.txt" >&2
    echo "ERROR: avoidance.py failed for $tranche (see $error_file)" >&2
    exit 1
fi
cat "$work/log.txt" >&2
echo "Done $tranche"
