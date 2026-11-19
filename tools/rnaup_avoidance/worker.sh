#!/bin/bash
# tools/rnaup_avoidance/worker.sh
# One window FASTA (a lists/tier_<n>/<seqname>.fa from 02_stratify.sh) ->
#   $RESULTS_DIR/results_<seqname>.csv
# Heavy lifting and column definitions: tools/rnaup_avoidance/avoidance.py.
#
# Calibration (AVOID_CALIB=1, set by calib_params in configs/tools/rnaup_avoidance.sh):
# runs the same work with a handful of shuffles and writes nothing under
# $TOOL_DIR, because lib/paths.sh points every run at the real results dir and a
# calibration result would otherwise be skipped as "done" by the real run.
#
# Usage: tools/rnaup_avoidance/worker.sh <fasta_file>

set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$(cd "$SRC_DIR/../.." && pwd)/lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

fasta_file="${1:?Usage: $0 <fasta_file>}"
[[ -r "$fasta_file" ]] || { echo "ERROR: cannot read $fasta_file" >&2; exit 1; }
command -v "$RNAUP_BIN" >/dev/null 2>&1 || { echo "ERROR: RNAup not found: $RNAUP_BIN" >&2; exit 1; }
command -v "$ESL_SHUFFLE_BIN" >/dev/null 2>&1 || { echo "ERROR: esl-shuffle not found: $ESL_SHUFFLE_BIN" >&2; exit 1; }
[[ -r "$NCRNA_FASTA" ]] || { echo "ERROR: ncRNA FASTA not found: $NCRNA_FASTA" >&2; exit 1; }

mkdir -p "$RESULTS_DIR" "$ERRORS_DIR" "$TMP_DIR"

seq_name=$(basename "$fasta_file" .fa)
output_file="$RESULTS_DIR/results_${seq_name}.csv"
error_file="$ERRORS_DIR/errors_${seq_name}.log"
calib="${AVOID_CALIB:-0}"

if [[ "$calib" != 1 && -s "$output_file" ]]; then
    echo "Skip $seq_name (done)"
    exit 0
fi

n_headers=$(grep -c '^>' "$fasta_file" || true)
if (( n_headers != 1 )); then
    echo "ERROR: expected exactly one FASTA record; found $n_headers in $fasta_file" >> "$error_file"
    exit 1
fi

# The manifest says which region (start_codon_region, 5UTR, ...) this window is.
region=$(awk -F'\t' -v s="$seq_name" '$1 == s { print $4; exit }' "$MANIFEST_TSV")
[[ -n "$region" ]] || { echo "ERROR: $seq_name not in $MANIFEST_TSV" >> "$error_file"; exit 1; }
export AVOID_REGION="$region"

# Job-local scratch: RNAup may drop files in its cwd despite --no_output_file.
scratch="${TMPDIR:-$TMP_DIR}"
mkdir -p "$scratch" 2>/dev/null || scratch="$TMP_DIR"
work=$(mktemp -d "$scratch/${seq_name}.XXXXXX")
trap 'rm -rf "$work"' EXIT
export TMPDIR="$work"

if [[ "$calib" == 1 ]]; then
    export RESULTS_DIR="$work/results" TOOL_DIR="$work" AVOID_CKPT_DIR="$work/ckpt"
else
    export AVOID_CKPT_DIR="$TOOL_DIR/ckpt/$seq_name"     # survives a restart
fi

# Progress goes to the SLURM log as it happens (tee); the errors file is only
# written on failure. avoidance.py logs to stderr and prints nothing on stdout.
if ! python3 "$SRC_DIR/avoidance.py" "$fasta_file" 2>&1 | tee "$work/log.txt" >&2; then
    cat "$work/log.txt" >> "$error_file"
    echo "ERROR: avoidance.py failed for $seq_name (see $error_file)" >&2
    exit 1
fi
echo "Done $seq_name"
