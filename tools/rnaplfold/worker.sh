#!/bin/bash
# tools/rnaplfold/worker.sh
# RNAplfold on one full mRNA, summarised over its start region.
#
#   $RESULTS_DIR/results_<seqname>.csv
#
# Column definitions: tools/rnaplfold/summarize.py. The whole mRNA is folded so
# the start window sees real flanking context; only the summary is windowed.
#
# Required env (via lib/paths.sh + configs/tools/rnaplfold.sh):
#   RESULTS_DIR, ERRORS_DIR, TMP_DIR, MANIFEST_TSV, RNAPLFOLD_BIN,
#   RNAPLFOLD_WINDOW, RNAPLFOLD_SPAN, RNAPLFOLD_U_VALUES,
#   RNAPLFOLD_START_UPSTREAM, RNAPLFOLD_START_DOWNSTREAM
#
# Usage: tools/rnaplfold/worker.sh <fasta_file>

set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$(cd "$SRC_DIR/../.." && pwd)/lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

fasta_file="${1:?Usage: $0 <fasta_file>}"
[[ -r "$fasta_file" ]] || { echo "ERROR: cannot read $fasta_file" >&2; exit 1; }
command -v "$RNAPLFOLD_BIN" >/dev/null 2>&1 || {
    echo "ERROR: RNAplfold executable not found: $RNAPLFOLD_BIN" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found" >&2; exit 1; }

mkdir -p "$RESULTS_DIR" "$ERRORS_DIR" "$TMP_DIR"

seq_name=$(basename "$fasta_file" .fa)
output_file="$RESULTS_DIR/results_${seq_name}.csv"
error_file="$ERRORS_DIR/errors_${seq_name}.log"

if [[ -s "$output_file" ]] && (( $(awk 'END{print NR}' "$output_file") > 1 )); then
    echo "Skip $seq_name (done)"
    exit 0
fi

n_headers=$(grep -c '^>' "$fasta_file" || true)
if (( n_headers != 1 )); then
    echo "ERROR: RNAplfold worker expects exactly one FASTA record; found $n_headers in $fasta_file" >> "$error_file"
    exit 1
fi

seq_length=$(awk '!/^>/ { gsub(/[[:space:]]/, ""); n += length($0) } END { print n+0 }' "$fasta_file")
if (( seq_length <= 0 )); then
    echo "ERROR: empty sequence for $seq_name" >> "$error_file"
    exit 1
fi

# RNAplfold -u must cover the AUG (3 nt) and every reported stretch length.
u_max=$(printf '%s\n' ${RNAPLFOLD_U_VALUES//,/ } 3 | sort -n | tail -n 1)

work=$(mktemp -d "$TMP_DIR/${seq_name}.XXXXXX")
trap 'rm -rf "$work"' EXIT

# -L is always passed: -W alone silently sets the span to the window size.
# RNAplfold writes <header>_lunp and <header>_dp.ps into its cwd.
( cd "$work" && "$RNAPLFOLD_BIN" -W "$RNAPLFOLD_WINDOW" -L "$RNAPLFOLD_SPAN" -u "$u_max" \
    < "$fasta_file" > /dev/null 2>>"$error_file" ) || {
    echo "ERROR: RNAplfold failed for $seq_name" >> "$error_file"
    exit 1
}

lunp=$(ls "$work"/*_lunp 2>/dev/null | head -n 1 || true)
[[ -n "$lunp" ]] || { echo "ERROR: no _lunp file from RNAplfold for $seq_name" >> "$error_file"; exit 1; }

python3 "$SRC_DIR/summarize.py" \
    --lunp "$lunp" --manifest "$MANIFEST_TSV" --seq-name "$seq_name" \
    --seq-length "$seq_length" \
    --upstream "$RNAPLFOLD_START_UPSTREAM" --downstream "$RNAPLFOLD_START_DOWNSTREAM" \
    --u-values "$RNAPLFOLD_U_VALUES" \
    --window "$RNAPLFOLD_WINDOW" --span "$RNAPLFOLD_SPAN" \
    > "$work/result.csv" 2>>"$error_file" || {
    echo "ERROR: summarize.py failed for $seq_name" >> "$error_file"
    exit 1
}

mv "$work/result.csv" "$output_file"
[[ -s "$error_file" ]] || rm -f "$error_file"
echo "Done $seq_name"
