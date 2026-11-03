#!/bin/bash
# tools/rnafold_structure/worker.sh
# Fold one 2-line FASTA with RNAfold and keep the structure.
#
#   results/<seqname>.fold      RNAfold stdout: header, sequence, "structure (mfe)"
#                               (with WITH_BPP=1 also the ensemble lines)
#   results/<seqname>.bpp.tsv   WITH_BPP=1: i, j, prob  (0-based, i<j, prob>=BPP_MIN_PROB)
#
# Flags mirror tools/rnafold/worker.sh (--noPS, no -C, default parameters) so
# the MFE equals the stored one. -p only adds the ensemble; the MFE is unchanged.
#
# Required env (via lib/paths.sh): RESULTS_DIR, ERRORS_DIR, TMP_DIR, RNAFOLD_BIN
# Optional env: WITH_BPP (0|1), BPP_MIN_PROB (default 0.01)
#
# Usage: tools/rnafold_structure/worker.sh <fasta_file>

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

fasta_file="${1:?Usage: $0 <fasta_file>}"
[[ -r "$fasta_file" ]] || { echo "ERROR: cannot read $fasta_file" >&2; exit 1; }
[[ -x "$RNAFOLD_BIN" ]] || { echo "ERROR: RNAfold not executable: $RNAFOLD_BIN" >&2; exit 1; }

WITH_BPP="${WITH_BPP:-0}"
BPP_MIN_PROB="${BPP_MIN_PROB:-0.01}"

mkdir -p "$RESULTS_DIR" "$ERRORS_DIR" "$TMP_DIR"

seq_name=$(basename "$fasta_file" .fa)
fold_file="$RESULTS_DIR/${seq_name}.fold"
bpp_file="$RESULTS_DIR/${seq_name}.bpp.tsv"
error_file="$ERRORS_DIR/errors_${seq_name}.log"

# A .fold file counts as done when it has a structure line; the bpp file is
# also needed in WITH_BPP mode. Results from calibration (MFE only) therefore
# do not block a later WITH_BPP=1 run.
if [[ -s "$fold_file" ]] && (( $(awk 'END{print NR}' "$fold_file") >= 3 )); then
    if [[ "$WITH_BPP" != 1 || -e "$bpp_file" ]]; then
        echo "Skip $seq_name (done)"
        exit 0
    fi
fi

work=$(mktemp -d "$TMP_DIR/${seq_name}.XXXXXX")
trap 'rm -rf "$work"' EXIT

flags=(--noPS)
[[ "$WITH_BPP" == 1 ]] && flags+=(-p)

# RNAfold writes <header>_dp.ps into its cwd; run it in the scratch dir.
( cd "$work" && "$RNAFOLD_BIN" "${flags[@]}" < "$fasta_file" > fold.out 2>>"$error_file" )

# Structure line = third line, ends with "(energy)".
if ! awk 'NR==3 && /\( *[-+]?[0-9.]+\)$/ {ok=1} END{exit !ok}' "$work/fold.out"; then
    echo "ERROR: no structure line from RNAfold for $seq_name" >> "$error_file"
    exit 1
fi

if [[ "$WITH_BPP" == 1 ]]; then
    dp=$(ls "$work"/*_dp.ps 2>/dev/null | head -n 1 || true)
    [[ -n "$dp" ]] || { echo "ERROR: no dot plot from RNAfold -p for $seq_name" >> "$error_file"; exit 1; }
    # ubox lines: "i j sqrt(prob) ubox" (1-based). Lower triangle holds lbox (MFE).
    awk -v floor="$BPP_MIN_PROB" 'BEGIN { print "i\tj\tprob" }
        $4 == "ubox" { p = $3 * $3; if (p >= floor) printf "%d\t%d\t%.5f\n", $1 - 1, $2 - 1, p }' \
        "$dp" > "$work/bpp.tsv"
    mv "$work/bpp.tsv" "$bpp_file"
fi

mv "$work/fold.out" "$fold_file"
[[ -s "$error_file" ]] || rm -f "$error_file"
echo "Done $seq_name (MFE $(sed -n 3p "$fold_file" | awk '{print $NF}' | tr -d '()'))"
