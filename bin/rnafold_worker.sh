#!/bin/bash
# rnafold_worker.sh
# Process a single FASTA: compute MFE, optionally generate shuffled MFEs and stats.
#
# Usage:  rnafold_worker.sh <fasta_file>
# Modes (set via N_SHUFFLES env var):
#   N_SHUFFLES=1000   full z-score pipeline
#   N_SHUFFLES=10     reduced shuffles (longer sequences)
#   N_SHUFFLES=0      MFE-only fast path (very long sequences)

set -euo pipefail

# Source config relative to this script
source "$(dirname "$0")/../config/config.main.sh"

fasta_file="${1:?Usage: $0 <fasta_file>}"
[[ -r "$fasta_file" ]] || { echo "ERROR: cannot read $fasta_file" >&2; exit 1; }

[[ -x "$RNAFOLD_BIN" ]] || { echo "ERROR: RNAfold not executable: $RNAFOLD_BIN" >&2; exit 1; }
if (( N_SHUFFLES > 0 )); then
    [[ -x "$ESL_SHUFFLE_BIN" ]] || {
        echo "ERROR: esl-shuffle not executable: $ESL_SHUFFLE_BIN" >&2; exit 1; }
fi

mkdir -p "$RESULTS_DIR" "$ERRORS_DIR" "$TMP_DIR"
(( N_SHUFFLES > 0 )) && mkdir -p "$RAW_SHUFFLES_DIR"

seq_name=$(basename "$fasta_file" .fa)
output_file="$RESULTS_DIR/results_${seq_name}.csv"
error_file="$ERRORS_DIR/errors_${seq_name}.log"

# --- Idempotent skip ---
if [[ -s "$output_file" && $(wc -l < "$output_file") -gt 1 ]]; then
    echo "Skip $seq_name (done)"
    exit 0
fi

# --- 1. Original MFE ---
extract_mfe() {
    awk 'match($0, /\(\s*([-+]?[0-9]*\.?[0-9]+)\s*\)$/, a) { print a[1]; exit }'
}
orig_mfe=$("$RNAFOLD_BIN" --noPS < "$fasta_file" 2>>"$error_file" | extract_mfe)

if [[ ! "$orig_mfe" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
    echo "ERROR: invalid orig MFE for $seq_name" >> "$error_file"
    exit 1
fi

# --- 2. MFE-only fast path ---
if (( N_SHUFFLES == 0 )); then
    {
        echo "seq_name,orig_mfe,median_mfe,pvalue,zscore,n_valid,n_missing,n_shuffles_requested"
        printf "%s,%.4f,NA,NA,NA,0,0,0\n" "$seq_name" "$orig_mfe"
    } > "$output_file"
    [[ -s "$error_file" ]] || rm -f "$error_file"
    echo "Done $seq_name (MFE-only: $orig_mfe)"
    exit 0
fi

# --- 3. Shuffle + stats ---
min_valid=$(( N_SHUFFLES * MIN_VALID_PERCENT / 100 ))
stats_file="$TMP_DIR/stats_${seq_name}.$$"
raw_file="$RAW_SHUFFLES_DIR/${seq_name}_raw_shuffles.csv"
raw_tmp="$TMP_DIR/raw_${seq_name}.$$"

echo "transcript_id,iteration,shuffled_mfe" > "$raw_tmp"

"$ESL_SHUFFLE_BIN" -d -N "$N_SHUFFLES" "$fasta_file" 2>>"$error_file" \
  | "$RNAFOLD_BIN" --noPS 2>>"$error_file" \
  | awk -v raw="$raw_tmp" -v s="$seq_name" '
      match($0, /\(\s*([-+]?[0-9]*\.?[0-9]+)\s*\)$/, a) {
          print a[1]
          print s "," ++iter "," a[1] >> raw
      }' \
  | sort -n > "$stats_file"

n_valid=$(wc -l < "$stats_file")
if (( n_valid < min_valid )); then
    echo "ERROR: insufficient shuffles for $seq_name ($n_valid/$N_SHUFFLES)" >> "$error_file"
    echo "Note: -d (dinucleotide) shuffle may fail to produce enough unique permutations for short sequences." >> "$error_file"
    rm -f "$stats_file" "$raw_tmp"
    exit 1
fi

mv "$raw_tmp" "$raw_file"

awk -v orig="$orig_mfe" -v expected="$N_SHUFFLES" -v s="$seq_name" -v out="$output_file" '
BEGIN { n = 0; more_stable = 0 }
{ values[++n] = $1; if ($1 <= orig) more_stable++ }
END {
    asort(values)
    median = (n % 2) ? values[int(n/2)+1] : (values[n/2] + values[n/2+1]) / 2
    for (i = 1; i <= n; i++) {
        d = values[i] - median
        mad_vals[i] = (d < 0) ? -d : d
    }
    asort(mad_vals)
    mad = (n % 2) ? mad_vals[int(n/2)+1] : (mad_vals[n/2] + mad_vals[n/2+1]) / 2
    missing = expected - n
    pvalue = (orig <= median) \
        ? (more_stable + missing + 1) / (expected + 1) \
        : (more_stable + 1) / (expected + 1)
    zfmt = (mad == 0) ? "NA" : sprintf("%.4f", (orig - median) / (mad * 1.4826))
    print "seq_name,orig_mfe,median_mfe,pvalue,zscore,n_valid,n_missing,n_shuffles_requested" > out
    printf "%s,%.4f,%.4f,%.4f,%s,%d,%d,%d\n", s, orig, median, pvalue, zfmt, n, missing, expected >> out
}' "$stats_file"

rm -f "$stats_file"
[[ -s "$error_file" ]] || rm -f "$error_file"
echo "Done $seq_name ($n_valid/$N_SHUFFLES shuffles)"
