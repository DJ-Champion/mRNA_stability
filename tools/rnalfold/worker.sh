#!/bin/bash
# tools/rnalfold/worker.sh
# Process one FASTA: minimum RNAlfold local MFE, shuffled null, p-value, z-score.
#
# Required env, normally set by lib/paths.sh after sourcing configs/tools/rnalfold.sh:
#   RESULTS_DIR, ERRORS_DIR, TMP_DIR, TOOL_DIR
#   RNALFOLD_BIN, ESL_SHUFFLE_BIN, RNALFOLD_WINDOW_LENGTH, MIN_VALID_PERCENT
#   N_SHUFFLES  (set by tier_params/calib_params; falls back to RNALFOLD_N_SHUFFLES)
#
# Output:
#   $RESULTS_DIR/results_<seqname>.csv
#   $TOOL_DIR/raw_shuffles/<seqname>_raw_shuffles.csv  when N_SHUFFLES > 0

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

fasta_file="${1:?Usage: $0 <fasta_file>}"
[[ -r "$fasta_file" ]] || { echo "ERROR: cannot read $fasta_file" >&2; exit 1; }

: "${N_SHUFFLES:=${RNALFOLD_N_SHUFFLES:-1000}}"
: "${RNALFOLD_WINDOW_LENGTH:=63}"
: "${MIN_VALID_PERCENT:=90}"

command -v "$RNALFOLD_BIN" >/dev/null 2>&1 || {
    echo "ERROR: RNALfold executable not found: $RNALFOLD_BIN" >&2; exit 1; }
if (( N_SHUFFLES > 0 )); then
    command -v "$ESL_SHUFFLE_BIN" >/dev/null 2>&1 || {
        echo "ERROR: esl-shuffle executable not found: $ESL_SHUFFLE_BIN" >&2; exit 1; }
fi

RAW_SHUFFLES_DIR="$TOOL_DIR/raw_shuffles"
mkdir -p "$RESULTS_DIR" "$ERRORS_DIR" "$TMP_DIR"
(( N_SHUFFLES > 0 )) && mkdir -p "$RAW_SHUFFLES_DIR"

seq_name=$(basename "$fasta_file" .fa)
output_file="$RESULTS_DIR/results_${seq_name}.csv"
error_file="$ERRORS_DIR/errors_${seq_name}.log"

# Idempotent skip.
if [[ -s "$output_file" ]] && (( $(awk 'END{print NR}' "$output_file") > 1 )); then
    echo "Skip $seq_name (done)"
    exit 0
fi

# Fail fast on unsupported input shapes. 02_stratify.sh should supply exactly
# one FASTA record per worker. UTR_pair-style 3-line constrained records are
# deliberately unsupported in the first RNAlfold pass and should be excluded
# via TOOL_REGIONS.
n_headers=$(grep -c '^>' "$fasta_file" || true)
if (( n_headers != 1 )); then
    echo "ERROR: RNAlfold worker expects exactly one FASTA record; found $n_headers in $fasta_file" >> "$error_file"
    exit 1
fi
if awk 'BEGIN{seen_seq=0} /^>/ {next} /^[.()<>|x]+$/ && seen_seq {found=1} NF {seen_seq=1} END{exit found ? 0 : 1}' "$fasta_file"; then
    echo "ERROR: constrained FASTA records are not supported by rnalfold first-pass worker: $fasta_file" >> "$error_file"
    exit 1
fi

# Extract the minimum local MFE from RNAlfold output.
# RNAlfold emits one local structure line per hit, each starting with dot-bracket
# notation, then echoes the sequence and a trailing summary line holding the summed
# energy of the whole non-overlapping decomposition. That summary is always at least
# as negative as any individual hit, so the match must be anchored to structure lines
# or it wins every comparison. We keep the most negative structure-line value.
extract_min_local_mfe() {
    awk '
        BEGIN { found = 0; lowest = 0 }
        /^[.()]+[[:space:]]+\(/ && match($0, /\([[:space:]]*[-+]?[0-9]*\.?[0-9]+[[:space:]]*\)/) {
            token = substr($0, RSTART, RLENGTH)
            gsub(/[()[:space:]]/, "", token)
            val = token + 0
            if (!found || val < lowest) {
                lowest = val
                found = 1
            }
        }
        END {
            if (found) printf "%.4f\n", lowest
            else exit 1
        }'
}

orig_min_local_mfe=$("$RNALFOLD_BIN" -i "$fasta_file" -L "$RNALFOLD_WINDOW_LENGTH" \
    2>>"$error_file" | extract_min_local_mfe) || {
        echo "ERROR: invalid original RNAlfold minimum local MFE for $seq_name" >> "$error_file"
        exit 1
    }

if [[ ! "$orig_min_local_mfe" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
    echo "ERROR: non-numeric original RNAlfold minimum local MFE for $seq_name: $orig_min_local_mfe" >> "$error_file"
    exit 1
fi

# Calibration fast path: original fold only.
if (( N_SHUFFLES == 0 )); then
    {
        echo "seq_name,min_local_mfe,median_shuffle_min_local_mfe,pvalue_min_local_mfe,zscore_min_local_mfe,n_valid,n_missing,n_shuffles_requested,rnalfold_window_length"
        printf "%s,%.4f,NA,NA,NA,0,0,0,%d\n" \
            "$seq_name" "$orig_min_local_mfe" "$RNALFOLD_WINDOW_LENGTH"
    } > "$output_file"
    [[ -s "$error_file" ]] || rm -f "$error_file"
    echo "Done $seq_name (RNAlfold original-only: $orig_min_local_mfe)"
    exit 0
fi

min_valid=$(( N_SHUFFLES * MIN_VALID_PERCENT / 100 ))
stats_file="$TMP_DIR/rnalfold_stats_${seq_name}.$$"
raw_tmp="$TMP_DIR/rnalfold_raw_${seq_name}.$$.csv"
raw_file="$RAW_SHUFFLES_DIR/${seq_name}_raw_shuffles.csv"

printf "seq_name,iteration,shuffle_min_local_mfe\n" > "$raw_tmp"

# Stream shuffled FASTA records through RNAlfold. For each shuffled record,
# keep the minimum local MFE within that record.
"$ESL_SHUFFLE_BIN" -d -N "$N_SHUFFLES" "$fasta_file" 2>>"$error_file" \
  | "$RNALFOLD_BIN" -L "$RNALFOLD_WINDOW_LENGTH" 2>>"$error_file" \
  | awk -v raw="$raw_tmp" -v s="$seq_name" '
        function flush_record() {
            if (seen_energy) {
                print lowest
                print s "," ++iter "," lowest >> raw
            }
            seen_record = 0
            seen_energy = 0
            lowest = 0
        }
        /^>/ {
            if (seen_record) flush_record()
            seen_record = 1
            next
        }
        # Anchored to structure lines; see extract_min_local_mfe for why the
        # trailing summary line must be excluded.
        /^[.()]+[[:space:]]+\(/ && match($0, /\([[:space:]]*[-+]?[0-9]*\.?[0-9]+[[:space:]]*\)/) {
            token = substr($0, RSTART, RLENGTH)
            gsub(/[()[:space:]]/, "", token)
            val = token + 0
            if (!seen_energy || val < lowest) {
                lowest = val
                seen_energy = 1
            }
        }
        END { if (seen_record) flush_record() }
    ' | sort -n > "$stats_file"

n_valid=$(awk 'END{print NR}' "$stats_file")
if (( n_valid < min_valid )); then
    echo "ERROR: insufficient valid RNAlfold shuffles for $seq_name ($n_valid/$N_SHUFFLES)" >> "$error_file"
    echo "Note: esl-shuffle -d can fail or produce unusable records for short/low-complexity sequences." >> "$error_file"
    rm -f "$stats_file" "$raw_tmp"
    exit 1
fi

mv "$raw_tmp" "$raw_file"

awk -v orig="$orig_min_local_mfe" \
    -v expected="$N_SHUFFLES" \
    -v s="$seq_name" \
    -v out="$output_file" \
    -v win="$RNALFOLD_WINDOW_LENGTH" '
function sort_numeric(a, n,    i, j, tmp) {
    for (i = 2; i <= n; i++) {
        tmp = a[i]
        j = i - 1
        while (j >= 1 && a[j] > tmp) {
            a[j + 1] = a[j]
            j--
        }
        a[j + 1] = tmp
    }
}
BEGIN { n = 0; more_stable = 0 }
{
    values[++n] = $1 + 0
    if (($1 + 0) <= orig) more_stable++
}
END {
    # Input is already sorted by sort -n, but sort defensively without gawk-only asort().
    sort_numeric(values, n)
    median = (n % 2) ? values[int(n/2)+1] : (values[n/2] + values[n/2+1]) / 2

    for (i = 1; i <= n; i++) {
        d = values[i] - median
        mad_vals[i] = (d < 0) ? -d : d
    }
    sort_numeric(mad_vals, n)
    mad = (n % 2) ? mad_vals[int(n/2)+1] : (mad_vals[n/2] + mad_vals[n/2+1]) / 2

    missing = expected - n
    pvalue = (orig <= median) \
        ? (more_stable + missing + 1) / (expected + 1) \
        : (more_stable + 1) / (expected + 1)
    zfmt = (mad == 0) ? "NA" : sprintf("%.4f", (orig - median) / (mad * 1.4826))

    print "seq_name,min_local_mfe,median_shuffle_min_local_mfe,pvalue_min_local_mfe,zscore_min_local_mfe,n_valid,n_missing,n_shuffles_requested,rnalfold_window_length" > out
    printf "%s,%.4f,%.4f,%.4f,%s,%d,%d,%d,%d\n", \
        s, orig, median, pvalue, zfmt, n, missing, expected, win >> out
}' "$stats_file"

rm -f "$stats_file"
[[ -s "$error_file" ]] || rm -f "$error_file"
echo "Done $seq_name ($n_valid/$N_SHUFFLES RNAlfold shuffles)"
