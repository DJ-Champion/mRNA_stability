#!/bin/bash
# tools/rnafold/worker.sh
# Process a single FASTA: compute MFE, optionally generate shuffled MFEs and stats.
#
# Two input shapes are supported:
#
#   (a) Standard 2-line FASTA (header + sequence)
#       Streaming pipe esl-shuffle | RNAfold for shuffles. Fast.
#
#   (b) Constrained 3-line FASTA (header + sequence + ViennaRNA hard-
#       constraint string). Detected automatically from a third line composed
#       only of [.()<>|x] characters.
#
#       Constraints currently produced by 01_extract.py for the UTR_pair
#       region encode cross-pair semantics:
#           <<<<<...xxxxxxx...>>>>>
#       i.e. 5UTR positions only pair with 3UTR positions, linker forced
#       unpaired. For shuffles, 5UTR and 3UTR halves are shuffled
#       INDEPENDENTLY then rejoined with the original linker — pooling them
#       would destroy the assay. The constraint string is reused as-is for
#       every shuffle (linker geometry is fixed by the extraction config).
#
# Required env (set by lib/paths.sh -> resolve_paths):
#   RESULTS_DIR, ERRORS_DIR, TMP_DIR, TOOL_DIR
#   RNAFOLD_BIN, ESL_SHUFFLE_BIN, MIN_VALID_PERCENT
#
# Required env (set by 04_submit.sh / sbatch / calibrate caller):
#   N_SHUFFLES   - 0 = MFE-only fast path; >0 = shuffle pipeline
#   DATASET, TOOL - for paths.sh to resolve correctly
#
# Usage: tools/rnafold/worker.sh <fasta_file>

set -euo pipefail

# Resolve paths from DATASET + TOOL env vars (set upstream by sbatch / submit).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib/paths.sh"
parse_pipeline_args "$@"   # no flags expected here; ensures globals init
resolve_paths

fasta_file="${1:?Usage: $0 <fasta_file>}"
[[ -r "$fasta_file" ]] || { echo "ERROR: cannot read $fasta_file" >&2; exit 1; }

[[ -x "$RNAFOLD_BIN" ]] || { echo "ERROR: RNAfold not executable: $RNAFOLD_BIN" >&2; exit 1; }
if (( N_SHUFFLES > 0 )); then
    [[ -x "$ESL_SHUFFLE_BIN" ]] || {
        echo "ERROR: esl-shuffle not executable: $ESL_SHUFFLE_BIN" >&2; exit 1; }
fi

# Tool-specific output dir, derived from TOOL_DIR
RAW_SHUFFLES_DIR="$TOOL_DIR/raw_shuffles"

mkdir -p "$RESULTS_DIR" "$ERRORS_DIR" "$TMP_DIR"
(( N_SHUFFLES > 0 )) && mkdir -p "$RAW_SHUFFLES_DIR"

seq_name=$(basename "$fasta_file" .fa)
output_file="$RESULTS_DIR/results_${seq_name}.csv"
error_file="$ERRORS_DIR/errors_${seq_name}.log"

# --- Idempotent skip ---
if [[ -s "$output_file" ]] && (( $(awk 'END{print NR}' "$output_file") > 1 )); then
    echo "Skip $seq_name (done)"
    exit 0
fi

# --- Detect constrained vs unconstrained input ---
# Constraint line = third record line (typically) consisting only of
# ViennaRNA hard-constraint characters.
read_record() {
    # Echo on stdout: hdr<TAB>seq<TAB>constraint   (constraint may be empty)
    awk '
        /^>/ { hdr = substr($0, 2); next }
        /^[.()<>|x]+$/ && seq != "" { cons = $0; next }
        { seq = seq $0 }
        END {
            gsub(/[ \t\r]/, "", seq)
            gsub(/[ \t\r]/, "", cons)
            printf "%s\t%s\t%s\n", hdr, seq, cons
        }
    ' "$1"
}

IFS=$'\t' read -r hdr full_seq constraint_str < <(read_record "$fasta_file")
is_constrained=0
rnafold_orig_flags=(--noPS)
if [[ -n "$constraint_str" ]]; then
    is_constrained=1
    rnafold_orig_flags=(--noPS -C)
fi

# --- 1. Original MFE ---
extract_mfe() {
    awk 'match($0, /\(\s*([-+]?[0-9]*\.?[0-9]+)\s*\)$/, a) { print a[1]; exit }'
}

orig_mfe=$("$RNAFOLD_BIN" "${rnafold_orig_flags[@]}" < "$fasta_file" \
            2>>"$error_file" | extract_mfe)

if [[ ! "$orig_mfe" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
    echo "ERROR: invalid orig MFE for $seq_name" >> "$error_file"
    exit 1
fi

# --- 2. MFE-only fast path (also used by calibration) ---
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

if (( is_constrained == 0 )); then
    # --- 3a. Unconstrained: streaming pipe (fast path, unchanged) ---
    "$ESL_SHUFFLE_BIN" -d -N "$N_SHUFFLES" "$fasta_file" 2>>"$error_file" \
      | "$RNAFOLD_BIN" --noPS 2>>"$error_file" \
      | awk -v raw="$raw_tmp" -v s="$seq_name" '
          match($0, /\(\s*([-+]?[0-9]*\.?[0-9]+)\s*\)$/, a) {
              print a[1]
              print s "," ++iter "," a[1] >> raw
          }' \
      | sort -n > "$stats_file"
else
    # --- 3b. Constrained: shuffle UTR halves independently, rejoin, fold ---
    # Parse 5UTR / linker / 3UTR boundaries from the constraint string.
    # Format guaranteed by 01_extract.py: '<' * n5 + 'x' * nlink + '>' * n3
    n5=$(awk -v c="$constraint_str" 'BEGIN{
        n=0; while (substr(c, n+1, 1) == "<") n++; print n
    }')
    nlink=$(awk -v c="$constraint_str" -v n5="$n5" 'BEGIN{
        n=0; while (substr(c, n5 + n + 1, 1) == "x") n++; print n
    }')
    n3=$(awk -v c="$constraint_str" -v skip="$((n5))" -v lk="$nlink" 'BEGIN{
        n=0; while (substr(c, skip + lk + n + 1, 1) == ">") n++; print n
    }')

    if (( n5 < 1 || n3 < 1 || (n5 + nlink + n3) != ${#constraint_str} )); then
        echo "ERROR: malformed UTR_pair constraint for $seq_name" >> "$error_file"
        echo "  constraint=$constraint_str" >> "$error_file"
        echo "  parsed n5=$n5 nlink=$nlink n3=$n3 (expected total ${#constraint_str})" >> "$error_file"
        rm -f "$raw_tmp"
        exit 1
    fi

    seq_5utr="${full_seq:0:n5}"
    linker="${full_seq:n5:nlink}"
    seq_3utr="${full_seq:n5+nlink:n3}"

    # FIFO-style temp fastas for esl-shuffle (it needs a file or '-')
    fa_5="$TMP_DIR/${seq_name}_5utr.$$.fa"
    fa_3="$TMP_DIR/${seq_name}_3utr.$$.fa"
    printf ">5utr\n%s\n" "$seq_5utr" > "$fa_5"
    printf ">3utr\n%s\n" "$seq_3utr" > "$fa_3"

    iter=0
    : > "$stats_file"
    for ((i = 1; i <= N_SHUFFLES; i++)); do
        s5=$("$ESL_SHUFFLE_BIN" -d "$fa_5" 2>>"$error_file" \
                | grep -v '^>' | tr -d '\n\r ')
        s3=$("$ESL_SHUFFLE_BIN" -d "$fa_3" 2>>"$error_file" \
                | grep -v '^>' | tr -d '\n\r ')

        # Skip iterations where shuffle returned nothing usable
        [[ -n "$s5" && -n "$s3" ]] || continue
        (( ${#s5} == n5 && ${#s3} == n3 )) || continue

        hybrid="${s5}${linker}${s3}"
        energy=$(printf ">shuf_%d\n%s\n%s\n" "$i" "$hybrid" "$constraint_str" \
                    | "$RNAFOLD_BIN" --noPS -C 2>>"$error_file" \
                    | extract_mfe)

        if [[ "$energy" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
            echo "$energy" >> "$stats_file"
            iter=$((iter + 1))
            echo "$seq_name,$iter,$energy" >> "$raw_tmp"
        fi
    done

    rm -f "$fa_5" "$fa_3"
    sort -n -o "$stats_file" "$stats_file"
fi

n_valid=$(awk 'END{print NR}' "$stats_file")
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
echo "Done $seq_name ($n_valid/$N_SHUFFLES shuffles${is_constrained:+, constrained})"