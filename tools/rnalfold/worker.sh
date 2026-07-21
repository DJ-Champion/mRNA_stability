#!/bin/bash
# tools/rnalfold/worker.sh
# Process one FASTA: local RNALfold structure features + dinucleotide-shuffle null.
#
# Per-sequence features (exact definitions in tools/rnalfold/summarize.py). All
# are computed over RNALfold's emitted structure lines only; the trailing
# summary line (summed energy of the non-overlapping decomposition) is excluded.
#   min_local_mfe                    strongest (most negative) local fold
#   median_local_mfe                 median energy across emitted local folds
#   hotspot_coverage_fraction        fraction of the region covered by local
#                                    folds below the hotspot MFE threshold
#                                    (overlapping hits merged); length-invariant
#   median_shuffle_min_local_mfe     null centre for min_local_mfe
#   pvalue_min_local_mfe             empirical p vs dinucleotide-shuffled null
#   zscore_min_local_mfe             robust (MAD) z of min vs null
#   median_shuffle_hotspot_coverage  null centre for coverage
#   hotspot_coverage_zscore          robust (MAD) z of coverage vs null
#   n_valid / n_missing / n_shuffles_requested       shuffle QC
#   rnalfold_window_length / hotspot_mfe_threshold   reproducibility params
#
# Required env, normally set by lib/paths.sh after sourcing configs/tools/rnalfold.sh:
#   RESULTS_DIR, ERRORS_DIR, TMP_DIR, TOOL_DIR
#   RNALFOLD_BIN, ESL_SHUFFLE_BIN, RNALFOLD_WINDOW_LENGTH, MIN_VALID_PERCENT
#   RNALFOLD_HOTSPOT_MFE
#   N_SHUFFLES  (set by tier_params/calib_params; falls back to RNALFOLD_N_SHUFFLES)
#
# Output:
#   $RESULTS_DIR/results_<seqname>.csv
#   $TOOL_DIR/raw_shuffles/<seqname>_raw_shuffles.csv  when N_SHUFFLES > 0

set -euo pipefail

RNALFOLD_SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$(cd "$RNALFOLD_SRC_DIR/../.." && pwd)/lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

SUMMARIZE="$RNALFOLD_SRC_DIR/summarize.py"

fasta_file="${1:?Usage: $0 <fasta_file>}"
[[ -r "$fasta_file" ]] || { echo "ERROR: cannot read $fasta_file" >&2; exit 1; }

: "${N_SHUFFLES:=${RNALFOLD_N_SHUFFLES:-1000}}"
: "${RNALFOLD_WINDOW_LENGTH:=63}"
: "${MIN_VALID_PERCENT:=90}"
: "${RNALFOLD_HOTSPOT_MFE:=-10.0}"

command -v "$RNALFOLD_BIN" >/dev/null 2>&1 || {
    echo "ERROR: RNALfold executable not found: $RNALFOLD_BIN" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 not found (required by rnalfold summarize.py)" >&2; exit 1; }
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

# Region length (spliced nt) for the coverage denominator. esl-shuffle -d
# preserves length, so the same value applies to every shuffled record.
seq_length=$(awk '!/^>/ { gsub(/[[:space:]]/, ""); n += length($0) } END { print n+0 }' "$fasta_file")
if (( seq_length <= 0 )); then
    echo "ERROR: empty sequence for $seq_name" >> "$error_file"
    exit 1
fi

HEADER="seq_name,min_local_mfe,median_local_mfe,median_shuffle_min_local_mfe,pvalue_min_local_mfe,zscore_min_local_mfe,hotspot_coverage_fraction,median_shuffle_hotspot_coverage,hotspot_coverage_zscore,n_valid,n_missing,n_shuffles_requested,rnalfold_window_length,hotspot_mfe_threshold"

# Original fold -> min_local_mfe, median_local_mfe, hotspot_coverage_fraction.
orig_stats=$("$RNALFOLD_BIN" -i "$fasta_file" -L "$RNALFOLD_WINDOW_LENGTH" 2>>"$error_file" \
    | python3 "$SUMMARIZE" original \
        --seq-length "$seq_length" --hotspot-mfe "$RNALFOLD_HOTSPOT_MFE" 2>>"$error_file") || {
        echo "ERROR: invalid original RNAlfold summary for $seq_name" >> "$error_file"
        exit 1
    }
IFS=$'\t' read -r orig_min orig_median orig_cov <<< "$orig_stats"

if [[ ! "$orig_min" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
    echo "ERROR: non-numeric original RNAlfold minimum local MFE for $seq_name: $orig_min" >> "$error_file"
    exit 1
fi

# Calibration fast path: original fold only (shuffle-derived columns are NA).
if (( N_SHUFFLES == 0 )); then
    {
        echo "$HEADER"
        printf "%s,%s,%s,NA,NA,NA,%s,NA,NA,0,0,0,%d,%s\n" \
            "$seq_name" "$orig_min" "$orig_median" "$orig_cov" \
            "$RNALFOLD_WINDOW_LENGTH" "$RNALFOLD_HOTSPOT_MFE"
    } > "$output_file"
    [[ -s "$error_file" ]] || rm -f "$error_file"
    echo "Done $seq_name (RNAlfold original-only: min=$orig_min median=$orig_median cov=$orig_cov)"
    exit 0
fi

min_valid=$(( N_SHUFFLES * MIN_VALID_PERCENT / 100 ))
raw_tmp="$TMP_DIR/rnalfold_raw_${seq_name}.$$.csv"
raw_file="$RAW_SHUFFLES_DIR/${seq_name}_raw_shuffles.csv"

# Stream shuffled FASTA records through RNAlfold; summarise the null in one
# pass. summarize.py writes the raw per-shuffle table to $raw_tmp and exits
# non-zero if too few shuffles folded (no output file is written in that case,
# so find_missing.sh will re-queue the sequence).
agg=$("$ESL_SHUFFLE_BIN" -d -N "$N_SHUFFLES" "$fasta_file" 2>>"$error_file" \
    | "$RNALFOLD_BIN" -L "$RNALFOLD_WINDOW_LENGTH" 2>>"$error_file" \
    | python3 "$SUMMARIZE" aggregate \
        --seq-length "$seq_length" --hotspot-mfe "$RNALFOLD_HOTSPOT_MFE" \
        --orig-min "$orig_min" --orig-coverage "$orig_cov" \
        --n-shuffles "$N_SHUFFLES" --min-valid "$min_valid" \
        --seq-name "$seq_name" --raw-out "$raw_tmp" 2>>"$error_file") || {
        echo "ERROR: RNAlfold shuffle aggregation failed or insufficient valid shuffles for $seq_name" >> "$error_file"
        echo "Note: esl-shuffle -d can fail or produce unusable records for short/low-complexity sequences." >> "$error_file"
        rm -f "$raw_tmp"
        exit 1
    }

IFS=$'\t' read -r med_shuf_min pval z_min med_shuf_cov z_cov n_valid n_missing <<< "$agg"

mv "$raw_tmp" "$raw_file"

{
    echo "$HEADER"
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%d,%d,%s\n" \
        "$seq_name" "$orig_min" "$orig_median" "$med_shuf_min" "$pval" "$z_min" \
        "$orig_cov" "$med_shuf_cov" "$z_cov" "$n_valid" "$n_missing" \
        "$N_SHUFFLES" "$RNALFOLD_WINDOW_LENGTH" "$RNALFOLD_HOTSPOT_MFE"
} > "$output_file"

[[ -s "$error_file" ]] || rm -f "$error_file"
echo "Done $seq_name ($n_valid/$N_SHUFFLES RNAlfold shuffles)"
