#!/bin/bash
# 03_calibrate.sh
# Time the longest sequences from each tier and emit SLURM resource recommendations.
# Run on a compute node:
#   srun --pty -p aoraki --time=01:00:00 -c 1 --mem=4G bash -c './bin/03_calibrate.sh'

set -euo pipefail
source "$(dirname "$0")/../config/config.main.sh"

[[ -s "$TIER_ROOT/.lengths.tsv" ]] || {
    echo "ERROR: lengths cache missing — run 02_stratify.sh first" >&2
    exit 1
}

read -ra tier_shuffles <<< "$TIER_SHUFFLES"
(( ${#tier_shuffles[@]} == 10 )) || {
    echo "ERROR: TIER_SHUFFLES must have 10 values (got ${#tier_shuffles[@]})" >&2
    exit 1
}

calib_root="${CALIBRATION_ROOT}/$(date +%Y%m%d_%H%M%S)"
mkdir -p "$calib_root"
report="$calib_root/recommendations.tsv"
detail="$calib_root/per_sample.tsv"

# --- Helpers ---
parse_elapsed() {
    awk -F: '
        NF==3 {print $1*3600 + $2*60 + $3}
        NF==2 {print $1*60 + $2}
        NF==1 {print $1}'
}

fmt_slurm_time() {
    local s=$1
    s=$(( ( (s + 59) / 60 ) * 60 ))    # round up to nearest minute
    printf "%02d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60))
}

# --- Headers ---
printf "tier\tn_shuffles\tn_samples\tmax_len\tmax_wall_s\tmax_rss_mb\trec_time\trec_mem\n" > "$report"
printf "tier\tseqname\tlength\tn_shuffles\twall_s\trss_mb\n" > "$detail"

echo "Calibration sandbox: $calib_root"
echo "Samples per tier: $SAMPLES_PER_TIER"
echo

# --- Per-tier loop ---
for tier in $(seq 1 10); do
    list="$TIER_ROOT/tier_${tier}.txt"
    [[ -s "$list" ]] || { echo "Tier $tier: empty — skipping"; continue; }

    n_shuf="${tier_shuffles[$((tier-1))]}"

    # Pick the SAMPLES_PER_TIER longest sequences from this tier
    # Lengths cache rows: <length>\t<seqname>\t<tier>
    # Tier list rows:     <path>
    samples=$(awk -F'\t' -v t="$tier" '$3 == t {print $1 "\t" $2}' \
                  "$TIER_ROOT/.lengths.tsv" \
              | sort -rn \
              | head -n "$SAMPLES_PER_TIER")

    [[ -n "$samples" ]] || { echo "Tier $tier: no samples found"; continue; }

    max_wall=0
    max_rss=0
    max_len=0
    n=0

    echo "=== Tier $tier (n_shuffles=$n_shuf) ==="

    while IFS=$'\t' read -r len seqname; do
        # Locate the .fa for this seqname in the tier dir
        fasta="$TIER_ROOT/tier_${tier}/${seqname}.fa"
        [[ -r "$fasta" ]] || { echo "  missing: $fasta" >&2; continue; }

        echo "  $seqname  (len=$len)"

        # Run worker under /usr/bin/time -v in an isolated sandbox so we don't
        # pollute real results.
        out=$(/usr/bin/time -v env \
                N_SHUFFLES="$n_shuf" \
                RESULTS_DIR="$calib_root/results" \
                RAW_SHUFFLES_DIR="$calib_root/raw" \
                ERRORS_DIR="$calib_root/errors" \
                TMP_DIR="$calib_root/tmp" \
                "$WORKER_SCRIPT" "$fasta" 2>&1) || {
            echo "    FAILED" >&2
            continue
        }

        wall=$(echo "$out" | awk '/Elapsed \(wall clock\)/ {print $NF}' | parse_elapsed)
        rss_kb=$(echo "$out" | awk '/Maximum resident set size/ {print $NF}')
        rss_mb=$(( rss_kb / 1024 ))
        wall_int=$(awk -v w="$wall" 'BEGIN{printf "%d", w + 0.5}')

        printf "    wall=%ss  rss=%dMB\n" "$wall_int" "$rss_mb"
        printf "%d\t%s\t%d\t%d\t%d\t%d\n" \
            "$tier" "$seqname" "$len" "$n_shuf" "$wall_int" "$rss_mb" >> "$detail"

        (( wall_int > max_wall )) && max_wall=$wall_int
        (( rss_mb   > max_rss  )) && max_rss=$rss_mb
        (( len      > max_len  )) && max_len=$len
        n=$((n+1))
    done <<< "$samples"

    if (( n == 0 )); then
        echo "Tier $tier: no successful samples — skipping recommendation"
        continue
    fi

    rec_time_s=$(( max_wall * SAFETY_FACTOR_TIME ))
    rec_mem=$((  max_rss  * SAFETY_FACTOR_MEM  ))
    (( rec_mem < 200 )) && rec_mem=200    # floor

    rec_time=$(fmt_slurm_time "$rec_time_s")
    printf "%d\t%d\t%d\t%d\t%d\t%d\t%s\t%dMB\n" \
        "$tier" "$n_shuf" "$n" "$max_len" "$max_wall" "$max_rss" \
        "$rec_time" "$rec_mem" >> "$report"
done

# --- Pretty print ---
echo
echo "Report: $report"
column -t -s $'\t' "$report"
echo
echo "Per-sample detail: $detail"

# --- Symlink "latest" for convenience ---
ln -sfn "$(basename "$calib_root")" "${CALIBRATION_ROOT}/latest"
echo
echo "Latest symlink: ${CALIBRATION_ROOT}/latest -> $(basename "$calib_root")"
