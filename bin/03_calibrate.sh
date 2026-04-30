#!/bin/bash
# 03_calibrate.sh
# Time representative samples per tier under the tool's calib_params (e.g.
# MFE-only for RNAfold), extrapolate to full work via the tool's
# predict_wall_s / predict_rss_mb hooks, and emit SLURM resource recommendations.
#
# Usage (run on a compute node):
#   srun --pty -p aoraki --time=01:00:00 -c 1 --mem=4G bash -c \
#       './bin/03_calibrate.sh --dataset human_liver --tool rnafold'
#
# With --verify: also runs ONE sample per tier with the actual tier_params
# (full work) and records measured-vs-predicted divergence in verify.tsv.

set -euo pipefail
source "$(dirname "$0")/../lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

[[ -n "$TOOL" ]] || { echo "ERROR: --tool required" >&2; exit 1; }

[[ -s "$LISTS_DIR/.lengths.tsv" ]] || {
    echo "ERROR: lengths cache missing — run 02_stratify.sh first" >&2
    exit 1
}

read -ra tier_shuf_check <<< "${TIER_SHUFFLES:-}"
# Tools without per-tier policy still work — predict_wall_s falls back to identity.

calib_root="$CALIBRATION_ROOT/$(date +%Y%m%d_%H%M%S)"
mkdir -p "$calib_root"
report="$calib_root/recommendations.tsv"
detail="$calib_root/per_sample.tsv"
verify_file="$calib_root/verify.tsv"

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

# Run worker under /usr/bin/time -v with given env-var overrides.
# Echoes "wall_s rss_mb" on success, returns nonzero on failure.
run_timed() {
    local fasta="$1" sandbox="$2" env_kvs="$3"
    local out
    local -a env_args=()
    if [[ -n "$env_kvs" ]]; then
        # env_kvs is comma-separated KEY=VALUE pairs
        IFS=',' read -ra kv_arr <<< "$env_kvs"
        for kv in "${kv_arr[@]}"; do
            env_args+=("$kv")
        done
    fi

    out=$(/usr/bin/time -v env \
            DATASET="$DATASET" TOOL="$TOOL" \
            RESULTS_DIR="$sandbox/results" \
            ERRORS_DIR="$sandbox/errors" \
            TMP_DIR="$sandbox/tmp" \
            TOOL_DIR="$sandbox" \
            "${env_args[@]}" \
            "$WORKER_SCRIPT" "$fasta" 2>&1) || return 1

    local wall rss_kb rss_mb wall_int
    wall=$(echo "$out" | awk '/Elapsed \(wall clock\)/ {print $NF}' | parse_elapsed)
    rss_kb=$(echo "$out" | awk '/Maximum resident set size/ {print $NF}')
    rss_mb=$(( rss_kb / 1024 ))
    wall_int=$(awk -v w="$wall" 'BEGIN{printf "%d", w + 0.5}')
    echo "$wall_int $rss_mb"
}

# --- Headers ---
printf "tier\tn_samples\tmax_len\tmeasured_wall_s\tmeasured_rss_mb\tpredicted_wall_s\tpredicted_rss_mb\trec_time\trec_mem\n" > "$report"
printf "tier\tseqname\tlength\tmeasured_wall_s\tmeasured_rss_mb\tpredicted_wall_s\tpredicted_rss_mb\n" > "$detail"
(( VERIFY )) && printf "tier\tseqname\tlength\tpredicted_wall_s\tmeasured_full_wall_s\tdivergence_pct\n" > "$verify_file"

calib_kvs=$(calib_params)
echo "Calibration sandbox: $calib_root"
echo "Tool: $TOOL  |  calib_params: ${calib_kvs:-<none>}"
echo "Samples per tier: $SAMPLES_PER_TIER"
echo

# --- Per-tier loop ---
for tier in $(seq 1 10); do
    list="$LISTS_DIR/tier_${tier}.txt"
    [[ -s "$list" ]] || { echo "Tier $tier: empty — skipping"; continue; }

    # Pick the SAMPLES_PER_TIER longest sequences from this tier.
    # Lengths cache rows: <length>\t<seqname>\t<tier>
    samples=$(awk -F'\t' -v t="$tier" -v n="$SAMPLES_PER_TIER" '
        $3 == t {
            i = count
            while (i > 0 && len[i] < $1) {
                len[i+1] = len[i]; name[i+1] = name[i]; i--
            }
            len[i+1] = $1; name[i+1] = $2
            if (count < n) count++
        }
        END { for (i = 1; i <= count; i++) print len[i] "\t" name[i] }
    ' "$LISTS_DIR/.lengths.tsv")

    [[ -n "$samples" ]] || { echo "Tier $tier: no samples found"; continue; }

    max_meas_wall=0; max_meas_rss=0
    max_pred_wall=0; max_pred_rss=0
    max_len=0; n=0

    echo "=== Tier $tier ==="

    while IFS=$'\t' read -r len seqname; do
        fasta="$LISTS_DIR/tier_${tier}/${seqname}.fa"
        [[ -r "$fasta" ]] || { echo "  missing: $fasta" >&2; continue; }

        echo "  $seqname  (len=$len)"

        # Per-sample sandbox so concurrent runs don't collide
        sandbox="$calib_root/sample_${tier}_${seqname}"
        mkdir -p "$sandbox"

        if ! result=$(run_timed "$fasta" "$sandbox" "$calib_kvs"); then
            echo "    FAILED" >&2
            continue
        fi
        read -r wall_int rss_mb <<< "$result"

        # Extrapolate to full work
        pred_wall=$(predict_wall_s "$wall_int" "$tier")
        pred_rss=$(predict_rss_mb "$rss_mb" "$tier")

        printf "    measured wall=%ss rss=%dMB  ->  predicted wall=%ss rss=%dMB\n" \
            "$wall_int" "$rss_mb" "$pred_wall" "$pred_rss"
        printf "%d\t%s\t%d\t%d\t%d\t%d\t%d\n" \
            "$tier" "$seqname" "$len" \
            "$wall_int" "$rss_mb" "$pred_wall" "$pred_rss" >> "$detail"

        (( wall_int  > max_meas_wall )) && max_meas_wall=$wall_int
        (( rss_mb    > max_meas_rss  )) && max_meas_rss=$rss_mb
        (( pred_wall > max_pred_wall )) && max_pred_wall=$pred_wall
        (( pred_rss  > max_pred_rss  )) && max_pred_rss=$pred_rss
        (( len       > max_len       )) && max_len=$len
        n=$((n + 1))
    done <<< "$samples"

    if (( n == 0 )); then
        echo "Tier $tier: no successful samples — skipping recommendation"
        continue
    fi

    rec_time_s=$(( max_pred_wall * SAFETY_FACTOR_TIME ))
    rec_mem=$((    max_pred_rss  * SAFETY_FACTOR_MEM  ))
    (( rec_mem < 200 )) && rec_mem=200    # floor

    rec_time=$(fmt_slurm_time "$rec_time_s")
    printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%dMB\n" \
        "$tier" "$n" "$max_len" \
        "$max_meas_wall" "$max_meas_rss" "$max_pred_wall" "$max_pred_rss" \
        "$rec_time" "$rec_mem" >> "$report"
done

# --- Optional --verify pass: one full-work run per tier ---
if (( VERIFY )); then
    echo
    echo "=== --verify: running ONE full-work sample per tier ==="
    while IFS=$'\t' read -r tier n_samples max_len meas_wall meas_rss pred_wall pred_rss rec_time rec_mem; do
        [[ "$tier" == "tier" ]] && continue

        # Pick the longest sample for this tier
        sample=$(awk -F'\t' -v t="$tier" '
            $3 == t && $1 > best { best=$1; bn=$2 }
            END { if (bn != "") print best "\t" bn }
        ' "$LISTS_DIR/.lengths.tsv")
        [[ -n "$sample" ]] || continue
        IFS=$'\t' read -r len seqname <<< "$sample"

        fasta="$LISTS_DIR/tier_${tier}/${seqname}.fa"
        [[ -r "$fasta" ]] || continue

        full_kvs=$(tier_params "$tier")
        echo "  Tier $tier: $seqname (len=$len) with full ${full_kvs:-<defaults>}"

        sandbox="$calib_root/verify_${tier}_${seqname}"
        mkdir -p "$sandbox"

        if ! result=$(run_timed "$fasta" "$sandbox" "$full_kvs"); then
            echo "    FAILED" >&2
            continue
        fi
        read -r full_wall full_rss <<< "$result"

        # Predicted from this sample's measured calib value would require
        # re-running calibration on this exact sequence; instead, compare
        # against the per_sample.tsv prediction for this same seqname.
        sample_pred=$(awk -F'\t' -v s="$seqname" '$2==s {print $6}' "$detail")
        sample_pred="${sample_pred:-NA}"

        if [[ "$sample_pred" != "NA" && "$sample_pred" -gt 0 ]]; then
            divergence=$(awk -v p="$sample_pred" -v m="$full_wall" \
                'BEGIN{printf "%.1f", (m - p) * 100.0 / p}')
        else
            divergence="NA"
        fi

        printf "%d\t%s\t%d\t%s\t%d\t%s\n" \
            "$tier" "$seqname" "$len" "$sample_pred" "$full_wall" "$divergence" \
            >> "$verify_file"
        printf "    predicted=%ss  measured_full=%ss  divergence=%s%%\n" \
            "$sample_pred" "$full_wall" "$divergence"
    done < "$report"
fi

# --- Pretty print ---
echo
echo "Report: $report"
column -t -s $'\t' "$report"
echo
echo "Per-sample detail: $detail"
if (( VERIFY )) && [[ -s "$verify_file" ]]; then
    echo
    echo "Verify: $verify_file"
    column -t -s $'\t' "$verify_file"
fi

# --- Symlink "latest" within this tool's calibration dir ---
ln -sfn "$(basename "$calib_root")" "$CALIBRATION_ROOT/latest"
echo
echo "Latest symlink: $CALIBRATION_ROOT/latest -> $(basename "$calib_root")"