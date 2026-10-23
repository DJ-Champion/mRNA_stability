#!/bin/bash
# 04_submit.sh
# Submit SLURM array jobs. Two modes:
#
# (1) Submit ALL tiers using calibration recommendations (default):
#       ./bin/04_submit.sh -d human_liver -t rnafold
#
# (2) Submit a specific list (overrides calibration), e.g. resubmissions:
#       ./bin/04_submit.sh -d human_liver -t rnafold --list lists/redo.txt \
#           SLURM_TIME=00:30:00 SLURM_MEM_PER_CPU=500MB
#       (or set those as env vars before the call)

set -euo pipefail
source "$(dirname "$0")/../lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

[[ -n "$TOOL" ]] || { echo "ERROR: --tool required" >&2; exit 1; }

mkdir -p "$SLURM_LOG_DIR"

# Aoraki caps submitted jobs at 5000 per user, and every array task counts.
# To stay under it, each array task processes PACK consecutive FASTAs (the
# worker is idempotent, so reruns skip finished ones) and its time limit is
# scaled by PACK.
time_to_s() {   # [D-]HH:MM:SS | MM:SS -> seconds
    local t="$1" d=0
    [[ "$t" == *-* ]] && { d=${t%%-*}; t=${t#*-}; }
    IFS=: read -ra p <<< "$t"
    case ${#p[@]} in
        3) echo $(( d*86400 + 10#${p[0]}*3600 + 10#${p[1]}*60 + 10#${p[2]} )) ;;
        2) echo $(( d*86400 + 10#${p[0]}*60 + 10#${p[1]} )) ;;
        *) echo $(( d*86400 + 10#${p[0]}*60 )) ;;
    esac
}
s_to_time() { printf "%02d:%02d:%02d" $(($1/3600)) $((($1%3600)/60)) $(($1%60)); }

# Free submit slots = cap minus whatever this user already has queued/running.
queued=$(squeue -u "$USER" -h -r 2>/dev/null | awk 'END{print NR+0}') || queued=0
SUBMIT_BUDGET=$(( MAX_SUBMITTED_JOBS - queued ))
(( SUBMIT_BUDGET > 0 )) || { echo "ERROR: no submit slots free ($queued jobs already queued, cap $MAX_SUBMITTED_JOBS)" >&2; exit 1; }

# $1 = number of FASTAs; $2 = PACK (FASTAs per task)
submit_array() {
    local fasta_list="$1" rec_time="$2" rec_mem="$3" tag="$4" tier="${5:-}" pack="${6:-1}"

    [[ -s "$fasta_list" ]] || { echo "Empty/missing list: $fasta_list — skipping"; return; }

    local total ntasks array_spec task_time
    total=$(awk 'END{print NR}' "$fasta_list")
    ntasks=$(( (total + pack - 1) / pack ))
    array_spec="0-$((ntasks - 1))%${CONCURRENT_LIMIT}"

    local secs=$(( $(time_to_s "$rec_time") * pack ))
    if (( secs > MAX_WALL_S )); then
        echo "ERROR: tier ${tag}: pack=$pack x $rec_time exceeds partition limit ($(s_to_time "$MAX_WALL_S"))." >&2
        echo "       Use a longer-limit partition (SLURM_PARTITION/MAX_WALL_S) or fewer, larger submissions." >&2
        exit 1
    fi
    task_time=$(s_to_time "$secs")

    # Tool's per-tier env contribution (e.g. N_SHUFFLES=1000 for rnafold).
    # Empty for custom-list mode (caller's env flows through via --export=ALL).
    local extra_export=""
    [[ -n "$tier" ]] && extra_export=$(tier_params "$tier")

    local export_str="ALL,DATASET=$DATASET,TOOL=$TOOL,FASTA_LIST=$fasta_list,PROJECT_ROOT=$PROJECT_ROOT,PACK=$pack"
    [[ -n "$extra_export" ]] && export_str="$export_str,$extra_export"

    echo "Submitting: list=$(basename "$fasta_list")  n=$total  tasks=$ntasks (pack=$pack)  time=$task_time  mem=$rec_mem${extra_export:+  [$extra_export]}  [tag=$tag]"

    sbatch \
        --job-name="${SLURM_JOB_NAME}_${TOOL}_${tag}" \
        --partition="$SLURM_PARTITION" \
        --time="$task_time" \
        --mem-per-cpu="$rec_mem" \
        --cpus-per-task="$SLURM_CPUS_PER_TASK" \
        --array="$array_spec" \
        --output="$SLURM_LOG_DIR/${TOOL}_%A_%a.out" \
        --error="$SLURM_LOG_DIR/${TOOL}_%A_%a.err" \
        --export="$export_str" \
        "$PROJECT_ROOT/slurm/array.sbatch"
}

# --- Mode 1: explicit list ---
if [[ -n "$EXPLICIT_LIST" ]]; then
    [[ -s "$EXPLICIT_LIST" ]] || { echo "ERROR: list not found or empty: $EXPLICIT_LIST" >&2; exit 1; }
    n=$(awk 'END{print NR}' "$EXPLICIT_LIST")
    pack=$(( (n + SUBMIT_BUDGET - 1) / SUBMIT_BUDGET ))
    submit_array "$EXPLICIT_LIST" "$SLURM_TIME" "$SLURM_MEM_PER_CPU" "custom" "" "$pack"
    exit 0
fi

# --- Mode 2: submit all tiers from calibration ---
recs="$CALIBRATION_ROOT/latest/recommendations.tsv"
if [[ ! -s "$recs" ]]; then
    echo "ERROR: no calibration recommendations at $recs" >&2
    echo "Run 03_calibrate.sh first, or submit explicitly with --list <file>." >&2
    exit 1
fi

echo "Reading calibration recommendations from: $recs"
echo

# Helper: filter a tier list to seqnames whose region is in TOOL_REGIONS.
# Echoes the path of a filtered list (in $LISTS_DIR) — or the original list
# if no filter is set. Returns nonzero if the filtered list is empty.
filter_list_by_region() {
    local input_list="$1"
    if [[ -z "${TOOL_REGIONS:-}" ]]; then
        echo "$input_list"; return 0
    fi
    local filtered="$LISTS_DIR/.${TOOL}_$(basename "$input_list")"
    awk -F'\t' -v allowed="$TOOL_REGIONS" '
        BEGIN { m = split(allowed, arr, " "); for (i=1;i<=m;i++) ok[arr[i]] = 1 }
        # First file: manifest. Header has $1 == "seqname"; skip.
        NR == FNR {
            if ($1 != "seqname" && ($4 in ok)) keep[$1] = 1
            next
        }
        # Second file: tier list of FASTA paths. Match basename without .fa.
        {
            seq = $0; sub(/.*\//, "", seq); sub(/\.fa$/, "", seq)
            if (seq in keep) print
        }
    ' "$MANIFEST_TSV" "$input_list" > "$filtered"
    [[ -s "$filtered" ]] || return 1
    echo "$filtered"
}

# recommendations.tsv columns:
#   tier  n_samples  max_len  measured_wall_s  measured_rss_mb
#   predicted_wall_s  predicted_rss_mb  rec_time  rec_mem
# Pass 1: resolve each tier's list and count FASTAs, to size PACK globally.
declare -a T_TIER T_LIST T_TIME T_MEM
grand_total=0
while IFS=$'\t' read -r tier n_samples max_len meas_wall meas_rss pred_wall pred_rss rec_time rec_mem; do
    [[ "$tier" == "tier" ]] && continue
    list="$LISTS_DIR/tier_${tier}.txt"
    if ! list=$(filter_list_by_region "$list"); then
        echo "Tier $tier: no seqnames match TOOL_REGIONS=${TOOL_REGIONS:-<all>} — skipping"
        continue
    fi
    [[ -s "$list" ]] || continue
    T_TIER+=("$tier"); T_LIST+=("$list"); T_TIME+=("$rec_time"); T_MEM+=("$rec_mem")
    grand_total=$(( grand_total + $(awk 'END{print NR}' "$list") ))
done < "$recs"

# Smallest PACK such that sum_t ceil(n_t/PACK) <= SUBMIT_BUDGET.
pack=1
while :; do
    tasks=0
    for l in "${T_LIST[@]}"; do
        n=$(awk 'END{print NR}' "$l"); tasks=$(( tasks + (n + pack - 1) / pack ))
    done
    (( tasks <= SUBMIT_BUDGET )) && break
    pack=$(( pack + 1 ))
done
echo "Total FASTAs: $grand_total  submit budget: $SUBMIT_BUDGET  -> pack=$pack per task ($tasks tasks)"
echo

# Pass 2: submit.
for i in "${!T_TIER[@]}"; do
    submit_array "${T_LIST[$i]}" "${T_TIME[$i]}" "${T_MEM[$i]}" "t${T_TIER[$i]}" "${T_TIER[$i]}" "$pack"
done
echo
echo "All tier submissions queued. Monitor: squeue -u \$USER --name=${SLURM_JOB_NAME}_${TOOL}_*"