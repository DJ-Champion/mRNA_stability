#!/bin/bash
# 02_stratify.sh
# Bin sequences from multifastas into 10 length tiers. Tier assignment uses
# lengths already in manifest.tsv (no re-scanning).
#
# Outputs (under $LISTS_DIR = $RUN_DIR/lists):
#   tier_<n>/<seqname>.fa   individual FASTA per sequence (working copy);
#                           may be 2-line (header+seq) or 3-line
#                           (header+seq+constraint) for UTR_pair-style records
#   tier_<n>.txt            list of paths (consumed by 04_submit.sh)
#   .lengths.tsv            cached length info (used by 03_calibrate.sh)
#
# Constraint-aware: if a non-header line consists only of ViennaRNA hard-
# constraint characters (.()<>|x), it is preserved verbatim as a third line
# in the per-sequence FASTA and excluded from sequence-length accounting.
# Tools that don't understand constraints can ignore this third line; tools
# that do (e.g. RNAfold -C) consume it directly.
#
# Auto-invalidates: if manifest is newer than the existing tier dirs they
# are wiped and rebuilt. Override with --no-invalidate.
#
# Usage:
#   ./bin/02_stratify.sh --dataset human_liver
#   ./bin/02_stratify.sh -d human_liver --no-invalidate

set -euo pipefail
source "$(dirname "$0")/../lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

# Tool-agnostic step: TOOL not required.
NO_INVALIDATE=0
for arg in "${PASSTHROUGH_ARGS[@]:-}"; do
    case "$arg" in
        --no-invalidate) NO_INVALIDATE=1 ;;
        -h|--help)
            sed -n '2,/^set/p' "$0" | sed 's/^# \{0,1\}//;$d'
            exit 0 ;;
    esac
done

[[ -s "$MANIFEST_TSV" ]] || {
    echo "ERROR: manifest not found at $MANIFEST_TSV" >&2
    echo "Run 01_extract.py --dataset $DATASET first." >&2
    exit 1
}

# Sanity check: dataset YAML newer than manifest = extraction is stale.
# 02_stratify only consumes manifest.tsv, so re-tiering on top of a stale
# manifest would silently propagate the staleness. Skipped if the user
# has explicitly opted out with --no-invalidate.
if (( NO_INVALIDATE == 0 )) && [[ "$DATASET_YAML" -nt "$MANIFEST_TSV" ]]; then
    echo "ERROR: $DATASET_YAML is newer than $MANIFEST_TSV" >&2
    echo "Dataset config has changed since the last extraction." >&2
    echo "Re-run: ./bin/01_extract.py --dataset $DATASET" >&2
    echo "(or pass --no-invalidate to bypass this check)" >&2
    exit 1
fi

read -ra bounds <<< "$TIER_BOUNDS"
(( ${#bounds[@]} == 9 )) || {
    echo "ERROR: TIER_BOUNDS must have 9 values (got ${#bounds[@]})" >&2
    exit 1
}

mkdir -p "$LISTS_DIR"

# --- Auto-invalidate stale tier dirs ---
if (( NO_INVALIDATE == 0 )); then
    needs_rebuild=0
    if [[ ! -f "$LISTS_DIR/.lengths.tsv" ]]; then
        needs_rebuild=1
    elif [[ "$MANIFEST_TSV" -nt "$LISTS_DIR/.lengths.tsv" ]]; then
        echo "Manifest newer than tier cache — rebuilding tiers."
        needs_rebuild=1
    fi
    if (( needs_rebuild )); then
        for i in $(seq 1 10); do
            rm -rf "$LISTS_DIR/tier_${i}" "$LISTS_DIR/tier_${i}.txt"
        done
    fi
fi

# --- Prepare tier dirs and lists ---
for i in $(seq 1 10); do
    mkdir -p "$LISTS_DIR/tier_${i}"
    : > "$LISTS_DIR/tier_${i}.txt"
done
length_cache="$LISTS_DIR/.lengths.tsv"
: > "$length_cache"

# --- Build a region-of-seq lookup ---
declare -A REGION_OF
while IFS=$'\t' read -r seqname gene_id tx_id region length rest; do
    [[ "$seqname" == "seqname" ]] && continue
    REGION_OF["$seqname"]="$region"
done < "$MANIFEST_TSV"

echo "Loaded ${#REGION_OF[@]} sequences from manifest."

# --- Group seqnames by region so each multifasta is read exactly once ---
declare -A SEQS_BY_REGION
for seqname in "${!REGION_OF[@]}"; do
    region="${REGION_OF[$seqname]}"
    SEQS_BY_REGION[$region]+="$seqname"$'\n'
done

# --- Stream each multifasta and split records into tier dirs ---
# Records may be 2-line (header + seq) or 3-line (header + seq + constraint).
# Constraint detection: a line of only [.()<>|x] characters (ViennaRNA hard-
# constraint alphabet). Sequence and constraint alphabets don't overlap.
for region in "${!SEQS_BY_REGION[@]}"; do
    mfa="$EXTRACT_DIR/extracted_${region}.fa"
    [[ -s "$mfa" ]] || {
        echo "WARNING: multifasta missing for region '$region': $mfa" >&2
        continue
    }

    awk -v cache="$length_cache" -v root="$LISTS_DIR" -v bounds_str="$TIER_BOUNDS" '
        BEGIN { n = split(bounds_str, b, " ") }
        function assign_tier(len,    i) {
            for (i = 1; i <= 9; i++) if (len < b[i]) return i
            return 10
        }
        function flush(    tier, path, list) {
            if (name == "") return
            tier = assign_tier(len)
            path = root "/tier_" tier "/" name ".fa"
            list = root "/tier_" tier ".txt"
            if (constraint != "") {
                printf ">%s\n%s\n%s\n", name, seq, constraint > path
            } else {
                printf ">%s\n%s\n", name, seq > path
            }
            close(path)
            print path >> list
            print len "\t" name "\t" tier >> cache
        }
        /^>/ {
            flush()
            name = substr($1, 2); seq = ""; constraint = ""; len = 0; next
        }
        # Constraint line: composed only of ViennaRNA hard-constraint chars.
        # Must be checked BEFORE the sequence-accumulator rule below.
        /^[.()<>|x]+$/ {
            constraint = $0
            next
        }
        { seq = seq $0; len += length($0) }
        END { flush() }
    ' "$mfa"
done

# --- Report ---
echo
printf "%-5s %-16s %s\n" "tier" "range" "count"
prev=0
total=0
for i in $(seq 1 10); do
    if (( i < 10 )); then
        upper="${bounds[$((i-1))]}"
        range="${prev}-${upper}"
        prev="$upper"
    else
        range=">${prev}"
    fi
    count=$(awk 'END{print NR}' "$LISTS_DIR/tier_${i}.txt")
    total=$((total + count))
    printf "%-5d %-16s %d\n" "$i" "$range" "$count"
done
echo "----- total: $total"