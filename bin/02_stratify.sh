#!/bin/bash
# 02_stratify.sh
# Bin sequences from multifastas into 10 length tiers, using lengths
# already computed by 01_extract.py and stored in manifest.tsv.
#
# Outputs:
#   lists/tier_<n>/<seqname>.fa   individual FASTA per sequence (the working copy)
#   lists/tier_<n>.txt            list of paths (consumed by 04_submit.sh)
#   lists/.lengths.tsv            cached length info (used by 03_calibrate.sh)
#
# Auto-invalidates: if the manifest is newer than the existing tier_<n> dirs,
# they are wiped and rebuilt. Override with --no-invalidate.

set -euo pipefail
source "$(dirname "$0")/../config/config.main.sh"

NO_INVALIDATE=0
for arg in "$@"; do
    case "$arg" in
        --no-invalidate) NO_INVALIDATE=1 ;;
        -h|--help)
            sed -n '2,/^set/p' "$0" | sed 's/^# \{0,1\}//;$d'
            exit 0 ;;
    esac
done

[[ -s "$MANIFEST_TSV" ]] || {
    echo "ERROR: manifest not found at $MANIFEST_TSV" >&2
    echo "Run 01_extract.py first." >&2
    exit 1
}

read -ra bounds <<< "$TIER_BOUNDS"
(( ${#bounds[@]} == 9 )) || {
    echo "ERROR: TIER_BOUNDS must have 9 values (got ${#bounds[@]})" >&2
    exit 1
}

# --- Auto-invalidate stale tier dirs ---
if (( NO_INVALIDATE == 0 )); then
    needs_rebuild=0
    if [[ ! -f "$TIER_ROOT/.lengths.tsv" ]]; then
        needs_rebuild=1
    elif [[ "$MANIFEST_TSV" -nt "$TIER_ROOT/.lengths.tsv" ]]; then
        echo "Manifest newer than tier cache — rebuilding tiers."
        needs_rebuild=1
    fi
    if (( needs_rebuild )); then
        for i in $(seq 1 10); do
            rm -rf "$TIER_ROOT/tier_${i}" "$TIER_ROOT/tier_${i}.txt"
        done
    fi
fi

# --- Prepare tier dirs and lists ---
for i in $(seq 1 10); do
    mkdir -p "$TIER_ROOT/tier_${i}"
    : > "$TIER_ROOT/tier_${i}.txt"
done
length_cache="$TIER_ROOT/.lengths.tsv"
: > "$length_cache"

# --- Build a seqname -> region lookup so we can find the source multifasta ---
# manifest.tsv columns: seqname, gene_id, transcript_id, region, length, ...
declare -A REGION_OF
declare -A LEN_OF
while IFS=$'\t' read -r seqname gene_id tx_id region length rest; do
    [[ "$seqname" == "seqname" ]] && continue   # header
    REGION_OF["$seqname"]="$region"
    LEN_OF["$seqname"]="$length"
done < "$MANIFEST_TSV"

echo "Loaded ${#REGION_OF[@]} sequences from manifest."

# --- Stream each multifasta once, splitting records into tier dirs ---
# Group seqnames by region so each multifasta is read exactly once.
declare -A SEQS_BY_REGION
for seqname in "${!REGION_OF[@]}"; do
    region="${REGION_OF[$seqname]}"
    SEQS_BY_REGION[$region]+="$seqname"$'\n'
done

assign_tier() {
    local len=$1
    local i
    for i in $(seq 0 8); do
        if (( len < bounds[i] )); then
            echo $((i + 1))
            return
        fi
    done
    echo 10
}

for region in "${!SEQS_BY_REGION[@]}"; do
    mfa="$MULTIFASTA_DIR/extracted_${region}.fa"
    [[ -s "$mfa" ]] || {
        echo "WARNING: multifasta missing for region '$region': $mfa" >&2
        continue
    }

    awk -v cache="$length_cache" -v root="$TIER_ROOT" -v bounds_str="$TIER_BOUNDS" '
        BEGIN {
            n = split(bounds_str, b, " ")
        }
        function assign_tier(len,    i) {
            for (i = 1; i <= 9; i++) if (len < b[i]) return i
            return 10
        }
        function flush(    tier, path, list) {
            if (name == "") return
            tier = assign_tier(len)
            path = root "/tier_" tier "/" name ".fa"
            list = root "/tier_" tier ".txt"
            printf ">%s\n%s\n", name, seq > path
            close(path)
            print path >> list
            print len "\t" name "\t" tier >> cache
        }
        /^>/ {
            flush()
            name = substr($1, 2); seq = ""; len = 0; next
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
    count=$(wc -l < "$TIER_ROOT/tier_${i}.txt")
    total=$((total + count))
    printf "%-5d %-16s %d\n" "$i" "$range" "$count"
done
echo "----- total: $total"
