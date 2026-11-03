#!/bin/bash
# bin/fetch_rmsk.sh
# Standalone step: download the UCSC hg38 RepeatMasker table and convert it to
# a BED6+ file in the same chromosome naming as canonical.gff (chr1, ...).
#
# Output (dataset-independent, so not under a dataset dir):
#   $RUNS_ROOT/_resources/hg38_rmsk/rmsk.bed.gz
#   $RUNS_ROOT/_resources/hg38_rmsk/rmsk.provenance.txt
#
# BED columns (0-based half-open, as BED):
#   1 chrom  2 start  3 end  4 repName  5 swScore  6 strand
#   7 repClass  8 repFamily  9 pctDiv (milliDiv/10)
# Strand is the repeat's orientation on the genome strand; step
# repeats_on_transcripts converts to transcript-relative sense/antisense.
#
# Only primary chromosomes (chr1-22, X, Y, M) are kept; alt/random/Un contigs
# are dropped because the genome FASTA used is the primary assembly.
#
# Usage:  ./bin/fetch_rmsk.sh            (needs internet; run on the cluster)
#         RUNS_ROOT=/path/to/runs ./bin/fetch_rmsk.sh

set -euo pipefail

BASE_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
: "${RUNS_ROOT:=$PROJECT_ROOT/runs}"
OUT_DIR="$RUNS_ROOT/_resources/hg38_rmsk"
OUT="$OUT_DIR/rmsk.bed.gz"

if [[ -s "$OUT" && "${1:-}" != "--force" ]]; then
    echo "Already present: $OUT (use --force to refetch)"
    exit 0
fi

mkdir -p "$OUT_DIR"
WORK="$(mktemp -d "$OUT_DIR/.fetch.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

# UCSC publishes no md5sum.txt for the database/ directory, so verify with the
# server's Content-Length plus a gzip integrity test.
echo "Downloading rmsk.txt.gz ..."
expected="$(curl -sSIL --retry 3 "$BASE_URL/rmsk.txt.gz" \
            | tr -d '\r' | awk 'tolower($1)=="content-length:"{n=$2} END{print n}')"
curl -fSL -s --retry 3 -o "$WORK/rmsk.txt.gz" "$BASE_URL/rmsk.txt.gz"

echo "Verifying download ..."
actual="$(stat -c %s "$WORK/rmsk.txt.gz")"
[[ -n "$expected" && "$expected" == "$actual" ]] || {
    echo "ERROR: size mismatch (server '$expected', got '$actual')" >&2; exit 1; }
gzip -t "$WORK/rmsk.txt.gz" || { echo "ERROR: corrupt gzip" >&2; exit 1; }

# rmsk.txt columns (1-based):
#  1 bin 2 swScore 3 milliDiv 4 milliDel 5 milliIns 6 genoName 7 genoStart
#  8 genoEnd 9 genoLeft 10 strand 11 repName 12 repClass 13 repFamily ...
echo "Converting to BED ..."
zcat "$WORK/rmsk.txt.gz" \
  | awk -F'\t' -v OFS='\t' '$6 ~ /^chr([0-9]+|X|Y|M)$/ {
        print $6, $7, $8, $11, $2, $10, $12, $13, $3/10 }' \
  | sort -k1,1V -k2,2n -k3,3n \
  | gzip -9 > "$WORK/rmsk.bed.gz"

n=$(zcat "$WORK/rmsk.bed.gz" | wc -l)
(( n > 4000000 )) || { echo "ERROR: only $n rows; expected ~5M" >&2; exit 1; }

mv "$WORK/rmsk.bed.gz" "$OUT"
{
    echo "source:    $BASE_URL/rmsk.txt.gz"
    echo "bytes:     $actual (matched server Content-Length; gzip -t ok)"
    echo "fetched:   $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "rows:      $n"
    echo "filter:    primary chromosomes only (chr1-22,X,Y,M)"
} > "$OUT_DIR/rmsk.provenance.txt"
echo "Wrote $OUT ($n rows)"
