#!/bin/bash
# tools/rnafold_shuffles/worker.sh
# Wrapper around tools/rnafold/worker.sh. The shared worker skips any sequence
# whose results file exists, but 03_calibrate.sh (N_SHUFFLES=0) writes MFE-only
# results for its sample sequences into the same results dir. Those would make
# a real run skip the sequence without producing shuffles. So, when shuffles
# are requested and a results file has no raw_shuffles file beside it, drop the
# stale results file and let the shared worker redo the sequence.
#
# Usage: tools/rnafold_shuffles/worker.sh <fasta_file>
# Optional env: REAL_WORKER (default tools/rnafold/worker.sh; for tests)

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

fasta_file="${1:?Usage: $0 <fasta_file>}"
seq_name=$(basename "$fasta_file" .fa)
results_file="$RESULTS_DIR/results_${seq_name}.csv"
raw_file="$TOOL_DIR/raw_shuffles/${seq_name}_raw_shuffles.csv"

if (( N_SHUFFLES > 0 )) && [[ -e "$results_file" && ! -s "$raw_file" ]]; then
    echo "Removing stale MFE-only results for $seq_name (no raw shuffles): $results_file"
    rm -f "$results_file"
fi

exec "${REAL_WORKER:-$PROJECT_ROOT/tools/rnafold/worker.sh}" "$fasta_file"
