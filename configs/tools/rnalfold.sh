#!/bin/bash
# configs/tools/rnalfold.sh
# Tool config for RNAlfold: minimum local MFE + dinucleotide-shuffle null.
#
# Optional dataset YAML block:
#
# tools:
#   rnalfold:
#     regions: [mRNA, CDS, 5UTR, 3UTR]
#     window_length: 63
#     n_shuffles: 1000
#     min_valid_percent: 90
#
# Shell/env overrides still win, e.g. RNALFOLD_N_SHUFFLES=100.

TOOL_NAME=rnalfold
WORKER_SCRIPT="$PROJECT_ROOT/tools/rnalfold/worker.sh"
COLLATE_SCRIPT="$PROJECT_ROOT/tools/rnalfold/collate.py"

# Read one scalar/list setting from the dataset YAML's tools.rnalfold block.
# Usage: _rnalfold_yaml key default
_rnalfold_yaml() {
    local key="$1" default="$2"
    python3 - "$DATASET_YAML" "$key" "$default" <<'PY'
import sys
try:
    import yaml
except ImportError:
    print(sys.argv[3])
    raise SystemExit

path, key, default = sys.argv[1], sys.argv[2], sys.argv[3]
try:
    with open(path) as f:
        cfg = yaml.safe_load(f) or {}
    val = ((cfg.get('tools') or {}).get('rnalfold') or {}).get(key, default)
except Exception:
    val = default

if isinstance(val, (list, tuple)):
    print(' '.join(str(x) for x in val))
elif isinstance(val, bool):
    print('true' if val else 'false')
elif val is None:
    print(default)
else:
    print(val)
PY
}

: "${RNALFOLD_BIN:=RNALfold}"
: "${ESL_SHUFFLE_BIN:=esl-shuffle}"

: "${RNALFOLD_WINDOW_LENGTH:=$(_rnalfold_yaml window_length 63)}"
: "${RNALFOLD_N_SHUFFLES:=$(_rnalfold_yaml n_shuffles 1000)}"
: "${MIN_VALID_PERCENT:=$(_rnalfold_yaml min_valid_percent 90)}"
: "${TOOL_REGIONS:=$(_rnalfold_yaml regions 'mRNA CDS 5UTR 3UTR')}"

# Make binaries visible if absolute paths were supplied.
if [[ "$RNALFOLD_BIN" == */* ]]; then
    PATH="$(dirname "$RNALFOLD_BIN"):$PATH"
fi
if [[ "$ESL_SHUFFLE_BIN" == */* ]]; then
    PATH="$(dirname "$ESL_SHUFFLE_BIN"):$PATH"
fi

export TOOL_NAME WORKER_SCRIPT COLLATE_SCRIPT
export RNALFOLD_BIN ESL_SHUFFLE_BIN RNALFOLD_WINDOW_LENGTH RNALFOLD_N_SHUFFLES
export MIN_VALID_PERCENT TOOL_REGIONS PATH

# Full work by tier. Kept simple for first pass: every tier uses the same
# shuffle count from YAML/config.
tier_params() {
    echo "N_SHUFFLES=$RNALFOLD_N_SHUFFLES"
}

# Calibration fast path: original minimum local MFE only.
calib_params() {
    echo "N_SHUFFLES=0"
}

# RNAlfold full work is approximately one original fold plus N shuffled folds.
predict_wall_s() {
    local measured_s="$1"
    awk -v m="$measured_s" -v n="$RNALFOLD_N_SHUFFLES" \
        'BEGIN { printf "%d\n", int(m * (n + 1) + 0.5) }'
}

# Peak memory should be governed mainly by sequence/window size, not shuffle count.
predict_rss_mb() {
    echo "$1"
}
