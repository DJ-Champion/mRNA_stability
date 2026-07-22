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
#     hotspot_mfe: -10.0            # coverage threshold (kcal/mol); see below
#
# Shell/env overrides still win, e.g. RNALFOLD_N_SHUFFLES=100.
#
# hotspot_mfe is the energy cutoff defining a local-structure "hotspot" for
# hotspot_coverage_fraction: a structure counts toward coverage iff its MFE is
# strictly below this value. It is a scientific parameter — sweep it (e.g.
# -10 / -15 / -20) and check the coverage ranking is stable. A fixed cutoff is
# deliberate: the shuffle-normalised hotspot_coverage_zscore uses the same
# threshold for observed and null, so composition sensitivity is handled there.

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

: "${RNALFOLD_BIN:=/home/chado47p/Software/ViennaRNA/bin/RNALfold}"
: "${ESL_SHUFFLE_BIN:=/home/chado47p/Software/hmmer3.4/bin/esl-shuffle}"

: "${RNALFOLD_WINDOW_LENGTH:=$(_rnalfold_yaml window_length 63)}"
: "${RNALFOLD_N_SHUFFLES:=$(_rnalfold_yaml n_shuffles 1000)}"
: "${MIN_VALID_PERCENT:=$(_rnalfold_yaml min_valid_percent 90)}"
: "${RNALFOLD_HOTSPOT_MFE:=$(_rnalfold_yaml hotspot_mfe -10.0)}"
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
export MIN_VALID_PERCENT RNALFOLD_HOTSPOT_MFE TOOL_REGIONS PATH

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
        'BEGIN { printf "%d\n", int((m +1) * (n + 1)) }'
}

# Peak memory should be governed mainly by sequence/window size, not shuffle count.
predict_rss_mb() {
    echo "$1"
}
