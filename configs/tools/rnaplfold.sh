#!/bin/bash
# configs/tools/rnaplfold.sh
# Tool config for RNAplfold: local accessibility (unpaired probability) of the
# start region (the AUG and its flanks), read from a fold of the FULL mRNA.
#
# Optional dataset YAML block:
#
# tools:
#   rnaplfold:
#     window: 300             # -W, averaging window
#     span: 150               # -L, max base-pair span (always set it; -W alone
#                             #     silently caps the span at the window size)
#     u_values: [1, 4, 8, 16] # stretch lengths reported (RNAplfold -u = max)
#     start_upstream: 30      # nt before the AUG   (default: extraction.codon_flank_length, else 30)
#     start_downstream: 30    # nt after the AUG    (default: extraction.codon_flank_length, else 30)
#
# The start window is [AUG - upstream, AUG + 3 + downstream), clipped to the
# transcript, which is the same span as the extracted start_codon_region.
# Shell/env overrides still win, e.g. RNAPLFOLD_WINDOW=500 RNAPLFOLD_SPAN=300.
#
# Fold the full mRNA, never the extracted start_codon_region: a fragment has
# artificial ends and its accessibilities differ from the in-transcript values
# for ~150 nt from the cut. Report window and span with any result.

TOOL_NAME=rnaplfold
WORKER_SCRIPT="$PROJECT_ROOT/tools/rnaplfold/worker.sh"
COLLATE_SCRIPT="$PROJECT_ROOT/tools/rnaplfold/collate.py"
TOOL_REGIONS=mRNA

# Read a setting from the dataset YAML. Usage: _rnaplfold_yaml key default [fallback_section.key]
# Looks in tools.rnaplfold.<key>, then the optional "section.key" fallback, then default.
_rnaplfold_yaml() {
    local key="$1" default="$2" fallback="${3:-}"
    python3 - "$DATASET_YAML" "$key" "$default" "$fallback" <<'PY'
import sys
path, key, default, fallback = sys.argv[1:5]
try:
    import yaml
    with open(path) as f:
        cfg = yaml.safe_load(f) or {}
    val = ((cfg.get('tools') or {}).get('rnaplfold') or {}).get(key)
    if val is None and fallback:
        sec, k = fallback.split('.')
        val = (cfg.get(sec) or {}).get(k)
    if val is None:
        val = default
except Exception:
    val = default

if isinstance(val, (list, tuple)):
    print(' '.join(str(x) for x in val))
else:
    print(val)
PY
}

: "${RNAPLFOLD_BIN:=/home/chado47p/Software/ViennaRNA/bin/RNAplfold}"
: "${RNAPLFOLD_WINDOW:=$(_rnaplfold_yaml window 300)}"
: "${RNAPLFOLD_SPAN:=$(_rnaplfold_yaml span 150)}"
: "${RNAPLFOLD_U_VALUES:=$(_rnaplfold_yaml u_values '1 4 8 16')}"
: "${RNAPLFOLD_START_UPSTREAM:=$(_rnaplfold_yaml start_upstream 30 extraction.codon_flank_length)}"
: "${RNAPLFOLD_START_DOWNSTREAM:=$(_rnaplfold_yaml start_downstream 30 extraction.codon_flank_length)}"

if [[ "$RNAPLFOLD_BIN" == */* ]]; then
    PATH="$(dirname "$RNAPLFOLD_BIN"):$PATH"
fi

export TOOL_NAME WORKER_SCRIPT COLLATE_SCRIPT TOOL_REGIONS PATH
export RNAPLFOLD_BIN RNAPLFOLD_WINDOW RNAPLFOLD_SPAN RNAPLFOLD_U_VALUES
export RNAPLFOLD_START_UPSTREAM RNAPLFOLD_START_DOWNSTREAM

# Time and memory are not measured yet: run calibration (03_calibrate.sh) and
# check the first real run. Defaults from lib/paths.sh scale 1:1 with the
# measured calibration run, which is reasonable because RNAplfold is roughly
# linear in length for a fixed window and span.
