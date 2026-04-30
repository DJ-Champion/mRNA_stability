#!/bin/bash
# 05_collate.sh
# Dispatcher. Sources paths and invokes the tool-specific collate script
# (e.g. tools/rnafold/collate.py).
#
# Usage:
#   ./bin/05_collate.sh -d human_liver -t rnafold
#   ./bin/05_collate.sh -d human_liver -t rnafold --include-raw-shuffles

set -euo pipefail
source "$(dirname "$0")/../lib/paths.sh"
parse_pipeline_args "$@"
resolve_paths

[[ -n "$TOOL" ]] || { echo "ERROR: --tool required" >&2; exit 1; }
[[ -r "$COLLATE_SCRIPT" ]] || { echo "ERROR: collate script not found: $COLLATE_SCRIPT" >&2; exit 1; }

# Pass through any extra args (e.g. --include-raw-shuffles)
exec python3 "$COLLATE_SCRIPT" "${PASSTHROUGH_ARGS[@]:-}"