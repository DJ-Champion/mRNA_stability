# lib/paths.sh
# Shared by every shell script. Sources cluster + dataset + tool configs and
# resolves all paths for a given (dataset, tool) pair.
#
# Usage from a script:
#   source "$(dirname "$0")/../lib/paths.sh"
#   parse_pipeline_args "$@"
#   resolve_paths
#
# Recognized flags (consumed by parse_pipeline_args):
#   --dataset NAME / -d NAME    (or DATASET env var)
#   --tool    NAME / -t NAME    (or TOOL env var; not all scripts need it)
#   --list    PATH / -l PATH    (custom FASTA list — used by 04_submit / find_missing)
#   --verify  / -V              (calibration verify mode)
#
# Anything else is left in PASSTHROUGH_ARGS for the caller to handle.

# Globals set by parse_pipeline_args
DATASET="${DATASET:-}"
TOOL="${TOOL:-}"
EXPLICIT_LIST=""
VERIFY=0
PASSTHROUGH_ARGS=()

parse_pipeline_args() {
    PASSTHROUGH_ARGS=()
    while (( $# > 0 )); do
        case "$1" in
            --dataset|-d) DATASET="$2"; shift 2 ;;
            --tool|-t)    TOOL="$2"; shift 2 ;;
            --list|-l)    EXPLICIT_LIST="$2"; shift 2 ;;
            --verify|-V)  VERIFY=1; shift ;;
            *)            PASSTHROUGH_ARGS+=("$1"); shift ;;
        esac
    done
}

# Resolve all paths and source the relevant configs.
# Sets:
#   PROJECT_ROOT, RUNS_ROOT, RUN_DIR
#   DATASET_YAML, EXTRACT_DIR, MANIFEST_TSV, LISTS_DIR
# If TOOL is set, also:
#   TOOL_CONFIG, TOOL_DIR, RESULTS_DIR, ERRORS_DIR, TMP_DIR,
#   SLURM_LOG_DIR, CALIBRATION_ROOT
# Sources (in order): cluster.sh, then tool config (if TOOL set).
# Defines fallback hook functions if the tool config didn't override them.
resolve_paths() {
    [[ -n "$DATASET" ]] || {
        echo "ERROR: --dataset NAME (or DATASET env) required" >&2
        return 1
    }

    # PROJECT_ROOT = parent of lib/
    PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    : "${RUNS_ROOT:=$PROJECT_ROOT/runs}"

    DATASET_YAML="$PROJECT_ROOT/configs/datasets/${DATASET}.yaml"
    [[ -r "$DATASET_YAML" ]] || {
        echo "ERROR: dataset config not found: $DATASET_YAML" >&2
        return 1
    }

    RUN_DIR="$RUNS_ROOT/$DATASET"
    EXTRACT_DIR="$RUN_DIR/extracted_regions"
    MANIFEST_TSV="$EXTRACT_DIR/manifest.tsv"
    LISTS_DIR="$RUN_DIR/lists"

    # Always source cluster config
    source "$PROJECT_ROOT/configs/cluster.sh"

    if [[ -n "$TOOL" ]]; then
        TOOL_CONFIG="$PROJECT_ROOT/configs/tools/${TOOL}.sh"
        [[ -r "$TOOL_CONFIG" ]] || {
            echo "ERROR: tool config not found: $TOOL_CONFIG" >&2
            return 1
        }
        # PROJECT_ROOT must be set before sourcing tool config (tool config
        # references it for WORKER_SCRIPT etc.)
        export PROJECT_ROOT
        source "$TOOL_CONFIG"

        TOOL_DIR="$RUN_DIR/$TOOL"
        RESULTS_DIR="$TOOL_DIR/results"
        ERRORS_DIR="$TOOL_DIR/errors"
        TMP_DIR="$TOOL_DIR/tmp"
        SLURM_LOG_DIR="$TOOL_DIR/slurm_logs"
        CALIBRATION_ROOT="$TOOL_DIR/calibration"

        # Default hook functions if the tool config didn't define them.
        # Each one is overridable per tool.
        type tier_params >/dev/null 2>&1 || tier_params() { echo ""; }
        type calib_params >/dev/null 2>&1 || calib_params() { echo ""; }
        # predict_wall_s <measured_seconds> <tier> -> predicted_seconds
        type predict_wall_s >/dev/null 2>&1 || predict_wall_s() { echo "$1"; }
        # predict_rss_mb <measured_mb> <tier> -> predicted_mb
        type predict_rss_mb >/dev/null 2>&1 || predict_rss_mb() { echo "$1"; }
    fi

    export PROJECT_ROOT RUNS_ROOT RUN_DIR DATASET TOOL
    export DATASET_YAML EXTRACT_DIR MANIFEST_TSV LISTS_DIR
    if [[ -n "$TOOL" ]]; then
        export TOOL_CONFIG TOOL_DIR RESULTS_DIR ERRORS_DIR TMP_DIR
        export SLURM_LOG_DIR CALIBRATION_ROOT
    fi
}