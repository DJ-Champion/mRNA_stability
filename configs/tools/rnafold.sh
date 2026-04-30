#!/bin/bash
# configs/tools/rnafold.sh
# RNAfold tool config. Sourced by lib/paths.sh AFTER cluster.sh and AFTER
# PROJECT_ROOT is set.
#
# A tool config must define:
#   TOOL_NAME            (sanity check, must match filename stem)
#   WORKER_SCRIPT        (per-sequence executable)
#   COLLATE_SCRIPT       (per-tool combiner)
#   tier_params <tier>   (env vars to export to sbatch for this tier;
#                         comma-joined "K=V,K=V" or empty)
#   calib_params         (env vars for fast calibration runs;
#                         comma-joined or empty)
#   predict_wall_s <measured_s> <tier>   -> predicted seconds
#   predict_rss_mb  <measured_mb> <tier> -> predicted MB
#
# Default fallbacks (in lib/paths.sh) make tier_params/calib_params empty and
# predict_*_s/_mb the identity. Tools without an extrapolation trick can
# omit the latter two.

TOOL_NAME=rnafold
WORKER_SCRIPT="$PROJECT_ROOT/tools/rnafold/worker.sh"
COLLATE_SCRIPT="$PROJECT_ROOT/tools/rnafold/collate.py"

# --- Binaries (adjust once, here) ---
: "${RNAFOLD_BIN:=/home/chado47p/Software/ViennaRNA/bin/RNAfold}"
: "${ESL_SHUFFLE_BIN:=/home/chado47p/Software/hmmer3.4/bin/esl-shuffle}"

# Make tools visible to sub-shells without per-task conda activation
PATH="$(dirname "$RNAFOLD_BIN"):$(dirname "$ESL_SHUFFLE_BIN"):$PATH"

# --- Per-tier shuffle policy (10 values for 10 tiers) ---
# 0 = MFE-only fast path (no shuffles).
: "${TIER_SHUFFLES:=1000 1000 1000 1000 1000 1000 1000 100 10 0}"

# --- Worker parameters ---
: "${MIN_VALID_PERCENT:=90}"      # fraction of shuffles that must succeed
: "${N_SHUFFLES:=1000}"            # default for ad-hoc / --list runs

export RNAFOLD_BIN ESL_SHUFFLE_BIN PATH \
       TIER_SHUFFLES MIN_VALID_PERCENT N_SHUFFLES \
       WORKER_SCRIPT COLLATE_SCRIPT TOOL_NAME

# --- Calibration hooks ---

# Run the worker in MFE-only mode (~1000x faster than full shuffles).
# We extrapolate to full work via predict_wall_s.
calib_params() {
    echo "N_SHUFFLES=0"
}

# t_full ~= t_mfe * (N + 1) -- esl-shuffle is O(L) per shuffle, RNAfold's
# folding is O(L^3), so shuffle generation is rounding error. The +1
# accounts for the original sequence.
predict_wall_s() {
    local measured=$1 tier=$2
    read -ra arr <<< "$TIER_SHUFFLES"
    local n=${arr[$((tier-1))]}
    awk -v m="$measured" -v n="$n" 'BEGIN { printf "%d", (m + 1) * (n + 1) + 1 }'
}

# Streaming pipe — peak RSS is per-sequence, doesn't scale with N.
predict_rss_mb() {
    echo "$1"
}

# --- Per-tier sbatch export ---
# 04_submit.sh joins this comma-onto its --export list.
tier_params() {
    local tier=$1
    read -ra arr <<< "$TIER_SHUFFLES"
    echo "N_SHUFFLES=${arr[$((tier-1))]}"
}