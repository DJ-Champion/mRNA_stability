#!/bin/bash
# configs/tools/rnacofold.sh
# RNAcofold tool config. Operates ONLY on UTR_pair regions: the worker reads
# the UTR_pair 3-line FASTA, parses 5UTR/3UTR boundaries from the constraint
# string, and re-formats the input as '5UTR&3UTR' for RNAcofold.
#
# RNAcofold computes the MFE of a heterodimer between two RNA strands
# directly — no linker, no constraints, no forced pairing. This is the
# methodologically correct way to ask "do these two UTRs cross-pair?"
# (cf. the constrained-RNAfold approach in configs/tools/rnafold.sh, which
# is retained for comparison but forces every base in both UTRs to pair).

TOOL_NAME=rnacofold
WORKER_SCRIPT="$PROJECT_ROOT/tools/rnacofold/worker.sh"
COLLATE_SCRIPT="$PROJECT_ROOT/tools/rnacofold/collate.py"

# --- Region restriction ---
# RNAcofold only makes sense on paired-strand input. 03_calibrate.sh and
# 04_submit.sh consult TOOL_REGIONS and filter samples / submission lists.
TOOL_REGIONS="UTR_pair"

# --- Binaries ---
: "${RNACOFOLD_BIN:=/home/chado47p/Software/ViennaRNA/bin/RNAcofold}"
: "${ESL_SHUFFLE_BIN:=/home/chado47p/Software/hmmer3.4/bin/esl-shuffle}"

PATH="$(dirname "$RNACOFOLD_BIN"):$(dirname "$ESL_SHUFFLE_BIN"):$PATH"

# --- Per-tier shuffle policy ---
# Methodology choice: original constrained_zscore_for_slurm.sh used N=100.
# Cofold is computationally cheaper than constrained-fold (no constraint
# overhead, simpler thermodynamic model), so we can match the larger N
# used for unconstrained RNAfold. Adjust here if a different null is wanted.
: "${TIER_SHUFFLES:=1000 1000 1000 1000 1000 1000 1000 100 10 10}"

# --- Worker parameters ---
: "${MIN_VALID_PERCENT:=90}"
: "${N_SHUFFLES:=1000}"

export RNACOFOLD_BIN ESL_SHUFFLE_BIN PATH \
       TIER_SHUFFLES MIN_VALID_PERCENT N_SHUFFLES \
       WORKER_SCRIPT COLLATE_SCRIPT TOOL_NAME TOOL_REGIONS

# --- Calibration hooks (mirror RNAfold; same extrapolation rationale) ---

# MFE-only fast path for calibration: skip the shuffle null.
calib_params() {
    echo "N_SHUFFLES=0"
}

# t_full ~= t_mfe * (N + 1). Per-shuffle work is dominated by RNAcofold
# itself; esl-shuffle on each strand is rounding error.
predict_wall_s() {
    local measured=$1 tier=$2
    read -ra arr <<< "$TIER_SHUFFLES"
    local n=${arr[$((tier-1))]}
    awk -v m="$measured" -v n="$n" 'BEGIN { printf "%d", (m + 1) * (n + 1) + 1 }'
}

# Streaming over shuffles — peak RSS set by sequence length, not by N.
predict_rss_mb() {
    echo "$1"
}

# --- Per-tier sbatch export ---
tier_params() {
    local tier=$1
    read -ra arr <<< "$TIER_SHUFFLES"
    echo "N_SHUFFLES=${arr[$((tier-1))]}"
}
