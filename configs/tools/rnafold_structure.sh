#!/bin/bash
# configs/tools/rnafold_structure.sh
# RNAfold on the full sequences of a hand-picked set, keeping the MFE
# dot-bracket structure (and, with WITH_BPP=1, base-pair probabilities).
# No shuffles. Same binary and flags as the stored z-scores (tools/rnafold/
# worker.sh: RNAfold --noPS, default parameters), so the MFE must reproduce
# the stored one; bin/alu_structure_pairing.py checks that.
#
#   runs/<dataset>/rnafold_structure/results/<seqname>.fold      RNAfold stdout
#   runs/<dataset>/rnafold_structure/results/<seqname>.bpp.tsv   (WITH_BPP=1)
#
# Not tiered: submit with an explicit list, e.g.
#   ./bin/04_submit.sh -d ultrafold_shuffle_set -t rnafold_structure \
#       --list runs/ultrafold_shuffle_set/ultrafold/shuffle_set.fasta_list.txt
#   WITH_BPP=1 SLURM_TIME=12:00:00 SLURM_MEM=... ./bin/04_submit.sh ...   (optional)
#
# The worker skips a sequence only when its outputs for the requested mode
# exist, so calibration results (MFE only, same file) never block a real run.

TOOL_NAME=rnafold_structure
WORKER_SCRIPT="$PROJECT_ROOT/tools/rnafold_structure/worker.sh"
# Nothing to collate; 05_collate.sh requires a readable script, so point at a
# no-op placeholder.
COLLATE_SCRIPT="$PROJECT_ROOT/tools/rnafold_structure/collate.py"

: "${RNAFOLD_BIN:=/home/chado47p/Software/ViennaRNA/bin/RNAfold}"
PATH="$(dirname "$RNAFOLD_BIN"):$PATH"
: "${WITH_BPP:=0}"

export RNAFOLD_BIN PATH WITH_BPP WORKER_SCRIPT COLLATE_SCRIPT TOOL_NAME

tier_params()  { echo "WITH_BPP=$WITH_BPP"; }
calib_params() { echo "WITH_BPP=0"; }
# Calibration times the MFE fold. -p adds the partition function and the
# pair-probability matrix: the multipliers below are rough guesses
# (~3x time, ~4x memory), not measured. Check the first real run and adjust.
predict_wall_s() {
    awk -v m="$1" -v p="$WITH_BPP" 'BEGIN { printf "%d", (p ? 3 * m : m) + 1 }'
}
predict_rss_mb() {
    awk -v m="$1" -v p="$WITH_BPP" 'BEGIN { printf "%d", (p ? 4 * m : m) }'
}
