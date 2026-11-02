#!/bin/bash
# configs/tools/rnafold_shuffles.sh
# RNAfold with a fixed, larger shuffle count, for a hand-picked set of
# sequences (e.g. the ultrafolded transcripts in the shuffle-distribution
# figure). Reuses tools/rnafold/worker.sh unchanged; only the output location
# and the shuffle policy differ:
#
#   runs/<dataset>/rnafold_shuffles/results/results_<seqname>.csv
#   runs/<dataset>/rnafold_shuffles/raw_shuffles/<seqname>_raw_shuffles.csv
#
# so the original runs/<dataset>/rnafold/ results are left alone (the worker
# skips any sequence whose results file already exists).
#
# Not tiered: submit with an explicit list, e.g.
#   N_SHUFFLES=500 SLURM_TIME=08:00:00 ./bin/04_submit.sh -d human_all \
#       -t rnafold_shuffles --list runs/human_all/ultrafold/shuffle_set.fasta_list.txt
#
# Shuffles are NOT seeded (esl-shuffle default), so a rerun gives a different
# sample; keep the raw_shuffles files from the run you plot.

TOOL_NAME=rnafold_shuffles
WORKER_SCRIPT="$PROJECT_ROOT/tools/rnafold/worker.sh"
COLLATE_SCRIPT="$PROJECT_ROOT/tools/rnafold/collate.py"

: "${RNAFOLD_BIN:=/home/chado47p/Software/ViennaRNA/bin/RNAfold}"
: "${ESL_SHUFFLE_BIN:=/home/chado47p/Software/hmmer3.4/bin/esl-shuffle}"
PATH="$(dirname "$RNAFOLD_BIN"):$(dirname "$ESL_SHUFFLE_BIN"):$PATH"

: "${N_SHUFFLES:=500}"
: "${MIN_VALID_PERCENT:=90}"

export RNAFOLD_BIN ESL_SHUFFLE_BIN PATH N_SHUFFLES MIN_VALID_PERCENT \
       WORKER_SCRIPT COLLATE_SCRIPT TOOL_NAME

tier_params()  { echo "N_SHUFFLES=$N_SHUFFLES"; }
calib_params() { echo "N_SHUFFLES=0"; }
# t_full ~= t_mfe * (N + 1)
predict_wall_s() { awk -v m="$1" -v n="$N_SHUFFLES" 'BEGIN { printf "%d", (m + 1) * (n + 1) + 1 }'; }
predict_rss_mb() { echo "$1"; }
