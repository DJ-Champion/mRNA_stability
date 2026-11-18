#!/bin/bash
# configs/tools/rnaup_avoidance.sh
# RNAup interaction / avoidance of mRNA windows against a curated set of 34
# abundant human ncRNAs (6 families), with a dinucleotide-shuffle robust Z.
# Method and columns: RNAUP_AVOIDANCE.md.
#
# Not tiered and not calibrated: windows have a fixed length, and the unit of
# work is a TRANCHE (M genes, one multi-record FASTA), not one sequence. Build
# the tranches, then submit them with the usual list mode:
#
#   ./bin/rnaup_make_tranches.py -d human_all                  # start_codon_pm30
#   ./bin/rnaup_make_tranches.py -d human_all --regions start_codon_pm30 stop_codon_pm30 last_100nt
#   SLURM_TIME=06:00:00 SLURM_MEM_PER_CPU=2000MB \
#       ./bin/04_submit.sh -d human_all -t rnaup_avoidance \
#           --list runs/human_all/lists/rnaup_start_codon_pm30.txt
#   ./bin/05_collate.sh -d human_all -t rnaup_avoidance
#
# Time per tranche ~ M * (N+1) * 34 RNAup target evaluations. Run one tranche
# with N_SHUFFLES=5 first and scale from that; predict_wall_s is a placeholder.
#
# Shell/env overrides win, e.g. N_SHUFFLES=100, RNAUP_ARGS="... -w 40".

TOOL_NAME=rnaup_avoidance
WORKER_SCRIPT="$PROJECT_ROOT/tools/rnaup_avoidance/worker.sh"
COLLATE_SCRIPT="$PROJECT_ROOT/tools/rnaup_avoidance/collate.py"

: "${RNAUP_BIN:=/home/chado47p/Software/ViennaRNA/bin/RNAup}"
: "${ESL_SHUFFLE_BIN:=/home/chado47p/Software/hmmer3.4/bin/esl-shuffle}"
PATH="$(dirname "$RNAUP_BIN"):$(dirname "$ESL_SHUFFLE_BIN"):$PATH"

# The 34-ncRNA FASTA is an input, not shipped in the repo. Family map: configs/ncrna/families.tsv.
: "${NCRNA_FASTA:=$PROJECT_ROOT/data/ncrna/curated-interaction-database.fa}"
: "${NCRNA_FAMILIES:=$PROJECT_ROOT/configs/ncrna/families.tsv}"

: "${RNAUP_ARGS:=--interaction_first --no_output_file -b -d2 --noLP}"
: "${N_SHUFFLES:=1000}"
: "${MIN_VALID_PERCENT:=90}"
: "${AVOID_KEEP_CKPT:=0}"

export TOOL_NAME WORKER_SCRIPT COLLATE_SCRIPT PATH RNAUP_BIN ESL_SHUFFLE_BIN \
       NCRNA_FASTA NCRNA_FAMILIES RNAUP_ARGS N_SHUFFLES MIN_VALID_PERCENT AVOID_KEEP_CKPT

tier_params()  { echo "N_SHUFFLES=$N_SHUFFLES"; }
calib_params() { echo "N_SHUFFLES=0"; }
