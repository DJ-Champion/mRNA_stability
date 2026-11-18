#!/bin/bash
# configs/tools/rnaup_avoidance.sh
# RNAup interaction / avoidance of mRNA windows against a curated set of 34
# abundant human ncRNAs (6 families), with a dinucleotide-shuffle robust Z.
# Method and columns: RNAUP_AVOIDANCE.md.
#
# Per-sequence like the other tools: one window FASTA per task, the standard
#   02_stratify -> 03_calibrate -> 04_submit -> 05_collate
# flow. Choose the window(s) with TOOL_REGIONS (a manifest region, default the
# primary start-codon window); several are fine:
#
#   TOOL_REGIONS="start_codon_region stop_codon_region" ./bin/04_submit.sh ...
#
#   manifest region      output window      kept only if (TOOL_REQUIRE_LEN)
#   start_codon_region   start_codon_pm30   exactly 63 nt  (5'UTR >= 30)
#   stop_codon_region    stop_codon_pm30    exactly 63 nt  (3'UTR >= 30)
#   tail_region          last_100nt         exactly 100 nt
#   5UTR                 utr5               >= 20 nt (whole 5'UTR, length varies)
#
# 04_submit applies TOOL_REQUIRE_LEN to the task list and the worker re-checks it, so
# no clipped window reaches RNAup. Shuffles are NOT seeded (as in the RNAfold tools);
# the worker persists them per gene so a restarted task reuses the same ones.
#
# A task works through (ncRNA, chunk of AVOID_CHUNK=100 targets) units on SLURM_CPUS_PER_TASK
# cores, longest ncRNA first, so any core count stays busy (8 is a good start; 34 is not needed).
# Set SLURM_CPUS_PER_TASK to the same value for BOTH 03_calibrate and 04_submit. RNAup cost grows
# about linearly with ncRNA length: ~0.06 s per 63-nt target for a 72-nt ncRNA, ~0.37 s for the
# 341-nt one, roughly 3-6 s for all 34 (RNAup 2.4.7, one core), so ~1-1.6 core-hours per gene at N=1000.
# Longer windows cost more (200 nt ~4x per target). Cap very long 5'UTRs with
# TOOL_REQUIRE_LEN="... 5UTR:20-1000".
#
# Shell/env overrides win, e.g. N_SHUFFLES=100, RNAUP_ARGS="... -w 40".

TOOL_NAME=rnaup_avoidance
WORKER_SCRIPT="$PROJECT_ROOT/tools/rnaup_avoidance/worker.sh"
COLLATE_SCRIPT="$PROJECT_ROOT/tools/rnaup_avoidance/collate.py"

: "${RNAUP_BIN:=/home/chado47p/Software/ViennaRNA/bin/RNAup}"
: "${ESL_SHUFFLE_BIN:=/home/chado47p/Software/hmmer3.4/bin/esl-shuffle}"
PATH="$(dirname "$RNAUP_BIN"):$(dirname "$ESL_SHUFFLE_BIN"):$PATH"

# The 34 curated ncRNAs and their family map. Family map: configs/ncrna/families.tsv.
: "${NCRNA_FASTA:=$PROJECT_ROOT/data/ncrna/curated-interaction-database.fa}"
: "${NCRNA_FAMILIES:=$PROJECT_ROOT/configs/ncrna/families.tsv}"

: "${RNAUP_ARGS:=--interaction_first --no_output_file -b -d2 --noLP}"
: "${N_SHUFFLES:=1000}"
: "${MIN_VALID_PERCENT:=90}"
: "${AVOID_KEEP_CKPT:=0}"
: "${AVOID_CHUNK:=100}"

export TOOL_NAME WORKER_SCRIPT COLLATE_SCRIPT PATH RNAUP_BIN ESL_SHUFFLE_BIN \
       NCRNA_FASTA NCRNA_FAMILIES RNAUP_ARGS N_SHUFFLES MIN_VALID_PERCENT AVOID_KEEP_CKPT AVOID_CHUNK

: "${TOOL_REGIONS:=start_codon_region}"
: "${TOOL_REQUIRE_LEN:=start_codon_region:63 stop_codon_region:63 tail_region:100 5UTR:>=20}"
: "${AVOID_CALIB_SHUFFLES:=2}"

export TOOL_REGIONS TOOL_REQUIRE_LEN AVOID_CALIB_SHUFFLES

tier_params()  { echo "N_SHUFFLES=$N_SHUFFLES"; }
# Calibration runs the full code path (all 34 ncRNAs) on the native window plus a few shuffles.
calib_params() { echo "N_SHUFFLES=$AVOID_CALIB_SHUFFLES,AVOID_CALIB=1"; }
# t_full ~= t_calib * (N+1) / (calib+1); the fixed part (python, one ncRNA fold each) is
# scaled too, so this errs high. Calibrate with the same SLURM_CPUS_PER_TASK as the submit.
predict_wall_s() { awk -v m="$1" -v n="$N_SHUFFLES" -v c="$AVOID_CALIB_SHUFFLES" 'BEGIN { printf "%d", m * (n + 1) / (c + 1) + 1 }'; }
# /usr/bin/time reports the largest single process, which is what --mem-per-cpu is for.
predict_rss_mb() { echo "$1"; }
