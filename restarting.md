# Cancellation safety  

The output CSV is only written after the shuffle pipeline completes and passes the min_valid threshold. Walk the worker:

1. Run shuffles → stats_file + raw_tmp in $TMP_DIR
2. Count valid shuffles; abort if < min_valid (no results CSV written)
3. mv raw_tmp raw_file to raw_shuffles/
4. awk summarises into $output_file

If scancel lands at step 1 or 2, no results CSV exists → next run does it cleanly. If it lands at step 4, the awk's first action is print "header" > out so a partially-killed step-4 leaves a single-line file, which fails the > 1 check → also gets re-run cleanly.
The narrow danger window is "killed between the header write and the data write inside the awk END block" — the data write is a single printf, so practically zero. I wouldn't worry about it.
Stale tmp/raw cruft: cancelled runs can leave files in $TMP_DIR (stats_*.$$, raw_*.$$) and incomplete files in raw_shuffles/. They don't break anything (raw_shuffles gets overwritten on re-run, tmp uses PID-suffixed names), but it's untidy.  

## The real gotcha: changed N_SHUFFLES

The skip check looks at file presence, not at what parameters were used. The output CSV does record n_shuffles_requested in its last column, but the worker doesn't read that.
So if you previously ran tier 8 with N_SHUFFLES=100, and you've since edited TIER_SHUFFLES in rnafold.sh to use 1000 for tier 8, then resubmitting will silently skip every tier-8 sequence that already has a 100-shuffle result. You'll end up with a mixed dataset where some rows are 100-shuffle and others are 1000-shuffle, with nothing flagging it.  

If that's a risk, this one-liner finds them:

``` awk -F, 'NR==1 {next} $8 != EXPECTED' EXPECTED=1000 \
    runs/<dataset>/rnafold/results/*.csv
```

Or for a full audit:
``` awk -F, 'FNR==2 {print $8, FILENAME}' runs/<dataset>/rnafold/results/*.csv \
    | sort | uniq -c
```  

### Cancel and clean  

First, kill what's queued:
``` scancel -u $USER --name=mrna_rnafold_t1
# ...or nuke everything matching the pattern:
scancel -u $USER --jobname=$(squeue -u $USER -h -o "%j" | grep mrna_rnafold | sort -u | paste -sd,)
# simpler if you don't have other jobs running:
scancel -u $USER
```  

Then clean. Three levels depending on how thorough you want to be:
Re-run this tool from scratch, keep extraction + tiering + calibration

``` rm -rf runs/<dataset>/rnafold/{results,errors,tmp,raw_shuffles,slurm_logs}
rm -f  runs/<dataset>/rnafold/combined.tsv
```  

This is the usual one. Tiering (lists/) and extracted FASTAs are tool-agnostic; calibration is expensive to redo and rarely changes meaningfully.

Re-run including recalibration (e.g. you upgraded ViennaRNA)
```
rm -rf runs/<dataset>/rnafold/
```  

Full reset (e.g. you changed gene list, re-extracted, re-tiered)
```  
rm -rf runs/<dataset>/
```  
The middle ground you may actually want
If only some sequences need redoing — the N_SHUFFLES mismatch case, or sequences whose error logs suggest borderline results — don't nuke everything. Use the resume pattern from the README:
``` 
# Generate list of seqnames whose CSVs aren't at full N_SHUFFLES
awk -F, 'FNR==2 && $8 < 1000 {print $1}' runs/<dataset>/rnafold/results/*.csv \
    | while read seq; do
        # find the FASTA path for this seqname in any tier list
        grep -lh "/${seq}\.fa$" runs/<dataset>/lists/tier_*.txt \
            | xargs -I{} grep "/${seq}\.fa$" {}
      done > /tmp/redo.txt

# Delete the stale results so the skip check doesn't trigger
while read fasta; do
    seq=$(basename "$fasta" .fa)
    rm -f "runs/<dataset>/rnafold/results/results_${seq}.csv"
done < /tmp/redo.txt

# Resubmit only those
./bin/04_submit.sh -d <dataset> -t rnafold --list /tmp/redo.txt \
    SLURM_TIME=02:00:00 SLURM_MEM_PER_CPU=500MB
```  

Faster than redoing everything, surgical, and --list mode bypasses the calibration loop entirely.