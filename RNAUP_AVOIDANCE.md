# RNAup mRNA:ncRNA interaction / avoidance

Per mRNA window: RNAup binding energies against 34 abundant human ncRNAs (6 families), summed, with a
dinucleotide-shuffle robust Z. Tool `rnaup_avoidance`; code in `tools/rnaup_avoidance/`,
`configs/tools/rnaup_avoidance.sh`.

## Run

Same flow as the other tools, one window FASTA per task:

    ./bin/01_extract.py -d human_all           # needs start_codon_region (the only window scored)
    ./bin/02_stratify.sh -d human_all
    export SLURM_CPUS_PER_TASK=8               # same value for calibrate and submit; any count works
    ./bin/03_calibrate.sh -d human_all -t rnaup_avoidance      # on a compute node
    ./bin/04_submit.sh -d human_all -t rnaup_avoidance
    ./bin/05_collate.sh -d human_all -t rnaup_avoidance

Timing test on 10 genes: `configs/datasets/human_10.yaml` (put the IDs in `data/human/gene_ids_human_10.txt`), then the
same steps with `-d human_10`. `N_SHUFFLES` defaults to 100 (`N_SHUFFLES=1000` for a finer null). Calibration runs the real code with
`AVOID_CALIB_SHUFFLES` (default 20) shuffles and extrapolates by (N+1)/(20+1); that includes the fixed per-run cost, so it
errs high. 03_calibrate rounds the measured wall to whole seconds before extrapolating, so a run of only a few seconds
(as with 2 shuffles on 8 cores) can be off by +-50%; 20 keeps that to a few percent. It writes nothing to `results/`.

Local timing (RNAup 2.4.7, one core, random 63-nt windows): cost per target grows about linearly with ncRNA length,
~0.06 s for a 72-nt ncRNA and ~0.37 s for the 341-nt one (20 of the 34 are 72-84 nt), roughly 3-6 s for all 34. So one gene at
N=100 is roughly 5-10 core-minutes (1-1.6 core-hours at N=1000; 200-nt windows cost ~4x more per target). Treat these as estimates until the 10-gene run.

## Window

Only the start codon +-30 nt is scored: manifest region `start_codon_region`, output `window` `start_codon_pm30`,
30 nt upstream + AUG + 30 nt downstream = 63 nt. A transcript with a 5'UTR shorter than 30 nt gives a clipped window
and is dropped. `TOOL_REGIONS=start_codon_region` and `TOOL_REQUIRE_LEN=start_codon_region:63` are fixed in the tool
config (env overrides are ignored); `04_submit.sh` applies them to the manifest when it builds each tier's task list,
and the worker refuses any other region or length, so nothing else reaches RNAup or esl-shuffle.

## Method

A gene's native window and N=100 (`N_SHUFFLES`) `esl-shuffle -d` shuffles (whole window, nothing held fixed) are cut
into chunks of `AVOID_CHUNK` (101) targets, so at N=100 each ncRNA sees all 101 in one RNAup stream. A work
unit is one ncRNA against one chunk: a stream with the ncRNA first, then the chunk, `RNAup --interaction_first --no_output_file -b -d2 --noLP`
(override with `RNAUP_ARGS`). From `(G_total = G_int + G_u_target + G_u_query)` we keep `G_total` and `G_int`;
a pair with no interaction line is 0/0 and stays in all sums. A target header missing from RNAup's output,
or an energy line that does not parse, fails the gene (truncation / unexpected format) instead of scoring as 0.

Units run `SLURM_CPUS_PER_TASK` at a time, longest ncRNA first, so cores stay busy despite the 72-341 nt spread
(34 units at N=100, so more than 34 cores is wasted).
Checkpoints in `<tool_dir>/ckpt/<seqname>/`: the shuffled sequences (written once) and one gz per unit. A restarted
task reuses both, so all 34 ncRNAs see the same shuffles. If the stored shuffles do not match the input (window
re-extracted, `N_SHUFFLES` changed) the checkpoint is discarded and rebuilt. `AVOID_KEEP_CKPT=1` keeps them.
Progress (`k/n units done`) streams to the SLURM log while the task runs.
The results file is written last; its existence means the gene is done.

The avoidance score is `avoidance_G_int_*`: `G_int`, RNAup's duplex (hybridisation) term of the predicted
interaction, without the opening penalties. `rnaup_G_total_*` (duplex + opening) is kept only for reference.
Energies in kcal/mol; more negative = stronger binding. Robust Z =
`(native - median(shuffles)) / (1.4826 * MAD(shuffles))`, `NA` when MAD = 0 (same as the RNAfold worker). The null is
shuffles of the window itself, not of the whole mRNA as in the RNAfold run, so it asks "more/less binding than
expected for this window's dinucleotide composition". Only the mRNA windows are shuffled; the ncRNAs are used as given. Single-ncRNA families (rRNA, SRP) often have MAD 0 -> `NA`.

## Output

`runs/<dataset>/rnaup_avoidance/`: `results/results_<seqname>.csv`, `combined.tsv` (joined with
the manifest on `seqname`, via the rnalfold collator), `avoidance_<window>.tsv` (one per window),
`raw_shuffles/<seqname>_shuffle_sums.csv.gz` (per-shuffle sums over the 34).

| columns | meaning |
|---|---|
| `window`, `window_length`, `n_shuffles`, `n_ncrna` | geometry (`region`/`length` are the manifest's source region) |
| `rnaup_G_total_{sum,mean,min}`, `avoidance_G_int_{sum,mean,min}` | native, over the 34 ncRNAs |
| `rnaup_G_total_sum_<F>`, `avoidance_G_int_sum_<F>` | per family `F` in tRNA rRNA SRP snRNA Y_RNA snoRNA |
| `*_sum_robust_z`, `*_sum_<F>_robust_z` | robust Z of the sum |
| `*_sum_shuf_median`, `*_sum_shuf_mad` | shuffle median and raw MAD (not scaled) of the total sum |
| `gc_content_mRNA_window`, `gc_ncRNA_mean`, `gc_diff_mRNA_vs_ncRNA` | GC covariate; ncRNA mean is the mean of the 34 per-ncRNA GCs |

## RNAup behaviour seen with 2.4.7

Positive `G_total` is reported as is (an unfavourable interaction, e.g. `0.25 = -37.50 + 20.73 + 17.02`); only
RNAup's `100.00 = ... + inf` "cannot score" lines are set to 0/0. `-w` defaults to 25 (maximum interaction length).
