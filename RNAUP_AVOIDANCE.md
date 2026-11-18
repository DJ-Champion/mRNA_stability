# RNAup mRNA:ncRNA interaction / avoidance

Per mRNA window: RNAup binding energies against 34 abundant human ncRNAs (6 families), summed, with a
dinucleotide-shuffle robust Z. Tool `rnaup_avoidance`; code in `tools/rnaup_avoidance/`,
`configs/tools/rnaup_avoidance.sh`.

## Run

Same flow as the other tools, one window FASTA per task:

    ./bin/01_extract.py -d human_all           # needs start_codon_region (+ stop_codon_region, tail_region, 5UTR if used)
    ./bin/02_stratify.sh -d human_all
    export SLURM_CPUS_PER_TASK=34              # same value for calibrate and submit: a task runs the 34 ncRNAs in parallel
    ./bin/03_calibrate.sh -d human_all -t rnaup_avoidance      # on a compute node
    ./bin/04_submit.sh -d human_all -t rnaup_avoidance
    ./bin/05_collate.sh -d human_all -t rnaup_avoidance

Pick windows with `TOOL_REGIONS` (default `start_codon_region`), set the same way for calibrate and submit:

    TOOL_REGIONS="start_codon_region stop_codon_region tail_region 5UTR" ./bin/04_submit.sh -d human_all -t rnaup_avoidance

Timing test on 10 genes: `configs/datasets/human_10.yaml` (put the IDs in `data/human/gene_ids_human_10.txt`), then the
same steps with `-d human_10`; `N_SHUFFLES=100` for a first look. Calibration runs the real code with
`AVOID_CALIB_SHUFFLES` (default 2) shuffles and extrapolates by (N+1)/(2+1); that includes the fixed per-run cost, so it
errs high. It writes nothing to `results/`.

Local timing (RNAup 2.4.7, 341-nt ncRNA, random windows): about 0.4 s per 63-nt target and 1.6 s per 200-nt target per
core, so one 63-nt gene at N=1000 is roughly 3.9 core-hours (34 ncRNAs x 1001 sequences).

## Windows

| manifest region (`TOOL_REGIONS`) | output `window` | kept only if (`TOOL_REQUIRE_LEN`) |
|---|---|---|
| `start_codon_region` | `start_codon_pm30` | exactly 63 nt (5'UTR >= 30) |
| `stop_codon_region` | `stop_codon_pm30` | exactly 63 nt (3'UTR >= 30) |
| `tail_region` | `last_100nt` | exactly 100 nt |
| `5UTR` | `utr5` | >= 20 nt; whole 5'UTR, length varies |

`04_submit.sh` applies `TOOL_REQUIRE_LEN` to the manifest when it builds each tier's task list, and the worker
re-checks it, so nothing clipped reaches RNAup. Rule syntax: `region:63` exact, `region:>=20` minimum,
`region:20-1000` range. Cost grows steeply with window length; consider capping 5UTR.

## Method

Per gene, per ncRNA (34 RNAup processes, `SLURM_CPUS_PER_TASK` at a time): one stream, ncRNA first, then every
gene's native window and N=1000 `esl-shuffle -d` shuffles, `RNAup --interaction_first --no_output_file -b -d2 --noLP`
(override with `RNAUP_ARGS`). From `(G_total = G_int + G_u_target + G_u_query)` we keep `G_total` and `G_int`;
a pair with no interaction line is 0/0 and stays in all sums. A target header missing from RNAup's output
fails the gene (truncation / unexpected format) instead of scoring as 0.

Checkpoints in `<tool_dir>/ckpt/<seqname>/`: the shuffled sequences (written once) and one gz per ncRNA. A restarted
task reuses both, so all 34 ncRNAs see the same shuffles. `AVOID_KEEP_CKPT=1` keeps them.
The results file is written last; its existence means the gene is done.

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
