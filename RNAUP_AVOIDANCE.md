# RNAup mRNA:ncRNA interaction / avoidance

Per mRNA window: RNAup binding energies against 34 abundant human ncRNAs (6 families), summed, with a
dinucleotide-shuffle robust Z. Tool `rnaup_avoidance`; code in `tools/rnaup_avoidance/`,
`bin/rnaup_make_tranches.py`, `configs/tools/rnaup_avoidance.sh`.

## Run

    # inputs: data/ncrna/curated-interaction-database.fa (not in the repo), configs/ncrna/families.tsv (in the repo)
    ./bin/rnaup_make_tranches.py -d human_all                       # start_codon_pm30; -M genes per tranche (10)
    ./bin/rnaup_make_tranches.py -d human_all --regions stop_codon_pm30 last_100nt utr5
    SLURM_TIME=06:00:00 SLURM_MEM_PER_CPU=2000MB \
        ./bin/04_submit.sh -d human_all -t rnaup_avoidance --list runs/human_all/lists/rnaup_start_codon_pm30.txt
    ./bin/05_collate.sh -d human_all -t rnaup_avoidance

Cluster timing test on 10 genes (`configs/datasets/human_10.yaml`; put the gene IDs in `data/human/gene_ids_human_10.txt`):

    ./bin/01_extract.py -d human_10
    ./bin/rnaup_make_tranches.py -d human_10 --regions start_codon_pm30 utr5 -M 10
    SLURM_CPUS_PER_TASK=34 SLURM_TIME=24:00:00 SLURM_MEM_PER_CPU=1000MB \
        ./bin/04_submit.sh -d human_10 -t rnaup_avoidance --list runs/human_10/lists/rnaup_start_codon_pm30.txt

Local timing (RNAup 2.4.7, 341-nt ncRNA, random windows): about 0.4 s per 63-nt target and 1.6 s per 200-nt target
per core, so one 10-gene tranche at N=1000 is roughly 39 core-hours for 63-nt windows. Use `N_SHUFFLES=100` for a first look.

Try one tranche first with `N_SHUFFLES=5` and read the timing and the first rows before the full run.
Needs `start_codon_region` (and `stop_codon_region`, `tail_region`, `5UTR` for the other windows) in the
dataset's `regions_to_extract`. Calibration/tiers are not used (fixed-length windows).

## Windows

| window | source region | length | keeps only |
|---|---|---|---|
| `start_codon_pm30` | `start_codon_region` | 63 | 5'UTR >= 30 |
| `stop_codon_pm30` | `stop_codon_region` | 63 | 3'UTR >= 30 |
| `last_100nt` | `tail_region` | 100 | mRNA >= 100 |
| `utr5` | `5UTR`, whole | varies | 5'UTR >= 20 (`--min-utr5-length`) |

Fixed windows are filtered on exact length, so nothing truncated reaches RNAup. `utr5` is variable-length: tranches are sorted by length and also capped at 1500 nt total (`--max-tranche-nt`), and cost grows with UTR length. Assumes `codon_flank_length: 30`,
`tail_length: 100` in the dataset YAML.

## Method

Per tranche, per ncRNA (34 RNAup processes, `SLURM_CPUS_PER_TASK` at a time): one stream, ncRNA first, then every
gene's native window and N=1000 `esl-shuffle -d` shuffles, `RNAup --interaction_first --no_output_file -b -d2 --noLP`
(override with `RNAUP_ARGS`). From `(G_total = G_int + G_u_target + G_u_query)` we keep `G_total` and `G_int`;
a pair with no interaction line is 0/0 and stays in all sums. A target header missing from RNAup's output
fails the tranche (truncation / unexpected format) instead of scoring as 0.

Checkpoints in `<tool_dir>/ckpt/<tranche>/`: the shuffled sequences (written once) and one gz per ncRNA. A restarted
task reuses both, so all ncRNAs in a tranche see the same shuffles. `AVOID_KEEP_CKPT=1` keeps them.
The results file is written last; its existence means the tranche is done.

Energies in kcal/mol; more negative = stronger binding. Robust Z =
`(native - median(shuffles)) / (1.4826 * MAD(shuffles))`, `NA` when MAD = 0 (same as the RNAfold worker). The null is
shuffles of the window itself, not of the whole mRNA as in the RNAfold run, so it asks "more/less binding than
expected for this window's dinucleotide composition". Only the mRNA windows are shuffled; the ncRNAs are used as given. Single-ncRNA families (rRNA, SRP) often have MAD 0 -> `NA`.

## Output

`runs/<dataset>/rnaup_avoidance/`: `results/results_<tranche>.csv` (one row per gene), `combined.tsv` (joined with
the manifest on `seqname`, via the rnalfold collator), `avoidance_<window>.tsv` (one per window),
`raw_shuffles/<tranche>_shuffle_sums.csv.gz` (per-shuffle sums over the 34).

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
