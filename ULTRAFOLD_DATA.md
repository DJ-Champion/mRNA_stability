# Ultrafolded-transcript data

Standalone data-gathering steps for studying transcripts with extreme negative
RNAfold z-scores ("ultrafolded", default z < −30 on the spliced mRNA). They are
not part of the numbered `01–05` pipeline; each reads existing run outputs and
writes under `runs/<dataset>/ultrafold/`. Plotting happens elsewhere
(`mRNA_stability_R`).

All transcript-space coordinates are **0-based, half-open, 5′→3′ on the spliced
transcript**: the same frame as `extracted_mRNA.fa`. Join on `transcript_id`.

Run order (`RUNS_ROOT` defaults to `<repo>/runs`):

    ./bin/ultrafold_candidates.py -d human_all --scores <rnafold_results.csv>
    ./bin/fetch_rmsk.sh                          # once; dataset-independent
    ./bin/repeats_on_transcripts.py -d human_all
    ./bin/icshape_on_transcripts.py -d human_all
    ./bin/ultrafold_gene_lists.py -d human_all --scores <rnafold_results.csv>

## `ultrafold_candidates.py` → `candidates.tsv`

One row per transcript whose `--region` (default `MRNA`) z-score is below
`--threshold` (default −30), most negative first.

| column | meaning |
|---|---|
| `transcript_id`, `gene_id`, `gene_name`, `chrom`, `strand` | identifiers / locus |
| `mrna_len`, `utr5_len`, `cds_len`, `utr3_len` | region lengths from `manifest.tsv` (the region boundaries every later step uses) |
| `mrna_score`, `mrna_zscore` | RNAfold MFE and z-score, spliced mRNA |
| `utr3_score`, `utr3_zscore` | same for the 3′UTR alone (`NA` if absent) — context, not a filter |

The scores CSV labels the full transcript `MRNA` (not `mRNA`).

## `fetch_rmsk.sh` → `runs/_resources/hg38_rmsk/rmsk.bed.gz`

UCSC hg38 `rmsk` table as BED6+: `chrom start end repName swScore strand
repClass repFamily pctDiv`. Primary chromosomes only. `rmsk.provenance.txt`
records source, size check and fetch date. UCSC publishes no checksum for this
directory, so integrity is the server `Content-Length` plus `gzip -t`.

## `repeats_on_transcripts.py`

Keeps only **exonic** repeat bases. A repeat spanning a splice junction becomes
one row per exon (shared `repeat_id`); pieces shorter than `--min-piece`
(default 10) are dropped.

`repeats.tsv` — one row per exonic repeat piece, all classes. Key columns:
`tx_start`, `tx_end`, `region` (5UTR/CDS/3UTR by majority overlap),
`orientation` (`sense`/`antisense` = repeat strand vs transcript strand),
`rep_name`, `rep_class`, `rep_family`, `pct_div`, `repeat_id`, `piece`,
`n_pieces`, genomic coordinates, `repeat_len`, `piece_len`.

`alu_pairs.tsv` — every oppositely oriented pair of Alus (`rep_family == Alu`)
lying in `--pair-region` (default 3′UTR); `a` is 5′ of `b`; `gap` is the
distance between them (0 = overlapping). Candidate base-pairing partners by
orientation only; it does not say they pair in the predicted structure.

`alu_summary.tsv` — per candidate: `n_alu`, `n_sense`, `n_antisense`, `alu_bp`,
`alu_frac`, `n_opposite_pairs`, `min_pair_gap`, `n_alu_elsewhere`
(Alus outside the pair region).

Checked against real data: for 443 full-length Alu hits, no sense element
contains the reverse-complement Alu motif and no antisense element contains the
forward one.

## `icshape_on_transcripts.py`

Projects the genome-coordinate icSHAPE score BEDs (HEK293, cytoplasm and
nucleoplasm; default location `~/Gardner_Lab/Gini-from-icSHAPE/human/`) onto the
spliced transcripts, strand-specifically.

`icshape_tracks.tsv` — one row per transcript position: `transcript_id`, `pos`,
`base`, `cytoplasm`, `nucleoplasm`. `NA` = not measured (never 0).

`icshape_summary.tsv` — per transcript × compartment × region (5UTR / CDS /
3UTR / mRNA): `length`, `n_valid`, `coverage`, `mean`, `gini` (NA when
`n_valid` <= 15). Validated against Gini-from-icSHAPE's `gini_long_memfix.tsv`:
all 41 comparable rows match on both `n_valid` and Gini.

**Coverage warning.** For the 45 candidates, icSHAPE covers at most 16% of any
3′UTR, and only 529 of 115,569 Alu nucleotides (0.5%); 2 transcripts have >= 50
valid nucleotides inside Alus. The score files carry no read depth, so coverage
here means "has a score".

## `ultrafold_gene_lists.py` → `gene_lists/`

Inputs for GO enrichment in R. `universe.tsv` is the background (every gene with
a z-score on the chosen region; use it as the enrichment universe, not the whole
genome). `ranked.tsv` is the same set most-ultrafolded first with
`rank_metric = -mrna_zscore` for ranked tests. `top_z<T>.txt` / `.symbols.txt`
list genes with z < T (T = -30, -25, -20 by default; file tags `m30` etc.).
Ultrafolded genes have very long 3′UTRs, so compare against a length-matched
control (`mrna_len`, `utr3_len` are in `universe.tsv`).
