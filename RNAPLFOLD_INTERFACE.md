# rnaplfold output: interface for downstream ingestion

Source: `runs/<dataset>/rnaplfold/combined.tsv` (written by `05_collate.sh -t rnaplfold`).
Producer: mRNA_stability repo, `tools/rnaplfold/`. Tab-separated, one header row, `NA` for missing.

## What it measures

Local RNA accessibility of the **translation start region**: the probability that nucleotides are
unpaired, from RNAplfold run on the **full mRNA** (spliced transcript, 5'UTR + CDS + 3'UTR), then
summarised over a window around the start codon (AUG). Higher = more accessible (less structured).
One row per transcript. Only `mRNA` records appear (no per-region rows).

## Keys

| Column | Notes |
|---|---|
| `seqname` | `{gene_id}_{transcript_id}_mRNA`; unique per row |
| `gene_id`, `transcript_id` | Ensembl IDs, transcript ID includes version (`ENST00000009589.8`). Join on these. |

## Columns

"Stretch value" = probability that `l` consecutive nt are all unpaired (`l` = 1 is the per-nucleotide
unpaired probability). `u<l>` in a column name is that `l` (default 1, 4, 8, 16; ingest by pattern, the set
is configurable). A stretch counts toward a window/region only if it lies wholly inside it.

**Metadata and geometry**

| Column | Type | Meaning |
|---|---|---|
| `region` | str | always `mRNA` |
| `length` | int | mRNA length (nt); equals `seq_length` |
| `selection_reason` | str | isoform-choice rule (e.g. `MANE_Select`) |
| `short_utrs` | `true`/`false` | extraction flag: 5' or 3' UTR shorter than the extraction's `min_utr_length` |
| `start_win_start`, `start_win_end`, `start_win_len` | int | primary window, **1-based inclusive**, transcript coordinates; 63 normally (30 nt + AUG + 30 nt), shorter when clipped at a transcript end |
| `seq_length` | int | folded length |
| `seq_shorter_than_window` | 0/1 | 1 if `seq_length < rnaplfold_window`: whole molecule in one window, end artefacts; treat as unreliable |
| `cds_offset` | int | 5'UTR length = number of nt before the AUG |
| `aug_unpaired` | float [0,1] | P(the 3 nt of the AUG are all unpaired) |
| `rnaplfold_window`, `rnaplfold_span`, `rnaplfold_u_max` | int | `-W` (300), `-L` (150), `-u` (16) |

**Primary start window** (30 nt upstream of the AUG to 30 nt downstream), for each `u<l>`

| Column | Meaning |
|---|---|
| `start_mean_unpaired_u<l>` | mean stretch value. Dominated by the most open positions |
| `start_meanlog_unpaired_u<l>` | mean of **log10** stretch value (values floored at 1e-12). Mean opening energy up to a constant; negative, more negative = more structured. Probably the better summary |
| `start_min_unpaired_u<l>` | bottleneck: the least accessible stretch in the window |
| `start_max_unpaired_u<l>` | best landing pad: the most accessible stretch. Near 1 for `u1`, so informative mainly for larger `l` |
| `start_median_unpaired_u<l>` | median stretch value |
| `start_openfrac_unpaired_u<l>` | fraction of stretches with value > 0.5 |
| `start_up_meanlog_unpaired_u<l>` | mean log10 over the part of the window **before** the AUG (AUG excluded); `NA` if no 5'UTR in window |
| `start_down_meanlog_unpaired_u<l>` | mean log10 over the part **after** the AUG |

**Flank windows** (symmetric: `f` nt each side of the AUG; default `f` = 15 and 60; the primary 30/30 is above; **not** RNAplfold's `-W`)

`flank<f>_mean_unpaired_u<l>`, `flank<f>_meanlog_unpaired_u<l>`

**Region means** of log10 stretch value, from the same whole-transcript fold

`utr5_meanlog_unpaired_u<l>`, `cds_meanlog_unpaired_u<l>`, `utr3_meanlog_unpaired_u<l>`, and
`cap50_meanlog_unpaired_u<l>` (transcript positions 1..50, i.e. cap-proximal). `NA` if the region is empty
(no 5'UTR) or too short for `l`. `cds_` and `utr3_` are `NA` when the manifest has no 3'UTR length.

Not in the file: the primary window's `start_upstream`/`start_downstream` (recoverable from `start_win_*`
and `cds_offset`) and the raw per-position tables (kept gzipped in `<tool_dir>/lunp/` only when the run sets
`keep_lunp: true`; row `i`, column `l` = P(the `l` nt ending at position `i` are unpaired)).

## Value behaviour and gotchas

- All probabilities are in [0, 1] and **span orders of magnitude** (`u16` often 1e-6 to 1e-2; written in
  scientific notation, e.g. `2.54353e-05`). Parse as float. Expect strong right skew; log10 or logit is
  usually appropriate for modelling.
- Means are non-increasing in `l`: `u1 >= u4 >= u8 >= u16`. `aug_unpaired` is normally below `u1`.
  The many statistics here are strongly correlated (different summaries of one fold, overlapping windows
  and stretches), not independent features. Decide the families to compare up front and adjust for GC
  content and length before reading anything into a correlation.
- `short_utrs = true` rows have a window clipped at transcript position 1 (`start_win_start = 1`,
  `start_win_len < 63`): no upstream context in the window. Downstream QC already removes short-UTR
  transcripts, so exclude/flag them with `short_utrs` or `cds_offset < 30`.
- Accessibility is **local** (pair span capped at `rnaplfold_span`). Long-range pairing (e.g. Alu-Alu
  across a long UTR) is invisible to it. Results depend on `W` and `L`; do not mix rows from runs with
  different `rnaplfold_window`/`rnaplfold_span`. Store these as metadata with each load.
- Constant within a run: `region`, `rnaplfold_*`. A run's parameters are the only provenance in the file.

## Suggested ingestion

- Table keyed by `(gene_id, transcript_id)`; one row per transcript, `NA` -> missing.
- Keep `rnaplfold_window`, `rnaplfold_span`, `rnaplfold_u_max` as per-load metadata (or columns), plus the
  dataset name and load date.
- Typical modelling inputs: `aug_unpaired`, `start_mean_unpaired_u1/u4/u8/u16` (means and probabilities log/logit transformed; `*_meanlog_*` are already log10),
  gated on `short_utrs == "false"` and `seq_shorter_than_window == 0`.
