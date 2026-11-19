# rnaup_avoidance output: interface for downstream ingestion

Source: `runs/<dataset>/rnaup_avoidance/combined.tsv` (written by `05_collate.sh -t rnaup_avoidance`).
Producer: mRNA_stability repo, `tools/rnaup_avoidance/` (method: `RNAUP_AVOIDANCE.md`). Tab-separated, one
header row, `NA` for missing. `avoidance_start_codon_pm30.tsv` in the same directory is the same table
(the collator splits by `window`, and only one window exists); ingest either, not both.

## What it measures

How strongly the **translation start window** (30 nt upstream of the AUG + AUG + 30 nt downstream = 63 nt)
is predicted to base-pair with 34 abundant human ncRNAs (20 tRNA, 1 rRNA, 1 SRP, 5 snRNA, 3 Y_RNA,
4 snoRNA). Each (window, ncRNA) pair is scored by ViennaRNA RNAup; the 34 energies are summed, overall and
per family. The same is done for 100 dinucleotide shuffles of the window (`esl-shuffle -d`), giving a
robust Z of the native sum against its shuffles. The hypothesis behind the name: start regions avoid
pairing with abundant ncRNAs, i.e. bind them **less** than composition predicts.

Only the 63-nt window is folded: no flanking mRNA context, and the ncRNAs are used as given (not shuffled).
One row per transcript.

## Keys

| Column | Notes |
|---|---|
| `seqname` | `{gene_id}_{transcript_id}_start_codon_region`; unique per row |
| `gene_id`, `transcript_id` | Ensembl IDs, transcript ID includes version (`ENST00000009589.8`). Join on these. |

## Coverage: which transcripts have a row

Only windows of exactly 63 nt are run, i.e. **5'UTR >= 30 nt**. A transcript with a shorter 5'UTR has no
row at all (not an `NA` row). For `human_all` that is 12,346 of 13,601 transcripts; the 1,255 absent ones
are exactly those the R cohort filter (`MIN_5UTR_LENGTH = 30`) removes anyway. A transcript that
failed or has not finished is also absent (no partial rows are ever written), so after the cohort
filter, missing values mean an incomplete run, not a property of the transcript.

## Columns

Energies are kcal/mol, **more negative = stronger predicted binding**. Two energy families, each with the
same statistics:

- `avoidance_G_int_*`: **the metric.** `G_int`, RNAup's duplex (hybridisation) energy of the best
  interaction, without the cost of opening either partner's structure.
- `rnaup_G_total_*`: `G_total = G_int + G_u_target + G_u_query` (duplex + opening penalties). Kept for
  reference only; can be positive for a single pair.

A pair with no predicted interaction contributes 0 to every sum (and 0 is a candidate for `_min`).

**Metadata (from the manifest join)**

| Column | Type | Meaning |
|---|---|---|
| `region` | str | always `start_codon_region` |
| `length` | int | manifest region length; always 63 |
| `selection_reason` | str | isoform-choice rule (e.g. `MANE_Select`) |
| `short_utrs` | `true`/`false` | extraction flag: 5' **or** 3' UTR below the extraction minimum. Here the 5'UTR is always >= 30, so `true` (44 rows in `human_all`) means a short **3'UTR**; the window itself is unaffected |

**Geometry and run size**

| Column | Type | Meaning |
|---|---|---|
| `window` | str | always `start_codon_pm30` |
| `window_length` | int | always 63 |
| `n_shuffles` | int | valid shuffles used for this row; normally 100, can be lower (the producer requires >= 90% of `N_SHUFFLES` to be valid, else the gene is dropped) |
| `n_ncrna` | int | always 34 |

**Native energies**, for `<E>` in `avoidance_G_int`, `rnaup_G_total`

| Column | Meaning |
|---|---|
| `<E>_sum` | sum over the 34 ncRNAs |
| `<E>_mean` | `<E>_sum / 34` (redundant with the sum; do not use both) |
| `<E>_min` | most negative single-ncRNA energy: the strongest individual partner |
| `<E>_sum_<F>` | sum over family `F` in `tRNA rRNA SRP snRNA Y_RNA snoRNA` (the six add up to `<E>_sum`) |

**Shuffle null**

| Column | Meaning |
|---|---|
| `<E>_sum_robust_z` | `(native - median(shuf)) / (1.4826 * MAD(shuf))` on the 34-ncRNA sum. `NA` when MAD = 0 |
| `<E>_sum_<F>_robust_z` | the same for each family sum |
| `<E>_sum_shuf_median` | median of the shuffled sums (total only; no per-family median/MAD) |
| `<E>_sum_shuf_mad` | raw MAD of the shuffled sums (**not** multiplied by 1.4826) |

**GC covariates**

| Column | Meaning |
|---|---|
| `gc_content_mRNA_window` | GC fraction of the 63-nt window |
| `gc_ncRNA_mean` | mean GC of the 34 ncRNAs; a constant |
| `gc_diff_mRNA_vs_ncRNA` | `gc_content_mRNA_window - gc_ncRNA_mean`; a shifted copy of the window GC |

Exact header, in file order (50 columns):

```
seqname gene_id transcript_id region length selection_reason short_utrs window window_length n_shuffles n_ncrna
rnaup_G_total_sum rnaup_G_total_mean rnaup_G_total_min avoidance_G_int_sum avoidance_G_int_mean avoidance_G_int_min
rnaup_G_total_sum_tRNA avoidance_G_int_sum_tRNA rnaup_G_total_sum_rRNA avoidance_G_int_sum_rRNA
rnaup_G_total_sum_SRP avoidance_G_int_sum_SRP rnaup_G_total_sum_snRNA avoidance_G_int_sum_snRNA
rnaup_G_total_sum_Y_RNA avoidance_G_int_sum_Y_RNA rnaup_G_total_sum_snoRNA avoidance_G_int_sum_snoRNA
rnaup_G_total_sum_robust_z avoidance_G_int_sum_robust_z
rnaup_G_total_sum_shuf_median rnaup_G_total_sum_shuf_mad avoidance_G_int_sum_shuf_median avoidance_G_int_sum_shuf_mad
rnaup_G_total_sum_tRNA_robust_z avoidance_G_int_sum_tRNA_robust_z rnaup_G_total_sum_rRNA_robust_z avoidance_G_int_sum_rRNA_robust_z
rnaup_G_total_sum_SRP_robust_z avoidance_G_int_sum_SRP_robust_z rnaup_G_total_sum_snRNA_robust_z avoidance_G_int_sum_snRNA_robust_z
rnaup_G_total_sum_Y_RNA_robust_z avoidance_G_int_sum_Y_RNA_robust_z rnaup_G_total_sum_snoRNA_robust_z avoidance_G_int_sum_snoRNA_robust_z
gc_content_mRNA_window gc_ncRNA_mean gc_diff_mRNA_vs_ncRNA
```

Ingest by name, not position. Mixed case (`G_int`, `tRNA`, `Y_RNA`, `mRNA`) is upstream's; the R loader's
`lowercase_names()` will fold it.

Not in this file: the per-shuffle sums, in `runs/<dataset>/rnaup_avoidance/raw_shuffles/<seqname>_shuffle_sums.csv.gz`
(one gzipped CSV per transcript, not collated; columns `seq_name,iteration,G_total_sum,G_int_sum`,
`iteration` 1..`n_shuffles`). Only needed to plot or re-derive the null; per-ncRNA energies are not kept.

## Value behaviour and gotchas

- **Sign of the Z is the whole point.** `avoidance_G_int_sum_robust_z > 0`: the native window binds the
  ncRNAs *less* than its shuffles do (avoidance). `< 0`: binds *more* than expected. The raw sums run the
  other way round in words: more negative = more binding. `avoidance_G_int_sum` is a binding energy,
  not an avoidance score, despite its name.
- Typical ranges (8 random `human_all` windows, real RNAup, native only; treat as a sanity check):
  `avoidance_G_int_sum` about -270 to -460, `rnaup_G_total_sum` about -130 to -225,
  `avoidance_G_int_min` about -17 to -30, `avoidance_G_int_sum_tRNA` about 60% of the total. Sums are
  essentially always negative; check the Z distribution on the first real load (centre, spread, `NA` rate).
- tRNAs are 20 of the 34 ncRNAs, so they dominate the total sum. The per-family sums and Zs are the way
  to see the others.
- `NA` Z: MAD of the shuffles is 0. Unlikely for the total, more likely for single-member families
  (`rRNA`, `SRP`) and small ones. Keep as missing; do not impute 0.
- The Z controls for the window's dinucleotide composition by construction; the raw sums do not, and
  are expected to track window GC. Adjust raw sums for window GC: the R repo already has
  `gc_content_start` for this same 63-nt window, so `gc_content_mRNA_window` duplicates it. `gc_diff_mRNA_vs_ncRNA` is the same variable
  shifted by a constant, and `gc_ncRNA_mean` is constant: do not use them as features.
- All these columns are different summaries of one set of 34 x 101 RNAup runs and are strongly
  correlated (`_sum`/`_mean` exactly; `G_int` with `G_total`; total with the tRNA family). Decide the
  families to compare up front.
- Values are written with 4 decimals (`-412.3700`), but energies are summed from RNAup's 2-decimal
  output, so the last two digits carry nothing. Z and GC are full-precision floats.
- The 63-nt window is folded on its own, without the surrounding mRNA; its own structure enters only
  through `G_total`'s opening term. This is not comparable to the RNAplfold start-window accessibility
  (folded in full-mRNA context), though the two are worth comparing.
- Provenance not in the file: RNAup arguments (default `--interaction_first --no_output_file -b -d2 --noLP`,
  max interaction length `-w` 25), RNAup version, and the ncRNA set (`data/ncrna/curated-interaction-database.fa`).
  The only run parameters in the file are `n_shuffles`, `n_ncrna`, `window`, `window_length`. Do not mix
  runs with different settings; record them as metadata with each load.
- Constant within a run: `region`, `length`, `window`, `window_length`, `n_ncrna`, `gc_ncRNA_mean`.
  `n_shuffles` is near-constant (see above).

## Suggested ingestion (mRNA_stability_R)

Follows the RNAplfold loader (`load_rnaplfold()` in `R/io/load_raw.R`) and the README's "new raw input
source" recipe.

- Copy the file to `data/raw/human/rnaup_avoidance_combined.tsv` (cf. `rnaplfold_combined.tsv`);
  read via `species_path(species, "rnaup_avoidance_combined.tsv")` with `read_if_exists()`, so a missing file
  returns `NULL` and the build carries on.
- New `load_rnaup_avoidance(species, all = FALSE)`: transcript-level, keyed on `transcript_id`, region token
  **last** and `start` (the `REGION_ALIASES` mapping of `start_codon_region`). Add it to the
  `transcript_level` list in `R/pipeline/build_dataset.R`.
- Checks (stop, as the RNAplfold loader does): required columns present; `transcript_id` unique and
  non-empty (a seqname missing from the manifest is joined with empty strings, not `NA`); `window`,
  `window_length` (63) and `n_ncrna` (34) constant. Message the `n_shuffles` range and the number of rows
  with `n_shuffles` below the mode. Return the run's `n_shuffles`/`n_ncrna`/`window` as an attribute (cf.
  `"rnaplfold_params"`) and attach it in `build_dataset()`.
- No row-level gating is needed: every row is a full 63-nt window. Do not `NA` the `short_utrs == "true"`
  rows (3'UTR flag; the window is fine; the cohort filter decides about them).
- Default (model) columns, kept small as with RNAplfold: `avoidance_G_int_sum_robust_z` (primary) and
  `avoidance_G_int_sum`. `all = TRUE` returns the rest (family sums and Zs, `_min`, the `G_total` set,
  shuffle median/MAD). Drop `_mean` and all three GC columns (`gc_content_start` already exists; at most assert it agrees with
  `gc_content_mRNA_window` as a join check).
- Naming suggestion, `<source>_<energy>_<stat>[_<family>]_start` in lower case:
  `rnaup_gint_sum_start`, `rnaup_gint_zscore_start` (cf. `rnafold_zscore_`), `rnaup_gint_min_start`,
  `rnaup_gint_sum_trna_start`, `rnaup_gint_zscore_trna_start`, `rnaup_gtotal_*`,
  `rnaup_gint_shuf_median_start`, `rnaup_gint_shuf_mad_start`. Check
  the result against `format_col_name()` in `R/utils/naming.R` and adjust to whatever it expects.
- Bump `CACHE_VERSION` (19 at the time of writing) with a changelog line in `R/config.R`, add rows to
  `R/feature_table.csv` for every new column pattern (state the sign: positive Z = avoidance), and run
  `Rscript scripts/check_feature_table.R`.
- Expect after the cohort filter: a value for essentially every transcript once the run is complete; fewer
  means unfinished or failed genes upstream (`runs/<dataset>/rnaup_avoidance/errors/`), not missing biology.
