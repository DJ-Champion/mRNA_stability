# rnaplfold features: quick reference

74 result columns in `combined.tsv` with the default settings (plus the key and manifest columns), but only **8 ideas** repeated across 4 stretch
lengths and several windows. Column-by-column detail: `RNAPLFOLD_INTERFACE.md`.

## Three terms

- **Not to be confused:** `rnaplfold_window` is RNAplfold's `-W` (300 nt, used while folding). The
  **start window** (63 nt around the AUG, `start_*` columns) and the **flank windows** (`flank15_*`,
  `flank60_*`) are only where scores are read from afterwards.

- **Stretch value** = probability that `l` consecutive nucleotides are all unpaired. `l = 1` is the plain
  per-nucleotide unpaired probability. Every feature summarises stretch values over some window.
- **`u<l>`** in a column name = the stretch length `l` (default 1, 4, 8, 16). Longer stretch = stricter
  test of "open", so smaller numbers.
- **Direction**: higher probability = more accessible = less structured. `*_meanlog_*` columns are
  log10 probabilities, so they are negative; **closer to 0 = more accessible**, more negative = more
  structured.

All features come from one RNAplfold run on the full mRNA (`-W 300 -L 150`), then are read off a window.
Fold is local (pairs span at most 150 nt), so long-range pairing is invisible.

## Where the numbers are taken from (windows and regions)

| Name | Positions | Notes |
|---|---|---|
| Primary start window | 30 nt before the AUG to 30 nt after it (63 nt) | `start_*` columns; clipped at transcript ends |
| Upstream / downstream | the part of the primary window before / after the AUG | AUG itself excluded |
| Flank windows | `f` nt either side of the AUG, `f` = 15 and 60 | `flank15_*`, `flank60_*` |
| 5'UTR, CDS, 3'UTR | whole regions | `utr5_*`, `cds_*`, `utr3_*` |
| Cap-proximal | transcript positions 1 to 50 | `cap50_*`; may extend into the CDS if the 5'UTR is short |

## Feature families

| Family | Columns (u = 1, 4, 8, 16) | Count | What it asks | How to read it |
|---|---|---|---|---|
| **Mean** | `start_mean_unpaired_u<l>` | 4 | Average accessibility of the start window | Dominated by the few most open positions; skewed |
| **Mean log** | `start_meanlog_unpaired_u<l>` | 4 | Same, but on the energy-like scale | Mean opening energy up to a constant. Probably the better "overall structure" summary |
| **Min** | `start_min_unpaired_u<l>` | 4 | How bad is the worst spot? | Bottleneck: a single stable helix gives a low value |
| **Max** | `start_max_unpaired_u<l>` | 4 | How good is the best spot? | Best landing pad. Near 1 for `u1` (little signal); useful for `u8`, `u16` |
| **Median** | `start_median_unpaired_u<l>` | 4 | Typical position, ignoring extremes | Robust version of the mean |
| **Open fraction** | `start_openfrac_unpaired_u<l>` | 4 | What share of the window is open? | 0 to 1; share of stretches with value > 0.5 |
| **Up / down split** | `start_up_meanlog_unpaired_u<l>`, `start_down_meanlog_unpaired_u<l>` | 8 | Is the structure before or after the AUG? | Upstream = what a scanning ribosome meets first. `NA` if the window has no 5'UTR part |
| **Flank windows** | `flank15_` and `flank60_` `{mean,meanlog}_unpaired_u<l>` | 16 | Does the answer depend on how much flank is included? | Compare against the primary 30/30 window. `u16` is often `NA` in `flank15` (31-nt window) |
| **Regions** | `utr5_`, `cds_`, `utr3_`, `cap50_` + `meanlog_unpaired_u<l>` | 16 | Is overall structure in a region related to the outcome? | Same scale as Mean log. `NA` if the region is empty or shorter than `l` |
| **AUG** | `aug_unpaired` | 1 | Is the start codon itself open? | 3-nt stretch value; one position, so noisy |

That is 4 + 4 + 4 + 4 + 4 + 4 + 8 + 16 + 16 + 1 = 65 feature columns, plus the 9 below = 74.

## Columns that are not features (9)

| Column | Use |
|---|---|
| `start_win_start`, `start_win_end`, `start_win_len` | The window actually used (normally 63 long). Short windows mean clipping at a transcript end |
| `seq_length` | Folded length |
| `seq_shorter_than_window` | 1 = mRNA shorter than 300 nt; unreliable. Use as a **filter** |
| `cds_offset` | 5'UTR length. Use for filtering short 5'UTRs |
| `rnaplfold_window`, `rnaplfold_span`, `rnaplfold_u_max` | Parameters, constant within a run. Record as metadata |

## Which ones overlap (don't count them as independent)

- All `u` lengths within a family are correlated: `u1` > `u4` > `u8` > `u16`, same window.
- Mean, median, mean log and open fraction are all "typical level" summaries of the same values.
- The 15, 30 and 60 windows overlap heavily (each contains the one before it).
- `start_up_` and `start_down_` split the primary window, so they combine into `start_meanlog_`.
- `aug_unpaired` sits inside every start window.

Independent-ish *ideas*: overall level (mean log), worst spot (min), best spot (max), upstream vs
downstream, flank size, and whole-region structure (5'UTR, CDS, 3'UTR).

## Suggested starting set (a suggestion, not tested)

To avoid fishing across 65 correlated features, decide a short list first, for example:
`start_meanlog_unpaired_u4`, `start_min_unpaired_u8`, `start_max_unpaired_u8`,
`start_up_meanlog_unpaired_u4`, `utr5_meanlog_unpaired_u4`, `utr3_meanlog_unpaired_u4`,
`cap50_meanlog_unpaired_u4`, then use the rest to check robustness. Adjust for GC content and length
before reading anything into a correlation.
