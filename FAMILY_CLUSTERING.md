# Family clustering for blocked cross-validation

Specification for assigning every gene in a cohort to a **sequence family**, so
that downstream modelling can hold all members of a family inside a single CV
fold. Two independently implementable halves:

| half | language | input | output |
|---|---|---|---|
| **Part 1 — pipeline** | Python + `seqkit` + `mmseqs` | `extracted_CDS.fa` per dataset | `family.tsv` (one row per gene) |
| **Part 2 — downstream** | R | `family.tsv` | outer/inner fold assignments |

The seam is `family.tsv`. Part 2 never reads the alignment table.

**Part 1 is implemented** — `bin/01c_family.py` plus `lib/family.py`:

```bash
./bin/01c_family.py --cohort human_only
./bin/01c_family.py -c human_only --recluster   # re-threshold, reuse the search
```

Part 2 is the specification below, for the downstream modelling pipeline.

---

## 0. What this is and is not

**Not phylogenetics.** The families do not need to be biologically correct.
The single guarantee required is: *no sequence in a test fold has a close
homologue in its training fold.* Single-linkage (connected components) delivers
exactly that guarantee by construction, because leakage is a pairwise property
and connected components are the transitive closure of pairwise edges.

The error modes are asymmetric:

* **Over-merging** — families larger than biological truth. Costs statistical
  power. Safe.
* **Under-merging** — a true homologue pair split across folds. Inflates the
  performance estimate. Unsafe, and undetectable from the model output.
* **Chaining to a giant component** — the degenerate case. Not "extra safe":
  if one block exceeds `1/k` of the corpus you cannot build balanced folds at
  all, and most of the data lands in one fold.

Design follows from that: bias toward over-merging, but instrument the
component-size distribution so chaining is visible before it reaches the model.

**Multi-species is in scope.** A cohort may contain 1..n species. Orthologues
merging into one family across species is intended — training on mouse
`Rpl13a` and testing on human `RPL13A` is leakage. All species are pooled into
a single search.

---

## Part 1 — Pipeline side

### 1.1 Cohort configuration

New config axis: `configs/cohorts/<name>.yaml`. A cohort names the datasets
pooled into one clustering run. A cohort of one dataset is the single-species
case; nothing special-cases it.

```yaml
# configs/cohorts/human_only.yaml
members:
  - dataset: human_test
    species: human
    transl_table: 1        # NCBI genetic code; 1 = standard

search:
  sensitivity: 7.5         # mmseqs -s; 7.5 = maximum
  max_seqs: 2000           # prefilter depth; see 1.4
  evalue: 1.0e-3           # search ceiling; MUST be >= loosest level below
  min_cov: 0.4             # mmseqs -c; MUST be <= loosest level below
  cov_mode: 0              # bidirectional coverage
  threads: 12

# Clustering levels. All are emitted as columns; downstream picks one.
# Every level must be >= the search settings above -- filtering can only
# tighten, never loosen, without re-running the search.
levels:
  strict: { min_nbs: 0.60, min_cov: 0.80, max_evalue: 1.0e-10, min_fident: 0.50 }
  medium: { min_nbs: 0.40, min_cov: 0.65, max_evalue: 1.0e-5,  min_fident: 0.30 }
  loose:  { min_nbs: 0.25, min_cov: 0.50, max_evalue: 1.0e-3,  min_fident: 0.20 }

selection:
  k_folds: 5               # ceiling is derived from k; see 1.7
  max_fold_fraction: 0.25  # largest family as a share of one fold
  # max_family_pct: 0.05   # explicit override, bypasses the derivation

translate:
  min_protein_len: 30
  max_internal_stop_frac: 0.05   # abort above this

binaries: {}                     # falls back to $SEQKIT_BIN / $MMSEQS_BIN, then PATH
  # seqkit: /opt/seqkit/2.13.0/seqkit
  # mmseqs: /opt/mmseqs/bin/mmseqs
```

Levels looser than `search:` abort at config load rather than silently
evaluating against edges that were never searched for. Levels are ordered by
their thresholds, not by YAML order, so "loosest" in the reports is well
defined regardless of how the block is written.

Outputs are namespaced by search parameters, because a different *search* needs
a new alignment table while different *levels* are only extra columns:

```
runs/_cohorts/<cohort>/family/<search_hash>/
    proteins.faa          pooled, namespaced, translated
    proteins.tsv          per-protein metadata + the gene universe
    proteins.ids          gene universe, one ID per line (for the R appendix)
    hits.tsv              alignment table (kept artefact)
    hits_rescued.tsv      --mask 0 re-search of low-complexity sequences,
                          written only when some sequence lacks a self-hit
    family.tsv            PRIMARY OUTPUT -- the seam
    locus_overlap_pairs.tsv  gene pairs sharing exonic sequence (see 1.6b)
    family_qc.tsv         per-level size distribution
    family_members.tsv    per-family composition (multi-species diagnostics)
    params.yaml           resolved config + tool versions
    translate_qc.tsv      per-dataset translation warnings
```

`<search_hash>` is a short hash of the resolved `search:` block plus the sorted
member list. Changing a level does **not** change the hash; neither does
`threads`, which must not affect the identity of a result.

---

### 1.2 Translation — do not use `seqkit translate --trim`

Verified against seqkit v2.13.0. The `--trim` help text claims it removes
trailing `X` and `*`. It actually **truncates the sequence at the first `X` or
`*`**:

| input CDS | raw translation | `--trim` output | correct |
|---|---|---|---|
| `ATGGCTTGGAAATAA` | `MAWK*` | `MAWK` | `MAWK` |
| `ATGGCTTGAAAAGCTTAA` | `MA*KA*` | `MA` | `MA*KA` |
| `ATGGCTNNNAAATAA` | `MAXK*` | `MA` | `MAXK` |

Consequence if used: selenoproteins (UGA read through as Sec, ~25 human genes)
and any CDS containing an assembly `N` are silently reduced to N-terminal
stubs. A stub fails the coverage filters, becomes a spurious singleton, and its
true homologues are then free to land in a different fold. That is
under-merging — the unsafe direction — and it is invisible downstream. The
problem scales with assembly quality, so it worsens exactly when the cohort
gains non-model species.

**Correct procedure.** Translate raw, post-process in Python:

```bash
seqkit translate -T <transl_table> -f 1 <in.fa> -o <raw.faa>
```

Then, per record:

1. Strip **one** trailing `*` if present (the terminal stop). Do not loop.
2. Replace every remaining internal `*` with `X`. Count these — mmseqs treats
   `X` as unmatchable, so internal stops weaken alignments rather than
   corrupting them.
3. Mark records shorter than 30 aa after step 1 as **unsearched**. Do not drop
   them. They are excluded from `proteins.faa` (too short to align
   meaningfully) but stay in the gene universe as flagged singletons, so a
   `left_join` from `manifest.tsv` remains lossless. Dropping them puts an
   unexplained NA in the blocking table — the same class of silent hole as the
   masking bug in 1.5.
4. Emit per-dataset counters to `translate_qc.tsv`: `dataset, species, n_cds,
   n_translated, n_searched, n_internal_stop, n_no_terminal_stop,
   n_too_short, n_with_X`.

`n_internal_stop` above a percent or so means the wrong `transl_table`, a frame
problem upstream, or a genuinely poor annotation. Fail loudly above 5%.
(Human MANE measures 19/13,601 = 0.14%.)

seqkit is doing real work here — it handles ambiguous codons properly
(`ACN`→`T`, `MGR`→`R`, etc.) rather than emitting `X`. Keep it; just don't let
it trim.

---

### 1.3 Pooled FASTA and ID namespacing

Concatenate all members into one `proteins.faa`. Header format:

```
><species>|<gene_id>|<transcript_id>
```

Pipe-delimited, exactly three fields, no spaces. `gene_id` is normalised per
the conventions in `METRICS.md` (namespace prefix and version suffix stripped);
`transcript_id` is verbatim from `manifest.tsv`.

The species field is mandatory even for single-species cohorts. It prevents ID
collisions when symbol-based gene lists are used, and it is what makes the
multi-species QC in 1.7 possible.

Assert IDs are unique across the pooled file. A duplicate is a config error
(same dataset listed twice, or two datasets sharing a species label) and must
abort, not deduplicate silently.

---

### 1.4 Search

One all-vs-all search, run permissively. All threshold decisions happen
downstream in 1.5 so that they live in exactly one auditable place.

```bash
mmseqs easy-search \
    proteins.faa proteins.faa hits.tsv "$TMPDIR" \
    -s 7.5 \
    -e 1e-3 \
    -c 0.4 --cov-mode 0 \
    --max-seqs 2000 \
    --threads "$THREADS" \
    -v 1 \
    --format-output "query,target,fident,alnlen,evalue,bits,qcov,tcov,qlen,tlen"
```

Notes:

* `-v 1` because mmseqs at its default verbosity emits several hundred lines
  per run and buries the pipeline's own logging. Errors and warnings still
  surface; `01c_family.py --verbose` restores `-v 3`.

* **Do not pass `--min-seq-id` here.** Identity filtering is applied
  downstream from the `fident` column. One filtering site, one place to audit.
* `--max-seqs 2000` raises the prefilter depth from the default 300. Large
  families (olfactory receptors ~400 members, C2H2 zinc fingers ~700) exceed
  300 and would have their hit lists truncated. Single linkage needs only
  enough edges to connect a component, so truncation is usually survivable,
  but *which* edges survive becomes order-dependent and therefore
  irreproducible. Raise it.
* The search settings are the outer bound of everything downstream. The
  `loose` level must sit at or inside `-e 1e-3` and `-c 0.4`. Validate this at
  config load and abort on violation — otherwise a level silently evaluates
  against edges that were never searched for.
* Record `mmseqs version` and `seqkit version` into `params.yaml`.

---

### 1.5 Edge construction and normalised bitscore

**Harvest self-bits before dropping self-hits.** This is order-critical.

```
self_bits[s] = bits of the row where query == target == s
```

Then drop self-hit rows.

**Low-complexity sequences need a rescue pass.** mmseqs masks low-complexity
regions with tantan in the prefilter (`--mask 1`, the default). A sequence
that is low-complexity along essentially its whole length retains no k-mers to
seed on as a *database* entry, so nothing finds it — including itself.

Two things are lost, not one:

1. **The denominator.** No self-hit means no `self_bits`, so `nbs` is
   undefined for every edge touching that sequence.
2. **The edges.** If the sequence's homologues are also masked, the edges
   between them disappear from the search entirely.

Both were observed. On human MANE, a ~33%-Cys keratin-associated protein had
no self-hit while appearing as a *query* with unambiguous paralogues
(e-values 1e-71 and 1e-56, full coverage both ways, 54–59% identity). On a
synthetic Cys/Gly/Ser-rich pair the failure was total: **zero rows** in the
main search despite the two being 99.5% identical.

Do **not** fall back to `self_bits[target]`. Where the query is the longer
sequence, `max()` would have selected `self_bits[query]`, so the smaller
denominator inflates `nbs` and admits edges that should have been rejected.

Do **not** rerun the whole search with `--mask 0`. Unmasked all-vs-all is
exactly what chains Cys-rich, Gly-rich and Ser-rich regions into one giant
component.

Instead, re-search only the affected sequences as queries against the *full*
database with `--mask 0`, every other parameter matching the main run. This
recovers their self-bitscores and their edges together, on the same scoring
scale — masking decides whether an alignment is *found*, not how it is
*scored* — and the recovered rows face the same level filters as everything
else. Restricting the unmasked queries to the handful that need it keeps the
chaining blast radius small.

Append the rescued rows to the edge stream and re-harvest self-bits from them.
If a protein still has no self-hit afterwards, abort with the offending IDs
listed rather than letting it through as a silent singleton.

**Normalised bitscore:**

```
nbs(q,t) = bits(q,t) / max(self_bits[q], self_bits[t])
```

Clamp to `[0, 1]`.

`max` in the denominator, not `min` and not the geometric mean, and not the
classic BSR convention of dividing by the query self-score. This is the primary
anti-chaining lever. Worked example — a 100 aa protein that fully matches one
domain inside a 900 aa multi-domain protein:

| denominator | value | nbs | edge kept? |
|---|---|---|---|
| `self_bits[q]` (~200) | 200 | ~1.0 | yes — chains |
| `min(...)` (~200) | 200 | ~1.0 | yes — chains |
| **`max(...)`** (~1800) | 1800 | **~0.11** | **no** |

Only `max` requires similarity across the *longer* sequence's full length,
which is what stops shared domains from welding unrelated families together.
It is also more robust across a phylogeny than raw identity, which matters once
the cohort spans more than one species: orthologue identity runs ~85%
human–mouse but ~65% human–zebrafish, so a fixed `min_fident` blocks leakage
tightly within mammals and loosely across vertebrates. Keep `min_fident` low
and let `nbs` and coverage carry the decision.

**Per level**, keep an edge `(q,t)` iff all hold:

```
nbs(q,t)   >= level.min_nbs
qcov       >= level.min_cov
tcov       >= level.min_cov       # both, not either -- this is cov_mode 0
evalue     <= level.max_evalue
fident     >= level.min_fident
```

Edges are undirected. mmseqs reports most pairs twice (`q→t` and `t→q`), so
deduplicate to an unordered pair before counting edges for QC or `n_edges`
will be roughly doubled. Do not shortcut this by counting only rows where
`query < target`: a masked target can be reported in **one direction only**
(the human KRTAP above appeared solely as a query), and that shortcut drops
such edges from the count entirely. Deduplicate on the unordered pair itself.
Direction never matters for the components — only for the QC total.

---

### 1.6 Clustering

Connected components over the filtered edge set, per level.

**The vertex set is the full gene universe, not the vertices present in the
edge list.** Genes with no surviving edge are singleton families and must
appear in the output. Building the graph only from edges silently drops them —
this is the single most likely bug in the whole procedure, and it produces a
`family.tsv` that looks fine, joins fine, and quietly omits most of the corpus.

Assert at the end: `nrow(family.tsv) == n genes in pooled FASTA`.

**Canonical family IDs.** mmseqs is threaded and names clusters after
representative sequences, so raw output is not stable across runs. Assign IDs
deterministically from the membership itself:

1. Sort families by size, descending.
2. Break ties by the lexicographically smallest `<species>|<gene_id>` member.
3. Label `fam_<level>_00001`, `fam_<level>_00002`, ...

Re-running on identical input must produce identical labels. Assert this in a
test.

Suggested implementation: `scipy.sparse.csgraph.connected_components` on a COO
matrix, or a plain union-find. Both are linear-ish; neither is a bottleneck
next to the search.

### 1.6b Locus overlap — the second edge source

Protein similarity is not the only reason two genes' features can be
non-independent. Genes whose **exons overlap in the genome** are transcribed
from the same DNA, so their extracted regions share sequence — and the genes
need not be homologous at all, which is precisely why protein clustering
cannot see it.

Found empirically via the audit in 1.8: ~67 human gene pairs showing 100%
3'UTR identity with no protein relationship, mostly **antisense**, all
confirmed as overlapping annotations. One pair shared 3,197 nt.

Detecting them from annotation coordinates rather than sequence identity
catches more — **427 pairs** on human MANE, 92% antisense — because the audit
only sees pairs whose shared segment dominates ≥50% of *both* regions, while a
coordinate overlap of any size counts. Overlap lengths run median 531 / p95
2841 / max 5899 bp, with only 19% below 200 bp, so the 100 bp floor is not
scraping up noise. Adding these edges cut cross-family 3'UTR identity at ≥90%
by 10× (142 genes → 14) at an unchanged blocking level.

```yaml
locus_overlap:
  enabled: true
  min_overlap_bp: 100
```

Three decisions worth stating:

* **Exons, not gene spans.** A gene nested inside another's intron overlaps in
  span while sharing no transcribed sequence; its features are independent and
  it must not be linked. Only exonic overlap means shared sequence.
* **Both strands count.** An antisense overlap still means the same DNA. The
  two extracted regions are reverse complements, so length and GC are identical
  by construction and structure is correlated. Orientation does not rescue
  independence.
* **Not a similarity threshold.** Overlap is binary, so the same edges are
  applied at every level.

Edges feed the same union-find as the sequence edges, so families are the
connected components of the union — a pair is linked if it is protein-similar
**or** locus-overlapping. `locus_overlap_pairs.tsv` records every pair with
its overlap length, strands and which levels it actually merged;
`family_qc.tsv` reports `n_locus_overlap_merges` per level.

Note that `canonical.gff` carries no `gene` features — `write_filtered_gff` in
`bin/01_extract.py` keeps only lines linked to a retained transcript, and
gene rows have no such link. Gene coordinates therefore come from the
`gene_id` attribute on exon rows.

---

### 1.7 Outputs

**`family.tsv`** — the seam. One row per gene, wide:

| column | description |
|---|---|
| `species` | from the pooled header |
| `gene_id` | normalised; join key |
| `transcript_id` | verbatim from `manifest.tsv` |
| `dataset` | source dataset name |
| `family_id_strict` | canonical ID at the strict level |
| `family_size_strict` | member count of that family |
| `family_id_medium` / `family_size_medium` | as above |
| `family_id_loose` / `family_size_loose` | as above |
| `protein_len` | aa, post-processing |
| `had_internal_stop` | `true` / `false` |
| `searched` | `false` for proteins excluded from the search (below `min_protein_len`); they are forced singletons |

Every gene gets a family at every level, singletons included, so a `left_join`
from `manifest.tsv` is lossless. `NA` never appears in a `family_id` column.
That holds for unsearchable proteins too — they carry `searched = false` and a
real singleton family ID rather than being omitted.

**`family_qc.tsv`** — one row per level, and the basis for choosing one:

| column | description |
|---|---|
| `level` | strict / medium / loose |
| `n_edges` | deduplicated undirected edges |
| `n_locus_overlap_merges` | families joined by locus overlap that sequence similarity had not already joined (see 1.6b) |
| `n_genes` | corpus size |
| `n_families` | components |
| `n_singletons` | families of size 1 |
| `pct_genes_in_multigene_families` | |
| `median_family_size` | over non-singletons |
| `max_family_size` | |
| **`max_family_pct`** | `max_family_size / n_genes` — the decision variable |
| `n_species_in_largest` | |

**Selection rule: take the loosest level whose `max_family_pct` stays under a
ceiling derived from k.** Loosest, because over-merging costs power while
under-merging inflates the estimate. The ceiling exists because a family must
pack into a split: it may occupy at most `max_fold_fraction` (default 0.25) of
one split, and a split is `1/k` of the corpus. At k=5 that is 5%; at k=10, 2.5%.

For a holdout design, `k_folds` is set from the **smallest** split, not the
number of folds — an 80-10-10 split takes `k_folds: 10`.

Deriving it from k rather than fixing it matters — a fixed percentage silently
becomes wrong when k changes, and it is a packing heuristic rather than a hard
boundary. A level that misses narrowly while its largest family still sits
well inside one split is often the better choice, so the run warns rather
than silently falling through to a much stricter level.

**Measured on human MANE (13,601 genes, k=10, with locus-overlap edges):**

| level | families | singletons | max family | max % |
|---|---|---|---|---|
| loose | 9,126 | 6,893 (50.7%) | 368 | 2.71% |
| **medium** | **10,605** | 8,797 (64.7%) | 282 | **2.07%** |
| strict | 12,158 | 11,074 (81.4%) | 24 | 0.18% |

**`medium` is the chosen blocking level.** The size distribution is only half
the argument; the audit in 1.8 supplies the other half. Cross-family
nucleotide identity at ≥80% — the band where similarity is genuine leakage:

| region | loose | medium |
|---|---|---|
| CDS | 2 genes | 5 genes |
| 3'UTR | 0.397% | 0.412% |
| 5'UTR | 0.169% | 0.243% |

Three genes separate them in CDS. Counting *any* detectable match instead
gives 2 versus 67, but that gap is distant paralogues whose nucleotide
sequences have diverged past recognition while their proteins remain
alignable — detectable as protein, not as DNA, and therefore not leakage.
Reading the "any hit" column as if it measured leakage overstates the case for
`loose`; the identity-resolved numbers are the ones to use.

`medium` also yields more clusters (10,605 vs 9,126), which marginally helps
cluster-robust standard errors, and it splits the Ras superfamily along known
subfamily lines — RAB (27), RHO (17), RAS (14) — where `loose` merges the
whole superfamily (94). Both are defensible; since `family.tsv` carries every
level as a column, refitting at `loose` is a one-line sensitivity check.

`strict` is unusable despite the reassuring 0.18%: its top three families are
*all* ZNF fragments (23, 15, 10), meaning it shreds the largest real gene
family in the genome. 81% singletons is the tell. A small `max_family_pct` is
not evidence of a good level — always confirm against `family_members.tsv`
and the audit.

An earlier fixed 2% ceiling picked `strict` here, because `medium` at 2.04%
missed by six genes and the rule fell straight through. That knife-edge is
what motivated both the k derivation and the near-miss warning.

Expressing the rule as a percentage rather than a count is what makes it
survive the move to multi-species: adding species grows the corpus and each
family roughly in proportion, so the threshold does not need retuning. That
holds for comparably-annotated genomes; it breaks if the cohort mixes very
different annotation depths (e.g. human plus yeast), in which case inspect
`family_members.tsv` directly.

**`family_members.tsv`** — one row per family per level:
`level, family_id, n_members, n_species, members_per_species, member_gene_ids`.

The multi-species diagnostic is the members-to-species ratio. `5 species ×
1 gene` is a clean orthologue group. `5 species × 40 genes` is either a real
large family or chaining, and the member list tells you which within a minute.

---

### 1.8 Auditing what the protein basis misses

Protein clustering captures CDS-level relatedness well and says nothing
directly about UTRs — yet UTR-derived features (uORFs, 5'/3' structure, 3'UTR
composition) carry much of the signal in a stability model. `bin/01d_family_audit.py`
measures that residual:

```bash
./bin/01d_family_audit.py --cohort human_only --level medium
./bin/01d_family_audit.py -c human_only -l medium --regions 3UTR,5UTR,CDS
```

For every gene and region it finds the most similar gene in a **different**
family. That is the leakage candidate — genes in different families may land
in different splits, so high nucleotide identity between them is similarity
the blocking does not prevent. Within-family similarity is reported alongside
as a control, since it is blocked by construction.

The audit is **split-agnostic**: it audits the family partition itself, not a
particular train/test assignment, so one run covers k-fold and holdout alike.

Outputs under `<family_dir>/audit/<level>/`:

| file | contents |
|---|---|
| `audit_summary.tsv` | per region: % of genes with any cross-family hit, identity percentiles, counts at ≥0.70/0.80/0.90/0.95 |
| `audit_per_gene.tsv` | per gene per region: max cross-family identity, the partner gene, coverages, max within-family identity |
| `audit_pairs.tsv` | every cross-family pair above `--report-identity`, sorted by identity, with a `strand` column |

**Both strands are searched**, pinned with `--strand 2`. An antisense match
means the two genes are transcribed from the same DNA, which is exactly the
locus-overlap leakage of 1.6b, so it must be found. Note that mmseqs already
behaves this way by default even though its `--strand` help reports a default
of `1` (forward only) — verified empirically; `--strand 1` really is forward
only, but `easy-search` defaults to `2`. Pinning keeps the intent visible and
survives a corrected default. mmseqs signals a reverse-complement hit by
reporting query coordinates in descending order, which is how the `strand`
column is derived.

**Reading it.** If cross-family identity has no meaningful tail above ~0.70–0.80
at real coverage, protein-based blocking is sufficient and the CDS-only basis
is defensible — with a figure rather than an assumption. If there is a tail,
`audit_pairs.tsv` names the responsible pairs.

**Use the identity-banded columns, not `n_with_cross_family_hit`.** The "any
hit" count includes weak matches that are not leakage, and comparing levels on
it is misleading. On human MANE, CDS "any hit" gives 2 genes at `loose` versus
67 at `medium` — a 33× gap that looks decisive — while at ≥80% identity the
same comparison is 2 versus 5. The difference is distant paralogues whose
nucleotide sequences have diverged past recognition but whose proteins remain
alignable: detectable, but not similar enough to leak. `pct_ge_80` and
`pct_ge_90` are the decision variables.

Shared-repeat matches (an Alu in two unrelated 3'UTRs) are largely excluded by
`--min-cov` (default 0.5, both query and target): 300 bp of Alu inside a 2 kb
UTR cannot reach 50% coverage. This is the same coverage discipline that stops
shared protein domains welding families together in 1.5. Short UTRs where the
repeat *is* most of the sequence will still surface, and there it is genuinely
arguable whether that is leakage or signal.

**If the tail is real**, the fix is not to re-label anything. The union-find in
1.6 does not care where edges come from, so nucleotide edges can feed the same
graph and families become components of the union — a pair is linked if it is
protein-similar *or* UTR-similar.

---

## Part 2 — Downstream (R)

Consumes `family.tsv` only.

### 2.0 The unit of observation is the gene

Worth stating explicitly, because it dissolves a question that otherwise looks
hard. The response (half-life) is measured per transcript, and the join recipe
in `METRICS.md` produces **one row per transcript** with region-prefixed
feature columns.

So a region is not an observation — it is a *column* on the gene's row. There
is no separate family label for the 5'UTR, and no separate decision about
which split the 3'UTR lands in: a gene's regions travel with the gene. One row,
one family label, one split assignment.

This changes only if a region-level model is ever built (predicting decay from
3'UTR features alone, with UTRs as observations). Then the unit changes and the
blocking has to follow it.

### 2.1 Ingest and choose a blocking level

```r
library(data.table)

fam <- fread("runs/_cohorts/<cohort>/family/<hash>/family.tsv",
             na.strings = "NA")

BLOCK_LEVEL <- "medium"   # measured choice for human MANE; see 1.7.
                          # Refit with "loose" as a sensitivity check —
                          # every level is already a column.
fam[, family_id := get(paste0("family_id_", BLOCK_LEVEL))]

stopifnot(
  !anyNA(fam$family_id),                     # every gene has a family
  !anyDuplicated(fam[, .(species, gene_id)]) # one row per gene
)
```

Join to the feature table on `gene_id` (and `species` for multi-species
cohorts), matching the conventions in `METRICS.md`.

### 2.2 Outer folds — Blocked Group K-Fold

> **This project uses an 80-10-10 holdout, not k-fold.** Read 2.2a first; the
> packer below is the general case that 2.2a specialises. The *blocking*
> requirement is identical either way — a family must never span two splits —
> so everything about `family.tsv` is unchanged.

The requirement is that every family lands entirely in one fold. The
complication is that family sizes are wildly unequal — one family of 400
alongside 15,000 singletons — so assigning families to folds at random gives
badly unbalanced folds.

**`caret::groupKFold` is not adequate here.** It partitions groups without
balancing on group size, so a fold that happens to receive the large families
can end up several times the size of another.

Use greedy longest-processing-time (LPT) bin-packing: sort families by size
descending, assign each to the currently-smallest fold. Same algorithm
`sklearn.model_selection.GroupKFold` uses, and it gets within a few percent of
optimal balance.

```r
#' Assign families to k folds, balancing on gene count.
#' Returns data.table(family_id, fold).
assign_folds <- function(fam, k, seed = 1L) {
  sizes <- fam[, .(n = .N), by = family_id]

  set.seed(seed)
  sizes <- sizes[sample(.N)]        # randomise, so equal-sized families
  sizes <- sizes[order(-n)]         # tie-break differently per seed

  load <- numeric(k)
  out  <- integer(nrow(sizes))
  for (i in seq_len(nrow(sizes))) {
    j       <- which.min(load)
    out[i]  <- j
    load[j] <- load[j] + sizes$n[i]
  }
  sizes[, fold := out][, .(family_id, fold)]
}

folds <- assign_folds(fam, k = 5L, seed = 42L)
fam   <- merge(fam, folds, by = "family_id", all.x = TRUE)
```

`sample(.N)` before `order(-n)` matters: without it, ties break by whatever
order the table happened to arrive in, and every seed produces the same
partition.

### 2.2a The 80-10-10 holdout used by this project

Two changes from 2.2, both small.

**Unequal bins.** The packer above assigns each family to the *emptiest* bin,
which is only correct when bins are the same size. For unequal targets, assign
to whichever bin is furthest below its quota: `which.max(target - load)`
instead of `which.min(load)`. That is the whole generalisation.

**Large families are pinned to training by design.** With 13,601 genes, a 10%
test split is ~1,360 genes, and blocking keeps families intact — so the
largest family (282 members at `medium`) would be ~21% of the test set
wherever it landed. A single holdout gives one draw and no error bar, so
rather than hope, the biggest families go to training deliberately.

State the criterion rather than eyeballing it. "Any family larger than 5% of
the smallest split" is reproducible and defensible; "the top few" is not. The
consequence must be reported: the test set is depleted of large families, so
its performance is not an unbiased estimate for a randomly chosen gene.

```r
#' Assign families to named splits with unequal targets, pinning large
#' families to the first split. Returns data.table(family_id, split).
assign_holdout <- function(fam,
                           props = c(train = 0.8, val = 0.1, test = 0.1),
                           pin_frac = 0.05, seed = 1L) {
  sizes <- fam[, .(n = .N), by = family_id]
  n_tot <- sum(sizes$n)

  # Families too large to sit in the smallest split without dominating it.
  pin_above <- pin_frac * min(props) * n_tot
  pinned <- sizes[n >  pin_above]
  rest   <- sizes[n <= pin_above]

  target <- props * n_tot
  load   <- setNames(numeric(length(props)), names(props))
  load[1] <- sum(pinned$n)                 # pinned genes already in split 1

  set.seed(seed)
  rest <- rest[sample(.N)][order(-n)]      # randomise, then LPT
  out  <- character(nrow(rest))
  for (i in seq_len(nrow(rest))) {
    j       <- which.max(target - load)    # furthest below quota, not emptiest
    out[i]  <- names(load)[j]
    load[j] <- load[j] + rest$n[i]
  }

  rbind(data.table(family_id = pinned$family_id, split = names(props)[1]),
        data.table(family_id = rest$family_id,   split = out))
}

splits <- assign_holdout(fam, seed = 42L)
fam    <- merge(fam, splits, by = "family_id", all.x = TRUE)
```

**Write the split to disk and read it back everywhere.** Generate it once,
save `gene_id → split`, and have every downstream script consume that file.
Regenerating it inside the modelling code means a rebuilt `family.tsv`, a
different R version, or even a changed row order can silently reshuffle which
genes are in test — results stop being reproducible with nothing looking
wrong. Record the `<search_hash>` of the family run alongside it so the split
is traceable to the clustering that produced it.

The assertions in 2.4 apply unchanged, substituting `split` for `fold`.

**2.3 does not apply.** A three-way split already separates model selection
(`val`) from the final estimate (`test`), which is what the nested inner loop
existed to provide.

### 2.3 Inner folds — nested, and drawn only from outer-train

> Not used by this project — the `val` split of 2.2a already serves this
> purpose. Retained for a k-fold design.

The inner split must respect the same family blocks, and must be built from
the outer training set alone. Re-running the same packer on the retained
families is sufficient:

```r
#' For outer fold `o`, return inner fold assignments over the training part.
inner_folds <- function(fam, outer_fold, k_inner = 5L, seed = 1L) {
  train <- fam[fold != outer_fold]
  assign_folds(train, k = k_inner, seed = seed + outer_fold)
}
```

Vary the seed per outer fold, otherwise every outer iteration reuses an
identical inner partition and the inner variance estimate is correlated across
outer folds.

### 2.4 Assertions to run before modelling

These are cheap and they catch the failure this whole exercise exists to
prevent:

```r
# 1. No family spans two outer folds -- the core guarantee.
stopifnot(fam[, uniqueN(fold), by = family_id][, max(V1)] == 1L)

# 2. Every gene assigned.
stopifnot(!anyNA(fam$fold))

# 3. Folds are reasonably balanced (within 20% of even).
sz <- fam[, .N, by = fold]$N
stopifnot(max(sz) / min(sz) < 1.2)

# 4. No gene_id appears in more than one fold.
stopifnot(fam[, uniqueN(fold), by = .(species, gene_id)][, max(V1)] == 1L)
```

Assertion 3 failing means the blocking level produced a family too large to
pack — go back to `family_qc.tsv` and pick a stricter level.

If the response variable is continuous (half-life), also check that its
distribution is comparable across folds — `fam[, summary(half_life), by =
fold]`. LPT balances counts, not the target. A large family with an atypical
half-life (histones, ribosomal proteins) sits entirely in one fold by
construction, which is correct for the leakage guarantee but can skew that
fold's target distribution. Worth seeing before it surprises you in the
results.

### 2.5 What the design tests

Random-over-families with all orthologues held together tests **generalisation
to novel gene families**. Report it as such. It is a different and stronger
claim than random-over-genes, and a different claim from leave-one-species-out
(which would test cross-species transfer). Stating which one is being measured
in the write-up avoids the usual ambiguity.

Note that pinning the largest families to training (2.2a) weakens this
slightly: the test set is depleted of large families, so it measures
generalisation to *small and singleton* families rather than to all of them.
That is a reasonable trade for a representative test set on a single draw,
but it should be stated rather than left implicit.

### 2.6 If the goal is inference rather than prediction

For a question of the form "does secondary structure affect half-life?" the
family label's main job is **not** defining the split. It is the **clustering
unit for uncertainty estimates**.

A 282-member family is not 282 independent observations — its members share
ancestry, regulation and sequence. Fit a model treating them as independent
and the confidence interval on the structure coefficient comes out too narrow,
so a null result cannot be trusted either. Family should therefore enter the
model as a random effect, or as the cluster for cluster-robust standard
errors. The held-out split then answers the separate, secondary question of
whether the relationship predicts out of sample.

Two things follow:

* **Fit on the training portion and report the coefficient with clustered
  uncertainty.** With ~10,900 genes and ~8,500 clusters in `train`, there is
  ample power; the 10% test split is for the predictive check, not the effect
  estimate.
* **Use the shuffle-based z-scores, not raw MFE.** Raw MFE is strongly
  confounded with length and GC — longer sequences fold to more negative MFE
  almost mechanically — so a naive regression will recover a length effect and
  call it structure. The pipeline already folds 1,000 dinucleotide-preserving
  shuffles per sequence for exactly this reason. It is the single most
  important control for the central claim.

Since every level is a column in `family.tsv`, refitting with
`family_id_loose` is a one-line sensitivity check, and "conclusions were
unchanged under a looser grouping" is worth more than the level choice itself.

---

## Appendix — clustering in R instead

Recommended approach is Part 1: cluster in Python, hand R a small `family.tsv`.
The alignment table is large (millions of rows multi-species), the family table
is one row per gene, and precomputing keeps a single artefact of record so two
analyses cannot silently diverge.

If clustering must happen in R anyway, read `hits.tsv` and:

```r
library(data.table); library(igraph)

COLS <- c("query","target","fident","alnlen","evalue","bits",
          "qcov","tcov","qlen","tlen")

# Read BOTH tables. hits_rescued.tsv holds the --mask 0 re-search of
# low-complexity sequences (see 1.5) and exists only when some sequence
# lacked a self-hit. Skipping it loses those sequences' self-bitscores AND
# their edges, silently turning real families into singletons.
hits <- fread("hits.tsv", col.names = COLS)
if (file.exists("hits_rescued.tsv")) {
  hits <- rbind(hits, fread("hits_rescued.tsv", col.names = COLS))
}

# 1. Harvest self-bits BEFORE dropping self-hits. Order-critical.
self_bits <- hits[query == target, .(id = query, sb = bits)]
setkey(self_bits, id)

hits <- hits[query != target]

# 2. Normalised bitscore, max() in the denominator (see 1.5).
hits[, sb_q := self_bits[.(query),  sb]]
hits[, sb_t := self_bits[.(target), sb]]
stopifnot(!anyNA(hits$sb_q), !anyNA(hits$sb_t))   # missing self-hit
hits[, nbs := pmin(1, bits / pmax(sb_q, sb_t))]

# 3. Filter at one level.
lvl  <- list(min_nbs = 0.40, min_cov = 0.65, max_evalue = 1e-5, min_fident = 0.30)
keep <- hits[nbs    >= lvl$min_nbs &
             qcov   >= lvl$min_cov & tcov >= lvl$min_cov &
             evalue <= lvl$max_evalue &
             fident >= lvl$min_fident]

# 4. Components -- SINGLETON-SAFE.
#    Vertices come from the full gene universe, NOT from the edge list.
#    Omitting `vertices=` silently drops every gene with no surviving edge,
#    which is most of the corpus.
all_ids <- unique(readLines("proteins.ids"))   # one ID per protein, all of them
g <- graph_from_data_frame(keep[, .(query, target)],
                           directed = FALSE,
                           vertices = data.frame(name = all_ids))

comp <- components(g, mode = "weak")

memb <- data.table(id = names(comp$membership),
                   raw_component = as.integer(comp$membership))
memb[, n := .N, by = raw_component]

# 5. Canonical IDs: size desc, tie-break on smallest member ID.
ord <- memb[, .(n = .N, min_id = min(id)), by = raw_component][order(-n, min_id)]
ord[, family_id := sprintf("fam_%05d", .I)]
memb <- merge(memb, ord[, .(raw_component, family_id)], by = "raw_component")

stopifnot(nrow(memb) == length(all_ids))   # nothing dropped
```

`proteins.ids` must be generated alongside `proteins.faa` in step 1.3 — the
full ID list is not recoverable from `hits.tsv`.
