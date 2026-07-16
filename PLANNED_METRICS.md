# Planned metrics

Per-transcript metrics under consideration for the pipeline. Not yet
implemented; design notes and priorities below. For the contract every
metric must follow, see `METRICS.md` and `metrics/__init__.py`.

## Candidates, in rough priority order

| plugin | shape | what it would measure |
|---|---|---|
| `gc_by_position` | wide | GC content at codon positions 1, 2, 3 across the CDS (GC1, GC2, GC3). GC3 is strongly associated with codon optimality and stability, mostly independent of CAI. |
| `kozak` | wide | Kozak-context strength at the main ATG; scalar 0–1 score from the canonical PWM (`gccRccATGG`, −3 and +4 strongest). Translation initiation efficiency proxy. |
| `stop_context` | wide | stop codon identity (TAA / TAG / TGA), the +4 nucleotide, and a binary `is_strong_stop` flag. Readthrough propensity correlate. |
| `utr3_motifs` | wide | densities of stability-modulating 3'UTR motifs: AU-rich elements (AUUUA, WWAUUUAWW), Pumilio sites (UGUANAUA), m6A consensus (DRACH), and polyadenylation signals (AAUAAA / AUUAAA) in the last ~50 nt. |

All four are sequence-only (no external reference data). `gc_by_position`,
`kozak`, and `stop_context` reuse the FASTA-parse + manifest-index pattern
established by `sequence_basic`, `codon_aa_counts`, and `uorf`.
`utr3_motifs` is the same pattern with hardcoded motif regexes.

## Explicitly excluded (with rationale)

- **CSC** — Fitted directly from mRNA half-life data, so using it as a
  feature in a half-life predictor would be methodologically circular.
- **CAI** — Implemented as `metrics/cai.py`. Reference weights are built
  from CoCoPUTs ribosomal protein genes (RPL*/RPS*) via
  `bin/build_cai_weights.py`.
- **miRNA seed sites** — Real biological signal but the reference-data
  burden (which miRNAs, which seed conservation cutoff) is high. Deferred
  until a specific hypothesis needs them.
- **Effective number of codons (Nc)** — Strongly correlated with CAI and
  GC3; adds little independent signal once GC3 is captured.
