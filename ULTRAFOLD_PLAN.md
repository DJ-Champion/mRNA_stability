# Ultrafolded transcripts: figure plan and questions for Paul

## 1. Where we are

"Ultrafolded" = RNAfold z-score below about −30 on the spliced mRNA (MFE compared with dinucleotide shuffles).

- 13,423 genes have a z-score; **45 are below −30** (83 below −25, 149 below −20).
- All 45 have very long 3′UTRs (about 3–11 kb), each with **4–23 Alu elements**, in both orientations relative to the transcript. Every one has at least 4 oppositely oriented Alu pairs by sequence (up to 112).
- Alu orientation is taken from RepeatMasker (UCSC hg38 `rmsk`), converted to sense/antisense relative to the transcript, and checked against the transcript sequences.
- "Opposite pair" is a sequence-orientation count only. It does not yet say whether two Alus pair in the predicted structure.

## 2. Plan by plot

**Plot A: z-score distribution, and native vs shuffled MFE**
- A1: z-score distribution for all genes, with the ultrafolded ones marked. Data exists.
- A2: four panels, one per transcript: histogram of 500 shuffled MFEs with the native MFE marked (kcal/mol on the x-axis). Shuffles are being rerun and kept.
- Caveat: for the longest transcripts the original z-scores came from only 10–100 shuffles, so the tail of A1 is noisy.
- Optional: an orientation-flip control (reverse-complement the antisense Alus so all face the same way, fold once). If orientation drives the extreme score, the MFE should rise sharply. Cheap (one fold per variant), but extra work.

**Plot B: architecture of one representative transcript**
- 5′UTR / CDS / 3′UTR model with RepeatMasker Alus overlaid.
- Alus marked as "oppositely oriented" from the **MFE structure** of the full mRNA (RNAfold dot-bracket), not from sequence orientation. The recomputed MFE will be checked against the stored score.
- Optional layer: base-pair probabilities (`RNAfold -p`) to show which Alu–Alu pairings are robust; Alu–Alu arcs as a layer.

**Plot C: icSHAPE vs prediction along the 3′UTR**
- icSHAPE reactivity (cytoplasm and nucleoplasm) against predicted unpaired probability.
- Predicted track: RNAplfold (`-W 500 -L 300 -u 1`). RNAlfold is not planned, since it gives energies, not per-nucleotide probabilities.
- Both are local predictors and cannot see the long-range Alu–Alu pairing. A global unpaired probability from `RNAfold -p` would include it.
- For "mismatch" to mean anything we also need the same comparison for ordinary, well-covered transcripts as a baseline.

## 3. The main problem: icSHAPE coverage

- Of the 45 ultrafolded transcripts, **icSHAPE scores 0.5% of Alu nucleotides in nucleoplasm and 0.03% in cytoplasm**. Only 2 have at least 50 scored nucleotides in their Alus.
- This is mostly a general depth problem: only 12% (cytoplasm) and 18% (nucleoplasm) of all genes have any 3′UTR score, and these candidates are poorer than average even outside Alus (0.2–1.5% against 5–9%).
- Alus lose about half more: in well-covered genes, Alu nucleotides are scored at about 0.5× the non-Alu rate. The cause is untested; multi-mapping is one possibility.
- Checking the 14 transcripts with z < −20 and a Gini value: 8 have **no** scored Alu nucleotides (a Gini value only means data exist somewhere in the transcript). Three have real coverage:

| gene | z | Alu nt in 3′UTR | Alu nt scored (cytoplasm) | Alus (opposite pairs) |
|---|---|---|---|---|
| SYAP1 | −27.8 | 1,903 | 44% | 7 (12) |
| PSMB2 | −24.8 | 1,226 | 44% | 4 (4) |
| DESI1 | −23.2 | 604 | 39% | 2 (1) |

- None is below −30. The usable examples are less extreme than the top of the list, and they are high-expression genes.

## 4. Questions for Paul

**Choice of transcript and threshold**
1. Is SYAP1 or PSMB2 (z about −25 to −28) acceptable as "the representative ultrafolded transcript", given that no transcript with z < −30 has usable Alu coverage? Should B and C use the same transcript?
2. If not, would Paul rather use the most extreme transcript for B and a different, covered transcript for C? Or show the coverage gap itself as a result?
3. Is there deeper or additional probing data (other icSHAPE or DMS-based datasets, or the raw reads behind the processed scores)?

**Plot A**
4. Include one ordinary, length-matched control transcript among the four panels?
5. Include the orientation-flip control?

**Plot B**
6. What cutoff defines an Alu as "paired with an opposite Alu" in the MFE structure: at least N base pairs, or at least X% of the Alu's length?

**Plot C**
7. Primary predicted track: global `RNAfold -p` (includes long-range pairing), RNAplfold, or both?
8. Is a baseline of icSHAPE-vs-prediction correlation over about 100 well-covered ordinary transcripts acceptable for judging "mismatch"?
9. Whole transcript or 3′UTR only? A transcriptome-wide Alu vs non-Alu reactivity summary alongside the single-transcript example?

**Other**
10. Is the GO enrichment still wanted? Gene lists and the 13,423-gene background are ready. Ultrafolded genes all have long 3′UTRs, so it would need a length-matched control.

## 5. Next steps once decided

1. MFE structures and base-pair probabilities for the candidate transcripts (cluster).
2. Shuffles (500 each) with MFEs kept; optionally the flipped-Alu folds.
3. RNAplfold, `RNAfold -p` unpaired track, and the baseline correlations.
4. Plotting in `mRNA_stability_R`.

Data produced so far: `runs/human_all/ultrafold/` (45 candidates) and `runs/human_all/ultrafold/top14/` (z < −20 with icSHAPE). Methods and column definitions: `ULTRAFOLD_DATA.md`.
