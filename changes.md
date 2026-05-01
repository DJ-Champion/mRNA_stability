# Log of changes

Changes made to pipeline. 

## 2026-05-01
New tool, three files:

configs/tools/rnacofold.sh — config with TOOL_REGIONS="UTR_pair" (the new region-restriction hook). Calibration hooks mirror RNAfold's: MFE-only fast path, t_full ≈ t_mfe × (N + 1), identity RSS predictor.
tools/rnacofold/worker.sh — reads the existing UTR_pair 3-line FASTA, parses 5UTR/3UTR boundaries from the constraint string (which is the only thing the constraint is used for here — discarded after parsing), reformats as 5UTR&3UTR, runs RNAcofold without constraints. Shuffles each strand independently with esl-shuffle -d, rejoins with &, refolds. Output schema is identical to the RNAfold worker's, so collation is uniform.
tools/rnacofold/collate.py — joins per-sequence results CSVs with manifest.tsv into combined.tsv. Schema includes manifest fields (gene_id, transcript_id, region, length, short_utrs) plus the result fields (orig_mfe, median_mfe, pvalue, zscore, etc.).

Infrastructure for region-restricted tools (additive — RNAfold unaffected):

lib/paths.sh — defaults TOOL_REGIONS="" after sourcing the tool config. Empty means "all regions" (existing behaviour).
bin/02_stratify.sh — .lengths.tsv gains a 4th column region, populated from the bash-side region loop. Existing column positions unchanged.
bin/03_calibrate.sh — sample-selection awk filters on column 4 against TOOL_REGIONS if set. Both the main per-tier picker and the --verify longest-sample picker. Tiers with zero applicable samples are skipped silently.
bin/04_submit.sh — new filter_list_by_region helper joins the tier list against the manifest by seqname, keeping only allowed regions. Applied in calibrated-mode only; --list mode is left untouched (escape hatch for "I know what I'm doing"). Tiers with zero matches log and skip.