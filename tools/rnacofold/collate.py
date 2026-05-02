#!/usr/bin/env python3
"""
tools/rnacofold/collate.py
Combine per-sequence results CSVs into a single TSV, joined with manifest
metadata (gene_id, transcript_id, region, length, etc.).

Required env (set by lib/paths.sh -> resolve_paths):
  RESULTS_DIR     directory containing results_<seqname>.csv files
  MANIFEST_TSV    extracted_regions/manifest.tsv
  TOOL_DIR        run dir for this tool — output lands here

Output:
  $TOOL_DIR/combined.tsv

Schema:
  seqname  gene_id  transcript_id  region  length  short_utrs
  selection_reason  orig_mfe  median_mfe  pvalue  zscore
  n_valid  n_missing  n_shuffles_requested

Records present in RESULTS_DIR but missing from MANIFEST_TSV are written
with empty manifest fields and a warning is logged. Records in MANIFEST_TSV
without results are silently omitted (use find_missing.sh to surface them).
"""
import csv
import glob
import logging
import os
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S'
)


REQUIRED_ENV = ("RESULTS_DIR", "MANIFEST_TSV", "TOOL_DIR")

OUT_KEYS = [
    "seqname", "gene_id", "transcript_id", "region", "length", "short_utrs",
    "selection_reason",
    "orig_mfe", "median_mfe", "pvalue", "zscore",
    "n_valid", "n_missing", "n_shuffles_requested",
]


def load_manifest(path):
    """Return {seqname -> manifest_row_dict}."""
    manifest = {}
    with open(path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            manifest[row["seqname"]] = row
    return manifest


def iter_results(results_dir):
    """Yield (seqname, results_row_dict) for every results_*.csv file."""
    pattern = os.path.join(results_dir, "results_*.csv")
    for path in sorted(glob.glob(pattern)):
        try:
            with open(path, newline='') as f:
                rows = list(csv.DictReader(f))
        except OSError as e:
            logging.warning(f"Cannot read {path}: {e}")
            continue
        if not rows:
            logging.warning(f"Empty results file: {path}")
            continue
        # Per-sequence files contain a single data row by convention.
        for row in rows:
            seqname = row.get("seq_name") or row.get("seqname")
            if not seqname:
                logging.warning(f"No seq_name in {path}")
                continue
            yield seqname, row


def main():
    missing = [v for v in REQUIRED_ENV if not os.environ.get(v)]
    if missing:
        logging.error(f"Missing required env vars: {', '.join(missing)}")
        sys.exit(1)

    results_dir = os.environ["RESULTS_DIR"]
    manifest_tsv = os.environ["MANIFEST_TSV"]
    tool_dir = os.environ["TOOL_DIR"]

    if not os.path.isfile(manifest_tsv):
        logging.error(f"Manifest not found: {manifest_tsv}")
        sys.exit(1)
    if not os.path.isdir(results_dir):
        logging.error(f"Results dir not found: {results_dir}")
        sys.exit(1)

    manifest = load_manifest(manifest_tsv)
    logging.info(f"Loaded {len(manifest)} manifest rows")

    out_path = os.path.join(tool_dir, "combined.tsv")
    n_written = 0
    n_orphan = 0

    with open(out_path, 'w', newline='') as out_f:
        writer = csv.DictWriter(out_f, fieldnames=OUT_KEYS, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()

        for seqname, r in iter_results(results_dir):
            m = manifest.get(seqname)
            if m is None:
                n_orphan += 1
                m = {}

            writer.writerow({
                "seqname": seqname,
                "gene_id": m.get("gene_id", ""),
                "transcript_id": m.get("transcript_id", ""),
                "region": m.get("region", ""),
                "length": m.get("length", ""),
                "short_utrs": m.get("short_utrs", ""),
                "selection_reason": m.get("selection_reason", ""),
                "orig_mfe": r.get("orig_mfe", ""),
                "median_mfe": r.get("median_mfe", ""),
                "pvalue": r.get("pvalue", ""),
                "zscore": r.get("zscore", ""),
                "n_valid": r.get("n_valid", ""),
                "n_missing": r.get("n_missing", ""),
                "n_shuffles_requested": r.get("n_shuffles_requested", ""),
            })
            n_written += 1

    if n_orphan:
        logging.warning(f"{n_orphan} result(s) have no matching manifest row "
                        f"(written with empty manifest fields)")
    logging.info(f"Wrote {n_written} records to {out_path}")


if __name__ == "__main__":
    main()