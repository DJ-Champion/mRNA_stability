#!/usr/bin/env python3
"""
05_collate.py
Combine per-sequence RNAfold result CSVs into one table per metric source.
Joins with manifest.tsv so the output has gene/transcript/region context.

Outputs:
  results/combined_rnafold.tsv          # one row per sequence, all RNAfold cols
  results/combined_rnafold_raw_shuffles.tsv  (optional, large)

Designed for join-at-analysis-time in R:
  manifest <- read_tsv("extracted_regions/manifest.tsv")
  rnafold  <- read_tsv("results/combined_rnafold.tsv")
  gc       <- read_tsv("results/combined_gc.tsv")        # future metrics
  full     <- manifest |> left_join(rnafold) |> left_join(gc)
"""
import argparse
import csv
import logging
import os
import sys
from glob import glob

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S'
)

def collate_rnafold(results_dir: str, manifest_path: str, out_path: str):
    """Concatenate per-sequence result CSVs, attach manifest metadata."""
    manifest = {}
    if os.path.exists(manifest_path):
        with open(manifest_path) as f:
            r = csv.DictReader(f, delimiter='\t')
            for row in r:
                manifest[row['seqname']] = row
        logging.info(f"Loaded manifest: {len(manifest)} records")
    else:
        logging.warning(f"Manifest not found: {manifest_path} — output will lack metadata")

    pattern = os.path.join(results_dir, "results_*.csv")
    files = sorted(glob(pattern))
    if not files:
        logging.error(f"No result files matched: {pattern}")
        sys.exit(1)
    logging.info(f"Found {len(files)} result files")

    # Build the combined table. Headers come from the first file; we add
    # manifest columns on the left for join convenience in R.
    manifest_cols = ['gene_id', 'transcript_id', 'region', 'length',
                     'selection_reason', 'short_utrs']

    written = 0
    with open(out_path, 'w', newline='') as fout:
        writer = None
        for fpath in files:
            with open(fpath) as fin:
                rdr = csv.DictReader(fin)
                for row in rdr:
                    seqname = row['seq_name']
                    meta = manifest.get(seqname, {})

                    # Compose the output row: manifest cols first, then results
                    out_row = {'seqname': seqname}
                    for c in manifest_cols:
                        out_row[c] = meta.get(c, '')
                    for k, v in row.items():
                        if k == 'seq_name':
                            continue
                        out_row[k] = v

                    if writer is None:
                        fieldnames = ['seqname'] + manifest_cols + \
                                     [k for k in row.keys() if k != 'seq_name']
                        writer = csv.DictWriter(fout, fieldnames=fieldnames,
                                                delimiter='\t')
                        writer.writeheader()
                    writer.writerow(out_row)
                    written += 1

    logging.info(f"Wrote {written} rows to {out_path}")

def collate_raw_shuffles(raw_dir: str, manifest_path: str, out_path: str):
    """Concatenate per-sequence raw shuffle CSVs. Large file — opt in."""
    manifest = {}
    if os.path.exists(manifest_path):
        with open(manifest_path) as f:
            r = csv.DictReader(f, delimiter='\t')
            for row in r:
                manifest[row['seqname']] = row

    pattern = os.path.join(raw_dir, "*_raw_shuffles.csv")
    files = sorted(glob(pattern))
    if not files:
        logging.error(f"No raw shuffle files matched: {pattern}")
        return
    logging.info(f"Found {len(files)} raw shuffle files")

    manifest_cols = ['gene_id', 'transcript_id', 'region']
    written = 0

    with open(out_path, 'w', newline='') as fout:
        writer = None
        for fpath in files:
            with open(fpath) as fin:
                rdr = csv.DictReader(fin)
                for row in rdr:
                    seqname = row['transcript_id']   # legacy field name
                    meta = manifest.get(seqname, {})
                    out_row = {'seqname': seqname}
                    for c in manifest_cols:
                        out_row[c] = meta.get(c, '')
                    out_row['iteration'] = row['iteration']
                    out_row['shuffled_mfe'] = row['shuffled_mfe']

                    if writer is None:
                        fieldnames = ['seqname'] + manifest_cols + \
                                     ['iteration', 'shuffled_mfe']
                        writer = csv.DictWriter(fout, fieldnames=fieldnames,
                                                delimiter='\t')
                        writer.writeheader()
                    writer.writerow(out_row)
                    written += 1

    logging.info(f"Wrote {written} rows to {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Collate RNAfold per-sequence outputs")
    parser.add_argument('--results-dir',
                        default=os.environ.get('RESULTS_DIR', 'results/rnafold'))
    parser.add_argument('--raw-dir',
                        default=os.environ.get('RAW_SHUFFLES_DIR',
                                                'results/rnafold_raw_shuffles'))
    parser.add_argument('--manifest',
                        default=os.environ.get('MANIFEST_TSV',
                                                'extracted_regions/manifest.tsv'))
    parser.add_argument('--out-dir', default='results',
                        help='Where to write combined tables')
    parser.add_argument('--include-raw-shuffles', action='store_true',
                        help='Also collate raw shuffle CSVs (large output)')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    collate_rnafold(args.results_dir, args.manifest,
                    os.path.join(args.out_dir, 'combined_rnafold.tsv'))

    if args.include_raw_shuffles:
        collate_raw_shuffles(args.raw_dir, args.manifest,
                             os.path.join(args.out_dir,
                                          'combined_rnafold_raw_shuffles.tsv'))

if __name__ == '__main__':
    main()
