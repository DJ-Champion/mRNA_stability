#!/usr/bin/env python3
"""
tools/rnafold/collate.py
Combine per-sequence RNAfold result CSVs into one TSV per dataset, joined
with the manifest so output has gene/transcript/region context.

Required env (set by bin/05_collate.sh -> resolve_paths):
  RESULTS_DIR, MANIFEST_TSV, TOOL_DIR
Optional flags (after pipeline args):
  --include-raw-shuffles   also concatenate raw shuffle CSVs (large)

Outputs:
  $TOOL_DIR/combined.tsv                   one row per sequence
  $TOOL_DIR/combined_raw_shuffles.tsv      optional, large
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


def load_manifest(path):
    if not os.path.exists(path):
        logging.warning(f"Manifest not found: {path} — output will lack metadata")
        return {}
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    logging.info(f"Loaded manifest: {len(rows)} records")
    return {row['seqname']: row for row in rows}


def collate_results(results_dir, manifest, out_path):
    pattern = os.path.join(results_dir, "results_*.csv")
    files = sorted(glob(pattern))
    if not files:
        logging.error(f"No result files matched: {pattern}")
        sys.exit(1)
    logging.info(f"Found {len(files)} result files")

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


def collate_raw_shuffles(raw_dir, manifest, out_path):
    pattern = os.path.join(raw_dir, "*_raw_shuffles.csv")
    files = sorted(glob(pattern))
    if not files:
        logging.warning(f"No raw shuffle files matched: {pattern}")
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
    parser.add_argument('--include-raw-shuffles', action='store_true',
                        help='Also collate raw shuffle CSVs (large output)')
    args, _ = parser.parse_known_args()

    results_dir = os.environ.get('RESULTS_DIR')
    manifest_path = os.environ.get('MANIFEST_TSV')
    tool_dir = os.environ.get('TOOL_DIR')
    if not (results_dir and manifest_path and tool_dir):
        logging.error("RESULTS_DIR, MANIFEST_TSV, and TOOL_DIR must be set "
                      "(invoke via bin/05_collate.sh).")
        sys.exit(1)

    os.makedirs(tool_dir, exist_ok=True)
    manifest = load_manifest(manifest_path)

    collate_results(results_dir, manifest,
                    os.path.join(tool_dir, 'combined.tsv'))

    if args.include_raw_shuffles:
        raw_dir = os.path.join(tool_dir, 'raw_shuffles')
        collate_raw_shuffles(raw_dir, manifest,
                             os.path.join(tool_dir, 'combined_raw_shuffles.tsv'))


if __name__ == '__main__':
    main()