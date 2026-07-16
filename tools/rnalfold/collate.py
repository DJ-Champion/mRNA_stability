#!/usr/bin/env python3
"""
tools/rnalfold/collate.py
Combine per-sequence RNAlfold result CSVs into one TSV per dataset, joined
with manifest metadata.

Required env, set by bin/05_collate.sh -> lib/paths.sh:
  RESULTS_DIR, MANIFEST_TSV, TOOL_DIR

Optional flag:
  --include-raw-shuffles

Outputs:
  $TOOL_DIR/combined.tsv
  $TOOL_DIR/combined_raw_shuffles.tsv      optional, large
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
from glob import glob

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S',
)

MANIFEST_COLS = ['gene_id', 'transcript_id', 'region', 'length',
                 'selection_reason', 'short_utrs']
RAW_MANIFEST_COLS = ['gene_id', 'transcript_id', 'region']


def load_manifest(path: str) -> dict[str, dict[str, str]]:
    if not os.path.exists(path):
        logging.warning("Manifest not found: %s — output will lack metadata", path)
        return {}
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    logging.info("Loaded manifest: %d records", len(rows))
    return {row['seqname']: row for row in rows}


def collate_results(results_dir: str, manifest: dict[str, dict[str, str]], out_path: str) -> None:
    pattern = os.path.join(results_dir, 'results_*.csv')
    files = sorted(glob(pattern))
    if not files:
        logging.error("No result files matched: %s", pattern)
        sys.exit(1)
    logging.info("Found %d result files", len(files))

    written = 0
    with open(out_path, 'w', newline='') as fout:
        writer = None
        for fpath in files:
            with open(fpath, newline='') as fin:
                reader = csv.DictReader(fin)
                if not reader.fieldnames or 'seq_name' not in reader.fieldnames:
                    logging.warning("Skipping malformed result file: %s", fpath)
                    continue
                result_cols = [c for c in reader.fieldnames if c != 'seq_name']
                for row in reader:
                    seqname = row['seq_name']
                    meta = manifest.get(seqname, {})
                    out_row = {'seqname': seqname}
                    for col in MANIFEST_COLS:
                        out_row[col] = meta.get(col, '')
                    for col in result_cols:
                        out_row[col] = row.get(col, '')

                    if writer is None:
                        fieldnames = ['seqname'] + MANIFEST_COLS + result_cols
                        writer = csv.DictWriter(fout, fieldnames=fieldnames,
                                                delimiter='\t', lineterminator='\n')
                        writer.writeheader()
                    writer.writerow(out_row)
                    written += 1

    if written == 0:
        logging.error("No rows written to %s", out_path)
        sys.exit(1)
    logging.info("Wrote %d rows to %s", written, out_path)


def collate_raw_shuffles(raw_dir: str, manifest: dict[str, dict[str, str]], out_path: str) -> None:
    pattern = os.path.join(raw_dir, '*_raw_shuffles.csv')
    files = sorted(glob(pattern))
    if not files:
        logging.warning("No raw shuffle files matched: %s", pattern)
        return
    logging.info("Found %d raw shuffle files", len(files))

    written = 0
    with open(out_path, 'w', newline='') as fout:
        writer = csv.DictWriter(
            fout,
            fieldnames=['seqname'] + RAW_MANIFEST_COLS + ['iteration', 'shuffle_min_local_mfe'],
            delimiter='\t',
            lineterminator='\n',
        )
        writer.writeheader()

        for fpath in files:
            with open(fpath, newline='') as fin:
                reader = csv.DictReader(fin)
                required = {'seq_name', 'iteration', 'shuffle_min_local_mfe'}
                if not reader.fieldnames or not required.issubset(reader.fieldnames):
                    logging.warning("Skipping malformed raw shuffle file: %s", fpath)
                    continue
                for row in reader:
                    seqname = row['seq_name']
                    meta = manifest.get(seqname, {})
                    out_row = {'seqname': seqname}
                    for col in RAW_MANIFEST_COLS:
                        out_row[col] = meta.get(col, '')
                    out_row['iteration'] = row.get('iteration', '')
                    out_row['shuffle_min_local_mfe'] = row.get('shuffle_min_local_mfe', '')
                    writer.writerow(out_row)
                    written += 1

    logging.info("Wrote %d rows to %s", written, out_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Collate RNAlfold per-sequence outputs")
    parser.add_argument('--include-raw-shuffles', action='store_true',
                        help='Also collate raw shuffle CSVs; this can be large')
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

    collate_results(results_dir, manifest, os.path.join(tool_dir, 'combined.tsv'))

    if args.include_raw_shuffles:
        collate_raw_shuffles(
            os.path.join(tool_dir, 'raw_shuffles'),
            manifest,
            os.path.join(tool_dir, 'combined_raw_shuffles.tsv'),
        )


if __name__ == '__main__':
    main()
