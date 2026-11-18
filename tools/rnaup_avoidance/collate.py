#!/usr/bin/env python3
"""
tools/rnaup_avoidance/collate.py
Combine per-tranche results into $TOOL_DIR/combined.tsv (joined with manifest
metadata by the shared rnalfold collator), then split it by window into
$TOOL_DIR/avoidance_<window>.tsv (e.g. avoidance_start_codon_pm30.tsv).

Required env, set by bin/05_collate.sh -> lib/paths.sh:
  RESULTS_DIR, MANIFEST_TSV, TOOL_DIR
"""
import csv
import logging
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rnalfold'))
from collate import collate_results, load_manifest  # noqa: E402


def main() -> None:
    results_dir = os.environ.get('RESULTS_DIR')
    manifest_path = os.environ.get('MANIFEST_TSV')
    tool_dir = os.environ.get('TOOL_DIR')
    if not (results_dir and manifest_path and tool_dir):
        logging.error("RESULTS_DIR, MANIFEST_TSV, and TOOL_DIR must be set "
                      "(invoke via bin/05_collate.sh).")
        sys.exit(1)
    os.makedirs(tool_dir, exist_ok=True)
    combined = os.path.join(tool_dir, 'combined.tsv')
    collate_results(results_dir, load_manifest(manifest_path), combined)

    outs, writers = {}, {}
    with open(combined, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            w = row['window']
            if w not in writers:
                outs[w] = open(os.path.join(tool_dir, f'avoidance_{w}.tsv'), 'w', newline='')
                writers[w] = csv.DictWriter(outs[w], fieldnames=reader.fieldnames,
                                            delimiter='\t', lineterminator='\n')
                writers[w].writeheader()
            writers[w].writerow(row)
    for fh in outs.values():
        fh.close()
    logging.info("Split into per-window files: %s", ', '.join(sorted(writers)))


if __name__ == '__main__':
    main()
