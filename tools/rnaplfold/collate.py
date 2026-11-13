#!/usr/bin/env python3
"""
tools/rnaplfold/collate.py
Combine per-sequence RNAplfold result CSVs into $TOOL_DIR/combined.tsv, joined
with manifest metadata. Same layout as the rnalfold collator, which it reuses.

Required env, set by bin/05_collate.sh -> lib/paths.sh:
  RESULTS_DIR, MANIFEST_TSV, TOOL_DIR
"""
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
    collate_results(results_dir, load_manifest(manifest_path),
                    os.path.join(tool_dir, 'combined.tsv'))


if __name__ == '__main__':
    main()
