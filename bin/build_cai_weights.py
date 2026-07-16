#!/usr/bin/env python3
"""
build_cai_weights.py — Build a CAI codon-weight table from a CoCoPUTs TSV.

CoCoPUTs data: https://hive.biochemistry.gwu.edu/review/codon
Select a species, download the CDS table (tab-delimited, all genes).

This script filters to ribosomal protein genes (RPL*/RPS*) as the reference
set, deduplicates isoforms, aggregates codon counts, and writes RSCU and
adaptation weights (w) for all 61 sense codons.

The output TSV is consumed by the 'cai' metrics plugin. Re-run this script
any time the CoCoPUTs source file changes.

Usage:
    ./bin/build_cai_weights.py \\
        --cocoputs  data/references/human/Human_CDS.tsv \\
        --output    data/references/human/cai_weights.tsv

    # Lower min-codons threshold if the reference set is small:
    ./bin/build_cai_weights.py \\
        --cocoputs  data/references/mouse/Mouse_CDS.tsv \\
        --output    data/references/mouse/cai_weights.tsv \\
        --min-codons 50
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

_THIS = Path(__file__).resolve()
_PROJECT_ROOT = _THIS.parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from metrics._cai_weights import build_weights  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('build_cai_weights')


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--cocoputs', '-c', required=True, metavar='TSV',
        help='Path to CoCoPUTs CDS table (tab-delimited).',
    )
    parser.add_argument(
        '--output', '-o', required=True, metavar='TSV',
        help='Path to write the weight table (codon, rscu, w).',
    )
    parser.add_argument(
        '--min-codons', type=int, default=100, metavar='N',
        help='Minimum # Codons per CoCoPUTs entry (default: 100).',
    )
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    cocoputs = Path(args.cocoputs)
    if not cocoputs.is_absolute():
        cocoputs = (_PROJECT_ROOT / cocoputs).resolve()
    if not cocoputs.exists():
        log.error("CoCoPUTs file not found: %s", cocoputs)
        sys.exit(1)

    output = Path(args.output)
    if not output.is_absolute():
        output = (_PROJECT_ROOT / output).resolve()

    try:
        build_weights(cocoputs, output, min_codons=args.min_codons)
    except Exception as e:
        log.error("Failed to build weights: %s", e)
        sys.exit(1)

    log.info("Done. Weight table written to %s", output)


if __name__ == '__main__':
    main()
