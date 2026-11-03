#!/usr/bin/env python3
"""
bin/ultrafold_candidates.py

Standalone step: list "ultrafolded" transcripts (extreme negative RNAfold
z-score) with the context later steps need.

Reads
  --scores   RNAfold results CSV (transcript_id, region, RNAfold_score,
             RNAfold_median, RNAfold_pval, RNAfold_zscore)
  runs/<dataset>/extracted_regions/{manifest.tsv, canonical.gff}

Writes
  runs/<dataset>/ultrafold/candidates.tsv      one row per candidate, sorted
                                               by z-score (most negative first)

Columns
  transcript_id, gene_id, gene_name, chrom, strand,
  mrna_len, utr5_len, cds_len, utr3_len,
  mrna_score, mrna_zscore,        RNAfold on the spliced transcript
  utr3_score, utr3_zscore         RNAfold on the 3'UTR alone (NA if absent)

Candidates are chosen on --region (default MRNA) alone; the 3'UTR columns are
context, not a filter. Stdlib only.

Usage
  ./bin/ultrafold_candidates.py -d human_all --scores path/to/rnafold_results.csv
  ./bin/ultrafold_candidates.py -d human_all --scores ... --threshold -25
  RUNS_ROOT=/path/to/runs ./bin/ultrafold_candidates.py ...
"""
import argparse
import csv
import logging
import os
import re
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger('ultrafold_candidates')

NA = 'NA'


def read_scores(path, region, threshold):
    """Return ({tid: row} for the chosen region below threshold,
               {tid: row} for 3UTR). Rows with a missing z-score are skipped."""
    hits, utr3 = {}, {}
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            reg = row['region']
            if reg not in (region, '3UTR'):
                continue
            try:
                z = float(row['RNAfold_zscore'])
            except (TypeError, ValueError):
                continue
            if reg == '3UTR':
                utr3[row['transcript_id']] = row
            if reg == region and z < threshold:
                hits[row['transcript_id']] = row
    return hits, utr3


def read_region_lengths(manifest_path, wanted):
    """{tid: {region: length}} for the wanted transcripts, from manifest.tsv."""
    out = {}
    with open(manifest_path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            tid = row['transcript_id']
            if tid in wanted:
                out.setdefault(tid, {})[row['region']] = row['length']
    return out


_ATTR_RE = re.compile(r'(\w+)=([^;]*)')


def read_transcript_loci(gff_path, wanted):
    """{tid: (gene_name, chrom, strand)} from 'transcript' lines of the GFF."""
    out = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9 or cols[2] != 'transcript':
                continue
            attrs = dict(_ATTR_RE.findall(cols[8]))
            tid = attrs.get('transcript_id') or attrs.get('ID', '')
            if tid in wanted:
                out[tid] = (attrs.get('gene_name', NA), cols[0], cols[6])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-d', '--dataset', required=True)
    ap.add_argument('--scores', required=True, help='RNAfold results CSV')
    ap.add_argument('--region', default='MRNA', help='region label to filter on (default MRNA)')
    ap.add_argument('--threshold', type=float, default=-30.0,
                    help='keep z-score strictly below this (default -30)')
    ap.add_argument('--runs-root', default=os.environ.get(
        'RUNS_ROOT', str(Path(__file__).resolve().parent.parent / 'runs')))
    args = ap.parse_args()

    run_dir = Path(args.runs_root) / args.dataset
    extract = run_dir / 'extracted_regions'
    manifest, gff = extract / 'manifest.tsv', extract / 'canonical.gff'
    for p in (manifest, gff):
        if not p.exists():
            sys.exit(f'ERROR: not found: {p}')

    hits, utr3 = read_scores(args.scores, args.region, args.threshold)
    log.info(f'{len(hits)} transcripts with {args.region} z < {args.threshold}')

    lengths = read_region_lengths(manifest, set(hits))
    loci = read_transcript_loci(gff, set(hits))

    missing = [t for t in hits if t not in lengths or t not in loci]
    if missing:
        log.warning(f'{len(missing)} candidates absent from manifest/GFF '
                    f'(dropped): {", ".join(missing[:5])}...')

    out_dir = run_dir / 'ultrafold'
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / 'candidates.tsv'
    cols = ['transcript_id', 'gene_id', 'gene_name', 'chrom', 'strand',
            'mrna_len', 'utr5_len', 'cds_len', 'utr3_len',
            'mrna_score', 'mrna_zscore', 'utr3_score', 'utr3_zscore']
    gene_of = {}
    with open(manifest, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['transcript_id'] in hits:
                gene_of[row['transcript_id']] = row['gene_id']

    n = 0
    with open(out, 'w') as fo:
        fo.write('\t'.join(cols) + '\n')
        for tid in sorted(hits, key=lambda t: float(hits[t]['RNAfold_zscore'])):
            if tid not in lengths or tid not in loci:
                continue
            ln, (gname, chrom, strand) = lengths[tid], loci[tid]
            u3 = utr3.get(tid, {})
            fo.write('\t'.join([
                tid, gene_of[tid], gname, chrom, strand,
                ln.get('mRNA', NA), ln.get('5UTR', NA), ln.get('CDS', NA), ln.get('3UTR', NA),
                hits[tid]['RNAfold_score'], hits[tid]['RNAfold_zscore'],
                u3.get('RNAfold_score', NA), u3.get('RNAfold_zscore', NA),
            ]) + '\n')
            n += 1
    log.info(f'wrote {n} rows -> {out}')


if __name__ == '__main__':
    main()
