#!/usr/bin/env python3
"""
bin/ultrafold_gene_lists.py

Standalone step: inputs for a GO enrichment of the ultrafolded genes, run in R
(clusterProfiler enrichGO / gseGO, or g:Profiler).

Reads
  --scores   RNAfold results CSV (as for ultrafold_candidates.py)
  runs/<dataset>/extracted_regions/{manifest.tsv, canonical.gff}

Writes (runs/<dataset>/ultrafold/gene_lists/)

  universe.tsv    every gene that was scored on --region: the BACKGROUND for
                  enrichment. Use this, not the whole genome, as the universe.
      gene_id (Ensembl, no version), gene_name, transcript_id,
      mrna_zscore, mrna_len, utr3_len
  ranked.tsv      universe sorted most-ultrafolded first, for ranked tests
      (adds rank_metric = -mrna_zscore, so larger = more ultrafolded)
  top_z<T>.txt          Ensembl gene IDs, one per line, z < T (no header)
  top_z<T>.symbols.txt  the same genes as symbols

Caveat for interpretation: ultrafolded genes have very long 3'UTRs, and long
transcripts are enriched for particular GO terms regardless of folding. Compare
against a length-matched control (universe.tsv carries mrna_len / utr3_len for
that) before reading an enrichment as being about folding.

Stdlib only.

Usage
  ./bin/ultrafold_gene_lists.py -d human_all --scores <rnafold_results.csv>
  ./bin/ultrafold_gene_lists.py -d human_all --scores ... --thresholds -40 -30
"""
import argparse
import csv
import logging
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ultrafold_candidates import read_transcript_loci  # noqa: E402

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger('ultrafold_gene_lists')


def tag(t):
    """-30.0 -> 'm30'; filename-safe threshold label."""
    return ('m' if t < 0 else 'p') + f'{abs(t):g}'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-d', '--dataset', required=True)
    ap.add_argument('--scores', required=True)
    ap.add_argument('--region', default='MRNA')
    ap.add_argument('--thresholds', type=float, nargs='+', default=[-30, -25, -20])
    ap.add_argument('--runs-root', default=os.environ.get(
        'RUNS_ROOT', str(Path(__file__).resolve().parent.parent / 'runs')))
    args = ap.parse_args()

    run_dir = Path(args.runs_root) / args.dataset
    manifest, gff = (run_dir / 'extracted_regions' / n for n in ('manifest.tsv', 'canonical.gff'))
    for p in (manifest, gff):
        if not p.exists():
            sys.exit(f'ERROR: not found: {p}')

    z = {}
    with open(args.scores, newline='') as f:
        for row in csv.DictReader(f):
            if row['region'] != args.region:
                continue
            try:
                z[row['transcript_id']] = float(row['RNAfold_zscore'])
            except (TypeError, ValueError):
                pass
    log.info(f'{len(z)} transcripts with a {args.region} z-score')

    gene_of, lens = {}, {}
    with open(manifest, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            t = row['transcript_id']
            if t in z:
                gene_of[t] = row['gene_id']
                lens.setdefault(t, {})[row['region']] = row['length']
    loci = read_transcript_loci(gff, set(gene_of))

    rows, seen = [], {}
    for t, zt in z.items():
        if t not in gene_of:
            continue
        g = gene_of[t]
        if g in seen:                       # canonical.gff is one transcript/gene; guard anyway
            log.warning(f'gene {g} appears twice ({seen[g]}, {t}); keeping the first')
            continue
        seen[g] = t
        rows.append(dict(gene_id=g, gene_name=loci.get(t, ('NA',))[0], transcript_id=t,
                         mrna_zscore=zt, mrna_len=lens[t].get('mRNA', 'NA'),
                         utr3_len=lens[t].get('3UTR', 'NA')))
    log.info(f'{len(rows)} genes in universe '
             f'({len(z) - len(rows)} scored transcripts absent from manifest/dropped)')
    rows.sort(key=lambda r: r['mrna_zscore'])

    out = run_dir / 'ultrafold' / 'gene_lists'
    out.mkdir(parents=True, exist_ok=True)
    cols = ['gene_id', 'gene_name', 'transcript_id', 'mrna_zscore', 'mrna_len', 'utr3_len']
    with open(out / 'universe.tsv', 'w') as fu, open(out / 'ranked.tsv', 'w') as fr:
        fu.write('\t'.join(cols) + '\n')
        fr.write('\t'.join(cols + ['rank_metric']) + '\n')
        # universe in gene_id order (no ranking implied); ranked in z order
        for r in sorted(rows, key=lambda r: r['gene_id']):
            fu.write('\t'.join(str(r[c]) for c in cols) + '\n')
        for r in rows:
            fr.write('\t'.join(str(r[c]) for c in cols) + f'\t{-r["mrna_zscore"]}\n')

    for thr in args.thresholds:
        sel = [r for r in rows if r['mrna_zscore'] < thr]
        base = out / f'top_z{tag(thr)}'
        base.with_suffix('.txt').write_text(''.join(r['gene_id'] + '\n' for r in sel))
        Path(str(base) + '.symbols.txt').write_text(
            ''.join(r['gene_name'] + '\n' for r in sel if r['gene_name'] != 'NA'))
        log.info(f'z < {thr:g}: {len(sel)} genes -> {base.name}.txt')
    log.info(f'wrote {out}')


if __name__ == '__main__':
    main()
