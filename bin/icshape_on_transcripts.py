#!/usr/bin/env python3
"""
bin/icshape_on_transcripts.py

Standalone step: project genome-coordinate icSHAPE scores onto the spliced
transcripts of a candidate list, giving a per-nucleotide reactivity track in
the same frame as extracted_mRNA.fa.

Reads
  runs/<dataset>/ultrafold/candidates.tsv
  runs/<dataset>/extracted_regions/{canonical.gff, extracted_mRNA.fa}
  --cytoplasm / --nucleoplasm   icSHAPE score BEDs (BED6: chrom start end NA
                                score strand; 0-based half-open; one row is a
                                run of identical scores). Strand-specific.

Writes (runs/<dataset>/ultrafold/)

  icshape_tracks.tsv   one row per transcript position
      transcript_id, pos (0-based, 5'->3'), base, cytoplasm, nucleoplasm
      Reactivity is NA where the BED has no score (no coverage). NA means "not
      measured", never zero.
  icshape_summary.tsv  one row per transcript x compartment x region
      transcript_id, compartment, region (5UTR/CDS/3UTR/mRNA), length,
      n_valid, coverage (n_valid/length), mean, gini
      gini is the population Gini of the valid scores, NA when n_valid <= 15
      (the cutoff of the Gini-from-icSHAPE pipeline), so it can be compared
      directly with that pipeline's gini_long.tsv.

Stdlib only.

Usage
  ./bin/icshape_on_transcripts.py -d human_all
"""
import argparse
import bisect
import csv
import logging
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.txcoords import load_transcripts  # noqa: E402

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger('icshape_on_transcripts')

NA = 'NA'
VALID_CUTOFF = 15          # matches Gini-from-icSHAPE: emit Gini only if n_valid > 15
GINI_DIR = Path.home() / 'Gardner_Lab' / 'Gini-from-icSHAPE'
COMPARTMENTS = ('cytoplasm', 'nucleoplasm')


def population_gini(scores):
    """Population Gini (no n-1 correction), identical to the Gini-from-icSHAPE
    pipeline's sorted-cumsum form. Needs a positive sum."""
    s = sorted(scores)
    n, cum, cc = len(s), 0.0, 0.0
    for x in s:
        cum += x
        cc += cum
    return (n + 1 - 2 * cc / cum) / n


def read_mrna_bases(fasta_path, wanted):
    """{tid: sequence} from extracted_mRNA.fa (header: GENE_TID_mRNA)."""
    seqs, name, buf = {}, None, []
    with open(fasta_path) as f:
        for line in list(f) + ['>']:
            line = line.strip()
            if line.startswith('>'):
                if name in wanted:
                    seqs[name] = ''.join(buf).upper()
                parts = line[1:].split('_')
                name = parts[1] if len(parts) >= 3 else None
                buf = []
            else:
                buf.append(line)
    return seqs


def merged_exons_by_chrom(txs):
    per = defaultdict(list)
    for tx in txs.values():
        per[tx.chrom].extend(tx.exons)
    out = {}
    for chrom, ivs in per.items():
        ivs.sort()
        merged = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[chrom] = ([m[0] for m in merged], [m[1] for m in merged])
    return out


def project_bed(bed_path, txs):
    """{tid: [score-or-None per tx position]} for one compartment."""
    idx = merged_exons_by_chrom(txs)
    by_chrom = defaultdict(list)
    for tx in txs.values():
        by_chrom[tx.chrom].append(tx)
    tracks = {t: [None] * tx.length for t, tx in txs.items()}
    n = hit = 0
    with open(bed_path) as f:
        for line in f:
            n += 1
            c = line.split('\t')
            chrom = c[0]
            if chrom not in idx:
                continue
            start, end = int(c[1]), int(c[2])
            starts, ends = idx[chrom]
            i = bisect.bisect_left(starts, end) - 1
            if i < 0 or ends[i] <= start:
                continue
            strand, score = c[5].strip(), float(c[4])
            for tx in by_chrom[chrom]:
                if tx.strand != strand:                   # strand-specific, as bedtools -s
                    continue
                for ts, te, _, _ in tx.map_interval(start, end):
                    tr = tracks[tx.tid]
                    for p in range(ts, te):
                        tr[p] = score
                    hit += 1
    log.info(f'{Path(bed_path).name}: scanned {n:,} rows, {hit:,} placed on candidates')
    return tracks


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-d', '--dataset', required=True)
    ap.add_argument('--runs-root', default=os.environ.get(
        'RUNS_ROOT', str(Path(__file__).resolve().parent.parent / 'runs')))
    ap.add_argument('--candidates', help='default: <run>/ultrafold/candidates.tsv')
    ap.add_argument('--cytoplasm', default=str(GINI_DIR / 'human' / 'icSHAPE_hg38_cy_invivo-score.bed'))
    ap.add_argument('--nucleoplasm', default=str(GINI_DIR / 'human' / 'icSHAPE_hg38_nucl_vivo-score.bed'))
    args = ap.parse_args()

    run_dir = Path(args.runs_root) / args.dataset
    out_dir = run_dir / 'ultrafold'
    cand_path = Path(args.candidates or out_dir / 'candidates.tsv')
    gff = run_dir / 'extracted_regions' / 'canonical.gff'
    fa = run_dir / 'extracted_regions' / 'extracted_mRNA.fa'
    beds = {'cytoplasm': Path(args.cytoplasm), 'nucleoplasm': Path(args.nucleoplasm)}
    for p in (cand_path, gff, fa, *beds.values()):
        if not p.exists():
            sys.exit(f'ERROR: not found: {p}')

    with open(cand_path, newline='') as f:
        cands = {r['transcript_id']: r for r in csv.DictReader(f, delimiter='\t')}
    txs = load_transcripts(gff, wanted=set(cands))
    bad = [t for t, tx in txs.items() if tx.length != int(cands[t]['mrna_len'])]
    if bad:
        sys.exit(f'ERROR: exon length != manifest mRNA length for {bad[:3]}...')
    bases = read_mrna_bases(fa, set(txs))
    for t in txs:
        if len(bases.get(t, '')) != txs[t].length:
            sys.exit(f'ERROR: mRNA FASTA length != exon length for {t}')

    proj = {comp: project_bed(p, txs) for comp, p in beds.items()}

    order = sorted(txs, key=lambda t: float(cands[t]['mrna_zscore']))
    fmt = lambda v: NA if v is None else f'{v:g}'
    with open(out_dir / 'icshape_tracks.tsv', 'w') as fo:
        fo.write('transcript_id\tpos\tbase\tcytoplasm\tnucleoplasm\n')
        for t in order:
            for p in range(txs[t].length):
                fo.write(f'{t}\t{p}\t{bases[t][p]}\t{fmt(proj["cytoplasm"][t][p])}'
                         f'\t{fmt(proj["nucleoplasm"][t][p])}\n')

    with open(out_dir / 'icshape_summary.tsv', 'w') as fo:
        fo.write('transcript_id\tcompartment\tregion\tlength\tn_valid\tcoverage\tmean\tgini\n')
        for t in order:
            c = cands[t]
            u5, cds, n = int(c['utr5_len']), int(c['cds_len']), int(c['mrna_len'])
            regions = {'5UTR': (0, u5), 'CDS': (u5, u5 + cds), '3UTR': (u5 + cds, n), 'mRNA': (0, n)}
            for comp in COMPARTMENTS:
                for reg, (a, b) in regions.items():
                    vals = [v for v in proj[comp][t][a:b] if v is not None]
                    L = b - a
                    mean = f'{sum(vals) / len(vals):.4f}' if vals else NA
                    gini = NA
                    if len(vals) > VALID_CUTOFF and sum(vals) > 0:
                        gini = f'{population_gini(vals):.6f}'
                    fo.write(f'{t}\t{comp}\t{reg}\t{L}\t{len(vals)}\t'
                             f'{len(vals) / L if L else 0:.4f}\t{mean}\t{gini}\n')
    log.info(f'wrote icshape_tracks.tsv and icshape_summary.tsv -> {out_dir}')


if __name__ == '__main__':
    main()
