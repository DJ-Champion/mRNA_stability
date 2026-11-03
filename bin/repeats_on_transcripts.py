#!/usr/bin/env python3
"""
bin/repeats_on_transcripts.py

Standalone step: place RepeatMasker annotations on the spliced transcripts of
a candidate list, tag each as sense / antisense relative to the transcript,
and find the oppositely oriented Alu pairs that can base-pair in cis.

Reads
  runs/<dataset>/ultrafold/candidates.tsv          (from ultrafold_candidates.py)
  runs/<dataset>/extracted_regions/canonical.gff   (exon structure)
  runs/_resources/hg38_rmsk/rmsk.bed.gz            (from fetch_rmsk.sh)

Writes (all under runs/<dataset>/ultrafold/)

  repeats.tsv      one row per repeat PIECE lying in exons (all repeat classes)
      transcript_id, repeat_id, piece, n_pieces,
      rep_name, rep_class, rep_family, pct_div, sw_score,
      tx_start, tx_end            0-based half-open, 5'->3' on the spliced
                                  transcript (same frame as extracted_mRNA.fa)
      region                      5UTR / CDS / 3UTR (majority overlap)
      orientation                 sense | antisense  (repeat strand vs transcript)
      g_chrom, g_start, g_end, g_strand   genomic, 0-based half-open
      repeat_len, piece_len       full repeat length / bases in this exon piece

  alu_pairs.tsv    one row per oppositely oriented Alu pair in --pair-region
      transcript_id, alu_a, alu_b (repeat_ids), name_a, name_b,
      a_start, a_end, b_start, b_end   (a is 5' of b)
      gap                         nt between the two (0 if they overlap)
      orient_a, orient_b          always differ; kept for convenience

  alu_summary.tsv  one row per candidate
      transcript_id, gene_name, region_len, n_alu, n_sense, n_antisense,
      alu_bp, alu_frac, n_opposite_pairs, min_pair_gap, n_alu_elsewhere

Notes
  * Only exonic bases are kept: a repeat wholly inside an intron is dropped,
    and one that spans a splice junction yields one row per exon (piece), all
    sharing repeat_id. Pieces shorter than --min-piece are dropped (default 10)
    to ignore repeats that only graze an exon edge.
  * Pairing counts treat a repeat as one element (its pieces' transcript span).
  * "Opposite pair" is a sequence-orientation criterion only. It does not
    assert that the pair actually pairs in the predicted structure.

Stdlib only.

Usage
  ./bin/repeats_on_transcripts.py -d human_all
  RUNS_ROOT=/path/to/runs ./bin/repeats_on_transcripts.py -d human_all
"""
import argparse
import bisect
import csv
import gzip
import logging
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.txcoords import load_transcripts, region_of  # noqa: E402

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger('repeats_on_transcripts')


def read_candidates(path):
    with open(path, newline='') as f:
        return {r['transcript_id']: r for r in csv.DictReader(f, delimiter='\t')}


def merged_exons_by_chrom(txs):
    """{chrom: (starts, ends)} - sorted, merged union of all exons, for a fast
    bisect pre-filter while streaming the (5M-row) repeat file."""
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


def overlaps_merged(idx, start, end):
    starts, ends = idx
    i = bisect.bisect_left(starts, end) - 1      # last block starting before `end`
    return i >= 0 and ends[i] > start


def collect_repeats(rmsk_path, txs, min_piece):
    """Stream the repeat BED; return a list of per-piece dicts (tx-space)."""
    idx = merged_exons_by_chrom(txs)
    by_chrom = defaultdict(list)
    for tx in txs.values():
        by_chrom[tx.chrom].append(tx)

    rows, n_seen = [], 0
    opener = gzip.open if str(rmsk_path).endswith('.gz') else open
    with opener(rmsk_path, 'rt') as f:
        for line in f:
            n_seen += 1
            c = line.rstrip('\n').split('\t')
            chrom = c[0]
            if chrom not in idx:
                continue
            start, end = int(c[1]), int(c[2])
            if not overlaps_merged(idx[chrom], start, end):
                continue
            rep_id = f'{chrom}:{start}-{end}:{c[3]}'
            for tx in by_chrom[chrom]:
                lo, hi = tx.span
                if end <= lo or start >= hi:
                    continue
                pieces = [p for p in tx.map_interval(start, end) if p[1] - p[0] >= min_piece]
                for k, (ts, te, gs, ge) in enumerate(pieces):
                    # repeat strand vs transcript strand -> orientation in tx frame
                    orient = 'sense' if c[5] == tx.strand else 'antisense'
                    rows.append(dict(
                        transcript_id=tx.tid, repeat_id=rep_id, piece=k + 1,
                        n_pieces=len(pieces), rep_name=c[3], rep_class=c[6],
                        rep_family=c[7], pct_div=c[8], sw_score=c[4],
                        tx_start=ts, tx_end=te, orientation=orient,
                        g_chrom=chrom, g_start=gs, g_end=ge, g_strand=c[5],
                        repeat_len=end - start, piece_len=te - ts))
    log.info(f'scanned {n_seen:,} repeat rows; kept {len(rows):,} exonic pieces')
    return rows


def write_tsv(path, cols, rows):
    with open(path, 'w') as fo:
        fo.write('\t'.join(cols) + '\n')
        for r in rows:
            fo.write('\t'.join(str(r[c]) for c in cols) + '\n')


def alu_elements(rows, family):
    """Collapse pieces to one element per Alu repeat_id: tx span + covered bp."""
    el = {}
    for r in rows:
        if r['rep_family'] != family:
            continue
        e = el.setdefault(r['repeat_id'], dict(
            repeat_id=r['repeat_id'], name=r['rep_name'], orientation=r['orientation'],
            region=r['region'], start=r['tx_start'], end=r['tx_end'], bp=0))
        e['start'] = min(e['start'], r['tx_start'])
        e['end'] = max(e['end'], r['tx_end'])
        e['bp'] += r['piece_len']
    return list(el.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-d', '--dataset', required=True)
    ap.add_argument('--runs-root', default=os.environ.get(
        'RUNS_ROOT', str(Path(__file__).resolve().parent.parent / 'runs')))
    ap.add_argument('--candidates', help='default: <run>/ultrafold/candidates.tsv')
    ap.add_argument('--rmsk', help='default: <runs>/_resources/hg38_rmsk/rmsk.bed.gz')
    ap.add_argument('--min-piece', type=int, default=10,
                    help='drop exonic repeat pieces shorter than this (default 10)')
    ap.add_argument('--alu-family', default='Alu', help='repFamily counted as Alu')
    ap.add_argument('--pair-region', default='3UTR', choices=['5UTR', 'CDS', '3UTR'],
                    help='region both Alus must lie in to be paired (default 3UTR)')
    args = ap.parse_args()

    runs = Path(args.runs_root)
    run_dir = runs / args.dataset
    out_dir = run_dir / 'ultrafold'
    cand_path = Path(args.candidates or out_dir / 'candidates.tsv')
    rmsk_path = Path(args.rmsk or runs / '_resources' / 'hg38_rmsk' / 'rmsk.bed.gz')
    gff_path = run_dir / 'extracted_regions' / 'canonical.gff'
    for p in (cand_path, rmsk_path, gff_path):
        if not p.exists():
            sys.exit(f'ERROR: not found: {p}')

    cands = read_candidates(cand_path)
    txs = load_transcripts(gff_path, wanted=set(cands))
    log.info(f'{len(txs)}/{len(cands)} candidate transcripts found in GFF')

    # Frame check: spliced exon length must equal the extracted mRNA length,
    # else the region boundaries from the manifest do not apply to our coordinates.
    bad = [t for t, tx in txs.items() if tx.length != int(cands[t]['mrna_len'])]
    if bad:
        sys.exit(f'ERROR: exon length != manifest mRNA length for {len(bad)} '
                 f'transcripts (e.g. {bad[:3]}); refusing to continue')

    rows = collect_repeats(rmsk_path, txs, args.min_piece)
    for r in rows:
        c = cands[r['transcript_id']]
        r['region'] = region_of(r['tx_start'], r['tx_end'], int(c['utr5_len']),
                                int(c['cds_len']), int(c['mrna_len']))
    rows.sort(key=lambda r: (r['transcript_id'], r['tx_start'], r['tx_end']))

    out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(out_dir / 'repeats.tsv',
              ['transcript_id', 'repeat_id', 'piece', 'n_pieces', 'rep_name', 'rep_class',
               'rep_family', 'pct_div', 'sw_score', 'tx_start', 'tx_end', 'region',
               'orientation', 'g_chrom', 'g_start', 'g_end', 'g_strand',
               'repeat_len', 'piece_len'], rows)

    by_tx = defaultdict(list)
    for r in rows:
        by_tx[r['transcript_id']].append(r)

    pair_rows, summary = [], []
    for tid, c in cands.items():
        if tid not in txs:
            continue
        els = alu_elements(by_tx.get(tid, []), args.alu_family)
        inreg = sorted((e for e in els if e['region'] == args.pair_region),
                       key=lambda e: (e['start'], e['end']))
        region_len = {'5UTR': c['utr5_len'], 'CDS': c['cds_len'], '3UTR': c['utr3_len']}[args.pair_region]
        pairs = []
        for i, a in enumerate(inreg):
            for b in inreg[i + 1:]:
                if a['orientation'] == b['orientation']:
                    continue
                gap = max(0, b['start'] - a['end'])
                pairs.append(dict(
                    transcript_id=tid, alu_a=a['repeat_id'], alu_b=b['repeat_id'],
                    name_a=a['name'], name_b=b['name'], a_start=a['start'], a_end=a['end'],
                    b_start=b['start'], b_end=b['end'], gap=gap,
                    orient_a=a['orientation'], orient_b=b['orientation']))
        pair_rows.extend(pairs)
        n_sense = sum(e['orientation'] == 'sense' for e in inreg)
        alu_bp = sum(e['bp'] for e in inreg)
        summary.append(dict(
            transcript_id=tid, gene_name=c['gene_name'], region_len=region_len,
            n_alu=len(inreg), n_sense=n_sense, n_antisense=len(inreg) - n_sense,
            alu_bp=alu_bp, alu_frac=f'{alu_bp / int(region_len):.4f}' if int(region_len) else 'NA',
            n_opposite_pairs=len(pairs),
            min_pair_gap=min((p['gap'] for p in pairs), default='NA'),
            n_alu_elsewhere=len(els) - len(inreg)))

    write_tsv(out_dir / 'alu_pairs.tsv',
              ['transcript_id', 'alu_a', 'alu_b', 'name_a', 'name_b', 'a_start', 'a_end',
               'b_start', 'b_end', 'gap', 'orient_a', 'orient_b'], pair_rows)
    summary.sort(key=lambda r: -r['n_opposite_pairs'])
    write_tsv(out_dir / 'alu_summary.tsv',
              ['transcript_id', 'gene_name', 'region_len', 'n_alu', 'n_sense', 'n_antisense',
               'alu_bp', 'alu_frac', 'n_opposite_pairs', 'min_pair_gap', 'n_alu_elsewhere'],
              summary)
    log.info(f'wrote repeats.tsv ({len(rows):,} rows), alu_pairs.tsv ({len(pair_rows):,}), '
             f'alu_summary.tsv ({len(summary)}) -> {out_dir}')


if __name__ == '__main__':
    main()
