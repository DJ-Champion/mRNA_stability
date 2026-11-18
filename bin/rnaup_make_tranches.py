#!/usr/bin/env python3
"""
bin/rnaup_make_tranches.py

Standalone step before `04_submit.sh -t rnaup_avoidance`: cut the extracted
region FASTAs into window FASTAs (strict length, no truncated windows) and pack
them into tranches of M genes, one multi-record FASTA each, for
tools/rnaup_avoidance. Uses what 01_extract.py already wrote; nothing is re-extracted.

Windows (--regions; default start_codon_pm30)
  start_codon_pm30  extracted_start_codon_region.fa  exactly 63 nt (needs 5'UTR >= 30)
  stop_codon_pm30   extracted_stop_codon_region.fa   exactly 63 nt (needs 3'UTR >= 30)
  last_100nt        extracted_tail_region.fa         exactly 100 nt
  utr5              whole 5'UTR, extracted_5UTR.fa; length varies, >= --min-utr5-length (default 20)
The pm30 / 100 nt names assume the dataset's codon_flank_length 30 and
tail_length 100 (the defaults). Windows of any other length are dropped; if none
survive for a window the script stops and says why. Windows with characters
other than A/C/G/T/U are dropped too.

Writes (under runs/<dataset>/)
  rnaup_avoidance/tranches/<window>/<window>_<NNNN>.fa
  lists/rnaup_<window>.txt          tranche paths, for 04_submit.sh --list

Stdlib only.

Usage
  ./bin/rnaup_make_tranches.py -d human_all
  ./bin/rnaup_make_tranches.py -d human_all --regions start_codon_pm30 utr5 -M 10
"""
import argparse
import csv
import logging
import os
import shutil
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger('rnaup_make_tranches')

# window -> (source extracted region, required window length, slice length or None)
WINDOWS = {
    'start_codon_pm30': ('start_codon_region', 63, None),
    'stop_codon_pm30': ('stop_codon_region', 63, None),
    'last_100nt': ('tail_region', 100, None),
    'utr5': ('5UTR', None, None),   # variable length
}


def iter_fasta(path):
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        yield name, ''.join(buf)


def window_records(extract_dir, window, min_len):
    src, need, _ = WINDOWS[window]
    path = extract_dir / f'extracted_{src}.fa'
    if not path.exists():
        sys.exit(f'ERROR: {path} not found; add {src} to regions_to_extract and rerun 01_extract.py')
    kept, short, odd, total = [], 0, 0, 0
    for name, seq in iter_fasta(path):
        total += 1
        seq = seq.upper()
        if (len(seq) != need) if need else (len(seq) < min_len):
            short += 1
        elif set(seq) - set('ACGTU'):
            odd += 1
        else:
            kept.append((name, seq))
    rule = f'not exactly {need} nt' if need else f'shorter than {min_len} nt'
    log.info('%s: %d records in %s; kept %d, %d %s, %d with non-ACGTU',
             window, total, path.name, len(kept), short, rule, odd)
    if not kept:
        sys.exit(f'ERROR: no usable {window} windows (check codon_flank_length / tail_length / --min-utr5-length)')
    return kept


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-d', '--dataset', required=True)
    ap.add_argument('--regions', nargs='+', default=['start_codon_pm30'], choices=sorted(WINDOWS))
    ap.add_argument('-M', '--genes-per-tranche', type=int, default=10)
    ap.add_argument('--min-utr5-length', type=int, default=20,
                    help='shortest whole 5\'UTR kept (default 20; tiny sequences have almost no distinct shuffles)')
    ap.add_argument('--max-tranche-nt', type=int, default=1500,
                    help='also close a tranche once its windows total this many nt (default 1500; '
                         'only matters for variable-length windows)')
    ap.add_argument('--runs-root', default=os.environ.get(
        'RUNS_ROOT', str(Path(__file__).resolve().parent.parent / 'runs')))
    args = ap.parse_args()

    run_dir = Path(args.runs_root) / args.dataset
    extract_dir = run_dir / 'extracted_regions'
    (run_dir / 'lists').mkdir(parents=True, exist_ok=True)
    for window in args.regions:
        recs = window_records(extract_dir, window, args.min_utr5_length)
        if WINDOWS[window][1] is None:
            recs.sort(key=lambda r: len(r[1]))    # similar lengths per tranche -> predictable run time

        out_dir = run_dir / 'rnaup_avoidance' / 'tranches' / window
        if out_dir.exists():
            shutil.rmtree(out_dir)       # tranche contents depend on M; never mix two layouts
        out_dir.mkdir(parents=True)
        paths = []
        tranches, cur, nt = [], [], 0
        for rec in recs:
            if cur and (len(cur) >= args.genes_per_tranche or nt + len(rec[1]) > args.max_tranche_nt):
                tranches.append(cur)
                cur, nt = [], 0
            cur.append(rec)
            nt += len(rec[1])
        tranches.append(cur)
        for t, chunk in enumerate(tranches, 1):
            p = out_dir / f'{window}_{t:04d}.fa'
            with open(p, 'w') as f:
                for name, seq in chunk:
                    f.write(f'>{name}\n{seq}\n')
            paths.append(str(p))
        lst = run_dir / 'lists' / f'rnaup_{window}.txt'
        lst.write_text('\n'.join(paths) + '\n')
        log.info('%s: %d genes -> %d tranches of <= %d; list %s',
                 window, len(recs), len(paths), args.genes_per_tranche, lst)


if __name__ == '__main__':
    main()
