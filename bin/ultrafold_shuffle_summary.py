#!/usr/bin/env python3
"""
bin/ultrafold_shuffle_summary.py

Standalone step: tables for the native-vs-shuffled MFE figure, from the output
of an rnafold run that kept its raw shuffles (tools/rnafold/worker.sh writes
them; use tool rnafold_shuffles, or the original rnafold run).

Reads (under runs/<dataset>/<--tool-dir>/)
  results/results_<seqname>.csv                 native MFE, shuffle stats
  raw_shuffles/<seqname>_raw_shuffles.csv       one row per shuffle
  optional --scores CSV                         stored RNAfold z-score, to compare

Writes (runs/<dataset>/ultrafold/)
  shuffle_distributions.tsv   one row per shuffle
      seqname, transcript_id, iteration, shuffled_mfe, shuffle_z
  shuffle_summary.tsv         one row per transcript
      seqname, transcript_id, native_mfe, n_shuffles, shuffle_median, shuffle_mad,
      shuffle_mean, shuffle_sd, native_z, empirical_p, shuffle_skew,
      shuffle_excess_kurtosis, stored_z (if --scores)

z-score definition (the same as the worker, so numbers are comparable):
    z(x) = (x - median) / (1.4826 * MAD)          median, MAD over the shuffles
native_z is z(native MFE). shuffle_z is z(shuffled MFE) over the same shuffles.
Each shuffle is standardised against its own distribution, so shuffle_z has
median 0 and robust SD 1 *by construction*: every transcript's shuffle curve
sits on the same centre and scale. What can differ is its SHAPE (skew,
kurtosis, tails), reported in the summary, and where the native value falls.
empirical_p = (shuffles <= native + 1) / (n + 1), floor 1/(n+1).

Stdlib only.

Usage
  ./bin/ultrafold_shuffle_summary.py -d human_all
  ./bin/ultrafold_shuffle_summary.py -d human_all --tool-dir rnafold --scores rnafold_results.csv
"""
import argparse
import csv
import logging
import math
import os
import statistics
import sys
from glob import glob
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger('ultrafold_shuffle_summary')

MAD_SCALE = 1.4826


def robust_stats(values):
    """(median, mad) of the shuffled MFEs; mad is the raw median absolute deviation."""
    med = statistics.median(values)
    return med, statistics.median(abs(v - med) for v in values)


def robust_z(x, med, mad):
    """Robust z as in tools/rnafold/worker.sh; None when MAD is 0."""
    return None if mad == 0 else (x - med) / (mad * MAD_SCALE)


def empirical_p(native, values):
    """Lower-tail empirical p with +1 correction (floor 1/(n+1))."""
    return (sum(v <= native for v in values) + 1) / (len(values) + 1)


def skew_kurt(values):
    """Population skewness and excess kurtosis (NaN-safe for tiny/constant input)."""
    n = len(values)
    m = sum(values) / n
    m2 = sum((v - m) ** 2 for v in values) / n
    if m2 == 0:
        return float('nan'), float('nan')
    m3 = sum((v - m) ** 3 for v in values) / n
    m4 = sum((v - m) ** 4 for v in values) / n
    return m3 / m2 ** 1.5, m4 / m2 ** 2 - 3


def fmt(v, nd=4):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return 'NA'
    return f'{v:.{nd}f}'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-d', '--dataset', required=True)
    ap.add_argument('--tool-dir', default='rnafold_shuffles',
                    help='directory under runs/<dataset>/ (default rnafold_shuffles)')
    ap.add_argument('--scores', help='RNAfold results CSV for the stored z-score (region MRNA)')
    ap.add_argument('--runs-root', default=os.environ.get(
        'RUNS_ROOT', str(Path(__file__).resolve().parent.parent / 'runs')))
    args = ap.parse_args()

    run_dir = Path(args.runs_root) / args.dataset
    tool_dir = run_dir / args.tool_dir
    res_files = sorted(glob(str(tool_dir / 'results' / 'results_*.csv')))
    if not res_files:
        sys.exit(f'ERROR: no results in {tool_dir / "results"}')

    stored = {}
    if args.scores:
        with open(args.scores, newline='') as f:
            for r in csv.DictReader(f):
                if r['region'] == 'MRNA':
                    stored[r['transcript_id']] = r['RNAfold_zscore']

    out_dir = run_dir / 'ultrafold'
    out_dir.mkdir(parents=True, exist_ok=True)
    n_done = 0
    with open(out_dir / 'shuffle_distributions.tsv', 'w') as fd, \
            open(out_dir / 'shuffle_summary.tsv', 'w') as fs:
        fd.write('seqname\ttranscript_id\titeration\tshuffled_mfe\tshuffle_z\n')
        fs.write('\t'.join(['seqname', 'transcript_id', 'native_mfe', 'n_shuffles',
                            'shuffle_median', 'shuffle_mad', 'shuffle_mean', 'shuffle_sd',
                            'native_z', 'empirical_p', 'shuffle_skew',
                            'shuffle_excess_kurtosis', 'stored_z']) + '\n')
        for rf in res_files:
            with open(rf, newline='') as f:
                row = next(csv.DictReader(f))
            seqname = row['seq_name']
            tid = seqname.split('_')[1] if seqname.count('_') >= 2 else seqname
            raw = tool_dir / 'raw_shuffles' / f'{seqname}_raw_shuffles.csv'
            if not raw.exists():
                log.warning(f'{seqname}: no raw shuffles file ({raw.name}); skipped')
                continue
            with open(raw, newline='') as f:
                rows = list(csv.DictReader(f))
            vals = [float(r['shuffled_mfe']) for r in rows]
            native = float(row['orig_mfe'])
            med, mad = robust_stats(vals)
            for r, v in zip(rows, vals):
                fd.write(f'{seqname}\t{tid}\t{r["iteration"]}\t{v:g}\t'
                         f'{fmt(robust_z(v, med, mad))}\n')
            sk, ku = skew_kurt(vals)
            sd = statistics.pstdev(vals)
            fs.write('\t'.join([
                seqname, tid, f'{native:g}', str(len(vals)), fmt(med), fmt(mad),
                fmt(sum(vals) / len(vals)), fmt(sd), fmt(robust_z(native, med, mad)),
                fmt(empirical_p(native, vals)), fmt(sk), fmt(ku),
                stored.get(tid, 'NA')]) + '\n')
            n_done += 1
            log.info(f'{seqname}: n={len(vals)} native z={fmt(robust_z(native, med, mad), 2)}')
    log.info(f'{n_done} transcripts -> {out_dir}')


if __name__ == '__main__':
    main()
