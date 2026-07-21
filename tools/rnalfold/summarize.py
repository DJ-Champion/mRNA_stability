#!/usr/bin/env python3
"""
tools/rnalfold/summarize.py

Summarise RNALfold output into per-sequence local-structure features, and
derive the shuffle-normalised features from a dinucleotide-shuffle null.

The worker pipes RNALfold's stdout straight into this script, so nothing is
buffered to disk between folding and summarising. Two modes:

  original   read the fold of ONE sequence; print
             "<min_local_mfe>\\t<median_local_mfe>\\t<hotspot_coverage_fraction>"

  aggregate  read the folds of N shuffled sequences (a multi-FASTA stream);
             stream per-shuffle statistics, write the raw per-shuffle table,
             and print the null-normalised summary:
             "<median_shuffle_min_local_mfe>\\t<pvalue_min_local_mfe>\\t
              <zscore_min_local_mfe>\\t<median_shuffle_hotspot_coverage>\\t
              <hotspot_coverage_zscore>\\t<n_valid>\\t<n_missing>"

All three per-sequence features are computed over RNALfold's emitted
locally-optimal *structure* lines only. RNALfold's trailing summary line — the
summed energy of the non-overlapping decomposition, which begins with
whitespace rather than dot-bracket and is always at least as negative as any
single hit — is excluded, matching the worker's anchored extractor.

Feature definitions
-------------------
min_local_mfe
    Most negative structure-line energy: the strongest single local fold.
    Reproduces the worker's previous `min_local_mfe` / `shuffle_min_local_mfe`
    exactly (same structure-line anchoring).

median_local_mfe
    Median energy across the emitted structure lines — the typical strength of
    a reported local fold. NA if no structures were emitted. The count of
    emitted structures is an RNALfold reporting property, so this is a typical
    strength, not a density.

hotspot_coverage_fraction
    Fraction of the region covered by the union of structures whose energy is
    strictly below the hotspot MFE threshold. Each structure occupies
    [start, start + len(structure) - 1] (1-based, inclusive); overlapping
    structures are merged, so clustered RNALfold hits are not double-counted
    (and a non-overlapping decomposition is handled identically). In [0, 1],
    and length-invariant — it does not inherit the region-length confound of
    min_local_mfe.

median_shuffle_* / *_zscore
    The observed statistic compared against the same statistic on each
    dinucleotide-shuffled sequence. The z-score is robust (median and MAD,
    scaled by 1.4826). For coverage a *positive* z means more structured than
    the composition-matched null; for min_local_mfe a *negative* z means the
    strongest fold is more stable than the null (sign follows the energy).
    NA when the null MAD is zero (all shuffles identical at that statistic).

pvalue_min_local_mfe
    Empirical one-sided p that a shuffle's strongest fold is at least as stable
    as the observed one, with a conservative correction for shuffles that
    failed to fold (identical formula to the previous worker).
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple

# RNALfold prints one locally-optimal structure per line, starting at column 0:
#     <dot-bracket>  ( <energy>)  <start>
# e.g. ".((((....)))).  ( -5.30)   127"
# A structure line therefore begins with dot-bracket then the parenthesised
# energy. This anchoring excludes the trailing summary line (" ( -17.80)",
# begins with whitespace) and the echoed sequence line. Records in a
# multi-FASTA stream are separated by the echoed ">header" lines.
_NUM = r'[-+]?[0-9]*\.?[0-9]+'
_STRUCT_ENERGY_RE = re.compile(r'^[().]+\s+\(\s*(' + _NUM + r')\s*\)')
# Full structure line including the 1-based start coordinate, needed to place
# the fold for coverage. Real RNALfold structure lines always carry it.
_STRUCT_FULL_RE = re.compile(r'^([().]+)\s+\(\s*' + _NUM + r'\s*\)\s+([0-9]+)\s*$')


def median(xs: Sequence[float]) -> Optional[float]:
    """Median of a sequence; None if empty."""
    n = len(xs)
    if n == 0:
        return None
    s = sorted(xs)
    mid = n // 2
    if n % 2:
        return s[mid]
    return (s[mid - 1] + s[mid]) / 2.0


def mad(xs: Sequence[float], center: float) -> float:
    """Median absolute deviation about `center` (unscaled). 0.0 if empty."""
    m = median([abs(x - center) for x in xs])
    return 0.0 if m is None else m


def merge_covered(intervals: List[Tuple[int, int]]) -> int:
    """Total nucleotides in the union of inclusive 1-based intervals."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    covered = 0
    cur_a, cur_b = intervals[0]
    for a, b in intervals[1:]:
        if a <= cur_b:                 # overlaps current run (start-sorted)
            if b > cur_b:
                cur_b = b
        else:
            covered += cur_b - cur_a + 1
            cur_a, cur_b = a, b
    covered += cur_b - cur_a + 1
    return covered


class RecordStats:
    """Per-sequence summary of one RNALfold record."""
    __slots__ = ('min_mfe', 'median_mfe', 'coverage', 'n_struct')

    def __init__(self, min_mfe: Optional[float], median_mfe: Optional[float],
                 coverage: Optional[float], n_struct: int) -> None:
        self.min_mfe = min_mfe
        self.median_mfe = median_mfe
        self.coverage = coverage
        self.n_struct = n_struct


def parse_record(lines: Iterable[str], hotspot_mfe: float,
                 seq_length: int) -> RecordStats:
    """Summarise the RNALfold output lines for a single sequence."""
    struct_energies: List[float] = []
    intervals: List[Tuple[int, int]] = []
    for line in lines:
        em = _STRUCT_ENERGY_RE.match(line)
        if not em:
            continue                    # summary line, sequence line, header
        energy = float(em.group(1))
        struct_energies.append(energy)
        if energy < hotspot_mfe:
            fm = _STRUCT_FULL_RE.match(line)
            if fm:                      # need the coordinate to place the fold
                struct, start = fm.group(1), int(fm.group(2))
                a = start
                b = start + len(struct) - 1
                if a < 1:
                    a = 1
                if seq_length > 0 and b > seq_length:
                    b = seq_length
                if b >= a:
                    intervals.append((a, b))
    min_mfe = min(struct_energies) if struct_energies else None
    med_mfe = median(struct_energies)
    coverage = merge_covered(intervals) / seq_length if seq_length > 0 else None
    return RecordStats(min_mfe, med_mfe, coverage, len(struct_energies))


def stream_records(line_iter: Iterable[str]) -> Iterator[List[str]]:
    """Split an RNALfold stream into per-sequence line groups on '>' headers."""
    rec: List[str] = []
    started = False
    for line in line_iter:
        if line.startswith('>'):
            if started:
                yield rec
            rec = []
            started = True
        else:
            rec.append(line)
    if started or rec:
        yield rec


def compute_original(line_iter: Iterable[str], hotspot_mfe: float,
                     seq_length: int) -> RecordStats:
    """Summarise the original (single-sequence) fold."""
    return parse_record(line_iter, hotspot_mfe, seq_length)


def aggregate_stats(shuf_mins: Sequence[float], shuf_covs: Sequence[float],
                    orig_min: float, orig_cov: float, n_shuffles: int) -> dict:
    """Null-normalised statistics for the observed min and coverage."""
    n_valid = len(shuf_mins)
    missing = n_shuffles - n_valid

    med_min = median(shuf_mins)
    mad_min = mad(shuf_mins, med_min) if med_min is not None else 0.0
    more_stable = sum(1 for x in shuf_mins if x <= orig_min)
    # Conservative missingness handling: if the observed min is at least as
    # stable as the null centre, count unfolded shuffles as "more stable".
    if med_min is not None and orig_min <= med_min:
        pvalue = (more_stable + missing + 1) / (n_shuffles + 1)
    else:
        pvalue = (more_stable + 1) / (n_shuffles + 1)
    z_min = None if mad_min == 0 else (orig_min - med_min) / (mad_min * 1.4826)

    med_cov = median(shuf_covs)
    mad_cov = mad(shuf_covs, med_cov) if med_cov is not None else 0.0
    z_cov = None if mad_cov == 0 else (orig_cov - med_cov) / (mad_cov * 1.4826)

    return {
        'median_shuffle_min_local_mfe': med_min,
        'pvalue_min_local_mfe': pvalue,
        'zscore_min_local_mfe': z_min,
        'median_shuffle_hotspot_coverage': med_cov,
        'hotspot_coverage_zscore': z_cov,
        'n_valid': n_valid,
        'n_missing': missing,
    }


def compute_aggregate(line_iter: Iterable[str], hotspot_mfe: float,
                      seq_length: int, orig_min: float, orig_cov: float,
                      n_shuffles: int, raw_writer=None) -> dict:
    """Stream shuffled folds, collect per-shuffle stats, return null summary.

    `raw_writer`, if given, is called as raw_writer(iteration, min, coverage)
    for each valid shuffle so the caller can persist the raw null. Memory is
    bounded by the number of shuffles (two scalars each), not by sequence
    length or by the number of emitted structures.
    """
    shuf_mins: List[float] = []
    shuf_covs: List[float] = []
    iteration = 0
    for rec in stream_records(line_iter):
        rs = parse_record(rec, hotspot_mfe, seq_length)
        if rs.min_mfe is None:          # unusable shuffle (no fold emitted)
            continue
        iteration += 1
        cov = rs.coverage if rs.coverage is not None else 0.0
        shuf_mins.append(rs.min_mfe)
        shuf_covs.append(cov)
        if raw_writer is not None:
            raw_writer(iteration, rs.min_mfe, cov)
    return aggregate_stats(shuf_mins, shuf_covs, orig_min, orig_cov, n_shuffles)


def _fmt(value: Optional[float], spec: str) -> str:
    return 'NA' if value is None else format(value, spec)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)

    po = sub.add_parser('original', help='summarise one sequence fold')
    po.add_argument('--seq-length', type=int, required=True)
    po.add_argument('--hotspot-mfe', type=float, required=True)

    pa = sub.add_parser('aggregate', help='summarise the shuffled null')
    pa.add_argument('--seq-length', type=int, required=True)
    pa.add_argument('--hotspot-mfe', type=float, required=True)
    pa.add_argument('--orig-min', type=float, required=True)
    pa.add_argument('--orig-coverage', type=float, required=True)
    pa.add_argument('--n-shuffles', type=int, required=True)
    pa.add_argument('--min-valid', type=int, required=True)
    pa.add_argument('--seq-name', required=True)
    pa.add_argument('--raw-out', required=True)

    args = parser.parse_args(argv)

    if args.mode == 'original':
        rs = compute_original(sys.stdin, args.hotspot_mfe, args.seq_length)
        if rs.min_mfe is None:
            sys.stderr.write("ERROR: no RNALfold structures found for original sequence\n")
            return 3
        sys.stdout.write("%s\t%s\t%s\n" % (
            format(rs.min_mfe, '.4f'),
            _fmt(rs.median_mfe, '.4f'),
            _fmt(rs.coverage, '.6f'),
        ))
        return 0

    # aggregate
    with open(args.raw_out, 'w', newline='') as raw_f:
        writer = csv.writer(raw_f, lineterminator='\n')
        writer.writerow(['seq_name', 'iteration', 'shuffle_min_local_mfe',
                         'shuffle_hotspot_coverage'])

        def raw_writer(iteration: int, mn: float, cov: float) -> None:
            writer.writerow([args.seq_name, iteration,
                             format(mn, '.4f'), format(cov, '.6f')])

        stats = compute_aggregate(sys.stdin, args.hotspot_mfe, args.seq_length,
                                  args.orig_min, args.orig_coverage,
                                  args.n_shuffles, raw_writer=raw_writer)

    if stats['n_valid'] < args.min_valid:
        sys.stderr.write(
            "ERROR: insufficient valid RNALfold shuffles: %d/%d (min %d)\n"
            % (stats['n_valid'], args.n_shuffles, args.min_valid))
        return 4

    sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\n" % (
        _fmt(stats['median_shuffle_min_local_mfe'], '.4f'),
        _fmt(stats['pvalue_min_local_mfe'], '.4f'),
        _fmt(stats['zscore_min_local_mfe'], '.4f'),
        _fmt(stats['median_shuffle_hotspot_coverage'], '.6f'),
        _fmt(stats['hotspot_coverage_zscore'], '.4f'),
        stats['n_valid'], stats['n_missing'],
    ))
    return 0


if __name__ == '__main__':
    sys.exit(main())
