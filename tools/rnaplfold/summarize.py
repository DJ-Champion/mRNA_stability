#!/usr/bin/env python3
"""
tools/rnaplfold/summarize.py
Summarise RNAplfold accessibility for one mRNA from its `_lunp` file (a fold of
the full mRNA). Prints a two-line CSV (header + row).

_lunp convention: row i, column l is the probability that the l nt ENDING at
position i (1-based) are all unpaired. A stretch of l nt wholly inside
[S, E] therefore ends at i in [S + l - 1, E]. NA entries are skipped. A
"stretch value" below is one such probability; for l = 1 it is the per-nt
unpaired probability.

Columns (u<l> = stretch length l, for each requested l):

  Window geometry
    start_win_start, start_win_end, start_win_len   primary window, 1-based
                                    inclusive, clipped to the transcript
    seq_length, seq_shorter_than_window   (1 if seq_length < RNAplfold window:
                                    whole molecule in one window, end artefacts)
    cds_offset                      5'UTR length (nt before the AUG)
    aug_unpaired                    P(the 3 nt of the AUG are all unpaired)

  Primary start window (AUG - upstream .. AUG + 3 + downstream)
    start_mean_unpaired_u<l>        mean of stretch values
    start_meanlog_unpaired_u<l>     mean of log10(stretch value); the mean
                                    opening energy up to a constant
    start_min_/start_max_/start_median_unpaired_u<l>
    start_openfrac_unpaired_u<l>    fraction of stretch values > 0.5
    start_up_meanlog_unpaired_u<l>  mean log10 over the part of the window
                                    before the AUG (NA if none)
    start_down_meanlog_unpaired_u<l>  ... after the AUG

  Window-size variants (symmetric flank f nt each side of the AUG)
    w<f>_mean_unpaired_u<l>, w<f>_meanlog_unpaired_u<l>

  Region means of log10 stretch value (whole-transcript fold)
    utr5_/cds_/utr3_meanlog_unpaired_u<l>   NA if the region is empty or its
                                    bounds are unknown (no three_utr_length)
    cap50_meanlog_unpaired_u<l>     transcript positions 1..50

  rnaplfold_window, rnaplfold_span, rnaplfold_u_max   parameters used

Stretch values below LOG_FLOOR are raised to it before taking logs (RNAplfold
can print 0 for vanishing probabilities).
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys

NA = 'NA'
LOG_FLOOR = 1e-12
CAP_LEN = 50


def read_lunp(path: str) -> dict[int, list[float | None]]:
    """position (1-based) -> [P(l=1), P(l=2), ...] with None for NA."""
    table: dict[int, list[float | None]] = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            table[int(parts[0])] = [None if x == 'NA' else float(x) for x in parts[1:]]
    return table


def utr_lengths(manifest_path: str, seqname: str) -> tuple[int, int | None]:
    """(5'UTR length, 3'UTR length or None) for seqname's transcript.

    Uses the manifest's five_utr_length / three_utr_length columns; manifests
    from before those columns fall back to the matching 5UTR / 3UTR row's
    length. The 5'UTR length is required (ValueError); a missing 3'UTR length
    gives None, which turns the CDS and 3'UTR region columns into NA.
    """
    with open(manifest_path, newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    me = next((r for r in rows if r['seqname'] == seqname), None)
    if me is None:
        raise ValueError(f"{seqname} not in manifest {manifest_path}")

    def lookup(col: str, region: str) -> int | None:
        if me.get(col, '') != '':
            return int(me[col])
        row = next((r for r in rows if r['gene_id'] == me['gene_id']
                    and r['transcript_id'] == me['transcript_id']
                    and r['region'] == region), None)
        return int(row['length']) if row is not None else None

    five = lookup('five_utr_length', '5UTR')
    if five is None:
        raise ValueError(
            f"no five_utr_length column and no 5UTR row for {seqname}; "
            "re-run 01_extract.py or add 5UTR to regions_to_extract")
    return five, lookup('three_utr_length', '3UTR')


def start_window(cds_offset: int, seq_len: int, upstream: int, downstream: int) -> tuple[int, int]:
    """1-based inclusive (S, E): AUG-upstream .. AUG+3+downstream, clipped."""
    start = max(1, cds_offset + 1 - upstream)
    end = min(seq_len, cds_offset + 3 + downstream)
    return start, end


def stretch_probs(lunp, l: int, first: int, last: int) -> list[float]:
    """Stretch values for l-stretches wholly inside [first, last] (1-based)."""
    out = []
    for i in range(first + l - 1, last + 1):
        row = lunp.get(i)
        if row is not None and len(row) >= l and row[l - 1] is not None:
            out.append(row[l - 1])
    return out


def mean(vals: list[float]) -> float | None:
    return sum(vals) / len(vals) if vals else None


def meanlog(vals: list[float]) -> float | None:
    return mean([math.log10(max(v, LOG_FLOOR)) for v in vals])


def summarize(lunp, seq_len: int, cds_offset: int, three_utr_len: int | None,
              upstream: int, downstream: int, u_values: list[int], window: int,
              flanks: list[int]) -> dict[str, object]:
    s, e = start_window(cds_offset, seq_len, upstream, downstream)
    out: dict[str, object] = {
        'start_win_start': s, 'start_win_end': e, 'start_win_len': e - s + 1,
        'seq_length': seq_len, 'seq_shorter_than_window': int(seq_len < window),
        'cds_offset': cds_offset,
    }
    aug = lunp.get(cds_offset + 3)
    out['aug_unpaired'] = aug[2] if aug is not None and len(aug) >= 3 else None

    # Part of the primary window before / after the AUG (AUG itself excluded).
    up_end = min(e, cds_offset)
    down_start = max(s, cds_offset + 4)

    for l in u_values:
        v = stretch_probs(lunp, l, s, e)
        out[f'start_mean_unpaired_u{l}'] = mean(v)
        out[f'start_meanlog_unpaired_u{l}'] = meanlog(v)
        out[f'start_min_unpaired_u{l}'] = min(v) if v else None
        out[f'start_max_unpaired_u{l}'] = max(v) if v else None
        out[f'start_median_unpaired_u{l}'] = statistics.median(v) if v else None
        out[f'start_openfrac_unpaired_u{l}'] = (sum(x > 0.5 for x in v) / len(v)) if v else None
        out[f'start_up_meanlog_unpaired_u{l}'] = meanlog(stretch_probs(lunp, l, s, up_end))
        out[f'start_down_meanlog_unpaired_u{l}'] = meanlog(stretch_probs(lunp, l, down_start, e))

    for f in flanks:
        fs, fe = start_window(cds_offset, seq_len, f, f)
        for l in u_values:
            v = stretch_probs(lunp, l, fs, fe)
            out[f'w{f}_mean_unpaired_u{l}'] = mean(v)
            out[f'w{f}_meanlog_unpaired_u{l}'] = meanlog(v)

    regions = {'utr5': (1, cds_offset), 'cap50': (1, min(CAP_LEN, seq_len))}
    if three_utr_len is not None:
        regions['cds'] = (cds_offset + 1, seq_len - three_utr_len)
        regions['utr3'] = (seq_len - three_utr_len + 1, seq_len)
    for name in ('utr5', 'cds', 'utr3', 'cap50'):
        first, last = regions.get(name, (1, 0))      # empty -> NA
        for l in u_values:
            out[f'{name}_meanlog_unpaired_u{l}'] = meanlog(stretch_probs(lunp, l, first, last))
    return out


def fmt(x) -> str:
    if x is None:
        return NA
    return f'{x:.6g}' if isinstance(x, float) else str(x)


def _ints(text: str) -> list[int]:
    return sorted({int(x) for x in text.replace(',', ' ').split()})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--lunp', required=True)
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--seq-name', required=True)
    ap.add_argument('--seq-length', type=int, required=True)
    ap.add_argument('--upstream', type=int, required=True)
    ap.add_argument('--downstream', type=int, required=True)
    ap.add_argument('--u-values', required=True, help='space/comma separated')
    ap.add_argument('--flanks', default='', help='extra symmetric flank sizes, space/comma separated')
    ap.add_argument('--window', type=int, required=True)
    ap.add_argument('--span', type=int, required=True)
    args = ap.parse_args()

    u_values = _ints(args.u_values)
    if not u_values or u_values[0] < 1:
        sys.exit('ERROR: u_values must be positive integers')
    flanks = _ints(args.flanks)
    u_max = max(u_values + [3])

    lunp = read_lunp(args.lunp)
    if not lunp:
        sys.exit(f'ERROR: no rows in {args.lunp}')
    five, three = utr_lengths(args.manifest, args.seq_name)
    row = summarize(lunp, args.seq_length, five, three, args.upstream, args.downstream,
                    u_values, args.window, flanks)
    row.update(rnaplfold_window=args.window, rnaplfold_span=args.span, rnaplfold_u_max=u_max)

    print(','.join(['seq_name'] + list(row)))
    print(','.join([args.seq_name] + [fmt(row[c]) for c in row]))


if __name__ == '__main__':
    main()
