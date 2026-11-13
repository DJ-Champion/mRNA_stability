#!/usr/bin/env python3
"""
tools/rnaplfold/summarize.py
Summarise RNAplfold accessibility over the start region of one mRNA.

Reads the RNAplfold `_lunp` file from a fold of the full mRNA and prints a
two-line CSV (header + row). Columns:

  start_win_start / start_win_end   1-based inclusive window actually used
                                    (clipped to the transcript)
  start_win_len
  cds_offset                        5'UTR length (nt before the AUG)
  aug_unpaired                      P(the 3 nt of the AUG are all unpaired)
  start_mean_unpaired_u<l>          for each requested l: mean over the window
                                    of P(l consecutive nt, wholly inside the
                                    window, all unpaired); NA if the window
                                    is shorter than l
  rnaplfold_window / rnaplfold_span / rnaplfold_u_max   reproducibility params

_lunp convention: row i, column l is the probability that the l nt ENDING at
position i (1-based) are all unpaired. A stretch of l nt inside window
[S, E] therefore ends at i in [S + l - 1, E]. NA entries are skipped.
"""
from __future__ import annotations

import argparse
import csv
import sys

NA = 'NA'


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


def five_utr_length(manifest_path: str, seqname: str) -> int:
    """5'UTR length for seqname's transcript.

    Uses the manifest's five_utr_length column; manifests from before that
    column existed fall back to the length of the matching 5UTR row.
    """
    with open(manifest_path, newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    me = next((r for r in rows if r['seqname'] == seqname), None)
    if me is None:
        raise ValueError(f"{seqname} not in manifest {manifest_path}")
    if me.get('five_utr_length', '') != '':
        return int(me['five_utr_length'])
    utr = next((r for r in rows if r['gene_id'] == me['gene_id']
                and r['transcript_id'] == me['transcript_id']
                and r['region'] == '5UTR'), None)
    if utr is None:
        raise ValueError(
            f"no five_utr_length column and no 5UTR row for {seqname}; "
            "re-run 01_extract.py or add 5UTR to regions_to_extract")
    return int(utr['length'])


def start_window(cds_offset: int, seq_len: int, upstream: int, downstream: int) -> tuple[int, int]:
    """1-based inclusive (S, E): AUG-upstream .. AUG+3+downstream, clipped."""
    start = max(1, cds_offset + 1 - upstream)
    end = min(seq_len, cds_offset + 3 + downstream)
    return start, end


def _mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def stretch_probs(lunp, l: int, first_end: int, last_end: int) -> list[float]:
    """P(l-stretch unpaired) for stretches ending at first_end..last_end."""
    out = []
    for i in range(first_end, last_end + 1):
        row = lunp.get(i)
        if row is not None and len(row) >= l and row[l - 1] is not None:
            out.append(row[l - 1])
    return out


def summarize(lunp, seq_len: int, cds_offset: int, upstream: int, downstream: int,
              u_values: list[int]) -> dict[str, object]:
    s, e = start_window(cds_offset, seq_len, upstream, downstream)
    out: dict[str, object] = {
        'start_win_start': s, 'start_win_end': e, 'start_win_len': e - s + 1,
        'cds_offset': cds_offset,
    }
    aug = lunp.get(cds_offset + 3)
    out['aug_unpaired'] = aug[2] if aug is not None and len(aug) >= 3 else None
    for l in u_values:
        out[f'start_mean_unpaired_u{l}'] = _mean(stretch_probs(lunp, l, s + l - 1, e))
    return out


def fmt(x) -> str:
    if x is None:
        return NA
    return f'{x:.6g}' if isinstance(x, float) else str(x)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--lunp', required=True)
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--seq-name', required=True)
    ap.add_argument('--seq-length', type=int, required=True)
    ap.add_argument('--upstream', type=int, required=True)
    ap.add_argument('--downstream', type=int, required=True)
    ap.add_argument('--u-values', required=True, help='space/comma separated')
    ap.add_argument('--window', type=int, required=True)
    ap.add_argument('--span', type=int, required=True)
    args = ap.parse_args()

    u_values = sorted({int(x) for x in args.u_values.replace(',', ' ').split()})
    if not u_values or u_values[0] < 1:
        sys.exit('ERROR: u_values must be positive integers')
    u_max = max(u_values + [3])

    lunp = read_lunp(args.lunp)
    if not lunp:
        sys.exit(f'ERROR: no rows in {args.lunp}')
    row = summarize(lunp, args.seq_length, five_utr_length(args.manifest, args.seq_name),
                    args.upstream, args.downstream, u_values)
    row.update(rnaplfold_window=args.window, rnaplfold_span=args.span, rnaplfold_u_max=u_max)

    cols = ['seq_name'] + list(row)
    print(','.join(cols))
    print(','.join([args.seq_name] + [fmt(row[c]) for c in row]))


if __name__ == '__main__':
    main()
