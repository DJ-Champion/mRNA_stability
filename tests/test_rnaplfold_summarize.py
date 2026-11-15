"""Synthetic checks for tools/rnaplfold/summarize.py.

Runs under pytest or standalone (``python tests/test_rnaplfold_summarize.py``).
The _lunp fixture is hand-built so every expected value is computable by hand,
and it pins the end-anchored indexing: row i, column l is the probability that
the l nt ENDING at i are unpaired.
"""
import math
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'tools' / 'rnaplfold'))

import summarize as S  # noqa: E402

# 12-nt transcript, AUG at 1-based 5..7 (cds_offset = 4). Column l=k is
# 0.05*i + 0.001*k for row i, NA where l > i (no stretch of that length ends there).
def _lunp_text(n=12, u=4):
    lines = ['#unpaired probabilities', '#i$\tl=1\t2\t3\t4']
    for i in range(1, n + 1):
        cells = [f'{0.05 * i + 0.001 * l:.3f}' if l <= i else 'NA' for l in range(1, u + 1)]
        lines.append(f'{i}\t' + '\t'.join(cells))
    return '\n'.join(lines) + '\n'


def _close(a, b, tol=1e-9):
    return a is not None and abs(a - b) <= tol


def _load():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / 'x_lunp'
        p.write_text(_lunp_text())
        return S.read_lunp(str(p))


def test_read_lunp_skips_comments_and_keeps_na():
    t = _load()
    assert sorted(t) == list(range(1, 13))
    assert t[1] == [0.051, None, None, None]
    assert _close(t[4][3], 0.204)


def test_start_window_and_clipping():
    # AUG at 5..7, 2 up / 2 down -> 3..9
    assert S.start_window(4, 12, 2, 2) == (3, 9)
    # clipped at both transcript ends
    assert S.start_window(4, 12, 30, 30) == (1, 12)


def _summ(t, up=2, down=2, u=(1, 3), three=2, flanks=(), window=300, offset=4, n=12):
    return S.summarize(t, n, offset, three, up, down, list(u), window, list(flanks))


def test_aug_and_stretch_means():
    t = _load()
    r = _summ(t)
    assert (r['start_win_start'], r['start_win_end'], r['start_win_len']) == (3, 9, 7)
    # AUG = 3 nt ending at 7, l=3 column
    assert _close(r['aug_unpaired'], 0.05 * 7 + 0.003)
    # u=1: ends 3..9
    assert _close(r['start_mean_unpaired_u1'], sum(0.05 * i + 0.001 for i in range(3, 10)) / 7)
    # u=3: stretch must lie inside the window, so ends 5..9
    assert _close(r['start_mean_unpaired_u3'], sum(0.05 * i + 0.003 for i in range(5, 10)) / 5)


def test_distribution_statistics_over_window():
    t = _load()
    r = _summ(t)
    v = [0.05 * i + 0.001 for i in range(3, 10)]          # u1, window 3..9
    assert _close(r['start_min_unpaired_u1'], min(v))
    assert _close(r['start_max_unpaired_u1'], max(v))
    assert _close(r['start_median_unpaired_u1'], 0.05 * 6 + 0.001)
    assert r['start_openfrac_unpaired_u1'] == 0.0           # nothing above 0.5
    assert _close(r['start_meanlog_unpaired_u1'], sum(math.log10(x) for x in v) / 7)


def test_openfrac_counts_values_above_half():
    t = _load()
    r = _summ(t, up=12, down=12)                             # whole 1..12
    # u1 values 0.051 .. 0.601; > 0.5 for i = 10, 11, 12
    assert _close(r['start_openfrac_unpaired_u1'], 3 / 12)


def test_up_and_down_segments_exclude_the_aug():
    t = _load()
    r = _summ(t)
    # window 3..9, AUG 5..7: upstream = 3..4, downstream = 8..9
    up = [math.log10(0.05 * i + 0.001) for i in (3, 4)]
    down = [math.log10(0.05 * i + 0.001) for i in (8, 9)]
    assert _close(r['start_up_meanlog_unpaired_u1'], sum(up) / 2)
    assert _close(r['start_down_meanlog_unpaired_u1'], sum(down) / 2)
    # no 5'UTR: upstream segment is empty -> None; downstream still defined
    r0 = _summ(t, offset=0)
    assert r0['start_up_meanlog_unpaired_u1'] is None
    assert r0['start_down_meanlog_unpaired_u1'] is not None


def test_flank_variants():
    t = _load()
    r = _summ(t, flanks=(1,))
    # flank 1 -> window 4..8
    assert _close(r['w1_mean_unpaired_u1'], sum(0.05 * i + 0.001 for i in range(4, 9)) / 5)
    assert 'w1_meanlog_unpaired_u1' in r


def test_region_means():
    t = _load()
    r = _summ(t, three=2)
    # 5'UTR 1..4, CDS 5..10, 3'UTR 11..12 (u1)
    f = lambda a, b: sum(math.log10(0.05 * i + 0.001) for i in range(a, b + 1)) / (b - a + 1)
    assert _close(r['utr5_meanlog_unpaired_u1'], f(1, 4))
    assert _close(r['cds_meanlog_unpaired_u1'], f(5, 10))
    assert _close(r['utr3_meanlog_unpaired_u1'], f(11, 12))
    assert _close(r['cap50_meanlog_unpaired_u1'], f(1, 12))   # capped at transcript end
    # unknown 3'UTR length -> CDS and 3'UTR NA, 5'UTR and cap kept
    r2 = _summ(t, three=None)
    assert r2['cds_meanlog_unpaired_u1'] is None and r2['utr3_meanlog_unpaired_u1'] is None
    assert r2['utr5_meanlog_unpaired_u1'] is not None
    # no 5'UTR -> utr5 NA
    assert _summ(t, offset=0)['utr5_meanlog_unpaired_u1'] is None


def test_log_floor_for_zero_probability():
    assert _close(S.meanlog([0.0, 1.0]), math.log10(S.LOG_FLOOR) / 2)


def test_seq_shorter_than_window_flag():
    t = _load()
    assert _summ(t, u=(1,), window=300)['seq_shorter_than_window'] == 1
    assert _summ(t, u=(1,), window=12)['seq_shorter_than_window'] == 0
    assert _summ(t, u=(1,), window=12)['seq_length'] == 12


def test_window_shorter_than_stretch_is_na():
    t = _load()
    r = _summ(t, up=0, down=0, u=(8,))   # window = AUG only (3 nt)
    assert r['start_mean_unpaired_u8'] is None
    assert r['start_min_unpaired_u8'] is None
    assert S.fmt(None) == 'NA'


def test_utr_lengths_columns_and_fallback():
    with tempfile.TemporaryDirectory() as d:
        new = Path(d) / 'new.tsv'
        new.write_text('seqname\tgene_id\ttranscript_id\tregion\tlength\tfive_utr_length\tthree_utr_length\n'
                       'G_T_mRNA\tG\tT\tmRNA\t100\t42\t17\n')
        assert S.utr_lengths(str(new), 'G_T_mRNA') == (42, 17)
        old = Path(d) / 'old.tsv'
        old.write_text('seqname\tgene_id\ttranscript_id\tregion\tlength\n'
                       'G_T_mRNA\tG\tT\tmRNA\t100\nG_T_5UTR\tG\tT\t5UTR\t37\n'
                       'G_T_3UTR\tG\tT\t3UTR\t11\n')
        assert S.utr_lengths(str(old), 'G_T_mRNA') == (37, 11)
        no3 = Path(d) / 'no3.tsv'
        no3.write_text('seqname\tgene_id\ttranscript_id\tregion\tlength\n'
                       'G_T_mRNA\tG\tT\tmRNA\t100\nG_T_5UTR\tG\tT\t5UTR\t37\n')
        assert S.utr_lengths(str(no3), 'G_T_mRNA') == (37, None)
        none = Path(d) / 'none.tsv'
        none.write_text('seqname\tgene_id\ttranscript_id\tregion\tlength\n'
                        'G_T_mRNA\tG\tT\tmRNA\t100\n')
        try:
            S.utr_lengths(str(none), 'G_T_mRNA')
        except ValueError:
            pass
        else:
            raise AssertionError('expected ValueError')


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn()
            print('ok', name)
