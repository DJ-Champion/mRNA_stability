"""Synthetic checks for tools/rnaplfold/summarize.py.

Runs under pytest or standalone (``python tests/test_rnaplfold_summarize.py``).
The _lunp fixture is hand-built so every expected value is computable by hand,
and it pins the end-anchored indexing: row i, column l is the probability that
the l nt ENDING at i are unpaired.
"""
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'tools' / 'rnaplfold'))

import summarize as S  # noqa: E402

# 12-nt transcript, AUG at 1-based 5..7 (cds_offset = 4). Column l=k is
# 0.1*i + 0.001*k for row i, NA where l > i (no stretch of that length ends there).
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


def test_aug_and_stretch_means():
    t = _load()
    r = S.summarize(t, 12, 4, 2, 2, [1, 3])
    assert (r['start_win_start'], r['start_win_end'], r['start_win_len']) == (3, 9, 7)
    # AUG = 3 nt ending at 7, l=3 column
    assert _close(r['aug_unpaired'], 0.05 * 7 + 0.003)
    # u=1: ends 3..9
    assert _close(r['start_mean_unpaired_u1'], sum(0.05 * i + 0.001 for i in range(3, 10)) / 7)
    # u=3: stretch must lie inside the window, so ends 5..9
    assert _close(r['start_mean_unpaired_u3'], sum(0.05 * i + 0.003 for i in range(5, 10)) / 5)


def test_window_shorter_than_stretch_is_na():
    t = _load()
    r = S.summarize(t, 12, 4, 0, 0, [8])   # window = AUG only (3 nt)
    assert r['start_mean_unpaired_u8'] is None
    assert S.fmt(None) == 'NA'


def test_five_utr_length_column_and_fallback():
    with tempfile.TemporaryDirectory() as d:
        new = Path(d) / 'new.tsv'
        new.write_text('seqname\tgene_id\ttranscript_id\tregion\tlength\tfive_utr_length\n'
                       'G_T_mRNA\tG\tT\tmRNA\t100\t42\n')
        assert S.five_utr_length(str(new), 'G_T_mRNA') == 42
        old = Path(d) / 'old.tsv'
        old.write_text('seqname\tgene_id\ttranscript_id\tregion\tlength\n'
                       'G_T_mRNA\tG\tT\tmRNA\t100\nG_T_5UTR\tG\tT\t5UTR\t37\n')
        assert S.five_utr_length(str(old), 'G_T_mRNA') == 37
        none = Path(d) / 'none.tsv'
        none.write_text('seqname\tgene_id\ttranscript_id\tregion\tlength\n'
                        'G_T_mRNA\tG\tT\tmRNA\t100\n')
        try:
            S.five_utr_length(str(none), 'G_T_mRNA')
        except ValueError:
            pass
        else:
            raise AssertionError('expected ValueError')


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn()
            print('ok', name)
