"""Synthetic checks for tools/rnalfold/summarize.py.

Runs under pytest (``python -m pytest tests/``) or standalone
(``python tests/test_rnalfold_summarize.py``), mirroring
``run_synthetic_check.py`` for environments without pytest.

The fixtures are hand-built RNALfold-format text with energies, coordinates,
and structure lengths chosen so every expected value is computable by hand.
They deliberately include the trailing summary line (" ( -17.80)") and a
sequence line so the test fails if either is ever counted as a local fold —
i.e. it guards the structure-line anchoring of the worker's parsing fix.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'tools' / 'rnalfold'))

import summarize as S  # noqa: E402

HOTSPOT = -5.0
SEQ_LEN = 20

# One sequence. Structures at [1,10] (-6.50) and [14,20] (-8.10) are hotspots;
# [5,10] (-3.20) is not. The " ( -17.80)" summary line and the sequence line
# must both be ignored.
ORIGINAL = """\
>orig
(((....))) ( -6.50)    1
((..)) ( -3.20)    5
((...)) ( -8.10)   14
AUGCAUGCAUGCAUGCAUGC
 ( -17.80)
"""

# Five shuffled records. shuf5 emits no structure line and must be skipped.
AGGREGATE = """\
>shuf1
(((....))) ( -7.00)    1
GCGCGCGCGCGCGCGCGCGC
 ( -12.00)
>shuf2
((..)) ( -4.00)    3
GCGCGCGCGCGCGCGCGCGC
>shuf3
((.)) ( -6.00)    1
((.)) ( -6.00)    3
GCGCGCGCGCGCGCGCGCGC
>shuf4
((................)) ( -5.50)    1
GCGCGCGCGCGCGCGCGCGC
>shuf5
NNNNNNNNNNNNNNNNNNNN
"""


def _close(a, b, tol=1e-4):
    return a is not None and b is not None and abs(a - b) <= tol


def test_merge_covered():
    assert S.merge_covered([]) == 0
    assert S.merge_covered([(1, 10)]) == 10
    assert S.merge_covered([(1, 10), (14, 20)]) == 17          # disjoint
    assert S.merge_covered([(1, 5), (3, 7)]) == 7              # overlap merged
    assert S.merge_covered([(3, 7), (1, 5)]) == 7              # order-independent
    assert S.merge_covered([(1, 5), (5, 8)]) == 8              # touching at 5


def test_median_and_mad():
    assert S.median([]) is None
    assert S.median([3.0]) == 3.0
    assert S.median([1.0, 3.0]) == 2.0
    assert S.median([1.0, 2.0, 3.0]) == 2.0
    assert S.mad([1.0, 1.0, 1.0], 1.0) == 0.0
    assert S.mad([-7.0, -6.0, -5.5, -4.0], -5.75) == 0.75


def test_original_features():
    rs = S.compute_original(ORIGINAL.splitlines(), HOTSPOT, SEQ_LEN)
    assert rs.n_struct == 3                       # summary + seq line excluded
    assert _close(rs.min_mfe, -8.10)              # not -17.80 (summary excluded)
    assert _close(rs.median_mfe, -6.50)
    assert _close(rs.coverage, 17 / 20)           # 0.85


def test_original_no_structures():
    rs = S.compute_original([">x", "AUGCAUGC", " ( -3.20)"], HOTSPOT, SEQ_LEN)
    assert rs.min_mfe is None                      # summary line is not a fold
    assert rs.median_mfe is None
    assert rs.coverage == 0.0


def test_aggregate_null():
    raw = []
    stats = S.compute_aggregate(
        AGGREGATE.splitlines(), HOTSPOT, SEQ_LEN,
        orig_min=-8.10, orig_cov=0.85, n_shuffles=5,
        raw_writer=lambda it, mn, cov: raw.append((it, mn, cov)),
    )

    assert stats['n_valid'] == 4                   # shuf5 skipped
    assert stats['n_missing'] == 1                 # 5 requested - 4 folded
    assert _close(stats['median_shuffle_min_local_mfe'], -5.75)
    assert _close(stats['pvalue_min_local_mfe'], 2 / 6)
    assert _close(stats['zscore_min_local_mfe'], -2.35 / (0.75 * 1.4826))
    assert _close(stats['median_shuffle_hotspot_coverage'], 0.425)
    assert _close(stats['hotspot_coverage_zscore'], 0.425 / (0.25 * 1.4826))

    # Raw null: one row per valid shuffle, in stream order.
    assert [it for it, _, _ in raw] == [1, 2, 3, 4]
    assert [round(mn, 4) for _, mn, _ in raw] == [-7.0, -4.0, -6.0, -5.5]
    covs = [cov for _, _, cov in raw]
    assert _close(covs[0], 0.5) and covs[1] == 0.0
    assert _close(covs[2], 7 / 20) and _close(covs[3], 1.0)


def test_aggregate_zero_mad_is_na():
    # All shuffles identical at both statistics -> MAD 0 -> z-scores NA.
    stats = S.aggregate_stats([-6.0, -6.0, -6.0], [0.5, 0.5, 0.5],
                              orig_min=-8.0, orig_cov=0.9, n_shuffles=3)
    assert stats['zscore_min_local_mfe'] is None
    assert stats['hotspot_coverage_zscore'] is None
    assert _close(stats['pvalue_min_local_mfe'], 1 / 4)


def _main():
    tests = [v for k, v in sorted(globals().items())
             if k.startswith('test_') and callable(v)]
    failed = 0
    for t in tests:
        try:
            t()
            print("PASS %s" % t.__name__)
        except AssertionError as e:
            failed += 1
            print("FAIL %s: %s" % (t.__name__, e))
    print("\n%d/%d passed" % (len(tests) - failed, len(tests)))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(_main())
