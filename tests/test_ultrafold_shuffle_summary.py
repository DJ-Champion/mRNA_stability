"""Tests for bin/ultrafold_shuffle_summary.py: z-score maths and end-to-end tables."""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))

from ultrafold_shuffle_summary import (MAD_SCALE, empirical_p, robust_stats,  # noqa: E402
                                       robust_z, skew_kurt)


def test_robust_z_matches_worker_formula():
    # worker.sh: (orig - median) / (MAD * 1.4826); median/MAD of the shuffles
    vals = [-100.0, -102.0, -98.0, -101.0, -99.0]       # median -100, deviations 0,2,2,1,1 -> MAD 1
    med, mad = robust_stats(vals)
    assert (med, mad) == (-100.0, 1.0)
    assert abs(robust_z(-250.0, med, mad) - (-150.0 / MAD_SCALE)) < 1e-12


def test_even_count_median_and_zero_mad():
    assert robust_stats([1.0, 2.0, 3.0, 4.0]) == (2.5, 1.0)
    assert robust_z(5.0, 3.0, 0.0) is None              # worker prints NA when MAD == 0


def test_empirical_p_has_floor():
    vals = [-10.0] * 99 + [-20.0]
    assert empirical_p(-50.0, vals) == 1 / 101           # beyond all shuffles
    assert abs(empirical_p(-15.0, vals) - 2 / 101) < 1e-12


def test_shuffle_z_has_zero_median_and_unit_robust_sd():
    # by construction: z of the shuffles against themselves
    vals = [-100 + (i % 7) * 1.5 + (i % 3) for i in range(101)]
    med, mad = robust_stats(vals)
    zs = [robust_z(v, med, mad) for v in vals]
    zmed, zmad = robust_stats(zs)
    assert abs(zmed) < 1e-9 and abs(zmad * MAD_SCALE - 1) < 1e-9


def test_skew_kurtosis_symmetric_and_skewed():
    sk, _ = skew_kurt([-2.0, -1.0, 0.0, 1.0, 2.0])
    assert abs(sk) < 1e-12
    sk2, _ = skew_kurt([0, 0, 0, 0, 10.0])
    assert sk2 > 1


def test_end_to_end_tables():
    root = Path(tempfile.mkdtemp())
    tool = root / 'ds' / 'rnafold_shuffles'
    (tool / 'results').mkdir(parents=True)
    (tool / 'raw_shuffles').mkdir()
    seq = 'ENSG1_ENST1.1_mRNA'
    (tool / 'results' / f'results_{seq}.csv').write_text(
        'seq_name,orig_mfe,median_mfe,pvalue,zscore,n_valid,n_missing,n_shuffles_requested\n'
        f'{seq},-250.0000,-100.0000,0.0099,-100.0000,5,0,5\n')
    vals = [-100.0, -102.0, -98.0, -101.0, -99.0]
    (tool / 'raw_shuffles' / f'{seq}_raw_shuffles.csv').write_text(
        'transcript_id,iteration,shuffled_mfe\n' +
        ''.join(f'{seq},{i + 1},{v}\n' for i, v in enumerate(vals)))
    subprocess.run([sys.executable, str(REPO / 'bin' / 'ultrafold_shuffle_summary.py'),
                    '-d', 'ds', '--runs-root', str(root)], check=True, capture_output=True)
    summ = list(csv.DictReader(open(root / 'ds' / 'ultrafold' / 'shuffle_summary.tsv'), delimiter='\t'))
    assert len(summ) == 1 and summ[0]['transcript_id'] == 'ENST1.1'
    assert abs(float(summ[0]['native_z']) - (-150 / MAD_SCALE)) < 1e-3
    dist = list(csv.DictReader(open(root / 'ds' / 'ultrafold' / 'shuffle_distributions.tsv'), delimiter='\t'))
    assert len(dist) == 5 and float(dist[0]['shuffle_z']) == 0.0
