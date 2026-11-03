"""Tests for bin/icshape_on_transcripts.py projection and Gini.

Transcript T on '-', exons (0-based half-open) [100,200) [300,400), length 200.
tx 0 == genomic 399; genomic 300 == tx 99; genomic 199 == tx 100.
"""
import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))

from icshape_on_transcripts import population_gini, project_bed  # noqa: E402
from lib.txcoords import Transcript  # noqa: E402


def _bed(rows):
    f = tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False)
    f.write(''.join('\t'.join(map(str, r)) + '\n' for r in rows))
    f.close()
    return f.name


def _project(rows):
    tx = Transcript('T', 'chr1', '-', exons=[(100, 200), (300, 400)]).finalize()
    path = _bed(rows)
    try:
        return project_bed(path, {'T': tx})['T']
    finally:
        os.unlink(path)


def test_gini_anchors():
    assert abs(population_gini([0, 1]) - 0.5) < 1e-12
    assert abs(population_gini([0, 0, 1]) - 2 / 3) < 1e-12
    assert population_gini([1, 1, 1]) == 0.0


def test_run_length_placed_on_minus_strand_in_tx_frame():
    tr = _project([('chr1', 396, 399, 'NA', 0.5, '-')])      # genomic 396..398 -> tx 1..3
    assert [i for i, v in enumerate(tr) if v is not None] == [1, 2, 3]
    assert set(v for v in tr if v is not None) == {0.5}


def test_strand_specific():
    tr = _project([('chr1', 396, 399, 'NA', 0.5, '+')])      # wrong strand -> ignored
    assert all(v is None for v in tr)


def test_junction_spanning_run_skips_intron():
    # genomic 198..302: 2 nt of exon1 (tx 100,101), 2 nt of exon2 (tx 97,98.. see below)
    tr = _project([('chr1', 198, 302, 'NA', 0.9, '-')])
    placed = [i for i, v in enumerate(tr) if v is not None]
    assert placed == [98, 99, 100, 101]                      # 300,301 -> tx 98,99 ; 198,199 -> tx 100,101


def test_intronic_and_other_chrom_ignored():
    tr = _project([('chr1', 220, 280, 'NA', 1.0, '-'), ('chr2', 100, 200, 'NA', 1.0, '-')])
    assert all(v is None for v in tr)
