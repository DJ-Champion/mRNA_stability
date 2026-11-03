"""Tests for lib.txcoords: genome -> spliced-transcript mapping.

Layout used throughout (genomic, 0-based half-open):
    exon A [100,200)   intron   exon B [300,400)   intron   exon C [500,600)
Transcript length 300. On '+', tx 0 == genomic 100; on '-', tx 0 == genomic 599.
"""
import os
import tempfile

from lib.txcoords import Transcript, load_transcripts, region_of

EXONS = [(100, 200), (300, 400), (500, 600)]


def _tx(strand):
    return Transcript('T', 'chr1', strand, exons=list(EXONS)).finalize()


def test_length():
    assert _tx('+').length == 300


def test_plus_within_exon():
    assert _tx('+').map_interval(110, 130) == [(10, 30, 110, 130)]
    assert _tx('+').map_interval(320, 340) == [(120, 140, 320, 340)]


def test_minus_within_exon():
    # last exon is first in the transcript on '-'; genomic 599 -> tx 0
    assert _tx('-').map_interval(580, 600) == [(0, 20, 580, 600)]
    assert _tx('-').map_interval(100, 110) == [(290, 300, 100, 110)]


def test_intronic_dropped():
    assert _tx('+').map_interval(210, 290) == []
    assert _tx('-').map_interval(210, 290) == []


def test_junction_spanning_splits_into_pieces():
    # 190..310 covers 10nt of A, 10nt of B (intron skipped)
    assert _tx('+').map_interval(190, 310) == [(90, 100, 190, 200), (100, 110, 300, 310)]
    # on '-' the same genomic interval is read B then A (5'->3')
    assert _tx('-').map_interval(190, 310) == [(190, 200, 300, 310), (200, 210, 190, 200)]


def test_exon_edge_clipping():
    assert _tx('+').map_interval(90, 110) == [(0, 10, 100, 110)]
    assert _tx('+').map_interval(590, 700) == [(290, 300, 590, 600)]


def test_region_of():
    assert region_of(0, 50, 100, 100, 300) == '5UTR'
    assert region_of(120, 180, 100, 100, 300) == 'CDS'
    assert region_of(250, 290, 100, 100, 300) == '3UTR'
    assert region_of(90, 130, 100, 100, 300) == 'CDS'      # 10 in 5UTR, 30 in CDS


def test_load_transcripts_gff_to_zero_based():
    gff = ('##gff-version 3\n'
           'chr1\tx\ttranscript\t101\t600\t.\t-\t.\tID=T;transcript_id=T;gene_id=G;gene_name=GN\n'
           'chr1\tx\texon\t101\t200\t.\t-\t.\ttranscript_id=T\n'
           'chr1\tx\texon\t301\t400\t.\t-\t.\ttranscript_id=T\n'
           'chr1\tx\texon\t501\t600\t.\t-\t.\ttranscript_id=T\n'
           'chr1\tx\texon\t1\t50\t.\t+\t.\ttranscript_id=OTHER\n')
    with tempfile.NamedTemporaryFile('w', suffix='.gff', delete=False) as f:
        f.write(gff)
        path = f.name
    try:
        txs = load_transcripts(path, wanted={'T'})
    finally:
        os.unlink(path)
    assert list(txs) == ['T']
    assert txs['T'].exons == EXONS
    assert (txs['T'].strand, txs['T'].gene_name) == ('-', 'GN')
