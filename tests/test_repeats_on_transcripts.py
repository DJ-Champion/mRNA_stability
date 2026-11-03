"""End-to-end test of bin/repeats_on_transcripts.py on a synthetic minus-strand
transcript.

Minus-strand transcript T, exons (genomic, 0-based half-open)
    [100,200) [300,400) [500,600)  -> length 300, tx 0 == genomic 599.
Region split from candidates.tsv: 5UTR 20, CDS 40, 3UTR 240 (3UTR = tx 60..300).

Repeats (BED, genomic strand given):
    R1 Alu  [520,580) '-'  -> tx 20..80, same strand as T  -> sense;    40 nt CDS + 20 nt 3UTR -> CDS by majority
    R2 Alu  [320,380) '+'  -> tx 120..180, opposite        -> antisense
    R3 Alu  [110,150) '+'  -> tx 250..290, opposite        -> antisense
    R4 Alu  [210,290) '+'  -> intronic                     -> dropped
    R5 L1   [330,360) '+'  -> tx 140..170, not Alu         -> kept in repeats, not paired
    R6 Alu  [190,310) '-'  -> spans A/B junction, 10nt+10nt -> two pieces, sense
"""
import csv
import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / 'bin' / 'repeats_on_transcripts.py'

GFF = ('chr1\tx\ttranscript\t101\t600\t.\t-\t.\tID=T;transcript_id=T;gene_id=G;gene_name=GN\n'
       'chr1\tx\texon\t101\t200\t.\t-\t.\ttranscript_id=T\n'
       'chr1\tx\texon\t301\t400\t.\t-\t.\ttranscript_id=T\n'
       'chr1\tx\texon\t501\t600\t.\t-\t.\ttranscript_id=T\n')
BED = [('chr1', 520, 580, 'AluA', 100, '-', 'SINE', 'Alu', 5.0),
       ('chr1', 320, 380, 'AluB', 100, '+', 'SINE', 'Alu', 5.0),
       ('chr1', 110, 150, 'AluC', 100, '+', 'SINE', 'Alu', 5.0),
       ('chr1', 210, 290, 'AluIntron', 100, '+', 'SINE', 'Alu', 5.0),
       ('chr1', 330, 360, 'L1X', 100, '+', 'LINE', 'L1', 5.0),
       ('chr1', 190, 310, 'AluJ', 100, '-', 'SINE', 'Alu', 5.0)]


def _run():
    root = Path(tempfile.mkdtemp())
    (root / 'ds' / 'extracted_regions').mkdir(parents=True)
    (root / 'ds' / 'ultrafold').mkdir()
    (root / '_resources' / 'hg38_rmsk').mkdir(parents=True)
    (root / 'ds' / 'extracted_regions' / 'canonical.gff').write_text(GFF)
    (root / 'ds' / 'ultrafold' / 'candidates.tsv').write_text(
        'transcript_id\tgene_id\tgene_name\tchrom\tstrand\tmrna_len\tutr5_len\tcds_len\tutr3_len\n'
        'T\tG\tGN\tchr1\t-\t300\t20\t40\t240\n')
    with gzip.open(root / '_resources' / 'hg38_rmsk' / 'rmsk.bed.gz', 'wt') as f:
        for r in BED:
            f.write('\t'.join(map(str, r)) + '\n')
    subprocess.run([sys.executable, str(SCRIPT), '-d', 'ds', '--runs-root', str(root)],
                   check=True, capture_output=True)
    rd = lambda n: list(csv.DictReader(open(root / 'ds' / 'ultrafold' / n), delimiter='\t'))
    return rd('repeats.tsv'), rd('alu_pairs.tsv'), rd('alu_summary.tsv')


def test_orientation_coordinates_and_filtering():
    reps, _, _ = _run()
    by = {(r['rep_name'], r['piece']): r for r in reps}
    assert 'AluIntron' not in {r['rep_name'] for r in reps}
    a = by[('AluA', '1')]
    assert (a['tx_start'], a['tx_end'], a['orientation']) == ('20', '80', 'sense')
    assert a['region'] == 'CDS'
    b = by[('AluB', '1')]
    assert (b['tx_start'], b['tx_end'], b['orientation'], b['region']) == ('120', '180', 'antisense', '3UTR')
    assert by[('L1X', '1')]['rep_class'] == 'LINE'


def test_junction_spanning_repeat_gives_two_pieces():
    reps, _, _ = _run()
    j = sorted((r for r in reps if r['rep_name'] == 'AluJ'), key=lambda r: int(r['piece']))
    assert len(j) == 2 and all(r['n_pieces'] == '2' for r in j)
    assert j[0]['piece_len'] == j[1]['piece_len'] == '10'
    assert j[0]['repeat_id'] == j[1]['repeat_id']


def test_pairs_only_opposite_orientation_in_3utr_and_non_alu_excluded():
    _, pairs, summ = _run()
    names = {frozenset((p['name_a'], p['name_b'])) for p in pairs}
    assert frozenset(('AluB', 'AluC')) not in names           # both antisense
    assert frozenset(('AluB', 'AluJ')) in names               # antisense x sense
    assert all(p['orient_a'] != p['orient_b'] for p in pairs)
    assert all('L1X' not in (p['name_a'], p['name_b']) for p in pairs)
    # 3UTR Alus: AluB (antisense), AluC (antisense), AluJ (sense); AluA is CDS
    assert {(p['name_a'], p['name_b']) for p in pairs} == {('AluB', 'AluJ'), ('AluJ', 'AluC')} \
        or {frozenset((p['name_a'], p['name_b'])) for p in pairs} == \
           {frozenset(('AluB', 'AluJ')), frozenset(('AluC', 'AluJ'))}
    s = summ[0]
    assert (s['n_alu'], s['n_sense'], s['n_antisense'], s['n_opposite_pairs'],
            s['n_alu_elsewhere']) == ('3', '1', '2', '2', '1')
