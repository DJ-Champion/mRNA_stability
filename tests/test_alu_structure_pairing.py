"""bin/alu_structure_pairing.py: dot-bracket parsing, Alu projection, cutoffs,
MFE validation, and a real RNAfold check that complementary Alus pair."""
import csv
import glob
import importlib.util
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / 'bin' / 'alu_structure_pairing.py'
spec = importlib.util.spec_from_file_location('alu_structure_pairing', SCRIPT)
asp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(asp)


def rep(rid, orient, pieces, family='Alu'):
    return [{'repeat_id': rid, 'rep_family': family, 'orientation': orient,
             'tx_start': str(s), 'tx_end': str(e)} for s, e in pieces]


# ---- parser

def test_parse_nested():
    assert asp.parse_dotbracket('((..))') == [5, 4, -1, -1, 1, 0]


def test_parse_unbalanced_raises():
    for s in ('(()', '())', ')('):
        with pytest.raises(ValueError):
            asp.parse_dotbracket(s)


def test_parse_other_bracket_raises():
    with pytest.raises(ValueError):
        asp.parse_dotbracket('[..]')


# ---- projection on a hand-built structure
# 0..9 Alu A (sense), 10..19 linker, 20..29 Alu B (antisense), 30..39 Alu C (sense)
# Structure: A[0:5] pairs B[29:24:-1]  (5 bp opposite); A[5:7] pairs C (2 bp same);
# B[20] pairs linker 15 (nonalu); rest unpaired.
def build():
    n = 40
    partner = [-1] * n
    def pair(i, j):
        partner[i], partner[j] = j, i
    for k in range(5):
        pair(k, 29 - k)       # A-B, nested
    pair(5, 39)               # A-C; the table is used as is (project() does not need nesting)
    pair(6, 38)
    pair(15, 20)
    return partner


def test_pair_counts_and_partner_types():
    partner = build()
    # nesting check is not needed by project(); it only reads the table
    alus = asp.build_alus(rep('A', 'sense', [(0, 10)]) + rep('B', 'antisense', [(20, 30)])
                          + rep('C', 'sense', [(30, 40)]) + rep('L1', 'sense', [(10, 20)], family='LINE'))
    assert [a['id'] for a in alus] == ['A', 'B', 'C']
    owner, ov = asp.base_to_alu(alus, 40)
    assert ov == 0
    pair_bp, bases = asp.project(partner, alus, owner)
    assert pair_bp == {(0, 1): 5, (0, 2): 2}
    a = bases[0]
    assert (a['opposite'], a['same'], a['self'], a['nonalu'], a['unpaired']) == (5, 2, 0, 0, 3)
    b = bases[1]
    assert (b['opposite'], b['nonalu'], b['unpaired']) == (5, 1, 4)
    assert sum(bases[0].values()) == alus[0]['len']


def test_self_pair_is_not_counted_between_alus():
    alus = asp.build_alus(rep('A', 'sense', [(0, 10)]))
    owner, _ = asp.base_to_alu(alus, 10)
    partner = asp.parse_dotbracket('((....))..')
    pair_bp, bases = asp.project(partner, alus, owner)
    assert pair_bp == {} and bases[0]['self'] == 4 and bases[0]['unpaired'] == 6


def test_split_exon_repeat_is_one_alu():
    alus = asp.build_alus(rep('S', 'sense', [(0, 10), (50, 60)]))
    assert len(alus) == 1 and alus[0]['len'] == 20 and (alus[0]['start'], alus[0]['end']) == (0, 60)
    owner, _ = asp.base_to_alu(alus, 70)
    assert owner[5] == 0 and owner[30] == -1 and owner[55] == 0


# ---- cutoffs: >= not >
def test_cutoff_boundaries():
    assert asp.called_paired(5, 100, 'nbp', 5)
    assert not asp.called_paired(4, 100, 'nbp', 5)
    assert asp.called_paired(10, 100, 'frac', 10)       # exactly 10%
    assert not asp.called_paired(9, 100, 'frac', 10)
    assert asp.called_paired(1, 300, 'frac', 0.3333) is True


# ---- stems
def test_collapse_stems():
    partner = asp.parse_dotbracket('(((.((...)).)))')
    stems = asp.collapse_stems(partner)
    # outer 3 bp (0,14),(1,13),(2,12); (4,10),(5,9) is a 2 bp stem
    assert stems == [(0, 3, 12, 15), (4, 6, 9, 11)]
    assert sum(e - s for s, e, _, _ in stems) == sum(p > i for i, p in enumerate(partner))


# ---- real RNAfold: a sequence-complementary Alu pair really pairs
def find_rnafold():
    cands = [os.environ.get('RNAFOLD_BIN'), shutil.which('RNAfold')]
    cands += glob.glob(os.path.expanduser('~/miniconda3/pkgs/viennarna-*/bin/RNAfold'))
    for c in cands:
        if c and os.access(c, os.X_OK):
            return c
    return None


def revcomp(s):
    return s[::-1].translate(str.maketrans('ACGU', 'UGCA'))


@pytest.mark.skipif(find_rnafold() is None, reason='RNAfold not available')
def test_complementary_alus_pair_in_real_structure():
    rnd = random.Random(1)
    alu = ''.join(rnd.choice('ACGU') for _ in range(80))
    spacer = ''.join(rnd.choice('AU') for _ in range(30))
    seq = alu + spacer + revcomp(alu)
    out = subprocess.run([find_rnafold(), '--noPS'], input=f'>t\n{seq}\n', text=True,
                         capture_output=True, check=True).stdout.splitlines()
    partner = asp.parse_dotbracket(out[2].split()[0])
    alus = asp.build_alus(rep('A', 'sense', [(0, 80)]) + rep('B', 'antisense', [(110, 190)]))
    owner, _ = asp.base_to_alu(alus, len(seq))
    pair_bp, bases = asp.project(partner, alus, owner)
    assert pair_bp.get((0, 1), 0) >= 60
    assert asp.called_paired(bases[0]['opposite'], 80, 'frac', 50)


# ---- end to end, including the MFE validation failure path
def write_inputs(root, stored_mfe, fold_mfe='-12.30', seq_len=20):
    ds = root / 'ds'
    (ds / 'ultrafold').mkdir(parents=True)
    (ds / 'rnafold_structure' / 'results').mkdir(parents=True)
    (ds / 'ultrafold' / 'shuffle_summary.tsv').write_text(
        f'seqname\ttranscript_id\tnative_mfe\nG_T1_mRNA\tT1\t{stored_mfe}\n')
    seq = 'A' * seq_len
    struct = '(((' + '.' * (seq_len - 6) + ')))'
    (ds / 'rnafold_structure' / 'results' / 'G_T1_mRNA.fold').write_text(
        f'>G_T1_mRNA\n{seq}\n{struct} ({fold_mfe})\n')
    (ds / 'ultrafold' / 'candidates.tsv').write_text(
        'transcript_id\tgene_name\tmrna_len\tutr5_len\tcds_len\tutr3_len\n'
        f'T1\tG\t{seq_len}\t2\t6\t{seq_len - 8}\n')
    (ds / 'ultrafold' / 'repeats.tsv').write_text(
        'transcript_id\trepeat_id\trep_family\torientation\ttx_start\ttx_end\n'
        'T1\tA1\tAlu\tsense\t0\t3\n'
        'T1\tA2\tAlu\tantisense\t17\t20\n')


def run_script(root):
    return subprocess.run(
        [sys.executable, str(SCRIPT), '-d', 'ds', '--repeats-dataset', 'ds', '--runs-root', str(root)],
        capture_output=True, text=True)


def test_end_to_end(tmp_path):
    write_inputs(tmp_path, stored_mfe='-12.3000')
    r = run_script(tmp_path)
    assert r.returncode == 0, r.stderr
    out = tmp_path / 'ds' / 'ultrafold' / 'structure'
    rows = list(csv.DictReader(open(out / 'alu_paired_by_cutoff.tsv'), delimiter='\t'))
    a1 = {(x['cutoff_type'], x['cutoff_value']): x['paired_with_opposite']
          for x in rows if x['alu'] == 'A1'}
    assert a1[('nbp', '1')] == 'true' and a1[('nbp', '5')] == 'false'
    assert a1[('frac', '50')] == 'true'     # 3 of 3 bases
    stems = list(csv.DictReader(open(out / 'arcs_collapsed.tsv'), delimiter='\t'))
    assert len(stems) == 1 and stems[0]['n_bp'] == '3' and stems[0]['arc_class'] == 'alu_opposite'
    model = list(csv.DictReader(open(out / 'model.tsv'), delimiter='\t'))[0]
    assert model['mfe_match'] == 'true'


def test_mfe_mismatch_stops_and_writes_nothing(tmp_path):
    write_inputs(tmp_path, stored_mfe='-13.5000')
    r = run_script(tmp_path)
    assert r.returncode == 1 and 'stored native_mfe' in r.stderr
    assert not (tmp_path / 'ds' / 'ultrafold' / 'structure').exists()


def test_length_mismatch_stops(tmp_path):
    write_inputs(tmp_path, stored_mfe='-12.3000')
    (tmp_path / 'ds' / 'ultrafold' / 'candidates.tsv').write_text(
        'transcript_id\tgene_name\tmrna_len\tutr5_len\tcds_len\tutr3_len\nT1\tG\t99\t2\t6\t91\n')
    r = run_script(tmp_path)
    assert r.returncode == 1 and 'mrna_len' in r.stderr
