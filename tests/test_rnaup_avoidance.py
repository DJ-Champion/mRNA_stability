"""tools/rnaup_avoidance + bin/rnaup_make_tranches.py, with a stub RNAup (the real binary is not needed)."""
import csv
import gzip
import os
import shutil
import subprocess
import sys
import random
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'tools' / 'rnaup_avoidance'))
sys.path.insert(0, str(REPO / 'bin'))

import avoidance  # noqa: E402
import rnaup_make_tranches  # noqa: E402

ESL = shutil.which('esl-shuffle') or '/opt/hmmer/bin/esl-shuffle'

# Stub with RNAup's summary-line format. Energy depends on (ncRNA, target) only,
# so the test can recompute every sum independently. Every 3rd target "does not interact".
STUB = '''#!/usr/bin/env python3
import sys
recs = []
for line in sys.stdin:
    line = line.strip()
    if line.startswith('>'): recs.append([line[1:], ''])
    elif line: recs[-1][1] += line
nc = recs[0]
print('>' + nc[0]); print(nc[1])
for name, seq in recs[1:]:
    print('>' + name); print(seq)
    k = seq.count('G') + len(nc[1]) % 5
    if k % 3:
        gi = -(k * 0.37); gt = gi + 1.9 + 1.3
        print('(((.&.))) 1,3 : 5,7 (%.2f = %.2f + 1.90 + 1.30)' % (gt, gi))
'''


def stub_energy(nc_seq, seq):
    k = seq.count('G') + len(nc_seq) % 5
    if not k % 3:
        return 0.0, 0.0
    gi = round(-(k * 0.37), 2)
    return round(gi + 1.9 + 1.3, 2), gi


def test_parse_rnaup_interacting_and_non_interacting():
    text = ('>nc\nACGU\n>s0\nGGGG\n(((.&.))) 1,3 : 5,7 (-5.20 = -8.40 + 1.90 + 1.30)\n'
            '>s1\nCCCC\n>s2\nAAAA\n(.&.) 1,2 : 3,4 (+0.40 = -1.00 + 0.50 + 0.90)\n')
    res = avoidance.parse_rnaup(text, ['s0', 's1', 's2'])
    assert res == {'s0': (-5.2, -8.4), 's1': (0.0, 0.0), 's2': (0.4, -1.0)}


def test_parse_rnaup_truncated_output_is_an_error_not_zeros():
    with pytest.raises(RuntimeError, match='lacks 1/2'):
        avoidance.parse_rnaup('>nc\nAC\n>s0\nGG\n', ['s0', 's1'])


def _rand(rng, n):
    return ''.join(rng.choice('ACGU') for _ in range(n))


@pytest.fixture
def env(tmp_path):
    rng = random.Random(7)
    stub = tmp_path / 'RNAup'
    stub.write_text(STUB)
    stub.chmod(0o755)
    fams = ['tRNA', 'rRNA', 'SRP', 'snRNA', 'Y_RNA', 'snoRNA', 'tRNA']
    ncs = {f'URS{i}_9606': (fam, _rand(rng, 70 + i)) for i, fam in enumerate(fams)}
    (tmp_path / 'nc.fa').write_text(''.join(f'>{k} desc\n{s}\n' for k, (_, s) in ncs.items()))
    (tmp_path / 'fam.tsv').write_text('# c\nncrna_id\tfamily\n' + ''.join(f'{k}\t{f}\n' for k, (f, _) in ncs.items()))
    genes = {f'G{i}_T{i}_start_codon_region': _rand(rng, 63) for i in range(3)}
    tr = tmp_path / 'tranches' / 'start_codon_pm30'
    tr.mkdir(parents=True)
    fa = tr / 'start_codon_pm30_0001.fa'
    fa.write_text(''.join(f'>{k}\n{s}\n' for k, s in genes.items()))
    tool = tmp_path / 'tool'
    e = dict(os.environ, RESULTS_DIR=str(tool / 'results'), TOOL_DIR=str(tool), TMPDIR=str(tmp_path / 'scratch'),
             RNAUP_BIN=str(stub), ESL_SHUFFLE_BIN=ESL, NCRNA_FASTA=str(tmp_path / 'nc.fa'),
             NCRNA_FAMILIES=str(tmp_path / 'fam.tsv'), N_SHUFFLES='20', MIN_VALID_PERCENT='90',
             AVOID_KEEP_CKPT='1', SLURM_CPUS_PER_TASK='2')
    # fam.tsv above lists tRNA twice, which is fine; ids are what must match the FASTA
    return dict(fa=fa, e=e, tool=tool, genes=genes, ncs=ncs)


def run(env):
    return subprocess.run([sys.executable, str(REPO / 'tools/rnaup_avoidance/avoidance.py'), str(env['fa'])],
                          env=env['e'], capture_output=True, text=True, check=True)


def read_results(env):
    with open(env['tool'] / 'results' / 'results_start_codon_pm30_0001.csv') as f:
        return list(csv.DictReader(f))


def test_end_to_end_sums_families_and_z(env):
    run(env)
    rows = read_results(env)
    assert [r['seq_name'] for r in rows] == list(env['genes'])
    with gzip.open(env['tool'] / 'raw_shuffles' / 'start_codon_pm30_0001_shuffle_sums.csv.gz', 'rt') as f:
        raw = list(csv.DictReader(f))
    assert len(raw) == 3 * 20
    r = rows[0]
    seq = env['genes'][r['seq_name']]
    tot = [stub_energy(s, seq) for _, (_, s) in env['ncs'].items()]
    assert float(r['rnaup_G_total_sum']) == pytest.approx(sum(t for t, _ in tot), abs=1e-3)
    assert float(r['avoidance_G_int_min']) == pytest.approx(min(i for _, i in tot), abs=1e-3)
    assert float(r['rnaup_G_total_mean']) == pytest.approx(sum(t for t, _ in tot) / 7, abs=1e-3)
    trna = [stub_energy(s, seq)[1] for _, (f, s) in env['ncs'].items() if f == 'tRNA']
    assert float(r['avoidance_G_int_sum_tRNA']) == pytest.approx(sum(trna), abs=1e-3)
    # Z recomputed from the persisted per-shuffle sums with the repo's formula
    shuf = [float(x['G_total_sum']) for x in raw if x['seq_name'] == r['seq_name']]
    med, mad = avoidance.robust_stats([round(x, 2) for x in shuf])
    z = avoidance.robust_z(round(sum(t for t, _ in tot), 2), med, mad)
    if z is None:                       # MAD == 0 -> NA, as in the RNAfold worker
        assert r['rnaup_G_total_sum_robust_z'] == 'NA'
    else:
        assert float(r['rnaup_G_total_sum_robust_z']) == pytest.approx(z, abs=1e-3)
    gc = (seq.count('G') + seq.count('C')) / 63
    assert float(r['gc_content_mRNA_window']) == pytest.approx(gc, abs=1e-4)
    assert float(r['gc_diff_mRNA_vs_ncRNA']) == pytest.approx(gc - float(r['gc_ncRNA_mean']), abs=1e-4)
    assert r['window'] == 'start_codon_pm30' and r['n_shuffles'] == '20' and r['n_ncrna'] == '7'


def test_restart_reuses_the_same_shuffles(env):
    run(env)
    first = (env['tool'] / 'raw_shuffles' / 'start_codon_pm30_0001_shuffle_sums.csv.gz')
    a = gzip.open(first, 'rt').read()
    (env['tool'] / 'results' / 'results_start_codon_pm30_0001.csv').unlink()
    first.unlink()
    # a lost ncRNA checkpoint is recomputed against the persisted shuffles, not new ones
    next((env['tool'] / 'ckpt' / 'start_codon_pm30_0001').glob('URS3_*.tsv.gz')).unlink()
    run(env)
    assert gzip.open(first, 'rt').read() == a


def test_ncrna_fasta_and_family_map_must_agree(env):
    (env['tool'].parent / 'fam.tsv').write_text('ncrna_id\tfamily\nURS0_9606\ttRNA\n')
    p = subprocess.run([sys.executable, str(REPO / 'tools/rnaup_avoidance/avoidance.py'), str(env['fa'])],
                       env=env['e'], capture_output=True, text=True)
    assert p.returncode != 0 and 'disagree' in p.stderr


def test_shipped_family_map_is_34_ncrnas_in_6_families():
    fam = avoidance.load_families(REPO / 'configs' / 'ncrna' / 'families.tsv')
    counts = {f: sum(v == f for v in fam.values()) for f in avoidance.FAMILIES}
    assert counts == {'tRNA': 20, 'rRNA': 1, 'SRP': 1, 'snRNA': 5, 'Y_RNA': 3, 'snoRNA': 4}
    assert len(fam) == 34


def test_tranche_builder_filters_short_windows(tmp_path):
    ex = tmp_path / 'runs' / 'd' / 'extracted_regions'
    ex.mkdir(parents=True)
    (ex / 'extracted_start_codon_region.fa').write_text(
        '>a\n' + 'A' * 63 + '\n>b\n' + 'C' * 40 + '\n>c\n' + 'G' * 62 + 'N\n>d\n' + 'U' * 63 + '\n>e\n' + 'G' * 63 + '\n')
    (ex / 'extracted_5UTR.fa').write_text('>a\n' + 'A' * 80 + '\n>b\n' + 'C' * 19 + '\n>c\n' + 'G' * 25 + '\n')
    subprocess.run([sys.executable, str(REPO / 'bin/rnaup_make_tranches.py'), '-d', 'd', '-M', '2',
                    '--regions', 'start_codon_pm30', 'utr5', '--runs-root', str(tmp_path / 'runs')],
                   check=True, capture_output=True)
    lst = (tmp_path / 'runs/d/lists/rnaup_start_codon_pm30.txt').read_text().split()
    assert [Path(p).name for p in lst] == ['start_codon_pm30_0001.fa', 'start_codon_pm30_0002.fa']
    assert [n for n, _ in rnaup_make_tranches.iter_fasta(lst[0])] == ['a', 'd']
    assert [n for n, _ in rnaup_make_tranches.iter_fasta(lst[1])] == ['e']
    u = (tmp_path / 'runs/d/lists/rnaup_utr5.txt').read_text().split()
    recs = list(rnaup_make_tranches.iter_fasta(u[0]))
    assert recs == [('c', 'G' * 25), ('a', 'A' * 80)]      # whole UTRs, >= 20 nt, sorted by length


def test_variable_length_tranche_end_to_end(env):
    env['fa'].write_text(''.join(f'>{k}\n{s[:30 + 10 * i]}\n' for i, (k, s) in enumerate(env['genes'].items())))
    run(env)
    assert [r['window_length'] for r in read_results(env)] == ['30', '40', '50']


REAL_2_4_7 = '''>nc
>t1
((((((((......((((((((&))))))))...))))))))   1,22  :   1,19  (0.25 = -37.50 + 20.73 + 17.02)
GGGCGCCCGGAUCCGGGCGCCC&GGGCGCCCAAAGGGCGCCC
>t2
WARNING: RNAduplex structure doesn't match any structure of RNAup structure ensemble
.&.   0,0   :   0,0   (100.00 = 100.00 + 2.31 + inf)
&
ERROR
diff inf
>t3
(((((((((((...((((((((&)))))))))))))))))))   1,22  :   1,19  (-2.48 = -40.24 + 20.73 + 17.02)
'''


def test_parse_real_rnaup_2_4_7_output():
    res = avoidance.parse_rnaup(REAL_2_4_7, ['t1', 't2', 't3'])
    assert res == {'t1': (0.25, -37.5), 't2': (0.0, 0.0), 't3': (-2.48, -40.24)}


@pytest.mark.skipif(not Path('/home/djchampion/miniforge3/envs/yaml_env/bin/RNAup').exists(), reason='needs real RNAup')
def test_real_rnaup_runs_without_stray_files(tmp_path):
    res = avoidance.run_rnaup('/home/djchampion/miniforge3/envs/yaml_env/bin/RNAup', avoidance.DEFAULT_RNAUP_ARGS.split(),
                              'nc', 'GGGCGCCCGGAUCCGGGCGCCC', [('s0', 'GGGCGCCCUUUGGGCGCCC'), ('s1', 'AAAAAAAAAAAAAA')], tmp_path)
    assert res['s0'][1] < -30 and res['s1'] == (0.0, 0.0)
    assert list(tmp_path.iterdir()) == []
