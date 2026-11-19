"""tools/rnaup_avoidance (avoidance.py and worker.sh), with a stub RNAup unless a test says real."""
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
    fa = tmp_path / 'start_codon_pm30_0001.fa'
    fa.write_text(''.join(f'>{k}\n{s}\n' for k, s in genes.items()))
    tool = tmp_path / 'tool'
    e = dict(os.environ, RESULTS_DIR=str(tool / 'results'), TOOL_DIR=str(tool), TMPDIR=str(tmp_path / 'scratch'),
             RNAUP_BIN=str(stub), ESL_SHUFFLE_BIN=ESL, NCRNA_FASTA=str(tmp_path / 'nc.fa'),
             NCRNA_FAMILIES=str(tmp_path / 'fam.tsv'), N_SHUFFLES='20', MIN_VALID_PERCENT='90',
             AVOID_KEEP_CKPT='1', AVOID_CHUNK='25',
             SLURM_CPUS_PER_TASK='2', AVOID_REGION='start_codon_region')
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
    # a lost work-unit checkpoint is recomputed against the persisted shuffles, not new ones
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


def test_length_ok_rules():
    spec = 'start_codon_region:63 5UTR:>=20 other:10-20'
    assert avoidance.length_ok('start_codon_region', 63, spec)
    assert not avoidance.length_ok('start_codon_region', 62, spec)
    assert avoidance.length_ok('5UTR', 20, spec) and not avoidance.length_ok('5UTR', 19, spec)
    assert avoidance.length_ok('other', 20, spec) and not avoidance.length_ok('other', 21, spec)
    assert avoidance.length_ok('tail_region', 7, spec)          # no rule -> any length


def test_violating_window_is_refused(env):
    # no TOOL_REQUIRE_LEN in env: the built-in start_codon_region:63 still applies (5'UTR < 30 -> clipped)
    env['fa'].write_text('>G_T_start_codon_region\nACGUACGUAC\n')
    p = subprocess.run([sys.executable, str(REPO / 'tools/rnaup_avoidance/avoidance.py'), str(env['fa'])],
                       env=env['e'], capture_output=True, text=True)
    assert p.returncode != 0 and 'violates TOOL_REQUIRE_LEN' in p.stderr


@pytest.mark.parametrize('region', ['5UTR', 'stop_codon_region', 'tail_region', ''])
def test_only_the_start_codon_window_is_scored(env, region):
    env['e']['AVOID_REGION'] = region
    p = subprocess.run([sys.executable, str(REPO / 'tools/rnaup_avoidance/avoidance.py'), str(env['fa'])],
                       env=env['e'], capture_output=True, text=True)
    assert p.returncode != 0 and 'is not scored' in p.stderr


def test_config_fixes_region_and_defaults():
    # env overrides of the region/length rules are ignored; N_SHUFFLES / AVOID_CHUNK keep their defaults
    e = {k: v for k, v in os.environ.items() if k not in ('N_SHUFFLES', 'AVOID_CHUNK')}
    e.update(PROJECT_ROOT=str(REPO), TOOL_REGIONS='5UTR tail_region', TOOL_REQUIRE_LEN='5UTR:>=20')
    out = subprocess.run(['bash', '-c', 'source configs/tools/rnaup_avoidance.sh; '
                          'echo "$TOOL_REGIONS|$TOOL_REQUIRE_LEN|$N_SHUFFLES|$AVOID_CHUNK"'],
                         cwd=REPO, env=e, capture_output=True, text=True, check=True).stdout.strip()
    assert out == 'start_codon_region|start_codon_region:63|100|101'


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


WORKER = REPO / 'tools' / 'rnaup_avoidance' / 'worker.sh'


def run_worker(env, tmp_path, calib, check=True):
    """worker.sh through lib/paths.sh, against a fabricated run dir for the human_10 dataset."""
    runs = tmp_path / 'runs'
    ex = runs / 'human_10' / 'extracted_regions'
    ex.mkdir(parents=True, exist_ok=True)
    name = next(iter(env['genes']))
    (ex / 'manifest.tsv').write_text(
        'seqname\tgene_id\ttranscript_id\tregion\tlength\n' + f'{name}\tG0\tT0\tstart_codon_region\t63\n')
    fa = tmp_path / f'{name}.fa'
    fa.write_text(f'>{name}\n{env["genes"][name]}\n')
    e = dict(env['e'], DATASET='human_10', TOOL='rnaup_avoidance', RUNS_ROOT=str(runs), N_SHUFFLES='6')
    e.pop('AVOID_REGION')
    e['AVOID_KEEP_CKPT'] = '0'
    for k in ('RESULTS_DIR', 'TOOL_DIR'):
        e.pop(k)
    if calib:
        e['AVOID_CALIB'] = '1'
    out = subprocess.run(['bash', str(WORKER), str(fa)], env=e, capture_output=True, text=True)
    if not check:
        return runs / 'human_10' / 'rnaup_avoidance', name, out
    assert out.returncode == 0, out.stderr
    return runs / 'human_10' / 'rnaup_avoidance', name, out.stdout


def test_worker_writes_results_and_skips_when_done(env, tmp_path):
    tool, name, out = run_worker(env, tmp_path, calib=False)
    res = tool / 'results' / f'results_{name}.csv'
    assert 'Done' in out and res.exists()
    assert 'start_codon_pm30' in res.read_text()          # region taken from the manifest
    assert not (tool / 'ckpt' / name).exists()            # checkpoints cleaned after success
    assert 'Skip' in run_worker(env, tmp_path, calib=False)[2]


def test_worker_failure_logs_to_errors_file_and_stderr(env, tmp_path):
    (env['tool'].parent / 'fam.tsv').write_text('ncrna_id\tfamily\nURS0_9606\ttRNA\n')
    tool, name, out = run_worker(env, tmp_path, calib=False, check=False)
    assert out.returncode != 0 and 'disagree' in out.stderr
    assert 'disagree' in (tool / 'errors' / f'errors_{name}.log').read_text()
    assert not (tool / 'results' / f'results_{name}.csv').exists()


def test_calibration_mode_leaves_no_results(env, tmp_path):
    tool, name, out = run_worker(env, tmp_path, calib=True)
    assert 'Done' in out
    assert not list((tool / 'results').glob('*')) and not (tool / 'raw_shuffles').exists()
    # and a real run afterwards is not mistaken for done
    assert 'Done' in run_worker(env, tmp_path, calib=False)[2]


def test_chunk_size_does_not_change_results(env):
    run(env)
    a = read_results(env)
    shutil.rmtree(env['tool'])
    env['e']['AVOID_CHUNK'] = '7'           # same numbers from different work units (shuffles differ -> compare natives)
    run(env)
    b = read_results(env)
    for ra, rb in zip(a, b):
        for c in ('rnaup_G_total_sum', 'avoidance_G_int_sum', 'rnaup_G_total_min', 'rnaup_G_total_sum_tRNA'):
            assert ra[c] == rb[c]


def test_parse_rnaup_unparseable_energy_line_is_an_error():
    text = '>nc\nACGU\n>s0\n(((.&.))) 1,3 : 5,7 (-5.2e0 = -8.40 + 1.90 + 1.30)\n'
    with pytest.raises(RuntimeError, match='cannot parse'):
        avoidance.parse_rnaup(text, ['s0'])


def test_parse_rnaup_inf_in_any_opening_term_is_no_interaction():
    text = '>nc\nACGU\n>s0\n.&.   0,0   :   0,0   (100.00 = 100.00 + inf + 2.31)\n'
    assert avoidance.parse_rnaup(text, ['s0']) == {'s0': (0.0, 0.0)}


def test_finalize_reads_units_in_chunk_order_past_1000_chunks(tmp_path):
    # 1002 sequences, chunk 1 -> 1002 unit files; c1.1000 sorts before c1.101 as text
    n_seq = 1002
    seqsets = [('G', ['ACGU'] * n_seq)]
    ncs = [('nc', 'ACGU', 'tRNA')]
    ckpt = tmp_path / 'ckpt'
    ckpt.mkdir()
    for ci in range(n_seq):
        with gzip.open(avoidance.unit_path(ckpt, 'nc', 1, ci), 'wt') as f:
            f.write(f'{-(ci + 1):.2f}\t{-(ci + 1):.2f}\n')
    avoidance.finalize('G', seqsets, ncs, ckpt, 1, n_seq, tmp_path / 'res', tmp_path / 'raw', 'w')
    with gzip.open(tmp_path / 'raw' / 'G_shuffle_sums.csv.gz', 'rt') as f:
        raw = list(csv.DictReader(f))
    assert [float(r['G_total_sum']) for r in raw] == [-(i + 1) for i in range(1, n_seq)]
    with open(tmp_path / 'res' / 'results_G.csv') as f:
        assert float(next(csv.DictReader(f))['rnaup_G_total_sum']) == -1


def test_stale_shuffle_checkpoint_is_discarded(env):
    run(env)
    ckpt = env['tool'] / 'ckpt' / 'start_codon_pm30_0001'
    (env['tool'] / 'results' / 'results_start_codon_pm30_0001.csv').unlink()
    # the window was re-extracted (same name, new sequence) and N_SHUFFLES changed
    name = next(iter(env['genes']))
    new = 'GC' * 31 + 'G'
    env['fa'].write_text(''.join(f'>{k}\n{new if k == name else s}\n' for k, s in env['genes'].items()))
    env['e']['N_SHUFFLES'] = '15'
    p = run(env)
    assert 'stale' in p.stderr
    rows = {r['seq_name']: r for r in read_results(env)}
    assert float(rows[name]['gc_content_mRNA_window']) == pytest.approx(1.0)
    assert {r['n_shuffles'] for r in rows.values()} == {'15'}
    natives = [s for rec, s in avoidance.read_fasta(ckpt / 'shuffles.fa.gz') if rec == f'{name}|0']
    assert natives == [new]
