"""tools/rnafold_shuffles/worker.sh: clears stale MFE-only results, keeps finished ones."""
import os
import subprocess
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WRAPPER = REPO / 'tools' / 'rnafold_shuffles' / 'worker.sh'


def run(n_shuffles, with_results, with_raw):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fa = tmp / 'seq1.fa'
        fa.write_text('>seq1\nACGU\n')
        stub = tmp / 'stub.sh'
        stub.write_text('#!/bin/bash\necho STUB-RAN\n')
        stub.chmod(0o755)
        tool = tmp / 'runs' / 'ultrafold_shuffle_set' / 'rnafold_shuffles'
        (tool / 'results').mkdir(parents=True)
        (tool / 'raw_shuffles').mkdir()
        res = tool / 'results' / 'results_seq1.csv'
        if with_results:
            res.write_text('seq_name,mfe\nseq1,-1\n')
        if with_raw:
            (tool / 'raw_shuffles' / 'seq1_raw_shuffles.csv').write_text('transcript_id,iteration,shuffled_mfe\nseq1,1,-1\n')
        env = dict(os.environ, DATASET='ultrafold_shuffle_set', TOOL='rnafold_shuffles', RUNS_ROOT=str(tmp / 'runs'),
                   N_SHUFFLES=str(n_shuffles), REAL_WORKER=str(stub))
        out = subprocess.run(['bash', str(WRAPPER), str(fa)], env=env, capture_output=True, text=True, check=True)
        return res.exists(), out.stdout


def test_stale_mfe_only_result_is_removed_when_shuffles_requested():
    exists, out = run(100, with_results=True, with_raw=False)
    assert not exists and 'Removing stale' in out and 'STUB-RAN' in out


def test_finished_result_with_raw_shuffles_is_kept():
    exists, out = run(100, with_results=True, with_raw=True)
    assert exists and 'Removing' not in out


def test_calibration_mode_never_removes():
    exists, out = run(0, with_results=True, with_raw=False)
    assert exists and 'Removing' not in out
