"""tools/rnafold_structure/worker.sh: keeps the structure, writes bpp only on request,
and a MFE-only result does not block a later WITH_BPP=1 run."""
import glob
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
WORKER = REPO / 'tools' / 'rnafold_structure' / 'worker.sh'


def find_rnafold():
    for c in [os.environ.get('RNAFOLD_BIN'), shutil.which('RNAfold'),
              *glob.glob(os.path.expanduser('~/miniconda3/pkgs/viennarna-*/bin/RNAfold'))]:
        if c and os.access(c, os.X_OK):
            return c


pytestmark = pytest.mark.skipif(find_rnafold() is None, reason='RNAfold not available')


def run(tmp, with_bpp):
    env = dict(os.environ, DATASET='ultrafold_shuffle_set', TOOL='rnafold_structure',
               RUNS_ROOT=str(tmp / 'runs'), RNAFOLD_BIN=find_rnafold(), WITH_BPP=str(with_bpp))
    return subprocess.run(['bash', str(WORKER), str(tmp / 'X.fa')], env=env,
                          capture_output=True, text=True, check=True).stdout


def test_structure_bpp_and_skip():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        (tmp / 'X.fa').write_text('>X\nGGGGAAAACCCC\n')
        res = tmp / 'runs' / 'ultrafold_shuffle_set' / 'rnafold_structure' / 'results'
        assert 'Done X (MFE -5.40)' in run(tmp, 0)
        assert (res / 'X.fold').read_text().splitlines()[2] == '((((....)))) ( -5.40)'
        assert not (res / 'X.bpp.tsv').exists()
        assert 'Skip' in run(tmp, 0)
        assert 'Done' in run(tmp, 1)                    # MFE-only result does not block -p
        bpp = (res / 'X.bpp.tsv').read_text().splitlines()
        assert bpp[0] == 'i\tj\tprob' and any(l.startswith('0\t11\t') for l in bpp[1:])
        assert 'Skip' in run(tmp, 1)
