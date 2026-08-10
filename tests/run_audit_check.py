#!/usr/bin/env python3
"""End-to-end check for bin/01d_family_audit.py.

Runs the family pipeline on the shared synthetic fixture, then audits it for
nucleotide-level similarity that protein-based blocking misses.

The fixture plants exactly one such case: S1 and S2 have unrelated proteins
(so they are singletons in different families) but ~92% identical 3'UTRs. That
pair is invisible to the protein search and must be surfaced by the audit —
it is the whole reason the audit exists. Every other 3'UTR is independent
random sequence, so a correct audit reports that pair and nothing else.

    python tests/run_audit_check.py
    python tests/run_audit_check.py --keep-outputs /tmp/audit_check -v

Requires seqkit and mmseqs; skips cleanly if either is absent.
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_family_check import COHORT, REPO_ROOT, build_fixture  # noqa: E402

LEVEL = 'medium'


def run(cmd: list[str], env: dict, verbose: bool) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, env=env, cwd=REPO_ROOT, text=True,
                          capture_output=not verbose)


def check(out_dir: Path, failures: list[str]) -> None:
    def expect(label: str, cond: bool) -> None:
        print(f"{'PASS' if cond else 'FAIL'}  {label}")
        if not cond:
            failures.append(label)

    with open(out_dir / 'audit_summary.tsv') as fh:
        summary = {r['region']: r for r in csv.DictReader(fh, delimiter='\t')}
    with open(out_dir / 'audit_pairs.tsv') as fh:
        pairs = list(csv.DictReader(fh, delimiter='\t'))
    with open(out_dir / 'audit_per_gene.tsv') as fh:
        per_gene = {r['gene']: r for r in csv.DictReader(fh, delimiter='\t')
                    if r['region'] == '3UTR'}

    expect("3UTR region audited", '3UTR' in summary)
    if '3UTR' not in summary:
        return

    def gene(name: str) -> str:
        return f'human|GENE{name}|TX{name}.1'

    # The planted cross-family pair must be found.
    s_pair = [p for p in pairs
              if {p['gene_a'], p['gene_b']} == {gene('S1'), gene('S2')}]
    expect("planted S1/S2 cross-family pair detected", len(s_pair) == 1)
    if s_pair:
        expect("detected at high identity (>=0.85)",
               float(s_pair[0]['identity']) >= 0.85)
        expect("detected at full coverage (>=0.9)",
               float(s_pair[0]['qcov']) >= 0.9 and float(s_pair[0]['tcov']) >= 0.9)
        expect("pair is genuinely cross-family",
               s_pair[0]['family_a'] != s_pair[0]['family_b'])

    # And nothing else, since every other 3'UTR is independent.
    expect("no spurious cross-family pairs", len(pairs) == 1)

    expect("S1 max cross-family identity is high",
           float(per_gene[gene('S1')]['max_cross_family_identity']) >= 0.85)
    expect("S2 names S1 as its partner",
           per_gene[gene('S2')]['cross_family_partner'] == gene('S1'))

    # Genes whose 3'UTRs are unrelated must show no cross-family similarity,
    # including same-family genes A1/A2 whose UTRs are independent random.
    for nm in ('A1', 'A2', 'B1', 'S3'):
        expect(f"{nm} has no cross-family 3UTR match",
               float(per_gene[gene(nm)]['max_cross_family_identity']) == 0.0)

    expect("summary counts exactly one gene pair at >=0.9 identity",
           int(summary['3UTR']['n_ge_90']) == 2)   # both members of the pair
    expect("summary reports a non-zero max",
           float(summary['3UTR']['max']) >= 0.85)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--keep-outputs', metavar='DIR')
    ap.add_argument('-v', '--verbose', action='store_true')
    args = ap.parse_args()

    for tool in ('seqkit', 'mmseqs'):
        if not shutil.which(tool):
            print(f"SKIP: {tool} not on PATH.")
            return 0

    tmp = Path(args.keep_outputs) if args.keep_outputs else \
        Path(tempfile.mkdtemp(prefix='audit_check_'))
    tmp.mkdir(parents=True, exist_ok=True)
    runs_root = tmp / 'runs'

    try:
        build_fixture(runs_root)
        env = dict(os.environ, RUNS_ROOT=str(runs_root))

        proc = run([sys.executable, str(REPO_ROOT / 'bin' / '01c_family.py'),
                    '--cohort', COHORT], env, args.verbose)
        if proc.returncode != 0:
            print("01c_family.py failed:")
            print(proc.stdout or ''); print(proc.stderr or '')
            return 1

        proc = run([sys.executable, str(REPO_ROOT / 'bin' / '01d_family_audit.py'),
                    '--cohort', COHORT, '--level', LEVEL]
                   + (['--verbose'] if args.verbose else []), env, args.verbose)
        if proc.returncode != 0:
            print("01d_family_audit.py failed:")
            print(proc.stdout or ''); print(proc.stderr or '')
            return 1

        hits = list((runs_root / '_cohorts' / COHORT / 'family')
                    .glob(f'*/audit/{LEVEL}/audit_summary.tsv'))
        if len(hits) != 1:
            print(f"expected one audit_summary.tsv, found {len(hits)}")
            return 1

        failures: list[str] = []
        check(hits[0].parent, failures)

        print()
        if failures:
            print(f"{len(failures)} FAILURE(S):")
            for f in failures:
                print("  " + f)
            return 1
        print("All checks passed.")
        return 0
    finally:
        if not args.keep_outputs:
            shutil.rmtree(tmp, ignore_errors=True)
        else:
            print(f"\nOutputs kept at {tmp}")


if __name__ == '__main__':
    sys.exit(main())
