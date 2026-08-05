#!/usr/bin/env python3
"""End-to-end check for bin/01c_family.py.

Builds a synthetic CDS set with a known family structure, runs the full
translate -> search -> cluster pipeline, and asserts the recovered families
match. Requires seqkit and mmseqs on PATH; skips cleanly if either is absent.

    python tests/run_family_check.py
    python tests/run_family_check.py --keep-outputs /tmp/family_check -v

Designed family structure:

    A1, A2, A3   ~95% / ~86% aa identity to A1   -> one family
    B1, B2       ~72% aa identity                -> one family
    S1, S2, S3   independent                     -> three singletons
    X1           A1 with a planted internal stop -> must join family A
    L1, L2       Cys/Gly/Ser-rich near-identical -> one family

Two load-bearing cases, both guarding silent under-merging:

X1 — `seqkit translate --trim` truncates at the first stop rather than
trimming the terminal one, which would reduce X1 to an ~80 aa stub, fail the
coverage floors, and leave it a spurious singleton.

L1/L2 — a keratin-associated-protein-like low-complexity pair. mmseqs masks
these with tantan in the prefilter (`--mask 1`, the default), so they vanish
from the main search entirely: no self-hit, and no edge to each other. Without
the `--mask 0` rescue pass they end up as two singletons despite being near
identical. Modelled on a real human MANE case (a ~33%-Cys KRTAP whose
paralogues hit at e-values ~1e-71 while it had no self-hit).
"""
from __future__ import annotations

import argparse
import csv
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
COHORT = 'synthetic_check'

BASES = 'ACGT'
STOPS = {'TAA', 'TAG', 'TGA'}
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES
          if a + b + c not in STOPS]
N_CODONS = 160

# One codon per residue; only used to turn the low-complexity protein fixture
# back into a CDS the pipeline can translate.
BACK_TRANSLATE = {
    'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT',
    'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTT', 'M': 'ATG', 'N': 'AAT',
    'P': 'CCT', 'Q': 'CAA', 'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT',
    'W': 'TGG', 'Y': 'TAT',
}

# Cys/Gly/Ser-rich motif in the style of a keratin-associated protein.
LOW_COMPLEXITY_CORE = (
    'MGCSGCSGGCGSSCGGCGSSCGGCGSGYGGCGSGCCVPVCCCKPVCCCVPACSCSSCGSCGGSKGVCGSCGG'
    'CKGGCGSCGGSKGGCGSSCCVPVCCSSSCGSCGGSKGVCGFRGGSKGGCGSCGCSQCSCYKPCCCSSGCGSS'
    'CCQSSCCKPSCSQSSCCKPCCSQSSCCKPCCCSSGCGSSCCQSSCCKPSCSQSSCCKPCCSQSSCC'
)


def back_translate(protein: str) -> str:
    return ''.join(BACK_TRANSLATE[aa] for aa in protein) + 'TAA'


def build_fixture(runs_root: Path) -> None:
    """Write extracted_CDS.fa for the `human_test` dataset under runs_root."""
    rng = random.Random(20260804)

    def cds(n):
        return [rng.choice(CODONS) for _ in range(n)]

    def mutate(codons, frac):
        out = list(codons)
        for i in rng.sample(range(len(out)), int(len(out) * frac)):
            out[i] = rng.choice(CODONS)
        return out

    def seq(codons):
        return 'ATG' + ''.join(codons) + 'TAA'

    a1 = cds(N_CODONS)
    b1 = cds(N_CODONS)
    x1 = list(a1)
    x1[80] = 'TGA'                       # planted internal stop

    lc1 = LOW_COMPLEXITY_CORE
    lc2 = LOW_COMPLEXITY_CORE.replace('GGCGSSCGG', 'GGCGSTCGG')

    records = {
        'A1': seq(a1), 'A2': seq(mutate(a1, 0.05)), 'A3': seq(mutate(a1, 0.14)),
        'B1': seq(b1), 'B2': seq(mutate(b1, 0.28)),
        'S1': seq(cds(N_CODONS)), 'S2': seq(cds(N_CODONS)), 'S3': seq(cds(N_CODONS)),
        'X1': seq(x1),
        'L1': back_translate(lc1), 'L2': back_translate(lc2),
    }

    out = runs_root / 'human_test' / 'extracted_regions'
    out.mkdir(parents=True, exist_ok=True)
    with open(out / 'extracted_CDS.fa', 'w') as fh:
        for name, s in records.items():
            fh.write(f">GENE{name}_TX{name}.1_CDS\n{s}\n")


def find_family_tsv(runs_root: Path) -> Path:
    hits = list((runs_root / '_cohorts' / COHORT / 'family').glob('*/family.tsv'))
    if len(hits) != 1:
        raise AssertionError(
            f"expected exactly one family.tsv, found {len(hits)}: {hits}")
    return hits[0]


def check(rows: list[dict], failures: list[str]) -> None:
    by_gene = {r['gene_id']: r for r in rows}

    def expect(label: str, cond: bool) -> None:
        print(f"{'PASS' if cond else 'FAIL'}  {label}")
        if not cond:
            failures.append(label)

    expect("all 11 genes present", len(rows) == 11)
    expect("no missing family assignment",
           all(r.get('family_id_medium') for r in rows))

    fam = {g: by_gene[g]['family_id_medium'] for g in by_gene}
    a_members = {'GENEA1', 'GENEA2', 'GENEA3', 'GENEX1'}
    expect("family A holds A1/A2/A3/X1",
           len({fam[g] for g in a_members}) == 1)
    expect("family B holds B1/B2", fam['GENEB1'] == fam['GENEB2'])
    expect("A and B are distinct families", fam['GENEA1'] != fam['GENEB1'])
    for s in ('GENES1', 'GENES2', 'GENES3'):
        expect(f"{s} is a singleton",
               by_gene[s]['family_size_medium'] == '1')
    expect("singletons are mutually distinct",
           len({fam['GENES1'], fam['GENES2'], fam['GENES3']}) == 3)

    # The --trim regression guard: X1 must retain full length and still cluster.
    expect("X1 flagged as having an internal stop",
           by_gene['GENEX1']['had_internal_stop'] == 'true')
    expect("X1 NOT truncated (same length as A1)",
           by_gene['GENEX1']['protein_len'] == by_gene['GENEA1']['protein_len'])
    expect("only X1 has an internal stop",
           sum(1 for r in rows if r['had_internal_stop'] == 'true') == 1)

    # The --mask 0 rescue guard: tantan removes this pair from the main search
    # entirely, so without the rescue pass they are two singletons.
    expect("low-complexity pair recovered into one family",
           fam['GENEL1'] == fam['GENEL2'])
    expect("low-complexity pair not merged into an unrelated family",
           by_gene['GENEL1']['family_size_medium'] == '2')

    expect("strict is at least as fine as medium",
           len({r['family_id_strict'] for r in rows})
           >= len({r['family_id_medium'] for r in rows}))
    expect("loose is at least as coarse as strict",
           len({r['family_id_loose'] for r in rows})
           <= len({r['family_id_strict'] for r in rows}))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--keep-outputs', metavar='DIR',
                        help="Write outputs here and leave them in place.")
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    for tool in ('seqkit', 'mmseqs'):
        if not shutil.which(tool):
            print(f"SKIP: {tool} not on PATH.")
            return 0

    tmp = Path(args.keep_outputs) if args.keep_outputs else \
        Path(tempfile.mkdtemp(prefix='family_check_'))
    tmp.mkdir(parents=True, exist_ok=True)
    runs_root = tmp / 'runs'

    try:
        build_fixture(runs_root)

        env = dict(os.environ, RUNS_ROOT=str(runs_root))
        cmd = [sys.executable, str(REPO_ROOT / 'bin' / '01c_family.py'),
               '--cohort', COHORT]
        if args.verbose:
            cmd.append('--verbose')
        proc = subprocess.run(cmd, env=env, cwd=REPO_ROOT,
                              capture_output=not args.verbose, text=True)
        if proc.returncode != 0:
            print("01c_family.py failed:")
            print(proc.stdout or '')
            print(proc.stderr or '')
            return 1

        with open(find_family_tsv(runs_root)) as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))

        failures: list[str] = []
        check(rows, failures)

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
