# Tests

Synthetic-fixture tests for the metric plugins. First instance is the
three-model NMD-fragility sanity check; the same pattern extends to the
other plugins as they land.

## Running

From the repo root:

    python -m pytest tests/

Or, for a quick non-pytest sanity check with diff output:

    python tests/run_synthetic_check.py
    python tests/run_synthetic_check.py --keep-outputs /tmp/nmd_check -v

The family-clustering end-to-end check is also pytest-free, and skips
cleanly if `seqkit` or `mmseqs` is not on PATH:

    python tests/run_family_check.py
    python tests/run_family_check.py --keep-outputs /tmp/family_check -v

The audit check runs the family pipeline then audits it, using the same
fixture (which plants one cross-family 3'UTR pair for it to find):

    python tests/run_audit_check.py
    python tests/run_audit_check.py --keep-outputs /tmp/audit_check -v

## Layout

    tests/
      conftest.py                          # puts repo root on sys.path
      run_synthetic_check.py               # CLI runner (no pytest required)
      run_family_check.py                  # family clustering, end-to-end
      run_audit_check.py                   # nucleotide leakage audit, end-to-end
      test_nmd_fragility_synthetic.py      # the three-model sanity check
      test_lib_family.py                   # lib/family.py units
      fixtures/
        _generate_and_verify.py            # generates the CDS FASTA + self-checks
        synthetic_5x200.gff3               # 5-exon × 200-nt transcript, + strand
        synthetic_5x200_minus.gff3         # same geometry, - strand twin
        synthetic_5x200_CDS.fa             # 999-nt CDS with targeted inserts (shared)
        manifest.tsv                       # minimal manifest for the synthetic transcript

The minus-strand GFF has identical genomic exon spans to the plus-strand
one, but `strand=-` and recomputed CDS phases. The FASTA represents the
transcribed sequence and is shared. Both GFFs parse to the same spliced
layout (junctions at 200/400/600/800, start_s=0); the
`_generate_and_verify.py` script confirms this.

## Why three assertion layers?

The NMD-fragility plugins are structurally similar but produce
different numbers, and most bugs would silently pass through a
single-model check. The fixture is deliberately constructed so that
each of three layers exposes a different fault class:

| Layer | What it tests                          | Bug class it catches                                    |
| ----- | -------------------------------------- | ------------------------------------------------------- |
| 1     | `zone_length` per model                | Threshold inequality, junction iteration, window arithmetic |
| 2     | `n_fragile_codons` / `n_alt_stops`     | In-frame vs out-of-frame iteration, codon-set membership    |
| 3     | `apply_nmd_rule=False` shifts the zone | `metric_config` not threaded through to the model fn       |

Each model produces a different expected value at every layer. A
discrepancy localises the bug to one model and one layer.

Every test is also parametrised over `strand ∈ {plus, minus}`, which
adds a fourth orthogonal axis: the strand-specific branches in
`_build_spliced_index`, `_genomic_to_spliced`, and
`_start_codon_genomic_pos`. Expected counts are identical between
strands by construction, so any divergence is a strand-handling bug.
Failure messages prefix the strand label (`[plus strand] ...` /
`[minus strand] ...`) for unambiguous diagnosis.

## Adding a new fixture

The synthetic-CDS pattern generalises cleanly:

1. **Choose a background**. GCC tandem is the cleanest neutral
   background for NMD-fragility — no fragile codon, no alt-stop in
   any frame, and trivially in-frame. For GC-content metrics, pick
   per-position bases that put each position in a single, known
   bucket (e.g. AAA codons for GC=0%, GCC for GC1=33%/GC2=33%/GC3=100%).

2. **Place targeted inserts**. Each insert should land on exactly one
   feature the plugin counts, with surrounding nt that don't add
   spurious hits. The current fixture's CGA and ATAACT inserts are
   the worked example.

3. **Compute expected counts by hand**, then add them to
   `fixtures/_generate_and_verify.py`'s `EXPECTED` dict. The script
   re-implements the plugin's logic and verifies the predictions —
   any drift between the test file and the generator is caught here
   before pytest sees it.

4. **Write the test** that asserts the expected counts exactly.
   Prefer parametrised tests so that adding a model or a config
   variant is a one-line change.

## Adding a new plugin

For each new metric plugin, add a `test_<plugin>_synthetic.py` that
mirrors the structure of `test_nmd_fragility_synthetic.py`:

- a `synthetic_paths` fixture that stages the relevant input files
  into `tmp_path`
- one or more parametrised tests over the configurations the plugin
  supports
- a `test_fixture_well_formed` smoke test

The cheap-win candidates that map directly onto the current fixture:

- `gc_by_position` — same FASTA, expected GC1/GC2/GC3 = 0/100/100
  for the GCC background; inserts at known positions shift exactly
  one of the three by a predictable delta.
- `kozak` / `stop_context` — read a short window around the start
  codon (positions 0–2) or stop codon (positions 996–998); the
  current fixture has both at known sequences.

## Regenerating the fixture

If the synthetic sequence or insert plan changes:

    cd tests/fixtures
    python _generate_and_verify.py

This rewrites `synthetic_5x200_CDS.fa` and self-checks the predictions
in `EXPECTED`. If the script fails, the test file's expected values
are out of sync with the fixture — fix one or the other before
committing.

## Dependencies

- `pytest` — `pip install pytest`
- `gffutils` — already required by the pipeline (used by the smoke
  test to validate the GFF fixture parses)

No project-level pytest config is needed. `conftest.py` injects the
repo root into `sys.path` so `from metrics._nmd_fragility_core import
run` resolves regardless of the working directory pytest is invoked
from. If a `pyproject.toml` with a proper package install lands later,
`conftest.py` can be slimmed or removed.
