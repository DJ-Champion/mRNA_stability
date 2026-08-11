#!/usr/bin/env python3
"""
01c_family.py
Cluster a cohort of datasets into sequence families for blocked cross-validation.

A cohort is a set of datasets (1..n species) pooled into one all-vs-all
protein search. Connected components of the filtered similarity graph are
single-linkage families — the transitive closure required so that no test
sequence has a close homologue in its training fold.

Three stages, each independently skippable when its output is current:

    translate  extracted_CDS.fa (per member) -> proteins.faa
    search     proteins.faa                  -> hits.tsv
    cluster    hits.tsv                      -> family.tsv (+ QC)

Outputs land at:
    $RUNS_ROOT/_cohorts/<cohort>/family/<search_hash>/

Usage:
  ./bin/01c_family.py --cohort human_only
  ./bin/01c_family.py -c human_only --recluster     # re-threshold, reuse search
  ./bin/01c_family.py -c human_only --force         # redo everything
  ./bin/01c_family.py --list-cohorts

See FAMILY_CLUSTERING.md for the specification.
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import statistics
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

# --- Path bootstrap -----------------------------------------------------------
_THIS = Path(__file__).resolve()
_PROJECT_ROOT = _THIS.parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from lib.family import (  # noqa: E402
    Level, SearchParams, UnionFind, canonical_families, find_overlapping_gene_pairs,
    gene_exons_from_gff, iter_fasta, make_protein_id, normalised_bitscore,
    post_process_translation, search_hash, split_protein_id,
    validate_level_within_search,
)
from lib.gff import normalise_id, split_composite_fasta_id  # noqa: E402
from lib.paths import resolve_paths  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%H:%M:%S'
)
log = logging.getLogger('family')

# hits.tsv column order — must match --format-output in _run_search().
H_QUERY, H_TARGET, H_FIDENT, H_ALNLEN, H_EVALUE, H_BITS, H_QCOV, H_TCOV, \
    H_QLEN, H_TLEN = range(10)

_FORMAT_OUTPUT = 'query,target,fident,alnlen,evalue,bits,qcov,tcov,qlen,tlen'

_DEFAULT_MIN_PROTEIN_LEN = 30
_DEFAULT_MAX_INTERNAL_STOP_FRAC = 0.05
_DEFAULT_MAX_FAMILY_PCT = 0.02
_DEFAULT_MIN_OVERLAP_BP = 100


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def _load_yaml(path: Path) -> dict:
    try:
        import yaml
    except ImportError:
        log.error("PyYAML not installed (pip install pyyaml).")
        sys.exit(1)
    with open(path, 'r') as f:
        return yaml.safe_load(f) or {}


class Cohort:
    """Resolved, validated cohort configuration."""

    def __init__(self, name: str, config: dict, project_root: Path):
        self.name = name
        self.project_root = project_root

        members = config.get('members') or []
        if not members:
            _die(f"Cohort '{name}' has no members.")
        self.members: list[dict] = []
        seen_species: set[str] = set()
        for m in members:
            dataset = m.get('dataset')
            species = m.get('species')
            if not dataset or not species:
                _die(f"Cohort '{name}': every member needs 'dataset' and 'species'.")
            if species in seen_species:
                _die(f"Cohort '{name}': species '{species}' listed twice. "
                     f"Protein IDs would collide.")
            seen_species.add(species)
            self.members.append({
                'dataset': dataset,
                'species': species,
                'transl_table': int(m.get('transl_table', 1)),
            })

        s = config.get('search') or {}
        self.search = SearchParams(
            sensitivity=float(s.get('sensitivity', 7.5)),
            max_seqs=int(s.get('max_seqs', 2000)),
            evalue=float(s.get('evalue', 1e-3)),
            min_cov=float(s.get('min_cov', 0.4)),
            cov_mode=int(s.get('cov_mode', 0)),
            threads=int(s.get('threads', os.cpu_count() or 4)),
        )

        levels_cfg = config.get('levels') or {}
        if not levels_cfg:
            _die(f"Cohort '{name}' defines no levels.")
        self.levels = [
            Level(name=lname,
                  min_nbs=float(lc['min_nbs']),
                  min_cov=float(lc['min_cov']),
                  max_evalue=float(lc['max_evalue']),
                  min_fident=float(lc['min_fident']))
            for lname, lc in levels_cfg.items()
        ]
        # Loosest first, by the thresholds themselves rather than YAML order.
        self.levels.sort(key=lambda L: (L.min_nbs, L.min_cov, L.min_fident))

        problems: list[str] = []
        for lv in self.levels:
            problems.extend(validate_level_within_search(lv, self.search))
        if problems:
            _die("Level(s) looser than the search that produced the edges:\n  "
                 + "\n  ".join(problems)
                 + "\n\nFiltering can only tighten. Either tighten the level or "
                   "loosen 'search:' and re-run the search.")

        # Third edge source. Genes whose exons overlap are transcribed from the
        # same DNA, so their region features are related whether or not they
        # are homologous — a gap protein clustering cannot see. Not a
        # similarity threshold, so the same edges apply at every level.
        lo = config.get('locus_overlap') or {}
        self.locus_overlap_enabled = bool(lo.get('enabled', True))
        self.min_overlap_bp = int(lo.get('min_overlap_bp', _DEFAULT_MIN_OVERLAP_BP))

        t = config.get('translate') or {}
        self.min_protein_len = int(t.get('min_protein_len', _DEFAULT_MIN_PROTEIN_LEN))
        self.max_internal_stop_frac = float(
            t.get('max_internal_stop_frac', _DEFAULT_MAX_INTERNAL_STOP_FRAC))

        # The real constraint on a family is that it must pack into a fold.
        # Deriving the ceiling from k keeps that explicit: a family may occupy
        # at most `max_fold_fraction` of one fold, and a fold is 1/k of the
        # corpus. A fixed percentage silently becomes wrong when k changes.
        sel = config.get('selection') or {}
        self.k_folds = int(sel.get('k_folds', 5))
        if self.k_folds < 2:
            _die(f"selection.k_folds must be >= 2 (got {self.k_folds}).")
        self.max_fold_fraction = float(sel.get('max_fold_fraction', 0.25))
        explicit = sel.get('max_family_pct')
        self.ceiling_is_explicit = explicit is not None
        self.max_family_pct = (float(explicit) if self.ceiling_is_explicit
                               else self.max_fold_fraction / self.k_folds)

        b = config.get('binaries') or {}
        self.seqkit_bin = b.get('seqkit') or os.environ.get('SEQKIT_BIN') or 'seqkit'
        self.mmseqs_bin = b.get('mmseqs') or os.environ.get('MMSEQS_BIN') or 'mmseqs'

    @property
    def hash(self) -> str:
        return search_hash(
            self.search,
            [(m['dataset'], m['species'], m['transl_table']) for m in self.members])


def _die(msg: str) -> None:
    log.error(msg)
    sys.exit(1)


def _require_binary(configured: str, label: str) -> str:
    resolved = shutil.which(configured)
    if not resolved:
        _die(f"{label} not found: {configured!r} is not on PATH. "
             f"Set it under 'binaries:' in the cohort YAML or via "
             f"{label.upper()}_BIN.")
    return resolved


def _tool_version(binary: str) -> str:
    try:
        out = subprocess.run([binary, 'version'], capture_output=True,
                             text=True, timeout=30)
        return (out.stdout or out.stderr).strip().splitlines()[0]
    except Exception:  # noqa: BLE001 — version reporting must never abort a run
        return 'unknown'


def _is_current(output: Path, inputs: Iterable[Path]) -> bool:
    if not output.exists():
        return False
    out_mtime = output.stat().st_mtime
    for ip in inputs:
        ip = Path(ip)
        if ip.exists() and ip.stat().st_mtime > out_mtime:
            return False
    return True


# ---------------------------------------------------------------------------
# Stage 1 — translate
# ---------------------------------------------------------------------------

def _member_cds_fasta(cohort: Cohort, member: dict) -> Path:
    paths = resolve_paths(member['dataset'], cohort.project_root,
                          species=member['species'])
    return paths.extract_dir / 'extracted_CDS.fa'


def _member_canonical_gff(cohort: Cohort, member: dict) -> Path:
    paths = resolve_paths(member['dataset'], cohort.project_root,
                          species=member['species'])
    return paths.canonical_gff


def stage_translate(cohort: Cohort, out_dir: Path, force: bool) -> None:
    faa = out_dir / 'proteins.faa'
    ptsv = out_dir / 'proteins.tsv'
    pids = out_dir / 'proteins.ids'
    qc = out_dir / 'translate_qc.tsv'

    inputs = [_member_cds_fasta(cohort, m) for m in cohort.members]
    missing = [str(p) for p in inputs if not p.exists()]
    if missing:
        _die("Missing extracted CDS FASTA(s). Run 01_extract.py first:\n  "
             + "\n  ".join(missing))

    if not force and all(p.exists() for p in (faa, ptsv, pids)) \
            and _is_current(faa, inputs):
        log.info("Translation current — skipping.")
        return

    seqkit = _require_binary(cohort.seqkit_bin, 'seqkit')
    tmp_dir = out_dir / 'tmp'
    tmp_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict] = []
    qc_rows: list[dict] = []
    seen_ids: dict[str, str] = {}

    for member in cohort.members:
        dataset, species = member['dataset'], member['species']
        cds_fa = _member_cds_fasta(cohort, member)
        raw_faa = tmp_dir / f'{species}.raw.faa'

        log.info(f"Translating {dataset} ({species}, table "
                 f"{member['transl_table']}) …")
        # No --trim: it truncates at the first '*' or 'X' rather than
        # trimming the terminal stop. See lib.family.post_process_translation.
        cmd = [seqkit, 'translate',
               '-T', str(member['transl_table']),
               '-f', '1',
               str(cds_fa), '-o', str(raw_faa)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            _die(f"seqkit translate failed for {dataset}:\n{proc.stderr}")

        n_cds = n_internal = n_no_term = n_short = n_with_x = n_bad_header = 0
        for rec_id, aa_raw in iter_fasta(raw_faa):
            n_cds += 1
            split = split_composite_fasta_id(rec_id, region='CDS')
            if split is None:
                n_bad_header += 1
                continue
            gene_raw, transcript_id, _ = split
            res = post_process_translation(aa_raw)
            if res.had_internal_stop:
                n_internal += 1
            if not res.had_terminal_stop:
                n_no_term += 1
            if 'X' in res.aa:
                n_with_x += 1
            # Too short to align meaningfully, but NOT dropped: it stays in
            # the gene universe as a flagged singleton so a left_join from
            # manifest.tsv remains lossless. Silently omitting it would put an
            # unexplained NA in the blocking table.
            too_short = len(res.aa) < cohort.min_protein_len
            if too_short:
                n_short += 1

            protein_id = make_protein_id(
                species, normalise_id(gene_raw), transcript_id)
            if protein_id in seen_ids:
                _die(f"Duplicate protein ID {protein_id!r} from datasets "
                     f"{seen_ids[protein_id]!r} and {dataset!r}. Check the "
                     f"cohort member list.")
            seen_ids[protein_id] = dataset

            records.append({
                'protein_id': protein_id,
                'species': species,
                'gene_id': normalise_id(gene_raw),
                'transcript_id': transcript_id,
                'dataset': dataset,
                'protein_len': len(res.aa),
                'had_internal_stop': 'true' if res.had_internal_stop else 'false',
                'searched': 'false' if too_short else 'true',
                'aa': res.aa,
            })

        if n_bad_header:
            _die(f"{dataset}: {n_bad_header} FASTA header(s) did not match the "
                 f"'<gene>_<transcript>_CDS' pattern. Silently skipping them "
                 f"would under-merge families; fix 01_extract.py output first.")

        n_translated = n_cds - n_bad_header
        frac_internal = (n_internal / n_cds) if n_cds else 0.0
        if frac_internal > cohort.max_internal_stop_frac:
            _die(f"{dataset}: {n_internal}/{n_cds} ({frac_internal:.1%}) CDS have "
                 f"internal stop codons, above the {cohort.max_internal_stop_frac:.0%} "
                 f"ceiling. Likely the wrong transl_table, a frame problem "
                 f"upstream, or a poor annotation.")

        qc_rows.append({
            'dataset': dataset, 'species': species,
            'n_cds': n_cds, 'n_translated': n_translated,
            'n_searched': n_translated - n_short,
            'n_internal_stop': n_internal, 'n_no_terminal_stop': n_no_term,
            'n_too_short': n_short, 'n_with_X': n_with_x,
        })
        log.info(f"  {n_translated} proteins ({n_internal} internal stop, "
                 f"{n_short} below {cohort.min_protein_len} aa kept as "
                 f"unsearched singletons)")

    if not records:
        _die("No proteins survived translation.")

    # Only searchable proteins go to mmseqs; the rest stay in the universe.
    n_searched = 0
    with open(faa, 'w') as fh:
        for r in records:
            if r['searched'] != 'true':
                continue
            n_searched += 1
            fh.write(f">{r['protein_id']}\n")
            seq = r['aa']
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + '\n')

    # proteins.ids is the gene universe — every protein, searched or not. The
    # cluster stage and the R recipe in FAMILY_CLUSTERING.md both need it:
    # singletons never appear in hits.tsv, so the universe is not recoverable
    # from the alignment table.
    with open(pids, 'w') as fh:
        for r in records:
            fh.write(r['protein_id'] + '\n')

    cols = ['protein_id', 'species', 'gene_id', 'transcript_id', 'dataset',
            'protein_len', 'had_internal_stop', 'searched']
    _write_tsv(ptsv, cols, [{k: r[k] for k in cols} for r in records])
    _write_tsv(qc, list(qc_rows[0].keys()), qc_rows)

    shutil.rmtree(tmp_dir, ignore_errors=True)
    n_unsearched = len(records) - n_searched
    log.info(f"Pooled {len(records)} proteins ({n_searched} searchable"
             + (f", {n_unsearched} kept as unsearched singletons"
                if n_unsearched else "") + f") → {faa}")


# ---------------------------------------------------------------------------
# Stage 2 — search
# ---------------------------------------------------------------------------

def stage_search(cohort: Cohort, out_dir: Path, force: bool,
                 verbose: bool = False) -> None:
    faa = out_dir / 'proteins.faa'
    hits = out_dir / 'hits.tsv'

    if not force and _is_current(hits, [faa]):
        log.info("Search current — skipping.")
        return

    mmseqs = _require_binary(cohort.mmseqs_bin, 'mmseqs')
    tmp_dir = out_dir / 'mmseqs_tmp'
    s = cohort.search

    # No --min-seq-id here: identity is filtered downstream from the `fident`
    # column so that every threshold decision lives in exactly one auditable
    # place. --max-seqs is raised well above the default 300 because large
    # families (olfactory receptors ~400, C2H2 zinc fingers ~700) would
    # otherwise have hit lists truncated non-deterministically.
    cmd = [mmseqs, 'easy-search', str(faa), str(faa), str(hits), str(tmp_dir),
           '-s', str(s.sensitivity),
           '-e', str(s.evalue),
           '-c', str(s.min_cov),
           '--cov-mode', str(s.cov_mode),
           '--max-seqs', str(s.max_seqs),
           '--threads', str(s.threads),
           '--format-output', _FORMAT_OUTPUT,
           # mmseqs at its default verbosity emits several hundred lines per
           # run and buries the pipeline's own logging. Errors and warnings
           # still surface at -v 1; use --verbose for the full trace.
           '-v', '3' if verbose else '1']

    log.info("Running mmseqs easy-search (all-vs-all) …")
    log.debug(' '.join(cmd))
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        _die(f"mmseqs easy-search failed (exit {proc.returncode}). "
             f"Re-run with --verbose for the full mmseqs trace.")

    shutil.rmtree(tmp_dir, ignore_errors=True)
    log.info(f"Alignment table → {hits} "
             f"({hits.stat().st_size / 1e6:.1f} MB)")


def rescue_masked(cohort: Cohort, out_dir: Path, missing: list[str]) -> Path:
    """Re-search sequences that the prefilter masked out, with masking off.

    mmseqs masks low-complexity regions with tantan in the prefilter
    (`--mask 1`, the default). A sequence that is low-complexity along
    essentially its whole length keeps no k-mers to seed on as a *database*
    entry, so nothing finds it — including itself. It can still act as a
    query, so it appears in hits.tsv with real, high-scoring homologues while
    having no self-hit row.

    Observed on human MANE: a ~33%-Cys keratin-associated protein with
    unambiguous paralogues (e-values ~1e-71, full coverage both ways) and no
    self-hit.

    Two things are lost, not one. The missing self-hit costs the normalised
    bitscore its denominator. But if *both* sides of a low-complexity family
    are masked, the edges between them disappear from the search entirely —
    verified on a synthetic Cys/Gly/Ser-rich pair, which produced no rows at
    all. Recovering only the self-bitscore would leave those paralogues as
    separate singletons: silent under-merging, the unsafe direction.

    So the affected sequences are re-searched as queries against the *full*
    database with `--mask 0`, recovering their self-bitscores and their edges
    in one pass. Every other search parameter matches the main run, so the
    recovered rows are on the same footing and face the same level filters.

    Rerunning the entire search unmasked is not an option — unmasked
    all-vs-all is exactly what chains Cys-rich, Gly-rich and Ser-rich regions
    into a single component. Restricting the unmasked queries to the handful
    that need it keeps that blast radius small, and any spurious low-complexity
    match still has to clear the nbs and coverage floors, which a partial
    match against a longer protein will not.

    Falling back to `self_bits[target]` instead would be unsafe: where the
    query is the longer sequence, max() would have selected `self_bits[query]`,
    so substituting the smaller denominator inflates nbs and admits edges that
    should have been rejected.
    """
    mmseqs = _require_binary(cohort.mmseqs_bin, 'mmseqs')
    faa = out_dir / 'proteins.faa'
    sub_faa = out_dir / 'tmp_masked_queries.faa'
    rescued = out_dir / 'hits_rescued.tsv'
    tmp_dir = out_dir / 'mmseqs_rescue_tmp'
    s = cohort.search

    want = set(missing)
    with open(sub_faa, 'w') as fh:
        for rid, seq in iter_fasta(faa):
            if rid in want:
                fh.write(f">{rid}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + '\n')

    cmd = [mmseqs, 'easy-search', str(sub_faa), str(faa),
           str(rescued), str(tmp_dir),
           '--mask', '0',                       # the whole point of this pass
           '-s', str(s.sensitivity),
           '-e', str(s.evalue),
           '-c', str(s.min_cov),
           '--cov-mode', str(s.cov_mode),
           '--max-seqs', str(s.max_seqs),
           '--threads', str(s.threads),
           '--format-output', _FORMAT_OUTPUT,
           '-v', '1']
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        _die(f"masked-sequence rescue search failed (exit {proc.returncode}):\n"
             f"{proc.stderr}")

    shutil.rmtree(tmp_dir, ignore_errors=True)
    sub_faa.unlink(missing_ok=True)
    return rescued


# ---------------------------------------------------------------------------
# Stage 3 — cluster
# ---------------------------------------------------------------------------

def stage_cluster(cohort: Cohort, out_dir: Path, force: bool) -> None:
    hits = out_dir / 'hits.tsv'
    ptsv = out_dir / 'proteins.tsv'
    family_tsv = out_dir / 'family.tsv'

    deps = [hits, ptsv, cohort_yaml_path(cohort)]
    if cohort.locus_overlap_enabled:
        deps += [_member_canonical_gff(cohort, m) for m in cohort.members]
    if not force and _is_current(family_tsv, deps):
        log.info("Clustering current — skipping.")
        return

    with open(ptsv) as fh:
        proteins = list(csv.DictReader(fh, delimiter='\t'))
    ids = [p['protein_id'] for p in proteins]
    index = {pid: i for i, pid in enumerate(ids)}
    n = len(ids)

    def harvest_self_bits(path: Path, into: dict[str, float]) -> None:
        """Self-bits must be read BEFORE self-hits are dropped. Order-critical:
        the normalised bitscore denominator comes from exactly these rows."""
        with open(path) as fh:
            for row in csv.reader(fh, delimiter='\t'):
                if row[H_QUERY] == row[H_TARGET]:
                    into[row[H_QUERY]] = float(row[H_BITS])

    log.info("Pass 1/2: harvesting self-bitscores …")
    hit_files = [hits]
    self_bits: dict[str, float] = {}
    harvest_self_bits(hits, self_bits)

    # Unsearched proteins are absent from proteins.faa by construction, so
    # they have no self-hit and never appear in any alignment row. They are
    # singletons by definition; only searched proteins need a denominator.
    searched = [p['protein_id'] for p in proteins if p.get('searched') == 'true']
    absent = [pid for pid in searched if pid not in self_bits]
    if absent:
        # tantan masked these out of the target database. They lose both their
        # denominator and, if their homologues were masked too, their edges.
        log.info(f"{len(absent)} protein(s) have no self-hit (low-complexity "
                 f"masking); re-searching them with --mask 0 …")
        rescued = rescue_masked(cohort, out_dir, absent)
        hit_files.append(rescued)
        harvest_self_bits(rescued, self_bits)

        still = [pid for pid in absent if pid not in self_bits]
        if still:
            shown = '\n  '.join(still[:10])
            _die(f"{len(still)} protein(s) still have no self-hit after the "
                 f"--mask 0 rescue pass, so their normalised bitscore is "
                 f"undefined. First few:\n  {shown}\n"
                 f"Exclude these sequences from the cohort; do not let them "
                 f"through as silent singletons.")

        n_rows = sum(1 for _ in open(rescued))
        log.info(f"Recovered {len(absent)} self-bitscore(s) and {n_rows} "
                 f"alignment row(s) → {rescued.name}")

    # Pass 2 — stream edges once, testing every level per row. Avoids holding
    # the (potentially multi-million-row) edge set in memory.
    log.info(f"Pass 2/2: building components for {len(cohort.levels)} level(s) …")
    ufs = {lv.name: UnionFind(n) for lv in cohort.levels}
    # Undirected edges, deduplicated by encoded index pair. A masked target
    # can be reported in one direction only, so counting rows where q < t
    # would silently drop those edges from the QC total. One int per accepted
    # edge; the alignment rows themselves are never held in memory.
    edge_sets: dict[str, set[int]] = {lv.name: set() for lv in cohort.levels}

    for hit_file in hit_files:
        with open(hit_file) as fh:
            for row in csv.reader(fh, delimiter='\t'):
                q, t = row[H_QUERY], row[H_TARGET]
                if q == t:
                    continue
                iq, it = index.get(q), index.get(t)
                if iq is None or it is None:
                    continue
                bits = float(row[H_BITS])
                nbs = normalised_bitscore(bits, self_bits[q], self_bits[t])
                qcov, tcov = float(row[H_QCOV]), float(row[H_TCOV])
                evalue, fident = float(row[H_EVALUE]), float(row[H_FIDENT])
                key = iq * n + it if iq < it else it * n + iq
                for lv in cohort.levels:
                    if lv.accepts(nbs, qcov, tcov, evalue, fident):
                        ufs[lv.name].union(iq, it)
                        edge_sets[lv.name].add(key)

    # Third edge source: exonic overlap in the genome. Applied after the
    # sequence edges so the merge counts report what locus overlap added
    # beyond homology.
    locus_merges = {lv.name: 0 for lv in cohort.levels}
    overlap_rows: list[dict] = []
    if cohort.locus_overlap_enabled:
        gene_index = {(p['species'], p['gene_id']): i
                      for i, p in enumerate(proteins)}
        for member in cohort.members:
            gff = _member_canonical_gff(cohort, member)
            if not gff.exists():
                log.warning(f"No canonical.gff for {member['dataset']} — "
                            f"skipping locus overlap for it.")
                continue
            species = member['species']
            pairs = find_overlapping_gene_pairs(
                gene_exons_from_gff(gff), cohort.min_overlap_bp)
            for ga, gb, chrom, bp, sa, sb in pairs:
                ia = gene_index.get((species, ga))
                ib = gene_index.get((species, gb))
                if ia is None or ib is None:
                    continue          # gene not in the cohort (filtered upstream)
                merged_at = []
                for lv in cohort.levels:
                    if ufs[lv.name].union(ia, ib):
                        locus_merges[lv.name] += 1
                        merged_at.append(lv.name)
                    edge_sets[lv.name].add(
                        ia * n + ib if ia < ib else ib * n + ia)
                overlap_rows.append({
                    'species': species, 'gene_a': ga, 'gene_b': gb,
                    'chrom': chrom, 'overlap_bp': bp,
                    'strand_a': sa, 'strand_b': sb,
                    'same_strand': 'true' if sa == sb else 'false',
                    'merged_levels': ','.join(merged_at) or 'none',
                })
        if overlap_rows:
            _write_tsv(out_dir / 'locus_overlap_pairs.tsv',
                       list(overlap_rows[0].keys()), overlap_rows)
            merges = ', '.join(f"{k}:{v}" for k, v in locus_merges.items())
            log.info(f"Locus overlap: {len(overlap_rows)} gene pair(s) sharing "
                     f">= {cohort.min_overlap_bp} bp of exonic sequence; "
                     f"families merged — {merges}")
        else:
            log.info("Locus overlap: no qualifying gene pairs.")

    family_of: dict[str, list[str]] = {}
    qc_rows: list[dict] = []
    member_rows: list[dict] = []

    for lv in cohort.levels:
        labels, groups = canonical_families(ufs[lv.name], ids, lv.name)
        family_of[lv.name] = labels

        sizes = [len(g) for g in groups]
        n_singletons = sum(1 for s in sizes if s == 1)
        multi = [s for s in sizes if s > 1]
        max_size = max(sizes)
        largest = groups[0]                      # sorted size-desc by construction
        qc_rows.append({
            'level': lv.name,
            'n_edges': len(edge_sets[lv.name]),
            'n_locus_overlap_merges': locus_merges[lv.name],
            'n_genes': n,
            'n_families': len(groups),
            'n_singletons': n_singletons,
            'pct_genes_in_multigene_families':
                f"{100.0 * (n - n_singletons) / n:.2f}",
            'median_family_size':
                f"{statistics.median(multi):.1f}" if multi else 'NA',
            'max_family_size': max_size,
            'max_family_pct': f"{100.0 * max_size / n:.3f}",
            'n_species_in_largest':
                len({split_protein_id(ids[i])[0] for i in largest}),
        })

        for fam_rank, members in enumerate(groups, start=1):
            label = labels[members[0]]
            per_species = Counter(split_protein_id(ids[i])[0] for i in members)
            member_rows.append({
                'level': lv.name,
                'family_id': label,
                'n_members': len(members),
                'n_species': len(per_species),
                'members_per_species': ';'.join(
                    f'{sp}:{c}' for sp, c in sorted(per_species.items())),
                'member_gene_ids': ';'.join(
                    sorted(split_protein_id(ids[i])[1] for i in members)),
            })

    # family.tsv — the seam. One row per gene, every gene present.
    level_names = [lv.name for lv in cohort.levels]
    cols = ['species', 'gene_id', 'transcript_id', 'dataset']
    for name in level_names:
        cols += [f'family_id_{name}', f'family_size_{name}']
    cols += ['protein_len', 'had_internal_stop', 'searched']

    size_of: dict[str, Counter] = {
        name: Counter(family_of[name]) for name in level_names}

    rows = []
    for i, p in enumerate(proteins):
        row = {
            'species': p['species'],
            'gene_id': p['gene_id'],
            'transcript_id': p['transcript_id'],
            'dataset': p['dataset'],
            'protein_len': p['protein_len'],
            'had_internal_stop': p['had_internal_stop'],
            'searched': p.get('searched', 'true'),
        }
        for name in level_names:
            fam = family_of[name][i]
            row[f'family_id_{name}'] = fam
            row[f'family_size_{name}'] = size_of[name][fam]
        rows.append(row)

    if len(rows) != n:
        _die(f"Internal error: {len(rows)} output rows for {n} proteins.")

    _write_tsv(family_tsv, cols, rows)
    _write_tsv(out_dir / 'family_qc.tsv', list(qc_rows[0].keys()), qc_rows)
    _write_tsv(out_dir / 'family_members.tsv',
               list(member_rows[0].keys()), member_rows)

    _report_selection(cohort, qc_rows)


def _report_selection(cohort: Cohort, qc_rows: list[dict]) -> None:
    """Log the size distribution and apply the selection rule.

    Rule: take the loosest level whose largest family stays under
    max_family_pct of the corpus. Loosest because over-merging costs power
    while under-merging inflates the performance estimate; bounded because a
    block larger than 1/k of the corpus cannot be packed into balanced folds.
    """
    log.info("Family size distribution (loosest first):")
    for r in qc_rows:
        log.info(f"  {r['level']:>8}  families={r['n_families']:>7}  "
                 f"singletons={r['n_singletons']:>7}  "
                 f"max={r['max_family_size']:>6} ({r['max_family_pct']}%)  "
                 f"species_in_largest={r['n_species_in_largest']}")

    ceiling = 100.0 * cohort.max_family_pct
    if cohort.ceiling_is_explicit:
        basis = "explicit selection.max_family_pct"
    else:
        basis = (f"{cohort.max_fold_fraction:g} of a fold at k={cohort.k_folds}")

    ok = [r for r in qc_rows if float(r['max_family_pct']) < ceiling]
    if ok:
        log.info(f"Recommended blocking level: '{ok[0]['level']}' "
                 f"(loosest with max_family_pct < {ceiling:g}%, {basis}).")
        runners_up = [r for r in qc_rows if r is not ok[0]
                      and float(r['max_family_pct']) >= ceiling]
        near = [r for r in runners_up
                if float(r['max_family_pct']) < ceiling * 1.5]
        if near:
            # A level that misses by a hair should not be dismissed silently:
            # the ceiling is a packing heuristic, not a hard boundary, and
            # over-merging is the safe direction.
            names = ', '.join(f"{r['level']} ({r['max_family_pct']}%)"
                              for r in near)
            log.warning(
                f"Level(s) just over the ceiling: {names}. The ceiling is a "
                f"packing heuristic, not a hard limit — a looser level whose "
                f"largest family is still well inside one fold "
                f"({100.0 / cohort.k_folds:.1f}% of the corpus) is usually the "
                f"better choice. Check family_members.tsv before deciding.")
    else:
        log.warning(
            f"No level keeps its largest family under {ceiling:g}% of the "
            f"corpus — every level has chained. Inspect family_members.tsv for "
            f"the largest family, then tighten min_nbs / min_cov.")


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def cohort_yaml_path(cohort: Cohort) -> Path:
    return cohort.project_root / 'configs' / 'cohorts' / f'{cohort.name}.yaml'


def write_params(cohort: Cohort, out_dir: Path) -> None:
    import yaml
    s = cohort.search
    payload = {
        'cohort': cohort.name,
        'search_hash': cohort.hash,
        'members': cohort.members,
        'search': {
            'sensitivity': s.sensitivity, 'max_seqs': s.max_seqs,
            'evalue': s.evalue, 'min_cov': s.min_cov,
            'cov_mode': s.cov_mode, 'threads': s.threads,
            'format_output': _FORMAT_OUTPUT,
        },
        'levels': {lv.name: {'min_nbs': lv.min_nbs, 'min_cov': lv.min_cov,
                             'max_evalue': lv.max_evalue,
                             'min_fident': lv.min_fident}
                   for lv in cohort.levels},
        'selection': {'max_family_pct': cohort.max_family_pct,
                      'k_folds': cohort.k_folds,
                      'max_fold_fraction': cohort.max_fold_fraction},
        'locus_overlap': {'enabled': cohort.locus_overlap_enabled,
                          'min_overlap_bp': cohort.min_overlap_bp},
        'translate': {'min_protein_len': cohort.min_protein_len,
                      'max_internal_stop_frac': cohort.max_internal_stop_frac},
        'tool_versions': {
            'seqkit': _tool_version(shutil.which(cohort.seqkit_bin) or cohort.seqkit_bin),
            'mmseqs': _tool_version(shutil.which(cohort.mmseqs_bin) or cohort.mmseqs_bin),
        },
    }
    with open(out_dir / 'params.yaml', 'w') as fh:
        yaml.dump(payload, fh, default_flow_style=False, sort_keys=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _list_cohorts(project_root: Path) -> list[str]:
    d = project_root / 'configs' / 'cohorts'
    if not d.is_dir():
        return []
    return sorted(p.stem for p in d.glob('*.yaml'))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cohort', '-c', default=os.environ.get('COHORT'),
                        help="Cohort name (configs/cohorts/<name>.yaml)")
    parser.add_argument('--force', action='store_true',
                        help="Redo every stage, including the search.")
    parser.add_argument('--recluster', action='store_true',
                        help="Redo clustering only, reusing hits.tsv. Use after "
                             "editing 'levels:'.")
    parser.add_argument('--list-cohorts', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.list_cohorts:
        cohorts = _list_cohorts(_PROJECT_ROOT)
        if not cohorts:
            print("(no cohorts found under configs/cohorts/)", file=sys.stderr)
            return
        print("Available cohorts:")
        for c in cohorts:
            print(f"  {c}")
        return

    if not args.cohort:
        _die("--cohort (or COHORT env var) required")

    yaml_path = _PROJECT_ROOT / 'configs' / 'cohorts' / f'{args.cohort}.yaml'
    if not yaml_path.is_file():
        _die(f"Cohort config not found: {yaml_path}")

    cohort = Cohort(args.cohort, _load_yaml(yaml_path), _PROJECT_ROOT)

    runs_root = Path(os.environ.get('RUNS_ROOT', _PROJECT_ROOT / 'runs'))
    out_dir = runs_root / '_cohorts' / cohort.name / 'family' / cohort.hash
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"Cohort '{cohort.name}': {len(cohort.members)} dataset(s), "
             f"{len(cohort.levels)} level(s), search {cohort.hash}")
    log.info(f"Output: {out_dir}")

    stage_translate(cohort, out_dir, force=args.force)
    stage_search(cohort, out_dir, force=args.force, verbose=args.verbose)
    stage_cluster(cohort, out_dir, force=args.force or args.recluster)
    write_params(cohort, out_dir)

    log.info(f"Done. Blocking table: {out_dir / 'family.tsv'}")


if __name__ == '__main__':
    main()
