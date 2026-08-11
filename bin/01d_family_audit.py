#!/usr/bin/env python3
"""
01d_family_audit.py
Audit the nucleotide-level similarity that protein-based family blocking misses.

01c_family.py groups genes by protein similarity. That captures CDS-level
relatedness well, but says nothing directly about UTRs — and UTR-derived
features (uORFs, 5'/3' structure, 3'UTR composition) carry much of the signal
in a stability model.

This measures the residual. For every gene and region it finds the most
similar gene in a DIFFERENT family:

    cross-family similarity = a leakage candidate.
        Genes in different families may land in different splits, so high
        nucleotide identity between them is similarity the blocking does not
        prevent.

    within-family similarity = the control.
        Already blocked by construction. Reported alongside so the
        cross-family numbers have a reference scale.

Split-agnostic by design: it audits the family partition itself, not a
particular train/test assignment, so the result holds for k-fold and holdout
alike.

Reading the result
------------------
If cross-family identity has no meaningful tail above ~0.70-0.80 at real
coverage, protein-based blocking is sufficient and the CDS-only basis is
defensible — with a figure to show for it. If there is a tail, `audit_pairs.tsv`
names exactly which gene pairs are responsible.

Note that shared-repeat matches (an Alu in two unrelated 3'UTRs) are largely
excluded by the coverage floor: 300 bp of Alu inside a 2 kb UTR cannot reach
50% coverage. Short UTRs where the repeat IS most of the sequence will still
appear, and are a judgement call — that is signal as much as leakage.

Usage:
  ./bin/01d_family_audit.py --cohort human_only --level loose
  ./bin/01d_family_audit.py -c human_only -l loose --regions 3UTR,5UTR,CDS
  ./bin/01d_family_audit.py -c human_only -l loose --min-cov 0.3 --force
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

_THIS = Path(__file__).resolve()
_PROJECT_ROOT = _THIS.parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from lib.family import iter_fasta, make_protein_id  # noqa: E402
from lib.gff import normalise_id, split_composite_fasta_id  # noqa: E402
from lib.paths import resolve_paths  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%H:%M:%S'
)
log = logging.getLogger('audit')

# hits column order — must match --format-output below.
(H_QUERY, H_TARGET, H_FIDENT, H_EVALUE, H_QCOV, H_TCOV, H_QLEN, H_TLEN,
 H_QSTART, H_QEND) = range(10)
_FORMAT_OUTPUT = 'query,target,fident,evalue,qcov,tcov,qlen,tlen,qstart,qend'

# UTR_pair is a 3-line format with a ViennaRNA constraint line, not FASTA.
_DEFAULT_SKIP = ('UTR_pair',)
_DEFAULT_REGIONS = ('3UTR', '5UTR')
_IDENTITY_BINS = (0.70, 0.80, 0.90, 0.95)


def _die(msg: str) -> None:
    log.error(msg)
    sys.exit(1)


def _load_yaml(path: Path) -> dict:
    try:
        import yaml
    except ImportError:
        _die("PyYAML not installed (pip install pyyaml).")
    with open(path) as f:
        return yaml.safe_load(f) or {}


def _write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def _percentile(sorted_vals: list[float], p: float) -> float:
    """Nearest-rank percentile. Empty input -> 0.0."""
    if not sorted_vals:
        return 0.0
    k = max(0, min(len(sorted_vals) - 1, int(round(p / 100.0 * len(sorted_vals) + 0.5)) - 1))
    return sorted_vals[k]


# ---------------------------------------------------------------------------
# Locating the family run
# ---------------------------------------------------------------------------

def find_family_dir(cohort: str, runs_root: Path, search_hash: str | None) -> Path:
    base = runs_root / '_cohorts' / cohort / 'family'
    if not base.is_dir():
        _die(f"No family output for cohort '{cohort}' at {base}. "
             f"Run bin/01c_family.py first.")
    if search_hash:
        d = base / search_hash
        if not d.is_dir():
            _die(f"No such search hash: {d}")
        return d
    candidates = sorted(p for p in base.iterdir()
                        if (p / 'family.tsv').exists())
    if not candidates:
        _die(f"No completed family run under {base} (no family.tsv).")
    if len(candidates) > 1:
        names = ', '.join(p.name for p in candidates)
        _die(f"Multiple search hashes under {base}: {names}. "
             f"Disambiguate with --search-hash.")
    return candidates[0]


# ---------------------------------------------------------------------------
# Region pooling
# ---------------------------------------------------------------------------

def pool_region(members: list[dict], region: str, out_path: Path) -> int:
    """Pool extracted_<region>.fa across members with namespaced IDs.

    Same '<species>|<gene_id>|<transcript_id>' scheme as the protein search, so
    IDs join directly to family.tsv.
    """
    n = 0
    with open(out_path, 'w') as out:
        for m in members:
            dataset, species = m['dataset'], m['species']
            paths = resolve_paths(dataset, _PROJECT_ROOT, species=species)
            fa = paths.extract_dir / f'extracted_{region}.fa'
            if not fa.exists():
                log.warning(f"  {dataset}: no extracted_{region}.fa — skipping")
                continue
            for rec_id, seq in iter_fasta(fa):
                split = split_composite_fasta_id(rec_id, region=region)
                if split is None or not seq:
                    continue
                gene_raw, transcript_id, _ = split
                pid = make_protein_id(species, normalise_id(gene_raw), transcript_id)
                out.write(f">{pid}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + '\n')
                n += 1
    return n


def search_region(mmseqs: str, fa: Path, out_tsv: Path, tmp: Path,
                  threads: int, min_cov: float, verbose: bool) -> None:
    """Nucleotide all-vs-all. --search-type 3 selects nucleotide-nucleotide.

    Deliberately permissive: identity is not filtered here, because the point
    of the audit is to see the identity *distribution*. The coverage floor is
    applied at search time only to keep the hit table bounded, and is set below
    the reporting floor so nothing near the decision boundary is lost.
    """
    cmd = [mmseqs, 'easy-search', str(fa), str(fa), str(out_tsv), str(tmp),
           '--search-type', '3',
           '-s', '7.5',
           '-e', '1e-3',
           '-c', str(max(0.1, min_cov - 0.2)),
           '--cov-mode', '0',
           '--max-seqs', '2000',
           '--threads', str(threads),
           # Both strands, pinned explicitly. An antisense match means the two
           # genes are transcribed from the same DNA, which is a leakage source
           # worth finding. mmseqs already behaves this way by default, but its
           # --strand help reports a default of 1 (forward only) while
           # easy-search actually behaves like 2 — verified empirically. Pinning
           # it keeps the intent visible and survives a future default change.
           '--strand', '2',
           '--format-output', _FORMAT_OUTPUT,
           '-v', '3' if verbose else '1']
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        _die(f"mmseqs nucleotide search failed (exit {proc.returncode}). "
             f"Re-run with --verbose for the full trace.")
    shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def analyse_region(region: str, hits: Path, family_of: dict[str, str],
                   present: set[str], min_cov: float):
    """Best cross-family and within-family match per gene.

    Returns (per_gene_rows, cross_pairs, summary_dict).
    """
    best_cross: dict[str, tuple[float, str, float, float]] = {}
    best_within: dict[str, float] = {}
    cross_pairs: list[dict] = []
    seen_pair: set[tuple[str, str]] = set()

    with open(hits) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            q, t = row[H_QUERY], row[H_TARGET]
            if q == t:
                continue
            fq, ft = family_of.get(q), family_of.get(t)
            if fq is None or ft is None:
                continue
            qcov, tcov = float(row[H_QCOV]), float(row[H_TCOV])
            if qcov < min_cov or tcov < min_cov:
                continue
            fident = float(row[H_FIDENT])

            if fq == ft:
                if fident > best_within.get(q, -1.0):
                    best_within[q] = fident
                continue

            # Cross-family: the leakage candidate.
            prev = best_cross.get(q)
            if prev is None or fident > prev[0]:
                best_cross[q] = (fident, t, qcov, tcov)

            key = (q, t) if q < t else (t, q)
            if key not in seen_pair:
                seen_pair.add(key)
                # mmseqs reports a reverse-complement match with the query
                # coordinates descending. An antisense hit means the pair share
                # the same DNA read in opposite directions — worth seeing as
                # such rather than looking like an unexplained perfect match.
                antisense = int(row[H_QSTART]) > int(row[H_QEND])
                cross_pairs.append({
                    'region': region,
                    'gene_a': key[0], 'gene_b': key[1],
                    'family_a': family_of[key[0]], 'family_b': family_of[key[1]],
                    'identity': f'{fident:.4f}',
                    'qcov': f'{qcov:.4f}', 'tcov': f'{tcov:.4f}',
                    'strand': 'antisense' if antisense else 'sense',
                    'evalue': row[H_EVALUE],
                })

    per_gene = []
    for gid in sorted(present):
        cross = best_cross.get(gid)
        per_gene.append({
            'region': region,
            'gene': gid,
            'family': family_of.get(gid, 'NA'),
            'max_cross_family_identity': f'{cross[0]:.4f}' if cross else '0.0000',
            'cross_family_partner': cross[1] if cross else 'NA',
            'cross_qcov': f'{cross[2]:.4f}' if cross else 'NA',
            'cross_tcov': f'{cross[3]:.4f}' if cross else 'NA',
            'max_within_family_identity':
                f'{best_within[gid]:.4f}' if gid in best_within else '0.0000',
        })

    # Percentiles over ALL genes with this region, no-hit counted as 0. That is
    # the honest denominator: "what fraction of genes have a similar gene in
    # another family", not "how similar are the ones that matched".
    vals = sorted(best_cross[g][0] if g in best_cross else 0.0 for g in present)
    n = len(present)
    summary = {
        'region': region,
        'n_genes': n,
        'n_with_cross_family_hit': len(best_cross),
        'pct_with_cross_family_hit': f'{100.0 * len(best_cross) / n:.2f}' if n else 'NA',
        'p50': f'{_percentile(vals, 50):.4f}',
        'p90': f'{_percentile(vals, 90):.4f}',
        'p95': f'{_percentile(vals, 95):.4f}',
        'p99': f'{_percentile(vals, 99):.4f}',
        'max': f'{vals[-1]:.4f}' if vals else 'NA',
    }
    for thr in _IDENTITY_BINS:
        k = sum(1 for v in vals if v >= thr)
        summary[f'n_ge_{int(thr * 100)}'] = k
        summary[f'pct_ge_{int(thr * 100)}'] = f'{100.0 * k / n:.3f}' if n else 'NA'
    return per_gene, cross_pairs, summary


def report(summaries: list[dict], min_cov: float) -> None:
    log.info(f"Cross-family nucleotide similarity (coverage floor {min_cov:g}, "
             f"both query and target):")
    log.info(f"  {'region':<8} {'genes':>7} {'any hit':>9} "
             f"{'p95':>7} {'p99':>7} {'>=0.8':>8} {'>=0.9':>8}")
    for s in summaries:
        log.info(f"  {s['region']:<8} {s['n_genes']:>7} "
                 f"{s['pct_with_cross_family_hit']:>8}% "
                 f"{s['p95']:>7} {s['p99']:>7} "
                 f"{s['pct_ge_80']:>7}% {s['pct_ge_90']:>7}%")

    worst = max((float(s['pct_ge_80']) for s in summaries
                 if s['pct_ge_80'] != 'NA'), default=0.0)
    if worst < 1.0:
        log.info(f"Verdict: under {worst:.2f}% of genes have a cross-family "
                 f"match at >=80% identity. Protein-based blocking looks "
                 f"sufficient — the CDS-only basis is defensible.")
    else:
        log.warning(
            f"Verdict: {worst:.2f}% of genes have a cross-family match at "
            f">=80% identity. That is nucleotide-level similarity the current "
            f"blocking does NOT prevent. Inspect audit_pairs.tsv; if the pairs "
            f"are genuine, add these edges as a second source into the family "
            f"graph rather than re-labelling anything.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cohort', '-c', default=os.environ.get('COHORT'),
                    help="Cohort name (configs/cohorts/<name>.yaml)")
    ap.add_argument('--level', '-l', required=True,
                    help="Blocking level to audit (the family_id_<level> column).")
    ap.add_argument('--regions', default=','.join(_DEFAULT_REGIONS),
                    help=f"Comma-separated regions. Default: {','.join(_DEFAULT_REGIONS)}")
    ap.add_argument('--min-cov', type=float, default=0.5,
                    help="Coverage floor (both query and target) for a pair to "
                         "count. Default 0.5 — this is what excludes "
                         "shared-repeat matches.")
    ap.add_argument('--report-identity', type=float, default=0.70,
                    help="List cross-family pairs at or above this identity. "
                         "Default 0.70")
    ap.add_argument('--search-hash', default=None,
                    help="Disambiguate when a cohort has several search runs.")
    ap.add_argument('--force', action='store_true',
                    help="Redo the nucleotide searches.")
    ap.add_argument('-v', '--verbose', action='store_true')
    args = ap.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    if not args.cohort:
        _die("--cohort (or COHORT env var) required")

    mmseqs = shutil.which(os.environ.get('MMSEQS_BIN') or 'mmseqs')
    if not mmseqs:
        _die("mmseqs not found on PATH (or via MMSEQS_BIN).")

    runs_root = Path(os.environ.get('RUNS_ROOT', _PROJECT_ROOT / 'runs'))
    family_dir = find_family_dir(args.cohort, runs_root, args.search_hash)
    family_tsv = family_dir / 'family.tsv'
    params = _load_yaml(family_dir / 'params.yaml')
    members = params.get('members') or []
    if not members:
        _die(f"No members recorded in {family_dir / 'params.yaml'}")
    threads = int((params.get('search') or {}).get('threads', os.cpu_count() or 4))

    # Family map, keyed the same way region FASTA IDs will be.
    col = f'family_id_{args.level}'
    family_of: dict[str, str] = {}
    with open(family_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if col not in (reader.fieldnames or []):
            avail = [c for c in (reader.fieldnames or []) if c.startswith('family_id_')]
            _die(f"No column '{col}' in {family_tsv}. Available: {', '.join(avail)}")
        for r in reader:
            pid = make_protein_id(r['species'], r['gene_id'], r['transcript_id'])
            family_of[pid] = r[col]

    log.info(f"Cohort '{args.cohort}' level '{args.level}': "
             f"{len(family_of)} genes, "
             f"{len(set(family_of.values()))} families")

    out_dir = family_dir / 'audit' / args.level
    out_dir.mkdir(parents=True, exist_ok=True)

    regions = [r.strip() for r in args.regions.split(',') if r.strip()]
    skipped = [r for r in regions if r in _DEFAULT_SKIP]
    if skipped:
        log.warning(f"Skipping non-FASTA region(s): {', '.join(skipped)}")
    regions = [r for r in regions if r not in _DEFAULT_SKIP]
    if not regions:
        _die("No auditable regions requested.")

    all_per_gene: list[dict] = []
    all_pairs: list[dict] = []
    summaries: list[dict] = []

    for region in regions:
        fa = out_dir / f'{region}.fa'
        hits = out_dir / f'{region}_hits.tsv'

        if args.force or not fa.exists():
            log.info(f"[{region}] pooling sequences …")
            n = pool_region(members, region, fa)
            if n == 0:
                log.warning(f"[{region}] no sequences — skipping region")
                continue
            log.info(f"[{region}] {n} sequences")

        if args.force or not hits.exists():
            log.info(f"[{region}] nucleotide all-vs-all …")
            search_region(mmseqs, fa, hits, out_dir / f'tmp_{region}',
                          threads, args.min_cov, args.verbose)

        present = {rid for rid, _ in iter_fasta(fa)} & set(family_of)
        per_gene, pairs, summary = analyse_region(
            region, hits, family_of, present, args.min_cov)
        all_per_gene += per_gene
        all_pairs += [p for p in pairs
                      if float(p['identity']) >= args.report_identity]
        summaries.append(summary)

    if not summaries:
        _die("No regions produced results.")

    _write_tsv(out_dir / 'audit_summary.tsv',
               list(summaries[0].keys()), summaries)
    _write_tsv(out_dir / 'audit_per_gene.tsv',
               list(all_per_gene[0].keys()), all_per_gene)
    all_pairs.sort(key=lambda p: -float(p['identity']))
    if all_pairs:
        _write_tsv(out_dir / 'audit_pairs.tsv',
                   list(all_pairs[0].keys()), all_pairs)
    else:
        (out_dir / 'audit_pairs.tsv').write_text(
            'region\tgene_a\tgene_b\tfamily_a\tfamily_b\t'
            'identity\tqcov\ttcov\tstrand\tevalue\n')

    report(summaries, args.min_cov)
    log.info(f"{len(all_pairs)} cross-family pair(s) at >= "
             f"{args.report_identity:g} identity → {out_dir / 'audit_pairs.tsv'}")
    log.info(f"Audit written to {out_dir}")


if __name__ == '__main__':
    main()
