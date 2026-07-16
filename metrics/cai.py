"""
cai — Codon Adaptation Index per transcript (Sharp & Li 1987).

Reads `extracted_CDS.fa` and scores each CDS as the geometric mean of
per-codon adaptation weights w_c. Stop codons, non-ACGT codons, and sense
codons with w == 0 in the reference are excluded from the mean.

Weight table resolution (first match wins):
  1. metric_config['weights_path']         — explicit path in dataset YAML
  2. paths.references_root / 'cai_weights.tsv'  — derived from species:

If the resolved file does not exist it is built on-the-fly from the CoCoPUTs
TSV pointed to by metric_config['cocoputs_tsv'] (see _cai_weights.py).

Dataset YAML (metrics.cai block):
  cocoputs_tsv : path to a CoCoPUTs *_CDS.tsv (relative to project root or
                 absolute). Required only when the weight table doesn't exist.
  weights_path : explicit weight-table path override (optional).
  min_codons   : minimum # Codons filter applied when building from CoCoPUTs
                 (default 100).

Output — wide, one row per manifest transcript:
  transcript_id  manifest ID
  cai            geometric mean of w_c over scored codons; NA if CDS absent
  n_codons_cai   count of codons included in the score; NA if CDS absent
"""
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Iterable, Optional

from lib.gff import (
    normalise_id,
    read_manifest_transcripts,
    split_composite_fasta_id,
    strip_prefix,
)
from metrics._cai_weights import build_weights, load_weights

log = logging.getLogger('metrics.cai')

OUTPUT_FILENAME = 'cai.tsv'
_NA = 'NA'
_REGION = 'CDS'
_DEFAULT_MIN_CODONS = 100


# ---------------------------------------------------------------------------
# Plugin contract
# ---------------------------------------------------------------------------

def get_input_paths(paths, metric_config) -> Iterable[Path]:
    """Return paths needed for staleness checking and pre-flight validation.

    If the weight table already exists it is listed as an input (so the
    output is considered stale if the weights are updated). If not, the
    CoCoPUTs source file is listed instead (must be present to build).
    """
    inputs = [
        paths.extract_dir / 'extracted_CDS.fa',
        paths.manifest_tsv,
    ]
    weights_path = _resolve_weights_path(paths, metric_config)
    if weights_path is not None and weights_path.exists():
        inputs.append(weights_path)
    else:
        cocoputs = _cocoputs_path(paths, metric_config)
        if cocoputs is not None:
            inputs.append(cocoputs)
    return inputs


def compute(paths, metric_config, output_path: Path) -> None:
    weights_path = _resolve_weights_path(paths, metric_config)
    if weights_path is None:
        raise ValueError(
            "Cannot resolve CAI weights path. Either set "
            "metrics.cai.weights_path in the dataset YAML, or set "
            "'species:' so that paths.references_root is populated."
        )

    if not weights_path.exists():
        cocoputs = _cocoputs_path(paths, metric_config)
        if cocoputs is None:
            raise ValueError(
                f"CAI weights not found at {weights_path} and no "
                "'cocoputs_tsv' is configured under metrics.cai to build it."
            )
        min_codons = int(metric_config.get('min_codons', _DEFAULT_MIN_CODONS))
        log.info("Building CAI weights from %s → %s", cocoputs, weights_path)
        build_weights(cocoputs, weights_path, min_codons=min_codons)

    log.info("Loading CAI weights from %s", weights_path)
    weights = load_weights(weights_path)
    log.info("Loaded weights for %d sense codons.", len(weights))

    cds_fa = paths.extract_dir / 'extracted_CDS.fa'
    transcript_order = read_manifest_transcripts(paths.manifest_tsv)
    manifest_index = _build_manifest_index(transcript_order)
    log.info("Manifest: %d transcripts.", len(transcript_order))

    rows: dict[str, dict] = {}
    n_records = n_unmatched = 0
    sample_unmatched: list[str] = []

    for record_id, seq in _iter_fasta(cds_fa):
        n_records += 1
        tid = _resolve_tid(record_id, manifest_index)
        if tid is None:
            n_unmatched += 1
            if len(sample_unmatched) < 5:
                sample_unmatched.append(record_id)
            continue
        if tid in rows:
            log.warning(
                "Duplicate match for manifest transcript %s (record %s); "
                "keeping first.",
                tid, record_id,
            )
            continue
        cai_score, n_codons = _score_cds(seq, weights, tid)
        rows[tid] = {'cai': cai_score, 'n_codons_cai': n_codons}

    if n_unmatched:
        log.warning(
            "%d/%d FASTA records did not match manifest; sample: %s",
            n_unmatched, n_records, sample_unmatched,
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as out:
        out.write('transcript_id\tcai\tn_codons_cai\n')
        for tid in transcript_order:
            row = rows.get(tid)
            if row is None:
                out.write(f'{tid}\t{_NA}\t{_NA}\n')
            else:
                cai = row['cai']
                nc = row['n_codons_cai']
                cai_str = f'{cai:.6f}' if cai is not None else _NA
                nc_str = str(nc) if nc is not None else _NA
                out.write(f'{tid}\t{cai_str}\t{nc_str}\n')

    n_with_data = sum(1 for tid in transcript_order if tid in rows)
    log.info(
        "Wrote %d rows to %s (%d with data, %d NA).",
        len(transcript_order), output_path.name,
        n_with_data, len(transcript_order) - n_with_data,
    )


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def _resolve_weights_path(paths, metric_config) -> Optional[Path]:
    explicit = metric_config.get('weights_path')
    if explicit:
        return _to_abs(explicit, paths.project_root)
    if paths.references_root is not None:
        return paths.references_root / 'cai_weights.tsv'
    return None


def _cocoputs_path(paths, metric_config) -> Optional[Path]:
    val = metric_config.get('cocoputs_tsv')
    if val is None:
        return None
    return _to_abs(val, paths.project_root)


def _to_abs(val: str, project_root: Path) -> Path:
    p = Path(val)
    return p if p.is_absolute() else (project_root / p).resolve()


# ---------------------------------------------------------------------------
# Manifest indexing — mirrors codon_aa_counts.py
# ---------------------------------------------------------------------------

def _build_manifest_index(transcript_ids: list[str]) -> dict[str, str]:
    index: dict[str, str] = {}
    for tid in transcript_ids:
        for key in (tid, strip_prefix(tid), normalise_id(tid)):
            if key:
                index.setdefault(key, tid)
    return index


def _resolve_tid(record_id: str, manifest_index: dict[str, str]) -> Optional[str]:
    parsed = split_composite_fasta_id(record_id, region=_REGION)
    candidates = []
    if parsed is not None:
        candidates.append(parsed[1])
    candidates.append(record_id)
    for cand in candidates:
        for key in (cand, strip_prefix(cand), normalise_id(cand)):
            if key and key in manifest_index:
                return manifest_index[key]
    return None


# ---------------------------------------------------------------------------
# FASTA reading — minimal, no biopython dependency
# ---------------------------------------------------------------------------

def _iter_fasta(fa_path: Path):
    header: Optional[str] = None
    chunks: list[str] = []
    with open(fa_path) as fh:
        for raw in fh:
            line = raw.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks)
                header = line[1:].split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)


# ---------------------------------------------------------------------------
# CAI scoring
# ---------------------------------------------------------------------------

def _score_cds(
    seq: str,
    weights: dict[str, float],
    transcript_id: str,
) -> tuple[Optional[float], Optional[int]]:
    """Return (cai, n_codons_cai) for one CDS.

    Included in the geometric mean: sense codons present in the weight table
    with w > 0. Excluded: stop codons, non-ACGT codons, sense codons with
    w == 0 (zero reference-set usage — logged once at weight-build time).
    Returns (None, None) when no codons are scoreable.
    """
    seq = seq.upper().replace('U', 'T')
    n = len(seq)
    trimmed = n - (n % 3)
    if n != trimmed:
        log.warning(
            "CDS length %d for %s is not a multiple of 3; dropping last %d nt.",
            n, transcript_id, n - trimmed,
        )

    log_w_sum = 0.0
    n_scored = 0
    for i in range(0, trimmed, 3):
        codon = seq[i:i + 3]
        w = weights.get(codon)
        if w is None or w <= 0:
            continue   # stop, non-ACGT, or zero-weight codon
        log_w_sum += math.log(w)
        n_scored += 1

    if n_scored == 0:
        return None, None

    return math.exp(log_w_sum / n_scored), n_scored
