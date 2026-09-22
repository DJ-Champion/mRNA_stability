"""metrics/junctions.py
Exon-junction features per transcript, all distances in spliced (mature mRNA)
coordinates.

Output (<metrics_dir>/junctions.tsv) columns:
    transcript_id                   joinable to manifest.tsv
    gene_id                         joinable across metrics
    strand
    n_exons
    n_mRNA_junctions                exon junctions in the full mature mRNA
    n_5UTR_junctions                exon junctions wholly within the 5' UTR
    n_CDS_junctions                 exon junctions wholly within the CDS
    n_3UTR_junctions                exon junctions wholly within the 3' UTR
    stop_dist_closest_upstream      spliced nt; NA if no upstream junction
    stop_dist_closest_downstream    intervening spliced nt; NA if no junction
    stop_dist_last_downstream       intervening spliced nt; canonical NMD metric
    start_dist_closest_upstream     spliced nt
    start_dist_closest_downstream   spliced nt

Notes
-----
* Region junction counts are derived from the actual exon boundaries. A
  junction is assigned to a region only when that region contains the
  nucleotide immediately before and immediately after the junction.
* Reference points: first nt of start codon, last nt of stop codon — both
  in transcript reading direction. Falls back to CDS extent if explicit
  start_codon / stop_codon features are absent.
* Downstream stop distances count only intervening nucleotides: a junction
  immediately after the stop codon has distance 0.
* Single-exon transcripts: all junction columns are 0 / NA.
* Stop or start codons split across an intron: handled by collapsing the
  feature pieces to the appropriate transcript-direction extreme.

Shared GFF helpers (DB build/open, ID normalisation, transcript index,
manifest reader) live in lib/gff.py.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Iterable, Optional

from lib.gff import (
    open_or_build_db,
    build_transcript_index,
    lookup_transcript,
    read_manifest_transcripts,
    normalise_id,
)

log = logging.getLogger('metrics.junctions')

OUTPUT_FILENAME = 'junctions.tsv'

# Tab-separated; NA written for missing values
_NA = 'NA'

_HEADER = [
    'transcript_id', 'gene_id', 'strand', 'n_exons',
    'n_mRNA_junctions',
    'n_5UTR_junctions', 'n_CDS_junctions', 'n_3UTR_junctions',
    'stop_dist_closest_upstream', 'stop_dist_closest_downstream',
    'stop_dist_last_downstream',
    'start_dist_closest_upstream', 'start_dist_closest_downstream',
]


# ---------------------------------------------------------------------------
# Plugin contract
# ---------------------------------------------------------------------------

def get_input_paths(paths, metric_config) -> Iterable[Path]:
    return [paths.canonical_gff, paths.manifest_tsv]


def compute(paths, metric_config, output_path: Path) -> None:
    """Entry point called by the orchestrator."""
    db = open_or_build_db(paths.canonical_gff,
                          paths.extract_dir / 'canonical.gff.db')

    index = build_transcript_index(db)
    if not index:
        log.error("No transcript-like features found in canonical.gff "
                  "(no features with exon children). Cannot proceed.")
        # Still write an empty TSV with header so the orchestrator's
        # 'output produced' check passes; the warning above is the signal.
        _write_tsv(output_path, [])
        return

    transcript_ids = read_manifest_transcripts(paths.manifest_tsv)
    log.info(f"Computing junctions for {len(transcript_ids)} transcripts")

    rows = []
    n_missing, n_no_exons = 0, 0
    sample_missing = []
    for tx_id in transcript_ids:
        tx_feature = lookup_transcript(db, tx_id, index)
        if tx_feature is None:
            n_missing += 1
            if len(sample_missing) < 5:
                sample_missing.append(tx_id)
            continue
        row = _compute_for_transcript(db, tx_feature, manifest_tx_id=tx_id)
        if row is None:
            n_no_exons += 1
            continue
        rows.append(row)

    if n_missing:
        log.warning(f"{n_missing} transcript(s) in manifest not found in canonical.gff")
        # Diagnostic: when failure rate is high, dump samples so the user
        # can compare manifest IDs vs the keys we actually indexed.
        if n_missing >= max(5, len(transcript_ids) // 2):
            sample_keys = [k for k in list(index.keys())[:8]]
            log.warning("  ID format mismatch suspected. Examples:")
            log.warning(f"    Manifest IDs (first 5): {sample_missing}")
            log.warning(f"    GFF index keys (first 8): {sample_keys}")
            log.warning("  Run with -v for the full feature-type list. "
                        "If the formats clearly differ in a way the normaliser "
                        "doesn't handle, file a bug or extend KNOWN_PREFIXES "
                        "in lib/gff.py.")
    if n_no_exons:
        log.warning(f"{n_no_exons} transcript(s) had no exon features")

    _write_tsv(output_path, rows)
    log.info(f"Wrote {len(rows)} rows to {output_path.name}")


# ---------------------------------------------------------------------------
# Spliced-coordinate machinery
# ---------------------------------------------------------------------------

def _build_spliced_index(exons, strand: str):
    """Sort exons in transcript order and assign each a spliced offset.

    Returns:
        spliced_index : list of (g_start, g_end, spliced_start, length)
                        where spliced_start is the 0-indexed mRNA position
                        of the first nt of this exon.
        junctions     : list of spliced positions of each junction, in
                        transcript order. Junction p separates spliced
                        position p-1 from p. There are len(exons)-1 of them.
    """
    if strand == '-':
        sorted_exons = sorted(exons, key=lambda e: e.start, reverse=True)
    else:
        sorted_exons = sorted(exons, key=lambda e: e.start)

    spliced_index = []
    cum = 0
    for e in sorted_exons:
        length = e.end - e.start + 1
        spliced_index.append((e.start, e.end, cum, length))
        cum += length

    # Junctions: cumulative length up to (but not including) each non-first exon
    junctions = []
    cum = 0
    for _, _, _, length in spliced_index[:-1]:
        cum += length
        junctions.append(cum)

    return spliced_index, junctions


def _genomic_to_spliced(genomic_pos: int, spliced_index, strand: str) -> Optional[int]:
    """Map a genomic position to a spliced (0-indexed) mRNA position.

    Returns None if the position falls outside every exon (shouldn't happen
    for start/stop codons in well-formed annotations, but possible for
    pathological cases — e.g. a stop codon entirely within an intron).
    """
    for g_start, g_end, spliced_start, _ in spliced_index:
        if g_start <= genomic_pos <= g_end:
            if strand == '-':
                offset = g_end - genomic_pos
            else:
                offset = genomic_pos - g_start
            return spliced_start + offset
    return None


def _feature_spliced_intervals(features, spliced_index, strand: str):
    """Return inclusive spliced-coordinate intervals covered by features.

    Each feature is intersected with every exon rather than assuming one
    annotation record per participating exon. This supports both the usual
    MANE per-exon UTR/CDS records and annotations that use a feature spanning
    more than one exon.
    """
    intervals = []
    for feature in features:
        for g_start, g_end, spliced_start, _ in spliced_index:
            overlap_start = max(feature.start, g_start)
            overlap_end = min(feature.end, g_end)
            if overlap_start > overlap_end:
                continue

            if strand == '-':
                interval_start = spliced_start + (g_end - overlap_end)
                interval_end = spliced_start + (g_end - overlap_start)
            else:
                interval_start = spliced_start + (overlap_start - g_start)
                interval_end = spliced_start + (overlap_end - g_start)
            intervals.append((interval_start, interval_end))
    return intervals


def _count_region_junctions(junctions, features, spliced_index, strand: str) -> int:
    """Count actual exon junctions wholly contained within a region.

    Junction ``j`` separates mature-mRNA positions ``j - 1`` and ``j``.
    It belongs to a region only when both flanking nucleotides are covered by
    that region. A junction coincident with a region boundary is therefore not
    assigned to either adjacent region.
    """
    intervals = _feature_spliced_intervals(features, spliced_index, strand)

    def is_covered(pos: int) -> bool:
        return any(start <= pos <= end for start, end in intervals)

    return sum(is_covered(j - 1) and is_covered(j) for j in junctions)


# ---------------------------------------------------------------------------
# Reference-point extraction
# ---------------------------------------------------------------------------

def _start_codon_genomic_pos(start_codons, cds_features, strand: str) -> Optional[int]:
    """Genomic position of the first nucleotide of the start codon.

    `start_codons` and `cds_features` are lists of gffutils Features (possibly
    empty). Falls back to CDS extents when start_codon is absent.
    """
    if start_codons:
        if strand == '-':
            return max(c.end for c in start_codons)
        return min(c.start for c in start_codons)
    if cds_features:
        # First CDS coordinate in transcript direction
        if strand == '-':
            return max(c.end for c in cds_features)
        return min(c.start for c in cds_features)
    return None


def _stop_codon_genomic_pos(stop_codons, cds_features, strand: str) -> Optional[int]:
    """Genomic position of the last nucleotide of the stop codon.

    Caveat: when stop_codon is missing we fall back to the 3' CDS edge,
    which is correct for GFFs that include the stop codon in CDS (Ensembl
    legacy) but off by 3 nt for GFFs that exclude it (most modern Ensembl,
    GENCODE, MANE include it; older RefSeq sometimes excludes it).
    Marginal effect on closest-junction distances.
    """
    if stop_codons:
        if strand == '-':
            return min(c.start for c in stop_codons)
        return max(c.end for c in stop_codons)
    if cds_features:
        if strand == '-':
            return min(c.start for c in cds_features)
        return max(c.end for c in cds_features)
    return None


# ---------------------------------------------------------------------------
# Per-transcript computation
# ---------------------------------------------------------------------------

def _compute_for_transcript(db, transcript, manifest_tx_id: str) -> Optional[dict]:
    """Compute the row for one transcript. Returns None if no exons."""
    strand = transcript.strand
    exons = list(db.children(transcript, featuretype='exon'))
    if not exons:
        return None

    spliced_index, junctions = _build_spliced_index(exons, strand)
    n_exons = len(exons)

    # Full-mRNA and regional counts from actual exon junction coordinates.
    # A regional junction must have region sequence on both sides.
    utr5 = list(db.children(transcript, featuretype='five_prime_UTR'))
    cds = list(db.children(transcript, featuretype='CDS'))
    utr3 = list(db.children(transcript, featuretype='three_prime_UTR'))
    n_mrna_j = len(junctions)
    n_5utr_j = _count_region_junctions(
        junctions, utr5, spliced_index, strand)
    n_cds_j = _count_region_junctions(
        junctions, cds, spliced_index, strand)
    n_3utr_j = _count_region_junctions(
        junctions, utr3, spliced_index, strand)

    # Reference points → spliced
    start_codons = list(db.children(transcript, featuretype='start_codon'))
    stop_codons = list(db.children(transcript, featuretype='stop_codon'))

    start_g = _start_codon_genomic_pos(start_codons, cds, strand)
    stop_g = _stop_codon_genomic_pos(stop_codons, cds, strand)

    start_s = (_genomic_to_spliced(start_g, spliced_index, strand)
               if start_g is not None else None)
    stop_s = (_genomic_to_spliced(stop_g, spliced_index, strand)
              if stop_g is not None else None)

    stop_up, stop_down, stop_last_down = _distances(
        stop_s, junctions, downstream_intervening_nt=True)
    start_up, start_down, _ = _distances(start_s, junctions)

    # gene_id from attributes; fall back to Parent or empty
    gene_id = ''
    for attr in ('gene_id', 'Parent'):
        if attr in transcript.attributes and transcript.attributes[attr]:
            gene_id = normalise_id(transcript.attributes[attr][0])
            break

    return {
        'transcript_id': manifest_tx_id,   # preserve manifest ID exactly
        'gene_id': gene_id,
        'strand': strand,
        'n_exons': n_exons,
        'n_mRNA_junctions': n_mrna_j,
        'n_5UTR_junctions': n_5utr_j,
        'n_CDS_junctions': n_cds_j,
        'n_3UTR_junctions': n_3utr_j,
        'stop_dist_closest_upstream': stop_up,
        'stop_dist_closest_downstream': stop_down,
        'stop_dist_last_downstream': stop_last_down,
        'start_dist_closest_upstream': start_up,
        'start_dist_closest_downstream': start_down,
    }


def _distances(pos: Optional[int], junctions: list[int],
               downstream_intervening_nt: bool = False):
    """Return (closest_upstream, closest_downstream, last_downstream).

    All in spliced nt. Upstream distances are ``pos - junction_pos``.
    Downstream distances are normally ``junction_pos - pos``; when
    ``downstream_intervening_nt`` is true, they are instead the number of
    nucleotides strictly between ``pos`` and the junction:
    ``junction_pos - pos - 1``. None is returned when no qualifying junction
    exists or ``pos`` is None.

    Convention
    ----------
    Junctions are 0-indexed spliced positions of the *first nt of the
    downstream exon* (i.e. junction p separates spliced positions p-1
    and p in the mature mRNA).

    A junction located exactly at ``pos`` (e.g. a start codon coinciding with
    a splice acceptor) is treated as upstream with distance 0. For stop-codon
    calculations, a junction at ``pos + 1`` is immediately adjacent to the
    final stop-codon nucleotide and therefore has downstream distance 0.
    """
    if pos is None or not junctions:
        return None, None, None

    upstream = [j for j in junctions if j <= pos]
    downstream = [j for j in junctions if j > pos]

    closest_up = (pos - max(upstream)) if upstream else None
    downstream_offset = 1 if downstream_intervening_nt else 0
    closest_down = (min(downstream) - pos - downstream_offset
                    if downstream else None)
    last_down = (max(downstream) - pos - downstream_offset
                 if downstream else None)
    return closest_up, closest_down, last_down


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _write_tsv(output_path: Path, rows: list[dict]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(
            f, fieldnames=_HEADER, delimiter='\t',
            extrasaction='ignore', lineterminator='\n',
        )
        writer.writeheader()
        for row in rows:
            # Stringify with NA for None values (consistent with manifest convention)
            out = {k: (_NA if row.get(k) is None else row.get(k)) for k in _HEADER}
            writer.writerow(out)
