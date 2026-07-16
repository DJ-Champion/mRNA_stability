"""
_cai_weights — build and load a CAI codon-weight table from a CoCoPUTs TSV.

CoCoPUTs (https://hive.biochemistry.gwu.edu/review/codon) provides per-CDS
raw codon counts from NCBI RefSeq. This module filters to ribosomal protein
genes (RPL*/RPS*, the classical highly-expressed reference set for CAI),
handles isoform deduplication, and computes RSCU → w (Sharp & Li 1987).

Exported:
    build_weights(cocoputs_tsv, output_path, min_codons=100) → None
    load_weights(weights_path) → dict[str, float]   # codon → w, sense only
"""
from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger('metrics._cai_weights')

# Standard genetic code, T form: codon → amino acid (None = stop codon)
_CODON_TABLE: dict[str, str | None] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": None, "TAG": None,
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": None, "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Synonymous codon families (stop codons excluded)
_AA_TO_CODONS: dict[str, list[str]] = {}
for _c, _a in _CODON_TABLE.items():
    if _a is not None:
        _AA_TO_CODONS.setdefault(_a, []).append(_c)

_ALL_CODONS_ALPHA: tuple[str, ...] = tuple(sorted(_CODON_TABLE.keys()))


def build_weights(
    cocoputs_tsv: Path,
    output_path: Path,
    min_codons: int = 100,
) -> None:
    """Build a CAI weight table from a CoCoPUTs CDS TSV.

    Reference set: ribosomal protein genes (Gene ID matching ^RPL or ^RPS,
    case-insensitive). These are the classical highly-expressed reference genes
    used by Sharp & Li (1987) for human CAI.

    Isoform deduplication: for each Gene ID, keep the Protein ID with the
    highest '# Codons' value. This prevents genes with many predicted isoforms
    from contributing multiple times to the aggregate counts.

    min_codons: entries with fewer than this many codons are excluded. The
    default (100) removes truncated / frameshifted predictions.

    Weight computation (Sharp & Li 1987):
        RSCU_c = count_c * n_syn / Σ(count_j  for j synonymous with c)
        w_c    = RSCU_c / max(RSCU_j  in same synonymous family)
    Single-codon families (Met: ATG, Trp: TGG) have RSCU = 1 and w = 1 by
    definition. If a full synonymous family has zero total usage in the
    reference set, uniform weights (rscu=1.0, w=1.0) are assigned and a
    warning is logged.

    Output TSV columns: codon, rscu, w. Stop codons are written as NA.
    """
    cocoputs_tsv = Path(cocoputs_tsv)
    output_path = Path(output_path)

    gene_best: dict[str, dict] = {}   # gene_id → {protein_id, n_codons, counts}

    with open(cocoputs_tsv, newline='') as fh:
        raw_header = fh.readline().rstrip('\r\n').split('\t')
        col = {name: idx for idx, name in enumerate(raw_header)}

        required_meta = ('Gene ID', 'Protein ID', '# Codons')
        missing_meta = [c for c in required_meta if c not in col]
        if missing_meta:
            raise ValueError(
                f"CoCoPUTs header is missing expected columns: {missing_meta}. "
                f"First columns seen: {raw_header[:6]}"
            )
        missing_codons = [c for c in _CODON_TABLE if c not in col]
        if missing_codons:
            raise ValueError(
                f"CoCoPUTs header is missing codon columns: {missing_codons[:10]}..."
            )

        codon_col = {codon: col[codon] for codon in _CODON_TABLE}
        gene_id_idx = col['Gene ID']
        protein_id_idx = col['Protein ID']
        n_codons_idx = col['# Codons']

        n_total = n_ribo = 0
        for line in fh:
            parts = line.rstrip('\r\n').split('\t')
            if len(parts) < len(raw_header):
                continue
            n_total += 1

            gene_id = parts[gene_id_idx].strip()
            gene_upper = gene_id.upper()
            if not (gene_upper.startswith('RPL') or gene_upper.startswith('RPS')):
                continue
            n_ribo += 1

            try:
                n_codons = int(parts[n_codons_idx])
            except ValueError:
                continue
            if n_codons < min_codons:
                continue

            # Keep only the longest isoform per gene
            if gene_id in gene_best and n_codons <= gene_best[gene_id]['n_codons']:
                continue

            counts: dict[str, int] = {}
            ok = True
            for codon, ci in codon_col.items():
                try:
                    counts[codon] = int(parts[ci])
                except (ValueError, IndexError):
                    ok = False
                    break
            if ok:
                gene_best[gene_id] = {
                    'protein_id': parts[protein_id_idx],
                    'n_codons': n_codons,
                    'counts': counts,
                }

    log.info(
        "%d total CoCoPUTs entries; %d ribosomal (RPL/RPS); "
        "%d genes after isoform deduplication + min_codons=%d filter.",
        n_total, n_ribo, len(gene_best), min_codons,
    )
    if not gene_best:
        raise ValueError(
            f"No ribosomal protein genes passed filters in {cocoputs_tsv}. "
            "Check that the file contains RPL/RPS genes and try lowering --min-codons."
        )

    # Aggregate codon counts across the selected reference genes
    aggregate: dict[str, int] = {c: 0 for c in _CODON_TABLE}
    for data in gene_best.values():
        for codon, cnt in data['counts'].items():
            aggregate[codon] += cnt

    rscu_w = _compute_rscu_w(aggregate)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as out:
        out.write('codon\trscu\tw\n')
        for codon in _ALL_CODONS_ALPHA:
            if codon in rscu_w:
                rscu, w = rscu_w[codon]
                out.write(f'{codon}\t{rscu:.6f}\t{w:.6f}\n')
            else:
                out.write(f'{codon}\tNA\tNA\n')   # stop codons

    log.info(
        "Wrote CAI weights for %d reference genes → %s",
        len(gene_best), output_path,
    )


def _compute_rscu_w(
    aggregate: dict[str, int],
) -> dict[str, tuple[float, float]]:
    """Return {codon: (rscu, w)} for all sense codons."""
    result: dict[str, tuple[float, float]] = {}
    for aa, codons in _AA_TO_CODONS.items():
        n_syn = len(codons)
        total = sum(aggregate.get(c, 0) for c in codons)

        if total == 0:
            log.warning(
                "AA %s: zero total usage across all reference genes; "
                "assigning uniform weights (rscu=1.0, w=1.0).",
                aa,
            )
            for c in codons:
                result[c] = (1.0, 1.0)
            continue

        rscus = {c: (aggregate.get(c, 0) * n_syn) / total for c in codons}
        max_rscu = max(rscus.values())
        for c, rscu in rscus.items():
            result[c] = (rscu, rscu / max_rscu)

    return result


def load_weights(weights_path: Path) -> dict[str, float]:
    """Load a weights TSV and return {codon: w} for sense codons only.

    Stop-codon rows (w == 'NA') are silently skipped. Sense codons with
    w == 0 are included; the caller decides whether to exclude them from
    scoring.
    """
    weights: dict[str, float] = {}
    with open(weights_path) as fh:
        next(fh)   # skip header
        for line in fh:
            parts = line.rstrip('\r\n').split('\t')
            if len(parts) < 3 or parts[2] == 'NA':
                continue
            codon = parts[0].strip().upper()
            try:
                weights[codon] = float(parts[2])
            except ValueError:
                continue
    return weights
