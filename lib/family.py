"""lib/family.py
Pure helpers for cohort-level sequence-family clustering.

Consumed by bin/01c_family.py. Everything here is side-effect free and
directly testable; all I/O, subprocess and orchestration lives in the
driver script.

See FAMILY_CLUSTERING.md for the specification these implement.
"""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from lib.gff import normalise_id

# Protein IDs are '<species>|<gene_id>|<transcript_id>'. Pipe is forbidden
# inside any component so the split is unambiguous.
ID_SEP = '|'

_SAFE_COMPONENT = re.compile(r'^[^|\s]+$')


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------

def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (record_id, sequence) from a FASTA file.

    Handles line-wrapped sequences (seqkit wraps at 60 by default).
    record_id is the header up to the first whitespace.
    """
    header: str | None = None
    chunks: list[str] = []
    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks)
                header = line[1:].split(None, 1)[0]
                chunks = []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, ''.join(chunks)


def make_protein_id(species: str, gene_id: str, transcript_id: str) -> str:
    """Build the namespaced protein ID, rejecting components that would
    make the ID ambiguous."""
    for label, part in (('species', species), ('gene_id', gene_id),
                        ('transcript_id', transcript_id)):
        if not _SAFE_COMPONENT.match(part):
            raise ValueError(
                f"{label} {part!r} contains whitespace or '{ID_SEP}'; "
                f"protein IDs would be ambiguous")
    return ID_SEP.join((species, gene_id, transcript_id))


def split_protein_id(protein_id: str) -> tuple[str, str, str]:
    """Inverse of make_protein_id."""
    parts = protein_id.split(ID_SEP)
    if len(parts) != 3:
        raise ValueError(f"Malformed protein ID: {protein_id!r}")
    return parts[0], parts[1], parts[2]


# ---------------------------------------------------------------------------
# Translation post-processing
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TranslationResult:
    aa: str
    had_internal_stop: bool
    had_terminal_stop: bool


def post_process_translation(aa: str) -> TranslationResult:
    """Clean a raw seqkit translation for search.

    Deliberately does NOT use `seqkit translate --trim`. Despite its help
    text ("remove all 'X' and '*' characters from the right end"), --trim
    truncates the sequence at the FIRST '*' or 'X' — verified against
    seqkit v2.13.0:

        MAWK*   -> MAWK    (correct)
        MA*KA*  -> MA      (should be MA*KA)
        MAXK*   -> MA      (should be MAXK)

    So a selenoprotein (UGA read through as Sec) or any CDS containing an
    assembly N is silently reduced to an N-terminal stub, which then fails
    the coverage filters and becomes a spurious singleton. That is
    under-merging — the unsafe direction for leakage — and it is invisible
    downstream.

    Instead: strip exactly one terminal stop, then convert any remaining
    internal stop to X. mmseqs treats X as unmatchable, so an internal stop
    weakens the alignment locally rather than corrupting the record.
    """
    had_terminal_stop = aa.endswith('*')
    if had_terminal_stop:
        aa = aa[:-1]                       # exactly one — do not loop
    had_internal_stop = '*' in aa
    if had_internal_stop:
        aa = aa.replace('*', 'X')
    return TranslationResult(aa, had_internal_stop, had_terminal_stop)


# ---------------------------------------------------------------------------
# Levels
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Level:
    """One clustering stringency. All five criteria must hold for an edge."""
    name: str
    min_nbs: float
    min_cov: float
    max_evalue: float
    min_fident: float

    def accepts(self, nbs: float, qcov: float, tcov: float,
                evalue: float, fident: float) -> bool:
        # Both coverages, not either: this is what --cov-mode 0 means, and
        # it is the main brake on domain-mediated chaining.
        return (nbs >= self.min_nbs
                and qcov >= self.min_cov
                and tcov >= self.min_cov
                and evalue <= self.max_evalue
                and fident >= self.min_fident)


@dataclass(frozen=True)
class SearchParams:
    sensitivity: float
    max_seqs: int
    evalue: float
    min_cov: float
    cov_mode: int
    threads: int


def validate_level_within_search(level: Level, search: SearchParams) -> list[str]:
    """Return human-readable violations if `level` is looser than `search`.

    Downstream filtering can only tighten. A level looser than the search
    would silently evaluate against edges that were never searched for, so
    this must abort rather than warn.
    """
    problems: list[str] = []
    if level.max_evalue > search.evalue:
        problems.append(
            f"level '{level.name}': max_evalue {level.max_evalue:g} is looser "
            f"than search evalue {search.evalue:g}")
    if level.min_cov < search.min_cov:
        problems.append(
            f"level '{level.name}': min_cov {level.min_cov:g} is looser "
            f"than search min_cov {search.min_cov:g}")
    return problems


# ---------------------------------------------------------------------------
# Normalised bitscore
# ---------------------------------------------------------------------------

def normalised_bitscore(bits: float, self_bits_q: float, self_bits_t: float) -> float:
    """bits(q,t) / max(self_bits[q], self_bits[t]), clamped to [0, 1].

    max() in the denominator — not min(), not the geometric mean, and not
    the classic BSR convention of dividing by the query self-score. This is
    the primary anti-chaining lever.

    A 100 aa protein fully matching one domain inside a 900 aa protein
    scores bits(q,t) ~ self_bits[q] ~ 200 against self_bits[t] ~ 1800:

        / self_bits[q]  -> ~1.00   edge kept, families weld together
        / min(...)      -> ~1.00   edge kept, families weld together
        / max(...)      -> ~0.11   edge rejected

    Only max() demands similarity across the longer sequence's full length.
    It also degrades more gracefully across a phylogeny than raw identity,
    which matters once a cohort spans several species.
    """
    denom = max(self_bits_q, self_bits_t)
    if denom <= 0:
        return 0.0
    return min(1.0, bits / denom)


# ---------------------------------------------------------------------------
# Union-find
# ---------------------------------------------------------------------------

class UnionFind:
    """Disjoint-set over [0, n), with path compression and union by rank.

    Connected components of the filtered similarity graph are exactly
    single-linkage clusters, which is the transitive closure required for
    the leakage guarantee.
    """
    __slots__ = ('_parent', '_rank', 'n')

    def __init__(self, n: int):
        self.n = n
        self._parent = list(range(n))
        self._rank = [0] * n

    def find(self, x: int) -> int:
        parent = self._parent
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:            # path compression
            parent[x], x = root, parent[x]
        return root

    def union(self, a: int, b: int) -> bool:
        """Merge the sets containing a and b. True if they were disjoint."""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        rank = self._rank
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        self._parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1
        return True


# ---------------------------------------------------------------------------
# Canonical family assignment
# ---------------------------------------------------------------------------

def canonical_families(uf: UnionFind, ids: Sequence[str],
                       level_name: str) -> tuple[list[str], list[list[int]]]:
    """Assign stable family IDs from membership alone.

    mmseqs is threaded and names clusters after representative sequences, so
    its raw labels are not reproducible. Ordering by (size desc, smallest
    member ID) depends only on the partition, so identical input yields
    identical labels.

    Returns (family_id_per_index, members_per_family) with both indexed
    consistently: members_per_family[k] are the indices whose family ID is
    the k-th assigned label.
    """
    groups: dict[int, list[int]] = {}
    for i in range(len(ids)):
        groups.setdefault(uf.find(i), []).append(i)

    # Precompute sort keys once — recomputing min() inside the comparator
    # would be O(n log n * family_size).
    keyed = [(-len(members), min(ids[i] for i in members), members)
             for members in groups.values()]
    keyed.sort(key=lambda t: (t[0], t[1]))

    family_of = [''] * len(ids)
    ordered_members: list[list[int]] = []
    width = max(5, len(str(len(keyed))))
    for rank, (_, _, members) in enumerate(keyed, start=1):
        label = f'fam_{level_name}_{rank:0{width}d}'
        for i in members:
            family_of[i] = label
        ordered_members.append(members)
    return family_of, ordered_members


# ---------------------------------------------------------------------------
# Search identity
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Locus overlap
# ---------------------------------------------------------------------------

_GENE_ID_RE = re.compile(r'(?:^|;)gene_id=([^;]+)')


def gene_exons_from_gff(path: Path) -> dict[str, tuple[str, str, list[tuple[int, int]]]]:
    """Merged exon intervals per gene: gene_id -> (chrom, strand, [(start, end)]).

    Two notes on why it is written this way.

    `canonical.gff` is produced by bin/01_extract.py write_filtered_gff, which
    keeps only lines linked to a retained transcript. Gene-level features have
    no such link, so the file contains NO `gene` rows — spans must come from
    the `gene_id` attribute on child features.

    Exons, not the gene span. Two genes nested one inside the other's intron
    share no transcribed sequence, so their extracted regions are independent
    and they must not be linked. Only exonic overlap means shared sequence.
    """
    raw: dict[str, tuple[str, str, list[list[int]]]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            m = _GENE_ID_RE.search(f[8])
            if not m:
                continue
            gid = normalise_id(m.group(1))
            entry = raw.get(gid)
            if entry is None:
                raw[gid] = (f[0], f[6], [[int(f[3]), int(f[4])]])
            else:
                entry[2].append([int(f[3]), int(f[4])])

    out: dict[str, tuple[str, str, list[tuple[int, int]]]] = {}
    for gid, (chrom, strand, ivs) in raw.items():
        ivs.sort()
        merged: list[list[int]] = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[gid] = (chrom, strand, [(s, e) for s, e in merged])
    return out


def find_overlapping_gene_pairs(
    gene_exons: dict[str, tuple[str, str, list[tuple[int, int]]]],
    min_overlap_bp: int = 100,
) -> list[tuple[str, str, str, int, str, str]]:
    """Gene pairs sharing at least `min_overlap_bp` of exonic sequence.

    Returns (gene_a, gene_b, chrom, overlap_bp, strand_a, strand_b) with
    gene_a < gene_b, sorted by descending overlap.

    Both strands count. An antisense overlap still means the two genes are
    transcribed from the same DNA, so their extracted regions are reverse
    complements — identical length and GC by construction, and correlated
    structure. That is non-independence regardless of orientation.
    """
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for gid, (chrom, _strand, ivs) in gene_exons.items():
        by_chrom.setdefault(chrom, []).extend((s, e, gid) for s, e in ivs)

    pair_bp: dict[tuple[str, str], int] = {}
    for ivs in by_chrom.values():
        ivs.sort()
        active: list[tuple[int, int, str]] = []      # (end, start, gene)
        for s, e, gid in ivs:
            if active:
                active = [a for a in active if a[0] > s]
            for ae, a_s, ag in active:
                if ag == gid:
                    continue
                ov = min(e, ae) - max(s, a_s)
                if ov > 0:
                    key = (ag, gid) if ag < gid else (gid, ag)
                    pair_bp[key] = pair_bp.get(key, 0) + ov
            active.append((e, s, gid))

    out = []
    for (ga, gb), bp in pair_bp.items():
        if bp < min_overlap_bp:
            continue
        chrom, sa, _ = gene_exons[ga]
        _, sb, _ = gene_exons[gb]
        out.append((ga, gb, chrom, bp, sa, sb))
    out.sort(key=lambda r: -r[3])
    return out


def search_hash(search: SearchParams, members: Iterable[tuple[str, str, int]]) -> str:
    """Short stable hash of everything that invalidates the alignment table.

    Covers the search parameters and the member list (dataset, species,
    genetic code). Deliberately excludes `levels` and `threads`: changing a
    level only adds columns, and thread count must not change the identity
    of the result.
    """
    payload = {
        'sensitivity': search.sensitivity,
        'max_seqs': search.max_seqs,
        'evalue': search.evalue,
        'min_cov': search.min_cov,
        'cov_mode': search.cov_mode,
        'members': sorted([list(m) for m in members]),
    }
    blob = json.dumps(payload, sort_keys=True, separators=(',', ':'))
    return hashlib.sha256(blob.encode()).hexdigest()[:8]
