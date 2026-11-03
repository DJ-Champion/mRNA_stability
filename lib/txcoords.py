"""
lib/txcoords.py
Genome <-> spliced-transcript coordinate mapping from a canonical.gff.

Transcript coordinates are 0-based, half-open, and run 5'->3' along the
spliced transcript, i.e. the same frame as the sequence in
extracted_mRNA.fa. Genomic coordinates passed to / returned from this module
are 0-based half-open (BED convention); GFF's 1-based inclusive coordinates
are converted on load.

Stdlib only.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

_ATTR_RE = re.compile(r'(\w+)=([^;]*)')


@dataclass
class Transcript:
    tid: str
    chrom: str
    strand: str                       # '+' or '-'
    gene_id: str = ''
    gene_name: str = ''
    exons: list = field(default_factory=list)   # ascending genomic (start0, end0)

    def finalize(self):
        self.exons.sort()
        return self

    @property
    def length(self) -> int:
        return sum(e - s for s, e in self.exons)

    @property
    def span(self):
        return self.exons[0][0], self.exons[-1][1]

    def map_interval(self, g_start0: int, g_end0: int):
        """Map a genomic interval onto the spliced transcript.

        Returns a list of (tx_start, tx_end, g_start, g_end) pieces, one per
        exon the interval overlaps, ordered 5'->3' along the transcript.
        Intronic stretches are dropped, so an interval spanning a splice
        junction yields >1 piece and an intronic one yields [].
        """
        pieces = []
        n = len(self.exons)
        # tx offset of each exon in transcript order
        if self.strand == '+':
            offs, acc = [0] * n, 0
            for i, (s, e) in enumerate(self.exons):
                offs[i] = acc
                acc += e - s
        else:
            offs, acc = [0] * n, 0
            for i in range(n - 1, -1, -1):
                offs[i] = acc
                acc += self.exons[i][1] - self.exons[i][0]

        for i, (es, ee) in enumerate(self.exons):
            os_, oe = max(g_start0, es), min(g_end0, ee)
            if os_ >= oe:
                continue
            if self.strand == '+':
                ts, te = offs[i] + (os_ - es), offs[i] + (oe - es)
            else:
                ts, te = offs[i] + (ee - oe), offs[i] + (ee - os_)
            pieces.append((ts, te, os_, oe))
        pieces.sort()
        return pieces


def load_transcripts(gff_path, wanted=None):
    """Read exon structure for transcripts from a GFF3 (canonical.gff).

    wanted: optional set of transcript IDs; others are skipped.
    Exons are attached to a transcript by their transcript_id attribute
    (falling back to Parent).
    """
    txs: dict[str, Transcript] = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9 or cols[2] not in ('exon', 'transcript'):
                continue
            attrs = dict(_ATTR_RE.findall(cols[8]))
            tid = attrs.get('transcript_id') or attrs.get('Parent', '')
            if cols[2] == 'transcript':
                tid = attrs.get('transcript_id') or attrs.get('ID', '')
            if wanted is not None and tid not in wanted:
                continue
            tx = txs.get(tid)
            if tx is None:
                tx = txs[tid] = Transcript(tid, cols[0], cols[6])
            if cols[2] == 'transcript':
                tx.gene_id = attrs.get('gene_id', tx.gene_id)
                tx.gene_name = attrs.get('gene_name', tx.gene_name)
            else:
                tx.exons.append((int(cols[3]) - 1, int(cols[4])))
    return {t: tx.finalize() for t, tx in txs.items() if tx.exons}


def region_of(tx_start, tx_end, utr5_len, cds_len, mrna_len):
    """Label a transcript-space interval as 5UTR / CDS / 3UTR by majority
    overlap; 'spans' is not used - the label is the region with most bases."""
    bounds = {'5UTR': (0, utr5_len),
              'CDS': (utr5_len, utr5_len + cds_len),
              '3UTR': (utr5_len + cds_len, mrna_len)}
    best, best_ov = 'NA', 0
    for name, (a, b) in bounds.items():
        ov = min(tx_end, b) - max(tx_start, a)
        if ov > best_ov:
            best, best_ov = name, ov
    return best
