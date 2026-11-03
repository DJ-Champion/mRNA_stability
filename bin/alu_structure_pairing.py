#!/usr/bin/env python3
"""
bin/alu_structure_pairing.py

Standalone step: project the full-mRNA RNAfold MFE structure onto the Alus and
say which Alus pair with an oppositely oriented Alu IN THE STRUCTURE (not just
by sequence orientation), at a sweep of cutoffs.

Reads
  runs/<dataset>/rnafold_structure/results/<seqname>.fold   (tools/rnafold_structure)
  runs/<dataset>/rnafold_structure/results/<seqname>.bpp.tsv   optional (WITH_BPP=1)
  runs/<dataset>/ultrafold/shuffle_summary.tsv     seqname -> transcript_id, native_mfe
  runs/<repeats-dataset>/ultrafold/repeats.tsv     (repeats_on_transcripts.py)
  runs/<repeats-dataset>/ultrafold/candidates.tsv  region lengths, gene names
  runs/<repeats-dataset>/ultrafold/alu_summary.tsv optional, orientation-only counts

Coordinates are 0-based half-open, 5'->3' on the spliced transcript: the frame
of extracted_mRNA.fa, the dot-bracket string and repeats.tsv.

Stops (exit 1) if, for any transcript, the folded sequence length differs from
candidates mrna_len, or the MFE differs from native_mfe in shuffle_summary.tsv
by more than --mfe-tol (default 0.01 kcal/mol). Nothing is written then.

Writes (runs/<dataset>/ultrafold/structure/)
  model.tsv                transcript_id, gene_name, mrna_len, utr5_len, cds_len, utr3_len,
                           mfe, stored_mfe, mfe_match
  alu_pair_counts.tsv      one row per Alu pair A,B (A 5' of B, A != B) with n_bp > 0
                           transcript_id, alu_a, alu_b, orient_a, orient_b, opposite, n_bp,
                           frac_a, frac_b, a_start, a_end, b_start, b_end
  alu_partners.tsv         one row per Alu: bases paired, split by partner type
                           transcript_id, alu, orientation, tx_start, tx_end, alu_len,
                           n_opposite, n_same, n_self, n_nonalu, n_unpaired
                           (counts are bases of this Alu; they sum to alu_len)
  alu_paired_by_cutoff.tsv one row per Alu x cutoff
                           transcript_id, alu, orientation, tx_start, tx_end, cutoff_type,
                           cutoff_value, paired_with_opposite, n_bp_opposite, n_bp_same,
                           n_bp_nonalu
  transcript_summary.tsv   one row per transcript x cutoff
                           transcript_id, gene_name, mrna_len, n_alu, n_sense, n_antisense,
                           cutoff_type, cutoff_value, n_paired_opposite,
                           n_opposite_pairs_orientation (alu_summary.tsv; sequence
                           orientation only, in its pair region), n_opposite_pairs_structure
                           (Alu pairs of opposite orientation with n_bp >= the cutoff, nbp only)
  arcs.tsv                 every base pair
                           transcript_id, i, j, alu_i, alu_j, orient_i, orient_j, arc_class
  arcs_collapsed.tsv       runs of stacked pairs (i,j),(i+1,j-1),... merged into one stem
                           transcript_id, i_start, i_end, j_start, j_end, n_bp, alu_i, alu_j,
                           orient_i, orient_j, arc_class     (arms half-open; the Alu of an
                           arm is the one holding most of its bases)
  bpp.tsv                  if bpp files exist: transcript_id, i, j, prob

arc_class: alu_opposite | alu_same | alu_self (both ends in one Alu) |
           alu_other (one end in an Alu, the other not) | other.

Cutoffs. "Alu A is paired with an opposite Alu" iff the number of A's bases paired
to an opposite-orientation Alu is >= N bases (nbp) or >= X% of A's exonic length
(frac). Comparison is >=. Cutoffs are inputs (--nbp, --frac), not constants.

Caveats
  * One MFE structure is one prediction; use bpp.tsv before calling a pair robust.
  * A repeat split over exons is one Alu (pieces share repeat_id); its length is
    the sum of its pieces.
  * If two Alu annotations overlap, a base belongs to the earlier-starting one.
  * Only () brackets are accepted (plain RNAfold); anything else raises.

Stdlib only.

Usage
  ./bin/alu_structure_pairing.py -d ultrafold_shuffle_set --repeats-dataset human_all
"""
import argparse
import csv
import logging
import os
import sys
from collections import defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger('alu_structure_pairing')

DEFAULT_NBP = (1, 5, 10, 15, 20, 30)
DEFAULT_FRAC = (5, 10, 25, 50)


# ---------------------------------------------------------------- parsing

def parse_dotbracket(struct):
    """Return partner[i] = j (0-based) or -1. Raises ValueError if not a balanced ()."""
    partner = [-1] * len(struct)
    stack = []
    for k, c in enumerate(struct):
        if c == '(':
            stack.append(k)
        elif c == ')':
            if not stack:
                raise ValueError(f'unbalanced: ")" at {k} has no partner')
            i = stack.pop()
            partner[i], partner[k] = k, i
        elif c != '.':
            raise ValueError(f'unsupported structure character {c!r} at {k}')
    if stack:
        raise ValueError(f'unbalanced: "(" at {stack[-1]} has no partner')
    return partner


def read_fold(path):
    """RNAfold stdout: (name, sequence, structure, mfe)."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 3 or not lines[0].startswith('>'):
        raise ValueError(f'{path}: not RNAfold output')
    name, seq = lines[0][1:].strip(), lines[1].strip()
    struct, _, rest = lines[2].partition(' ')
    mfe = float(rest.strip().strip('()').strip())
    if len(struct) != len(seq):
        raise ValueError(f'{path}: structure length {len(struct)} != sequence length {len(seq)}')
    return name, seq, struct, mfe


def read_tsv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


# ---------------------------------------------------------------- Alus

def build_alus(repeat_rows):
    """Alu rows of one transcript (repeats.tsv, any family) -> list of dicts, 5' first."""
    by_id = defaultdict(list)
    for r in repeat_rows:
        if r['rep_family'] == 'Alu':
            by_id[r['repeat_id']].append(r)
    alus = []
    for rid, pieces in by_id.items():
        ivs = sorted((int(p['tx_start']), int(p['tx_end'])) for p in pieces)
        alus.append({'id': rid, 'orientation': pieces[0]['orientation'],
                     'start': ivs[0][0], 'end': ivs[-1][1], 'pieces': ivs,
                     'len': sum(e - s for s, e in ivs)})
    alus.sort(key=lambda a: (a['start'], a['end'], a['id']))
    return alus


def base_to_alu(alus, n):
    """owner[k] = index into alus or -1. Overlaps: earlier-starting Alu wins."""
    owner = [-1] * n
    n_overlap = 0
    for idx, a in enumerate(alus):
        for s, e in a['pieces']:
            for k in range(s, min(e, n)):
                if owner[k] == -1:
                    owner[k] = idx
                else:
                    n_overlap += 1
    return owner, n_overlap


def pair_class(a, b, alus):
    """a, b: Alu indices or -1 for the two ends of a pair."""
    if a == -1 and b == -1:
        return 'other'
    if a == -1 or b == -1:
        return 'alu_other'
    if a == b:
        return 'alu_self'
    return 'alu_opposite' if alus[a]['orientation'] != alus[b]['orientation'] else 'alu_same'


def project(partner, alus, owner):
    """Pair counts and per-Alu partner-type base counts.

    pair_bp[(a,b)] with a<b: base pairs between Alu a and Alu b.
    bases[a] = dict(opposite, same, self, nonalu, unpaired) in bases of Alu a.
    """
    pair_bp = defaultdict(int)
    bases = [{'opposite': 0, 'same': 0, 'self': 0, 'nonalu': 0, 'unpaired': 0} for _ in alus]
    for i, j in enumerate(partner):
        if j > i:
            a, b = owner[i], owner[j]
            if a != -1 and b != -1 and a != b:
                pair_bp[(min(a, b), max(a, b))] += 1
    for k, a in enumerate(owner):
        if a == -1:
            continue
        j = partner[k]
        if j == -1:
            bases[a]['unpaired'] += 1
            continue
        b = owner[j]
        c = pair_class(a, b, alus)
        key = {'alu_opposite': 'opposite', 'alu_same': 'same',
               'alu_self': 'self', 'alu_other': 'nonalu'}[c]
        bases[a][key] += 1
    return dict(pair_bp), bases


def called_paired(n_opposite, alu_len, ctype, cval):
    """paired-with-opposite call; >= at the boundary."""
    if ctype == 'nbp':
        return n_opposite >= cval
    if ctype == 'frac':
        return n_opposite * 100 >= cval * alu_len     # integers: no float edge
    raise ValueError(ctype)


def collapse_stems(partner):
    """Runs of (i,j),(i+1,j-1),... -> list of (i_start, j_end_excl_of_inner...) as
    (i_start, i_end, j_start, j_end) half-open, i_end-i_start == j_end-j_start == n_bp.
    j arm: j_start is the innermost 3' position, j_end-1 the outermost."""
    stems = []
    k, n = 0, len(partner)
    while k < n:
        j = partner[k]
        if j > k:
            m = 1
            while k + m < n and partner[k + m] == j - m and j - m > k + m:
                m += 1
            stems.append((k, k + m, j - m + 1, j + 1))
            k += m
        else:
            k += 1
    return stems


def arm_alu(owner, s, e):
    """Alu index holding most bases of [s,e), or -1 if none (ties: lowest index)."""
    cnt = defaultdict(int)
    for k in range(s, e):
        if owner[k] != -1:
            cnt[owner[k]] += 1
    if not cnt:
        return -1
    return max(sorted(cnt), key=lambda a: cnt[a])


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('-d', '--dataset', default='ultrafold_shuffle_set')
    ap.add_argument('--repeats-dataset', default='human_all',
                    help='dataset holding ultrafold/repeats.tsv + candidates.tsv')
    ap.add_argument('--runs-root', default=os.environ.get('RUNS_ROOT', 'runs'))
    ap.add_argument('--fold-dir', help='default runs/<dataset>/rnafold_structure/results')
    ap.add_argument('--out-dir', help='default runs/<dataset>/ultrafold/structure')
    ap.add_argument('--nbp', type=int, nargs='+', default=list(DEFAULT_NBP))
    ap.add_argument('--frac', type=float, nargs='+', default=list(DEFAULT_FRAC),
                    help='percent of the Alu length')
    ap.add_argument('--mfe-tol', type=float, default=0.01)
    args = ap.parse_args()

    runs = Path(args.runs_root)
    ds_ultra = runs / args.dataset / 'ultrafold'
    rep_ultra = runs / args.repeats_dataset / 'ultrafold'
    fold_dir = Path(args.fold_dir) if args.fold_dir else runs / args.dataset / 'rnafold_structure' / 'results'
    out_dir = Path(args.out_dir) if args.out_dir else ds_ultra / 'structure'

    summary = read_tsv(ds_ultra / 'shuffle_summary.tsv')
    cand = {r['transcript_id']: r for r in read_tsv(rep_ultra / 'candidates.tsv')}
    reps = defaultdict(list)
    for r in read_tsv(rep_ultra / 'repeats.tsv'):
        reps[r['transcript_id']].append(r)
    alu_sum = {}
    if (rep_ultra / 'alu_summary.tsv').exists():
        alu_sum = {r['transcript_id']: r for r in read_tsv(rep_ultra / 'alu_summary.tsv')}

    model, pair_rows, partner_rows, cut_rows, sum_rows = [], [], [], [], []
    arc_rows, stem_rows, bpp_rows = [], [], []
    problems = []

    for srow in summary:
        seqname, tid = srow['seqname'], srow['transcript_id']
        fold = fold_dir / f'{seqname}.fold'
        if tid not in cand:
            problems.append(f'{tid}: not in candidates.tsv')
            continue
        if not fold.exists():
            problems.append(f'{seqname}: missing {fold}')
            continue
        c = cand[tid]
        _, seq, struct, mfe = read_fold(fold)
        n = len(seq)
        stored = float(srow['native_mfe'])
        if n != int(c['mrna_len']):
            problems.append(f'{seqname}: folded length {n} != candidates mrna_len {c["mrna_len"]}')
        match = abs(mfe - stored) <= args.mfe_tol
        if not match:
            problems.append(f'{seqname}: MFE {mfe} != stored native_mfe {stored}')
        model.append({'transcript_id': tid, 'gene_name': c['gene_name'], 'mrna_len': c['mrna_len'],
                      'utr5_len': c['utr5_len'], 'cds_len': c['cds_len'], 'utr3_len': c['utr3_len'],
                      'mfe': f'{mfe:.2f}', 'stored_mfe': f'{stored:.4f}', 'mfe_match': str(match).lower()})
        if problems:
            continue          # keep checking the rest, but write nothing

        partner = parse_dotbracket(struct)
        alus = build_alus(reps.get(tid, []))
        owner, n_overlap = base_to_alu(alus, n)
        if n_overlap:
            log.warning('%s: %d bases lie in overlapping Alu annotations', seqname, n_overlap)
        pair_bp, bases = project(partner, alus, owner)

        for (a, b), nbp in sorted(pair_bp.items()):
            A, B = alus[a], alus[b]
            pair_rows.append({
                'transcript_id': tid, 'alu_a': A['id'], 'alu_b': B['id'],
                'orient_a': A['orientation'], 'orient_b': B['orientation'],
                'opposite': str(A['orientation'] != B['orientation']).lower(), 'n_bp': nbp,
                'frac_a': f'{nbp / A["len"]:.4f}', 'frac_b': f'{nbp / B["len"]:.4f}',
                'a_start': A['start'], 'a_end': A['end'], 'b_start': B['start'], 'b_end': B['end']})

        cutoffs = [('nbp', v) for v in args.nbp] + [('frac', v) for v in args.frac]
        for idx, A in enumerate(alus):
            bs = bases[idx]
            partner_rows.append({
                'transcript_id': tid, 'alu': A['id'], 'orientation': A['orientation'],
                'tx_start': A['start'], 'tx_end': A['end'], 'alu_len': A['len'],
                'n_opposite': bs['opposite'], 'n_same': bs['same'], 'n_self': bs['self'],
                'n_nonalu': bs['nonalu'], 'n_unpaired': bs['unpaired']})
            for ctype, cval in cutoffs:
                cut_rows.append({
                    'transcript_id': tid, 'alu': A['id'], 'orientation': A['orientation'],
                    'tx_start': A['start'], 'tx_end': A['end'],
                    'cutoff_type': ctype, 'cutoff_value': f'{cval:g}',
                    'paired_with_opposite': str(called_paired(bs['opposite'], A['len'], ctype, cval)).lower(),
                    'n_bp_opposite': bs['opposite'], 'n_bp_same': bs['same'], 'n_bp_nonalu': bs['nonalu']})

        n_sense = sum(a['orientation'] == 'sense' for a in alus)
        asum = alu_sum.get(tid, {})
        for ctype, cval in cutoffs:
            n_par = sum(called_paired(bs['opposite'], A['len'], ctype, cval)
                        for A, bs in zip(alus, bases))
            if ctype == 'nbp':
                n_struct = sum(1 for (a, b), nbp in pair_bp.items()
                               if alus[a]['orientation'] != alus[b]['orientation'] and nbp >= cval)
            else:
                n_struct = 'NA'
            sum_rows.append({
                'transcript_id': tid, 'gene_name': c['gene_name'], 'mrna_len': n, 'n_alu': len(alus),
                'n_sense': n_sense, 'n_antisense': len(alus) - n_sense,
                'cutoff_type': ctype, 'cutoff_value': f'{cval:g}', 'n_paired_opposite': n_par,
                'n_opposite_pairs_orientation': asum.get('n_opposite_pairs', 'NA'),
                'n_opposite_pairs_structure': n_struct})

        for i, j in enumerate(partner):
            if j > i:
                a, b = owner[i], owner[j]
                arc_rows.append({
                    'transcript_id': tid, 'i': i, 'j': j,
                    'alu_i': alus[a]['id'] if a != -1 else '', 'alu_j': alus[b]['id'] if b != -1 else '',
                    'orient_i': alus[a]['orientation'] if a != -1 else '',
                    'orient_j': alus[b]['orientation'] if b != -1 else '',
                    'arc_class': pair_class(a, b, alus)})
        for i_s, i_e, j_s, j_e in collapse_stems(partner):
            a, b = arm_alu(owner, i_s, i_e), arm_alu(owner, j_s, j_e)
            stem_rows.append({
                'transcript_id': tid, 'i_start': i_s, 'i_end': i_e, 'j_start': j_s, 'j_end': j_e,
                'n_bp': i_e - i_s,
                'alu_i': alus[a]['id'] if a != -1 else '', 'alu_j': alus[b]['id'] if b != -1 else '',
                'orient_i': alus[a]['orientation'] if a != -1 else '',
                'orient_j': alus[b]['orientation'] if b != -1 else '',
                'arc_class': pair_class(a, b, alus)})

        bpp = fold_dir / f'{seqname}.bpp.tsv'
        if bpp.exists():
            for r in read_tsv(bpp):
                bpp_rows.append({'transcript_id': tid, 'i': r['i'], 'j': r['j'], 'prob': r['prob']})

    n_ok = len(model)
    for p in problems:
        log.error(p)
    if problems:
        log.error('stopping: nothing written (%d problem(s))', len(problems))
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)

    def write(name, rows, cols):
        with open(out_dir / name, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
        log.info('wrote %s (%d rows)', out_dir / name, len(rows))

    write('model.tsv', model, ['transcript_id', 'gene_name', 'mrna_len', 'utr5_len', 'cds_len',
                               'utr3_len', 'mfe', 'stored_mfe', 'mfe_match'])
    write('alu_pair_counts.tsv', pair_rows, ['transcript_id', 'alu_a', 'alu_b', 'orient_a', 'orient_b',
                                             'opposite', 'n_bp', 'frac_a', 'frac_b',
                                             'a_start', 'a_end', 'b_start', 'b_end'])
    write('alu_partners.tsv', partner_rows, ['transcript_id', 'alu', 'orientation', 'tx_start', 'tx_end',
                                             'alu_len', 'n_opposite', 'n_same', 'n_self', 'n_nonalu',
                                             'n_unpaired'])
    write('alu_paired_by_cutoff.tsv', cut_rows, ['transcript_id', 'alu', 'orientation', 'tx_start', 'tx_end',
                                                 'cutoff_type', 'cutoff_value', 'paired_with_opposite',
                                                 'n_bp_opposite', 'n_bp_same', 'n_bp_nonalu'])
    write('transcript_summary.tsv', sum_rows, ['transcript_id', 'gene_name', 'mrna_len', 'n_alu',
                                               'n_sense', 'n_antisense', 'cutoff_type', 'cutoff_value',
                                               'n_paired_opposite', 'n_opposite_pairs_orientation',
                                               'n_opposite_pairs_structure'])
    arc_cols = ['alu_i', 'alu_j', 'orient_i', 'orient_j', 'arc_class']
    write('arcs.tsv', arc_rows, ['transcript_id', 'i', 'j'] + arc_cols)
    write('arcs_collapsed.tsv', stem_rows, ['transcript_id', 'i_start', 'i_end', 'j_start', 'j_end',
                                            'n_bp'] + arc_cols)
    if bpp_rows:
        write('bpp.tsv', bpp_rows, ['transcript_id', 'i', 'j', 'prob'])
    log.info('%d transcript(s); MFE matches the stored value for all', n_ok)
    return 0


if __name__ == '__main__':
    sys.exit(main())
