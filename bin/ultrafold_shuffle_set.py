#!/usr/bin/env python3
"""
bin/ultrafold_shuffle_set.py

Standalone step: pick the N shortest ultrafolded transcripts (shortest = fastest
to fold) and write one single-sequence FASTA per transcript, plus the list file
that bin/04_submit.sh --list consumes.

Reads
  runs/<dataset>/ultrafold/candidates.tsv           (ultrafold_candidates.py)
  runs/<dataset>/extracted_regions/extracted_mRNA.fa
  or --ids FILE   one transcript_id per line, instead of "N shortest"

Writes (runs/<dataset>/ultrafold/)
  shuffle_set.tsv              the chosen transcripts (id, gene, mrna_len, z)
  shuffle_set/<seqname>.fa     2-line FASTA; <seqname> is the mRNA FASTA header
                               (GENE_TRANSCRIPT_mRNA), so rnafold collate and
                               manifest.tsv join on it
  shuffle_set.fasta_list.txt   absolute paths, one per line

Run this where the job will run (the list holds absolute paths).

Usage
  ./bin/ultrafold_shuffle_set.py -d human_all -n 5
  ./bin/ultrafold_shuffle_set.py -d human_all --ids my_ids.txt
Stdlib only.
"""
import argparse
import csv
import logging
import os
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger('ultrafold_shuffle_set')


def read_fasta_by_transcript(path, wanted):
    """{tid: (header, sequence)} for headers of the form GENE_TID_mRNA."""
    out, name, hdr, buf = {}, None, None, []
    with open(path) as f:
        for line in list(f) + ['>']:
            line = line.strip()
            if line.startswith('>'):
                if name in wanted:
                    out[name] = (hdr, ''.join(buf))
                hdr = line[1:]
                parts = hdr.split('_')
                name = parts[1] if len(parts) >= 3 else None
                buf = []
            else:
                buf.append(line)
    return out


def pick(candidates, n=None, ids=None):
    """Rows of the chosen transcripts. Shortest n by mrna_len unless ids given."""
    if ids is not None:
        by = {r['transcript_id']: r for r in candidates}
        missing = [i for i in ids if i not in by]
        if missing:
            raise ValueError(f'not in candidates.tsv: {missing}')
        return [by[i] for i in ids]
    return sorted(candidates, key=lambda r: int(r['mrna_len']))[:n]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-d', '--dataset', required=True)
    ap.add_argument('-n', type=int, default=5, help='how many shortest (default 5)')
    ap.add_argument('--ids', help='file of transcript_ids; overrides -n')
    ap.add_argument('--runs-root', default=os.environ.get(
        'RUNS_ROOT', str(Path(__file__).resolve().parent.parent / 'runs')))
    args = ap.parse_args()

    run_dir = Path(args.runs_root).resolve() / args.dataset
    out_dir = run_dir / 'ultrafold'
    cand_path = out_dir / 'candidates.tsv'
    fa_path = run_dir / 'extracted_regions' / 'extracted_mRNA.fa'
    for p in (cand_path, fa_path):
        if not p.exists():
            sys.exit(f'ERROR: not found: {p}')

    with open(cand_path, newline='') as f:
        cands = list(csv.DictReader(f, delimiter='\t'))
    ids = None
    if args.ids:
        ids = [l.strip() for l in open(args.ids) if l.strip()]
    try:
        chosen = pick(cands, args.n, ids)
    except ValueError as e:
        sys.exit(f'ERROR: {e}')

    seqs = read_fasta_by_transcript(fa_path, {r['transcript_id'] for r in chosen})
    fa_dir = out_dir / 'shuffle_set'
    fa_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    with open(out_dir / 'shuffle_set.tsv', 'w') as fo:
        fo.write('transcript_id\tseqname\tgene_name\tmrna_len\tmrna_zscore\n')
        for r in chosen:
            tid = r['transcript_id']
            if tid not in seqs:
                sys.exit(f'ERROR: {tid} not found in {fa_path}')
            hdr, seq = seqs[tid]
            if len(seq) != int(r['mrna_len']):
                sys.exit(f'ERROR: {tid} FASTA length {len(seq)} != manifest {r["mrna_len"]}')
            fp = fa_dir / f'{hdr}.fa'
            fp.write_text(f'>{hdr}\n{seq}\n')
            paths.append(str(fp))
            fo.write(f'{tid}\t{hdr}\t{r["gene_name"]}\t{r["mrna_len"]}\t{r["mrna_zscore"]}\n')
            log.info(f'{r["gene_name"]:8} {tid}  {r["mrna_len"]} nt  z={r["mrna_zscore"]}')
    (out_dir / 'shuffle_set.fasta_list.txt').write_text('\n'.join(paths) + '\n')
    log.info(f'wrote {len(paths)} FASTAs and shuffle_set.fasta_list.txt -> {out_dir}')


if __name__ == '__main__':
    main()
