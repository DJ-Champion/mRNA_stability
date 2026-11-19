#!/usr/bin/env python3
"""
tools/rnaup_avoidance/avoidance.py
RNAup mRNA:ncRNA interaction / avoidance for one mRNA window (plus its shuffles).

Input  : a FASTA of native start-codon windows (start codon +-30 nt, 63 nt),
         normally one record (the per-sequence file bin/02_stratify.sh writes).
         Several records also work. Only these windows are shuffled; the ncRNAs
         are used as given.
Output : $RESULTS_DIR/results_<seqname>.csv   one row per gene (seq_name key)
         $TOOL_DIR/raw_shuffles/<seqname>_shuffle_sums.csv.gz
                                              per-shuffle sums over all ncRNAs

Per gene: the native window and N dinucleotide shuffles (esl-shuffle -d, the
same method as tools/rnafold/worker.sh) go to RNAup in ONE stream per ncRNA,
ncRNA first (--interaction_first), so the ncRNA accessibility is computed once
per window and ncRNA, not once per shuffled mRNA. G_total and G_int are parsed from
"(G_total = G_int + G_u_target + G_u_query)"; pairs with no interaction count
as 0.0 / 0.0 and stay in every sum.

Work units and checkpointing: the task's cores are fed (ncRNA, chunk of
AVOID_CHUNK targets) units, longest ncRNA first. The shuffled sequences are
written once to $CKPT/shuffles.fa.gz and each unit's parsed energies to
$CKPT/<ncrna>.c<chunk>.<NNN>.tsv.gz.
A restarted task reuses both, so every ncRNA sees the SAME shuffles
(esl-shuffle is unseeded; regenerating them on restart would mix shuffle sets
inside one sum). If the stored shuffles no longer match the input (native
window changed, different N_SHUFFLES) the whole checkpoint is discarded and
rebuilt. Checkpoints are deleted after the results file is written unless
AVOID_KEEP_CKPT=1.

Robust Z, as tools/rnafold/worker.sh and bin/ultrafold_shuffle_summary.py:
    (native - median(shuf)) / (1.4826 * MAD(shuf));  NA when MAD == 0.
Energies are kcal/mol; more negative = stronger binding.

Env (set by configs/tools/rnaup_avoidance.sh and lib/paths.sh):
    RESULTS_DIR TOOL_DIR TMP_DIR ERRORS_DIR  RNAUP_BIN ESL_SHUFFLE_BIN
    NCRNA_FASTA NCRNA_FAMILIES  N_SHUFFLES MIN_VALID_PERCENT
    RNAUP_ARGS (default "--interaction_first --no_output_file -b -d2 --noLP")
    AVOID_REGION (manifest region of the window; must be start_codon_region), AVOID_CKPT_DIR,
    TOOL_REQUIRE_LEN (default start_codon_region:63), AVOID_CHUNK (default 101), AVOID_KEEP_CKPT, SLURM_CPUS_PER_TASK (parallel units per task)
"""
from __future__ import annotations

import csv
import gzip
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / 'bin'))
from ultrafold_shuffle_summary import MAD_SCALE, robust_stats, robust_z  # noqa: E402,F401

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger('rnaup_avoidance')

# manifest region -> window name used in the output files. Only the start-codon window is scored.
WINDOW_NAMES = {'start_codon_region': 'start_codon_pm30'}
DEFAULT_REQUIRE_LEN = 'start_codon_region:63'    # 30 + AUG + 30; shorter = 5'UTR < 30, clipped
FAMILIES = ['tRNA', 'rRNA', 'SRP', 'snRNA', 'Y_RNA', 'snoRNA']
DEFAULT_RNAUP_ARGS = '--interaction_first --no_output_file -b -d2 --noLP'

# "(((.&.))) 1,10 : 5,14 (-5.20 = -8.40 + 1.90 + 1.30)"
_NUM = r'\s*([-+]?(?:\d+\.?\d*|inf|nan))\s*'
ENERGY_RE = re.compile(rf'\({_NUM}={_NUM}\+{_NUM}\+{_NUM}\)')


def read_fasta(path):
    """[(name, seq)]; name = first header token, seq upper-cased, T -> U. Reads .gz too."""
    opener = gzip.open if str(path).endswith('.gz') else open
    out, name, buf = [], None, []
    with opener(path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(buf).upper().replace('T', 'U')))
                name, buf = (line[1:].split() or [''])[0], []
            elif line:
                buf.append(line)
    if name is not None:
        out.append((name, ''.join(buf).upper().replace('T', 'U')))
    return out


def length_ok(region, length, spec):
    """TOOL_REQUIRE_LEN ("region:63 region:>=20 region:20-1000") check; no entry = any length.
    bin/04_submit.sh applies the same rule to the manifest when it builds the task list."""
    for item in (spec or '').split():
        r, _, rule = item.partition(':')
        if r != region:
            continue
        if rule.startswith('>='):
            return length >= int(rule[2:])
        if '-' in rule:
            lo, hi = rule.split('-')
            return int(lo) <= length <= int(hi)
        return length == int(rule)
    return True


def gc_fraction(seq):
    return (seq.count('G') + seq.count('C')) / len(seq)


def load_families(path):
    fam = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip() or line.startswith('ncrna_id'):
                continue
            k, v = line.rstrip('\n').split('\t')[:2]
            fam[k] = v
    return fam


def load_ncrnas(fasta, fam_path):
    """[(id, seq, family)] in FASTA order; the FASTA and the family map must agree exactly."""
    ncs = read_fasta(fasta)
    fam = load_families(fam_path)
    ids = [n for n, _ in ncs]
    if len(set(ids)) != len(ids):
        raise SystemExit(f'duplicate ncRNA ids in {fasta}')
    if set(ids) != set(fam):
        raise SystemExit(f'ncRNA FASTA and {fam_path} disagree: only in FASTA {sorted(set(ids) - set(fam))}, '
                         f'only in map {sorted(set(fam) - set(ids))}')
    bad = sorted(set(fam.values()) - set(FAMILIES))
    if bad:
        raise SystemExit(f'unknown family names in {fam_path}: {bad}')
    return [(n, s, fam[n]) for n, s in ncs]


def parse_rnaup(text, target_names):
    """RNAup stdout -> {target_name: (G_total, G_int)}.

    A target gets the first energy line (contains '&' and '=') that follows its
    '>' header. Targets whose header appears with no such line are
    non-interacting (0.0, 0.0). If a header is missing altogether, or an energy
    line does not parse, the run was truncated or the output format is not what
    we expect, and we raise rather than score it as 0.
    """
    want = set(target_names)
    res, seen, cur = {}, set(), None
    for line in text.splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0] if line[1:].split() else ''
            cur = name if name in want else None
            if cur is not None:
                seen.add(cur)
            continue
        if cur is None or cur in res or '&' not in line or '=' not in line:
            continue
        m = ENERGY_RE.search(line)
        if not m:
            raise RuntimeError(f'cannot parse RNAup energy line for {cur}: {line.strip()[:200]}')
        gt, gi = float(m.group(1)), float(m.group(2))
        # RNAup 2.4.7 prints "(100.00 = 100.00 + x + inf)" (then ERROR) for a pair it
        # cannot score: no interaction, so 0/0 like a missing line.
        res[cur] = (0.0, 0.0) if (gt >= 100.0 or 'inf' in line or 'nan' in line) else (gt, gi)
    missing = want - seen
    if missing:
        raise RuntimeError(f'RNAup output lacks {len(missing)}/{len(want)} target headers '
                           f'(e.g. {sorted(missing)[:3]}); truncated run or unexpected output format')
    return {n: res.get(n, (0.0, 0.0)) for n in target_names}


def run_rnaup(rnaup_bin, rnaup_args, nc_name, nc_seq, targets, scratch):
    """One RNAup process: ncRNA first, then every target. Returns {target: (G_total, G_int)}."""
    stream = f'>{nc_name}\n{nc_seq}\n' + ''.join(f'>{n}\n{s}\n' for n, s in targets)
    cwd = tempfile.mkdtemp(dir=scratch, prefix='rnaup_')
    try:
        p = subprocess.run([rnaup_bin, *rnaup_args], input=stream, capture_output=True,
                           text=True, cwd=cwd)
        if p.returncode != 0:
            raise RuntimeError(f'RNAup exit {p.returncode} for {nc_name}: {p.stderr.strip()[:500]}')
        stray = os.listdir(cwd)
        if stray:
            log.warning('RNAup left %d file(s) in scratch for %s (e.g. %s); removed',
                        len(stray), nc_name, stray[0])
        return parse_rnaup(p.stdout, [n for n, _ in targets])
    finally:
        shutil.rmtree(cwd, ignore_errors=True)


def make_shuffles(esl_bin, name, seq, n):
    """N dinucleotide-preserving shuffles of one window (esl-shuffle -d -N, as rnafold/worker.sh)."""
    with tempfile.NamedTemporaryFile('w', suffix='.fa') as tf:
        tf.write(f'>{name}\n{seq}\n')
        tf.flush()
        p = subprocess.run([esl_bin, '-d', '-N', str(n), tf.name], capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f'esl-shuffle failed for {name}: {p.stderr.strip()[:300]}')
    seqs = []
    for line in p.stdout.splitlines():
        if line.startswith('>'):
            seqs.append([])
        elif line.strip() and seqs:
            seqs[-1].append(line.strip())
    return [''.join(s).upper().replace('T', 'U') for s in seqs]


def build_sequence_set(genes, esl_bin, n_shuf, min_valid, ckpt):
    """[(gene_name, [native, shuf1..shufN])]; shuffles persisted in ckpt so restarts reuse them."""
    path = ckpt / 'shuffles.fa.gz'
    if path.exists():
        grouped = {}
        for rec, s in read_fasta(path):
            grouped.setdefault(rec.rsplit('|', 1)[0], []).append(s)
        if ([(g, seqs[0]) for g, seqs in grouped.items()] == genes
                and all(len(seqs) == n_shuf + 1 for seqs in grouped.values())):
            return list(grouped.items())
        # Every unit checkpoint was computed against these shuffles, so they go too.
        log.warning('%s does not match the input (native window or N_SHUFFLES=%d changed); '
                    'stale checkpoint discarded', path, n_shuf)
        for p in ckpt.iterdir():
            p.unlink()
    out, dropped = [], []
    for g, seq in genes:
        shufs = make_shuffles(esl_bin, g, seq, n_shuf)
        shufs = [s for s in shufs if len(s) == len(seq)]
        if len(shufs) < min_valid:
            dropped.append((g, len(shufs)))
            continue
        out.append((g, [seq] + shufs))
    for g, k in dropped:
        log.error('%s: only %d/%d valid shuffles (min %d); gene skipped', g, k, n_shuf, min_valid)
    tmp = ckpt / 'shuffles.fa.gz.tmp'
    with gzip.open(tmp, 'wt') as f:
        for g, seqs in out:
            for i, s in enumerate(seqs):
                f.write(f'>{g}|{i}\n{s}\n')
    os.replace(tmp, path)
    return out


def unit_path(ckpt, name, chunk, ci):
    return ckpt / f'{name}.c{chunk}.{ci:03d}.tsv.gz'


def unit_checkpoint(ckpt, nc, chunk, ci, targets, rnaup_bin, rnaup_args, scratch):
    """One work unit: ncRNA `nc` against targets[ci*chunk:(ci+1)*chunk]. ncRNA first in the
    stream, so its accessibility is computed once per unit (~0.5 s, about 2% of a 100-target unit)."""
    name, seq, _ = nc
    path = unit_path(ckpt, name, chunk, ci)
    if path.exists():
        return path
    part = targets[ci * chunk:(ci + 1) * chunk]
    res = run_rnaup(rnaup_bin, rnaup_args, name, seq, part, scratch)
    tmp = ckpt / (path.name + '.tmp')
    with gzip.open(tmp, 'wt') as f:
        for t, _ in part:
            gt, gi = res[t]
            f.write(f'{gt:.2f}\t{gi:.2f}\n')
    os.replace(tmp, path)
    return path


def fmt(v):
    return 'NA' if v is None else f'{v:.4f}'


def finalize(tranche, seqsets, ncs, ckpt, chunk, n_chunks, results_dir, raw_dir, window):
    """Sum ncRNA checkpoints -> per-gene result rows + per-shuffle sums."""
    n_seq = sum(len(s) for _, s in seqsets)
    fam_idx = {f: [k for k, nc in enumerate(ncs) if nc[2] == f] for f in FAMILIES}
    # sums[metric][group] = list over sequences; metric 0 = G_total, 1 = G_int
    groups = ['all'] + FAMILIES
    sums = {(m, g): [0.0] * n_seq for m in (0, 1) for g in groups}
    native_min = {m: {} for m in (0, 1)}      # gene -> min over ncRNAs of the native value
    starts, i = {}, 0
    for g, s in seqsets:
        starts[g] = i
        i += len(s)
    native_pos = {starts[g]: g for g, _ in seqsets}
    for k, (name, _, family) in enumerate(ncs):
        j = -1
        for ci in range(n_chunks):      # chunk order by index, not by sorted file name
            with gzip.open(unit_path(ckpt, name, chunk, ci), 'rt') as f:
                for line in f:
                    j += 1
                    a, b = line.split('\t')
                    for m, v in ((0, float(a)), (1, float(b))):
                        sums[(m, 'all')][j] += v
                        sums[(m, family)][j] += v
                        if j in native_pos:
                            d = native_min[m]
                            gname = native_pos[j]
                            d[gname] = v if gname not in d else min(d[gname], v)
        if j + 1 != n_seq:
            raise RuntimeError(f'checkpoint {name} has {j + 1} rows, expected {n_seq}')

    nc_gc = sum(gc_fraction(s) for _, s, _ in ncs) / len(ncs)
    mname = {0: 'rnaup_G_total', 1: 'avoidance_G_int'}
    cols = []
    for m in (0, 1):
        cols += [f'{mname[m]}_sum', f'{mname[m]}_mean', f'{mname[m]}_min']
    cols += [f'{mname[m]}_sum_{f}' for f in FAMILIES for m in (0, 1)]
    cols += [f'{mname[m]}_sum_robust_z' for m in (0, 1)]
    cols += [f'{mname[m]}_sum_shuf_{x}' for m in (0, 1) for x in ('median', 'mad')]
    cols += [f'{mname[m]}_sum_{f}_robust_z' for f in FAMILIES for m in (0, 1)]
    head = (['seq_name', 'window', 'window_length', 'n_shuffles', 'n_ncrna']
            + cols + ['gc_content_mRNA_window', 'gc_ncRNA_mean', 'gc_diff_mRNA_vs_ncRNA'])

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(raw_dir, exist_ok=True)
    out_csv = Path(results_dir) / f'results_{tranche}.csv'
    tmp_csv = out_csv.with_suffix('.csv.tmp')
    raw_path = Path(raw_dir) / f'{tranche}_shuffle_sums.csv.gz'
    with open(tmp_csv, 'w', newline='') as fo, gzip.open(str(raw_path) + '.tmp', 'wt') as fr:
        w = csv.writer(fo, lineterminator='\n')
        w.writerow(head)
        fr.write('seq_name,iteration,G_total_sum,G_int_sum\n')
        for g, seqs in seqsets:
            a, n = starts[g], len(seqs) - 1
            row = {'seq_name': g, 'window': window, 'window_length': len(seqs[0]),
                   'n_shuffles': n, 'n_ncrna': len(ncs)}
            for m in (0, 1):
                tot = sums[(m, 'all')]
                nat = round(tot[a], 2)
                shuf = [round(x, 2) for x in tot[a + 1:a + 1 + n]]
                med, mad = robust_stats(shuf)
                row[f'{mname[m]}_sum'] = nat
                row[f'{mname[m]}_mean'] = nat / len(ncs)
                row[f'{mname[m]}_min'] = native_min[m][g]
                row[f'{mname[m]}_sum_robust_z'] = robust_z(nat, med, mad)
                row[f'{mname[m]}_sum_shuf_median'] = med
                row[f'{mname[m]}_sum_shuf_mad'] = mad
                for f in FAMILIES:
                    fs = sums[(m, f)]
                    fnat = round(fs[a], 2)
                    fmed, fmad = robust_stats([round(x, 2) for x in fs[a + 1:a + 1 + n]])
                    row[f'{mname[m]}_sum_{f}'] = fnat
                    row[f'{mname[m]}_sum_{f}_robust_z'] = robust_z(fnat, fmed, fmad)
            for it in range(1, n + 1):
                fr.write(f'{g},{it},{sums[(0, "all")][a + it]:.2f},{sums[(1, "all")][a + it]:.2f}\n')
            gc = gc_fraction(seqs[0])
            row.update(gc_content_mRNA_window=gc, gc_ncRNA_mean=nc_gc,
                       gc_diff_mRNA_vs_ncRNA=gc - nc_gc)
            w.writerow([row[c] if isinstance(row[c], (str, int)) else fmt(row[c]) for c in head])
    os.replace(str(raw_path) + '.tmp', raw_path)
    os.replace(tmp_csv, out_csv)   # results file last: its existence means "gene done"
    return out_csv


def main():
    if len(sys.argv) != 2:
        sys.exit('Usage: avoidance.py <window_fasta>')
    fasta = Path(sys.argv[1])
    tranche = fasta.stem         # = manifest seqname for the per-sequence flow
    env = os.environ
    results_dir, tool_dir = env['RESULTS_DIR'], Path(env['TOOL_DIR'])
    n_shuf = int(env.get('N_SHUFFLES', '100'))
    min_valid = n_shuf * int(env.get('MIN_VALID_PERCENT', '90')) // 100
    rnaup_args = shlex.split(env.get('RNAUP_ARGS') or DEFAULT_RNAUP_ARGS)
    region = env.get('AVOID_REGION', '')
    if region not in WINDOW_NAMES:
        sys.exit(f'{tranche}: region {region or "<unset>"} is not scored; only {sorted(WINDOW_NAMES)}')
    window = WINDOW_NAMES[region]

    ncs = load_ncrnas(env['NCRNA_FASTA'], env['NCRNA_FAMILIES'])
    genes = read_fasta(fasta)

    bad = [(g, len(q)) for g, q in genes
           if not length_ok(region, len(q), env.get('TOOL_REQUIRE_LEN', DEFAULT_REQUIRE_LEN))]
    if bad:
        sys.exit(f'{tranche}: window length violates TOOL_REQUIRE_LEN for {region}: {bad}')

    ckpt = Path(env.get('AVOID_CKPT_DIR') or tool_dir / 'ckpt' / tranche)
    ckpt.mkdir(parents=True, exist_ok=True)
    scratch = env.get('TMPDIR') or env.get('TMP_DIR')
    os.makedirs(scratch, exist_ok=True)

    seqsets = build_sequence_set(genes, env['ESL_SHUFFLE_BIN'], n_shuf, min_valid, ckpt)
    if not seqsets:
        sys.exit(f'{tranche}: no genes with enough valid shuffles')
    targets = [(f's{i}', s) for i, s in enumerate(q for _, seqs in seqsets for q in seqs)]
    log.info('%s: %d genes, %d sequences per ncRNA, %d ncRNAs', tranche, len(seqsets),
             len(targets), len(ncs))

    # Work units = (ncRNA, chunk of targets), longest ncRNA first. RNAup time grows about linearly
    # with ncRNA length (72 nt ~0.06 s/target, 341 nt ~0.37 s), so one process per ncRNA leaves most
    # cores idle behind the longest one; many small units keep any core count busy.
    chunk = int(env.get('AVOID_CHUNK', '101'))
    n_chunks = -(-len(targets) // chunk)
    units = [(nc, ci) for nc in sorted(ncs, key=lambda x: -len(x[1])) for ci in range(n_chunks)]
    workers = max(1, int(env.get('SLURM_CPUS_PER_TASK', '1')))
    log.info('%d work units (%d ncRNAs x %d chunks of <= %d) on %d worker(s)',
             len(units), len(ncs), n_chunks, chunk, workers)
    step = max(1, len(units) // 20)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(unit_checkpoint, ckpt, nc, chunk, ci, targets, env['RNAUP_BIN'], rnaup_args, scratch)
                for nc, ci in units]
        for k, fut in enumerate(as_completed(futs), 1):
            fut.result()
            if k % step == 0 or k == len(units):
                log.info('%d/%d units done', k, len(units))

    out = finalize(tranche, seqsets, ncs, ckpt, chunk, n_chunks, results_dir, tool_dir / 'raw_shuffles', window)
    if env.get('AVOID_KEEP_CKPT') != '1':
        shutil.rmtree(ckpt, ignore_errors=True)
    log.info('wrote %s', out)


if __name__ == '__main__':
    main()
