#!/usr/bin/env python3
"""
01_extract.py
Extract genomic regions (mRNA, CDS, 5UTR, 3UTR, start/stop codon regions, tail,
UTR_pair) from a GFF + genome FASTA pair, for a list of target gene IDs.

Outputs (under $RUNS_ROOT/<dataset>/extracted_regions/):
  - extracted_<region>.fa     (multifasta per region)
  - manifest.tsv              (canonical metadata table — used by all downstream steps)
  - extraction_summary.csv    (per-gene QC log)
  - run_manifest.yaml         (run-level reproducibility metadata)
  - utr_pair_geometry.tsv     (ONLY if UTR_pair is requested — per-record
                               5UTR/linker/3UTR lengths)

Usage:
  ./bin/01_extract.py --dataset human_liver
  ./bin/01_extract.py -d human_liver --force
"""
import sys
import os
import re
import csv
import subprocess
import shutil
import tempfile
import logging
import argparse
import datetime
from collections import defaultdict
from contextlib import ExitStack
from typing import NamedTuple

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S'
)

try:
    import yaml
    from Bio import SeqIO
except ImportError as e:
    logging.error(f"Missing Python dependency - {e}")
    logging.info(f"Active Python: {sys.executable}")
    logging.info("Ensure you have activated the correct conda environment before running.")
    sys.exit(1)


def check_dependencies():
    if shutil.which("gffread") is None:
        logging.error("'gffread' command not found.")
        logging.info("Please install via conda: conda install -c bioconda gffread")
        sys.exit(1)


class TranscriptSelection(NamedTuple):
    stripped_id: str
    orig_id: str
    reason: str


def normalise_gff_id(raw: str) -> str:
    """Strip namespace prefixes (e.g. 'gene:', 'transcript:') and version suffixes."""
    return raw.split(':')[-1].split('.')[0]


def load_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def load_gene_ids(path):
    with open(path, 'r') as f:
        return set(normalise_gff_id(line.strip()) for line in f if line.strip())


def is_extraction_current(out_dir, manifest_path, source_files, requested_regions):
    """Skip extraction if outputs exist and are newer than all source files."""
    if not os.path.exists(manifest_path):
        return False
    for region in requested_regions:
        if not os.path.exists(os.path.join(out_dir, f"extracted_{region}.fa")):
            return False
    # If UTR_pair was requested, sidecar must also be present
    if "UTR_pair" in requested_regions:
        if not os.path.exists(os.path.join(out_dir, "utr_pair_geometry.tsv")):
            return False
    manifest_mtime = os.path.getmtime(manifest_path)
    for src in source_files:
        if os.path.exists(src) and os.path.getmtime(src) > manifest_mtime:
            return False
    return True


def write_run_manifest(out_dir, dataset_name, dataset_yaml_path, config,
                       requested_genes, priority_tags):
    extr = config.get('extraction', {})
    utr_pair = extr.get('utr_pair', {})
    manifest = {
        "dataset_name": dataset_name,
        "dataset_yaml": dataset_yaml_path,
        "run_timestamp": datetime.datetime.now().isoformat(),
        "python_executable": sys.executable,
        "genome_fasta": config['inputs']['genome_fasta'],
        "annotation_gff": config['inputs']['annotation_gff'],
        "gene_list": config['inputs']['gene_list'],
        "n_requested_genes": len(requested_genes),
        "isoform_priority_tags": priority_tags,
        "output_mode": extr.get('output_mode', 'multifasta'),
        "min_utr_length": extr.get('min_utr_length', 30),
        "filter_short_utrs": extr.get('filter_short_utrs', False),
        "codon_flank_length": extr.get('codon_flank_length', 30),
        "tail_length": extr.get('tail_length', 100),
        "regions_to_extract": extr.get('regions_to_extract', []),
        "utr_pair_linker_length": utr_pair.get('linker_length', 7),
        "utr_pair_linker_char": utr_pair.get('linker_char', 'N'),
        "utr_pair_constraint_style": utr_pair.get('constraint_style', 'cross_pair'),
    }
    path = os.path.join(out_dir, "run_manifest.yaml")
    with open(path, 'w') as f:
        yaml.dump(manifest, f, default_flow_style=False)
    logging.info(f"Run manifest saved to: {path}")


def map_requested_transcripts(gff_path, requested_genes, priority_tags):
    logging.info("Scanning GFF for requested gene transcripts...")
    valid_txs = set()
    gene_to_tx = defaultdict(list)
    tx_data = {}

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue

            if "\ttranscript\t" in line or "\tmRNA\t" in line:
                t_match = re.search(r'transcript_id=([^; \n]+)', line)
                if not t_match:
                    continue

                # gene_id= first (MANE style), Parent= fallback (Ensembl style)
                g_match = re.search(r'gene_id=([^; \n]+)', line) \
                          or re.search(r'Parent=([^; \n]+)', line)
                if not g_match:
                    continue

                g_id = normalise_gff_id(g_match.group(1))
                if g_id not in requested_genes:
                    continue

                t_id = normalise_gff_id(t_match.group(1))
                orig_t_id = t_match.group(1).split(':')[-1]
                valid_txs.add(t_id)
                gene_to_tx[g_id].append(t_id)
                tx_data[t_id] = {'tags': set(), 'cds_len': 0, 'tx_id': orig_t_id}

                for tag in priority_tags:
                    if re.search(r'\b' + re.escape(tag) + r'\b', line):
                        tx_data[t_id]['tags'].add(tag)

    return valid_txs, gene_to_tx, tx_data


def update_cds_lengths_in_place(gff_path, valid_txs, tx_data):
    logging.info("Calculating CDS lengths for target transcripts...")
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            if "\tCDS\t" in line:
                parent_match = re.search(r'Parent=([^; \n]+)', line)
                if parent_match:
                    p_id = normalise_gff_id(parent_match.group(1))
                    if p_id in valid_txs:
                        parts = line.split('\t')
                        tx_data[p_id]['cds_len'] += int(parts[4]) - int(parts[3]) + 1


def select_best_transcripts(requested_genes, gene_to_tx, tx_data, priority_tags):
    logging.info("Selecting best transcript per gene...")
    selected = {}
    for g_id in requested_genes:
        transcripts = gene_to_tx.get(g_id, [])
        if not transcripts:
            continue

        best_tx = None
        reason = "longest_cds"
        for tag in priority_tags:
            tagged = [t for t in transcripts if tag in tx_data[t]['tags']]
            if tagged:
                best_tx = max(tagged, key=lambda t: tx_data[t]['cds_len'])
                reason = tag
                break
        if not best_tx:
            best_tx = max(transcripts, key=lambda t: tx_data[t]['cds_len'])

        selected[g_id] = TranscriptSelection(
            stripped_id=best_tx,
            orig_id=tx_data[best_tx]['tx_id'],
            reason=reason
        )
    return selected


def write_filtered_gff(input_gff, output_gff, selected_transcripts):
    logging.info("Writing filtered GFF...")
    keep_tx_ids = set(tx.stripped_id for tx in selected_transcripts.values())

    written = 0
    with open(input_gff, 'r') as f_in, open(output_gff, 'w') as f_out:
        for line in f_in:
            if line.startswith("#"):
                continue
            match = re.search(r'(?:transcript_id=|Parent=)([^; \n]+)', line)
            if match and normalise_gff_id(match.group(1)) in keep_tx_ids:
                f_out.write(line)
                written += 1
    logging.info(f"Wrote {written} matching lines to temporary GFF.")


def build_utr_pair_record(seq_5utr, seq_3utr, linker_char, linker_len):
    """Construct a UTR_pair hybrid sequence and ViennaRNA cross-pair constraint.

    Returns (hybrid_seq, constraint_str) where:
      hybrid_seq    = 5UTR + linker + 3UTR
      constraint_str = '<' * |5UTR| + 'x' * linker_len + '>' * |3UTR|

    The '<' / '>' constraint forces base pairs to form only between the two
    UTRs (RNAfold -C semantics); 'x' forces the linker positions unpaired.
    """
    linker = linker_char * linker_len
    hybrid = f"{seq_5utr}{linker}{seq_3utr}"
    constraint = ('<' * len(seq_5utr)
                  + 'x' * linker_len
                  + '>' * len(seq_3utr))
    return hybrid, constraint


def extract_and_slice_sequences(temp_fa, out_dir, extraction_cfg, selected_transcripts):
    """Parse sequences from gffread, slice regions, write multifasta + manifest rows.

    Returns (logs, manifest_rows, processed_genes, utr_pair_geometry_rows).
    """
    logging.info("Parsing sequences and slicing regions...")

    mode = extraction_cfg.get('output_mode', 'multifasta')
    min_utr = extraction_cfg.get('min_utr_length', 30)
    filter_short = extraction_cfg.get('filter_short_utrs', False)
    flank = extraction_cfg.get('codon_flank_length', 30)
    tail_len = extraction_cfg.get('tail_length', 100)
    requested_regions = list(extraction_cfg.get('regions_to_extract', []))
    want_utr_pair = "UTR_pair" in requested_regions

    utr_pair_cfg = extraction_cfg.get('utr_pair', {}) or {}
    utr_pair_linker_len = int(utr_pair_cfg.get('linker_length', 7))
    utr_pair_linker_char = str(utr_pair_cfg.get('linker_char', 'N'))[:1] or 'N'

    if want_utr_pair and utr_pair_linker_len < 1:
        logging.warning("UTR_pair: linker_length < 1; forcing to 1.")
        utr_pair_linker_len = 1

    os.makedirs(out_dir, exist_ok=True)
    logs, manifest_rows, processed_genes = [], [], set()
    utr_pair_geometry_rows = []

    tx_lookup = {tx.stripped_id: (g_id, tx) for g_id, tx in selected_transcripts.items()}

    with ExitStack() as stack:
        multifasta_files = {}
        if mode == "multifasta":
            for region in requested_regions:
                multifasta_files[region] = stack.enter_context(
                    open(os.path.join(out_dir, f"extracted_{region}.fa"), 'w')
                )

        for record in SeqIO.parse(temp_fa, "fasta"):
            tx_id_stripped = normalise_gff_id(record.id)
            sequence = str(record.seq)
            seq_len = len(sequence)

            if tx_id_stripped not in tx_lookup:
                continue

            g_id, tx_info = tx_lookup[tx_id_stripped]
            processed_genes.add(g_id)

            log_entry = {
                "Gene_ID": g_id, "Transcript_ID": tx_info.orig_id,
                "Selection_Reason": tx_info.reason,
                "Total_Length": seq_len, "CDS_Length": 0,
                "5UTR_Length": 0, "3UTR_Length": 0,
                "UTR_pair_Length": 0, "Status": "Success"
            }

            match = re.search(r'CDS=(\d+)-(\d+)', record.description)
            if not match:
                log_entry["Status"] = "No_CDS"
                logs.append(log_entry)
                continue

            cds_start = int(match.group(1))
            end_idx = int(match.group(2))
            start_idx = cds_start - 1

            log_entry["5UTR_Length"] = start_idx
            log_entry["3UTR_Length"] = seq_len - end_idx
            log_entry["CDS_Length"] = end_idx - start_idx

            short_utrs = (log_entry["5UTR_Length"] < min_utr or
                          log_entry["3UTR_Length"] < min_utr)
            if short_utrs:
                log_entry["Status"] = "Short_or_Missing_UTRs"
                if filter_short:
                    logs.append(log_entry)
                    continue

            seq_5utr = sequence[:start_idx]
            seq_3utr = sequence[end_idx:]

            regions = {
                "mRNA": sequence,
                "CDS": sequence[start_idx:end_idx],
                "5UTR": seq_5utr,
                "3UTR": seq_3utr,
                "start_codon_region": sequence[max(0, start_idx - flank):
                                                min(seq_len, start_idx + 3 + flank)],
                "stop_codon_region": sequence[max(0, end_idx - 3 - flank):
                                               min(seq_len, end_idx + flank)],
                "tail_region": sequence[-tail_len:] if seq_len >= tail_len else sequence
            }

            for region_name, seq_string in regions.items():
                if region_name not in requested_regions or not seq_string:
                    continue

                seqname = f"{g_id}_{tx_info.orig_id}_{region_name}"
                header = f">{seqname}"

                if mode == "multifasta":
                    multifasta_files[region_name].write(f"{header}\n{seq_string}\n")
                else:
                    with open(os.path.join(out_dir, f"{seqname}.fa"), 'w') as f:
                        f.write(f"{header}\n{seq_string}\n")

                manifest_rows.append({
                    "seqname": seqname,
                    "gene_id": g_id,
                    "transcript_id": tx_info.orig_id,
                    "region": region_name,
                    "length": len(seq_string),
                    "selection_reason": tx_info.reason,
                    "short_utrs": "true" if short_utrs else "false",
                })

            # --- UTR_pair: hybrid sequence with cross-pair constraint ---
            # Only emitted if both UTRs are non-empty. Short UTRs are still
            # included unless filter_short_utrs dropped this record above.
            if want_utr_pair and seq_5utr and seq_3utr:
                hybrid_seq, constraint_str = build_utr_pair_record(
                    seq_5utr, seq_3utr,
                    utr_pair_linker_char, utr_pair_linker_len
                )
                seqname = f"{g_id}_{tx_info.orig_id}_UTR_pair"
                # Three-line FASTA: header / sequence / constraint
                record_text = f">{seqname}\n{hybrid_seq}\n{constraint_str}\n"

                if mode == "multifasta":
                    multifasta_files["UTR_pair"].write(record_text)
                else:
                    with open(os.path.join(out_dir, f"{seqname}.fa"), 'w') as f:
                        f.write(record_text)

                manifest_rows.append({
                    "seqname": seqname,
                    "gene_id": g_id,
                    "transcript_id": tx_info.orig_id,
                    "region": "UTR_pair",
                    "length": len(hybrid_seq),
                    "selection_reason": tx_info.reason,
                    "short_utrs": "true" if short_utrs else "false",
                })
                utr_pair_geometry_rows.append({
                    "seqname": seqname,
                    "len_5utr": len(seq_5utr),
                    "linker_len": utr_pair_linker_len,
                    "len_3utr": len(seq_3utr),
                })
                log_entry["UTR_pair_Length"] = len(hybrid_seq)

            logs.append(log_entry)

    return logs, manifest_rows, processed_genes, utr_pair_geometry_rows


def write_manifest_tsv(manifest_rows, out_dir):
    if not manifest_rows:
        logging.warning("No manifest rows to write.")
        return
    keys = ["seqname", "gene_id", "transcript_id", "region",
            "length", "selection_reason", "short_utrs"]
    path = os.path.join(out_dir, "manifest.tsv")
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
        w.writeheader()
        w.writerows(manifest_rows)
    logging.info(f"Manifest written: {path} ({len(manifest_rows)} records)")


def write_utr_pair_geometry(rows, out_dir):
    """Sidecar file: per-UTR_pair-record component lengths.

    Schema:  seqname  len_5utr  linker_len  len_3utr
    """
    if not rows:
        return
    path = os.path.join(out_dir, "utr_pair_geometry.tsv")
    keys = ["seqname", "len_5utr", "linker_len", "len_3utr"]
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    logging.info(f"UTR_pair geometry sidecar written: {path} ({len(rows)} records)")


def write_summary_log(logs, requested_genes, selected, processed_genes, out_dir):
    logging.info("Compiling extraction summary...")
    log_path = os.path.join(out_dir, "extraction_summary.csv")

    missing_from_gff = requested_genes - set(selected.keys())
    for g_id in missing_from_gff:
        logs.append({"Gene_ID": g_id, "Transcript_ID": "None",
                     "Selection_Reason": "Not_Found_in_GFF",
                     "Total_Length": 0, "CDS_Length": 0,
                     "5UTR_Length": 0, "3UTR_Length": 0,
                     "UTR_pair_Length": 0, "Status": "Missing"})

    missing_cds = set(selected.keys()) - processed_genes
    for g_id in missing_cds:
        tx_info = selected[g_id]
        logs.append({"Gene_ID": g_id, "Transcript_ID": tx_info.orig_id,
                     "Selection_Reason": tx_info.reason,
                     "Total_Length": 0, "CDS_Length": 0,
                     "5UTR_Length": 0, "3UTR_Length": 0,
                     "UTR_pair_Length": 0, "Status": "No_CDS_in_FASTA"})

    keys = ["Gene_ID", "Transcript_ID", "Selection_Reason",
            "Total_Length", "CDS_Length", "5UTR_Length", "3UTR_Length",
            "UTR_pair_Length", "Status"]
    with open(log_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(logs)
    logging.info(f"Summary saved to: {log_path}")

    counts = defaultdict(int)
    for entry in logs:
        counts[entry["Status"]] += 1
    logging.info("Status breakdown: " +
                 ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))


def resolve_paths(dataset_name):
    """Pure-Python equivalent of lib/paths.sh's resolve_paths for extraction."""
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    runs_root = os.environ.get('RUNS_ROOT', os.path.join(project_root, 'runs'))

    yaml_path = os.path.join(project_root, 'configs', 'datasets',
                             f"{dataset_name}.yaml")
    if not os.path.isfile(yaml_path):
        logging.error(f"Dataset config not found: {yaml_path}")
        sys.exit(1)

    out_dir = os.path.join(runs_root, dataset_name, 'extracted_regions')
    return yaml_path, out_dir


def main():
    check_dependencies()

    parser = argparse.ArgumentParser(description="Extract genomic regions from GFF/FASTA.")
    parser.add_argument("--dataset", "-d",
                        default=os.environ.get('DATASET'),
                        help="Dataset name (config at configs/datasets/<name>.yaml)")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--force", action="store_true",
                        help="Re-extract even if outputs are current")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not args.dataset:
        logging.error("--dataset (or DATASET env var) required")
        sys.exit(1)

    yaml_path, out_dir = resolve_paths(args.dataset)
    config = load_config(yaml_path)

    REQUIRED = [('inputs', 'genome_fasta'), ('inputs', 'annotation_gff'),
                ('inputs', 'gene_list')]
    for section, key in REQUIRED:
        if key not in config.get(section, {}):
            logging.error(f"Missing required config key: '{section}.{key}'.")
            sys.exit(1)

    genome_fa = config['inputs']['genome_fasta']
    anno_gff = config['inputs']['annotation_gff']
    genes_file = config['inputs']['gene_list']

    extr_cfg = config.get('extraction', {})
    priority_tags = extr_cfg.get('isoform_priority', [])
    requested_regions = extr_cfg.get('regions_to_extract', [])

    os.makedirs(out_dir, exist_ok=True)
    manifest_path = os.path.join(out_dir, "manifest.tsv")

    if not args.force and is_extraction_current(
            out_dir, manifest_path,
            [genome_fa, anno_gff, genes_file, yaml_path],
            requested_regions):
        logging.info(f"Outputs in '{out_dir}' are current — skipping extraction.")
        logging.info("Use --force to re-extract.")
        return

    requested_genes = load_gene_ids(genes_file)
    logging.info(f"Dataset: {args.dataset}")
    logging.info(f"Loaded {len(requested_genes)} target gene IDs.")

    if not priority_tags:
        logging.warning("No isoform_priority tags defined. Falling back to longest CDS.")

    valid_txs, gene_to_tx, tx_data = map_requested_transcripts(
        anno_gff, requested_genes, priority_tags)
    update_cds_lengths_in_place(anno_gff, valid_txs, tx_data)
    selected = select_best_transcripts(
        requested_genes, gene_to_tx, tx_data, priority_tags)
    logging.info(f"Found suitable transcripts for {len(selected)} genes.")

    fd_gff, temp_gff = tempfile.mkstemp(suffix=".gff")
    fd_fa, temp_fa = tempfile.mkstemp(suffix=".fa")
    os.close(fd_gff); os.close(fd_fa)

    try:
        write_filtered_gff(anno_gff, temp_gff, selected)

        logging.info("Running gffread to extract base mRNAs...")
        try:
            subprocess.run(['gffread', '-w', temp_fa, '-g', genome_fa, temp_gff],
                           check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"gffread failed:\n{e.stderr}")
            sys.exit(1)

        logs, manifest_rows, processed_genes, utr_pair_geom = \
            extract_and_slice_sequences(temp_fa, out_dir, extr_cfg, selected)
        write_manifest_tsv(manifest_rows, out_dir)
        if "UTR_pair" in requested_regions:
            write_utr_pair_geometry(utr_pair_geom, out_dir)
        write_summary_log(logs, requested_genes, selected, processed_genes, out_dir)

    finally:
        for p in (temp_gff, temp_fa):
            if os.path.exists(p):
                os.remove(p)

    write_run_manifest(out_dir, args.dataset, yaml_path, config,
                       requested_genes, priority_tags)
    logging.info("Extraction complete.")


if __name__ == "__main__":
    main()