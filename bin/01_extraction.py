#!/usr/bin/env python3
"""
Extract genomic regions (mRNA, CDS, 5UTR, 3UTR, start/stop codon regions, tail)
from a GFF + genome FASTA pair, for a list of target gene IDs.

Outputs:
  - extracted_<region>.fa  (multifasta per region)
  - manifest.tsv           (canonical metadata table — used by all downstream steps)
  - extraction_summary.csv (per-gene QC log)
  - run_manifest.yaml      (run-level reproducibility metadata)
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

# --- Setup Default Logging ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S'
)

# --- Dependency Checks ---
try:
    import yaml
    from Bio import SeqIO
except ImportError as e:
    logging.error(f"Missing Python dependency - {e}")
    logging.info(f"Active Python: {sys.executable}")
    logging.info("Ensure you have activated the correct conda environment before running.")
    sys.exit(1)

def check_dependencies():
    """Ensure external command line tools are available."""
    if shutil.which("gffread") is None:
        logging.error("'gffread' command not found.")
        logging.info("Please install via conda: conda install -c bioconda gffread")
        sys.exit(1)

# --- Data Structures ---
class TranscriptSelection(NamedTuple):
    stripped_id: str
    orig_id: str
    reason: str

# --- Functions ---

def normalise_gff_id(raw: str) -> str:
    """Strip namespace prefixes (e.g. 'gene:', 'transcript:') and version suffixes."""
    return raw.split(':')[-1].split('.')[0]

def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def load_gene_ids(file_path):
    with open(file_path, 'r') as f:
        return set(normalise_gff_id(line.strip()) for line in f if line.strip())

def is_extraction_current(out_dir, manifest_path, source_files, requested_regions):
    """Skip extraction if outputs exist and are newer than all source files."""
    if not os.path.exists(manifest_path):
        return False
    for region in requested_regions:
        mfa = os.path.join(out_dir, f"extracted_{region}.fa")
        if not os.path.exists(mfa):
            return False
    manifest_mtime = os.path.getmtime(manifest_path)
    for src in source_files:
        if os.path.exists(src) and os.path.getmtime(src) > manifest_mtime:
            return False
    return True

def write_run_manifest(out_dir, config, requested_genes, priority_tags):
    manifest = {
        "run_timestamp": datetime.datetime.now().isoformat(),
        "python_executable": sys.executable,
        "genome_fasta": config['inputs']['genome_fasta'],
        "annotation_gff": config['inputs']['annotation_gff'],
        "gene_list": config['inputs']['gene_list'],
        "n_requested_genes": len(requested_genes),
        "isoform_priority_tags": priority_tags,
        "output_mode": config.get('outputs', {}).get('output_mode', 'multifasta'),
        "min_utr_length": config.get('outputs', {}).get('min_utr_length', 30),
        "filter_short_utrs": config.get('outputs', {}).get('filter_short_utrs', False),
        "codon_flank_length": config.get('outputs', {}).get('codon_flank_length', 30),
        "tail_length": config.get('outputs', {}).get('tail_length', 100),
        "regions_to_extract": config.get('processing', {}).get('regions_to_extract', []),
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

                # Try gene_id= first (human MANE style), fall back to Parent= (Ensembl style)
                g_match = re.search(r'gene_id=([^; \n]+)', line)
                if not g_match:
                    g_match = re.search(r'Parent=([^; \n]+)', line)
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
    """Pass 1b: Calculate CDS lengths exclusively for our valid transcripts."""
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
    selected_transcripts = {}

    for g_id in requested_genes:
        transcripts = gene_to_tx.get(g_id, [])
        if not transcripts:
            continue

        best_tx = None
        reason = "longest_cds"

        for tag in priority_tags:
            tagged_txs = [t for t in transcripts if tag in tx_data[t]['tags']]
            if tagged_txs:
                best_tx = max(tagged_txs, key=lambda t: tx_data[t]['cds_len'])
                reason = tag
                break

        if not best_tx:
            best_tx = max(transcripts, key=lambda t: tx_data[t]['cds_len'])

        selected_transcripts[g_id] = TranscriptSelection(
            stripped_id=best_tx,
            orig_id=tx_data[best_tx]['tx_id'],
            reason=reason
        )

    return selected_transcripts

def write_filtered_gff(input_gff, output_gff, selected_transcripts):
    logging.info("Writing filtered GFF...")
    keep_tx_ids = set(tx.stripped_id for tx in selected_transcripts.values())

    written = 0
    with open(input_gff, 'r') as f_in, open(output_gff, 'w') as f_out:
        for line in f_in:
            if line.startswith("#"):
                continue

            match = re.search(r'(?:transcript_id=|Parent=)([^; \n]+)', line)
            if match:
                t_id = normalise_gff_id(match.group(1))
                if t_id in keep_tx_ids:
                    f_out.write(line)
                    written += 1

    logging.info(f"Wrote {written} matching lines to temporary GFF.")

def extract_and_slice_sequences(temp_fa, config, selected_transcripts):
    """Parse sequences from gffread, slice regions, write multifasta + manifest rows."""
    logging.info("Parsing sequences and slicing regions...")

    outputs_cfg = config.get('outputs', {})
    proc_cfg = config.get('processing', {})

    out_dir = outputs_cfg.get('output_dir', 'extracted_regions')
    mode = outputs_cfg.get('output_mode', 'multifasta')
    min_utr = outputs_cfg.get('min_utr_length', 30)
    filter_short = outputs_cfg.get('filter_short_utrs', False)
    flank = outputs_cfg.get('codon_flank_length', 30)
    tail_len = outputs_cfg.get('tail_length', 100)
    requested_regions = list(proc_cfg.get('regions_to_extract', []))

    os.makedirs(out_dir, exist_ok=True)
    logs = []
    manifest_rows = []
    processed_genes = set()

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
                "5UTR_Length": 0, "3UTR_Length": 0, "Status": "Success"
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

            regions = {
                "mRNA": sequence,
                "CDS": sequence[start_idx:end_idx],
                "5UTR": sequence[:start_idx],
                "3UTR": sequence[end_idx:],
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

            logs.append(log_entry)

    return logs, manifest_rows, processed_genes

def write_manifest_tsv(manifest_rows, out_dir):
    """Write the canonical metadata table used by all downstream steps."""
    manifest_path = os.path.join(out_dir, "manifest.tsv")
    if not manifest_rows:
        logging.warning("No manifest rows to write.")
        return

    keys = ["seqname", "gene_id", "transcript_id", "region",
            "length", "selection_reason", "short_utrs"]
    with open(manifest_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
        w.writeheader()
        w.writerows(manifest_rows)
    logging.info(f"Manifest written: {manifest_path} ({len(manifest_rows)} records)")

def write_summary_log(logs, requested_genes, selected_transcripts, processed_genes, out_dir):
    logging.info("Compiling extraction summary...")
    log_path = os.path.join(out_dir, "extraction_summary.csv")

    missing_from_gff = requested_genes - set(selected_transcripts.keys())
    for g_id in missing_from_gff:
        logs.append({"Gene_ID": g_id, "Transcript_ID": "None",
                     "Selection_Reason": "Not_Found_in_GFF",
                     "Total_Length": 0, "CDS_Length": 0,
                     "5UTR_Length": 0, "3UTR_Length": 0, "Status": "Missing"})

    missing_cds = set(selected_transcripts.keys()) - processed_genes
    for g_id in missing_cds:
        tx_info = selected_transcripts[g_id]
        logs.append({"Gene_ID": g_id, "Transcript_ID": tx_info.orig_id,
                     "Selection_Reason": tx_info.reason,
                     "Total_Length": 0, "CDS_Length": 0,
                     "5UTR_Length": 0, "3UTR_Length": 0, "Status": "No_CDS_in_FASTA"})

    keys = ["Gene_ID", "Transcript_ID", "Selection_Reason",
            "Total_Length", "CDS_Length", "5UTR_Length", "3UTR_Length", "Status"]
    with open(log_path, 'w', newline='') as f:
        dict_writer = csv.DictWriter(f, fieldnames=keys)
        dict_writer.writeheader()
        dict_writer.writerows(logs)

    logging.info(f"Summary saved to: {log_path}")

    # Print quick QC summary to stdout
    status_counts = defaultdict(int)
    for entry in logs:
        status_counts[entry["Status"]] += 1
    logging.info("Status breakdown: " +
                 ", ".join(f"{k}={v}" for k, v in sorted(status_counts.items())))

def main():
    check_dependencies()

    parser = argparse.ArgumentParser(description="Extract genomic regions from GFF/FASTA.")
    parser.add_argument("config", help="Path to config.yaml")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable DEBUG logging")
    parser.add_argument("--force", action="store_true",
                        help="Re-extract even if outputs are current")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    config = load_config(args.config)

    REQUIRED_KEYS = [('inputs', 'genome_fasta'), ('inputs', 'annotation_gff'),
                     ('inputs', 'gene_list')]
    for section, key in REQUIRED_KEYS:
        if key not in config.get(section, {}):
            logging.error(f"Missing required config key: '{section}.{key}'.")
            sys.exit(1)

    genome_fa = config['inputs']['genome_fasta']
    anno_gff = config['inputs']['annotation_gff']
    genes_file = config['inputs']['gene_list']
    priority_tags = config.get('processing', {}).get('isoform_priority', [])
    requested_regions = config.get('processing', {}).get('regions_to_extract', [])

    out_dir = config.get('outputs', {}).get('output_dir', 'extracted_regions')
    os.makedirs(out_dir, exist_ok=True)

    manifest_path = os.path.join(out_dir, "manifest.tsv")

    # Idempotency check
    if not args.force and is_extraction_current(
            out_dir, manifest_path,
            [genome_fa, anno_gff, genes_file],
            requested_regions):
        logging.info(f"Outputs in '{out_dir}' are current — skipping extraction.")
        logging.info("Use --force to re-extract.")
        return

    requested_genes = load_gene_ids(genes_file)
    logging.info(f"Loaded {len(requested_genes)} target gene IDs.")

    if not priority_tags:
        logging.warning("No isoform_priority tags defined. Falling back to longest CDS.")

    valid_txs, gene_to_tx, tx_data = map_requested_transcripts(
        anno_gff, requested_genes, priority_tags)
    update_cds_lengths_in_place(anno_gff, valid_txs, tx_data)
    selected_txs = select_best_transcripts(
        requested_genes, gene_to_tx, tx_data, priority_tags)
    logging.info(f"Found suitable transcripts for {len(selected_txs)} genes.")

    fd_gff, temp_gff = tempfile.mkstemp(suffix=".gff")
    fd_fa, temp_fa = tempfile.mkstemp(suffix=".fa")
    os.close(fd_gff)
    os.close(fd_fa)

    try:
        write_filtered_gff(anno_gff, temp_gff, selected_txs)

        logging.info("Running gffread to extract base mRNAs...")
        try:
            subprocess.run(['gffread', '-w', temp_fa, '-g', genome_fa, temp_gff],
                           check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"gffread failed:\n{e.stderr}")
            sys.exit(1)

        logs, manifest_rows, processed_genes = extract_and_slice_sequences(
            temp_fa, config, selected_txs)
        write_manifest_tsv(manifest_rows, out_dir)
        write_summary_log(logs, requested_genes, selected_txs, processed_genes, out_dir)

    finally:
        if os.path.exists(temp_gff): os.remove(temp_gff)
        if os.path.exists(temp_fa): os.remove(temp_fa)

    write_run_manifest(out_dir, config, requested_genes, priority_tags)
    logging.info("Pipeline complete.")

if __name__ == "__main__":
    main()
