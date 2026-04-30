# mRNA_stability
Code for experiments and interpretation on mRNA stability


# Pipeline files

### 01_extraction.py  
Requires inputs:
    gene_list, a text file with one gene ID per line that will be extracted
    genome_fasta, the reference genome in FASTA format
    annotation_gff, the annotation file in GFF format

Recommend to output multifasta files, as next step in pipeline will split into individual.

### 02_stratify.sh


## Log

Added pipeline script, 01_extraction.py
Added yaml template for extraction, 01_extraction_template.yaml
