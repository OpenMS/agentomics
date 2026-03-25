# FASTA Taxonomy Splitter

Split a multi-organism FASTA file by taxonomy parsed from headers.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Split by OS= field (UniProt format)
python fasta_taxonomy_splitter.py --input combined.fasta --pattern "OS=([^=]+) OX=" --output-dir split/

# Custom pattern
python fasta_taxonomy_splitter.py --input combined.fasta --pattern "\[(.+?)\]$" --output-dir split/
```
