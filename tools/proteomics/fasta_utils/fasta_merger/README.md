# FASTA Merger

Merge multiple FASTA files with optional deduplication.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Simple merge
python fasta_merger.py --inputs db1.fasta db2.fasta --output merged.fasta

# Merge with deduplication by identifier
python fasta_merger.py --inputs db1.fasta db2.fasta --remove-duplicates --output merged.fasta

# Merge with deduplication by sequence
python fasta_merger.py --inputs db1.fasta db2.fasta --remove-duplicates --dedup-by sequence --output merged.fasta
```
