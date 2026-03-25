# FASTA Subset Extractor

Extract proteins from a FASTA database by accession list, keyword, or length range.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Extract by accession list
python fasta_subset_extractor.py --input db.fasta --accessions list.txt --output subset.fasta

# Extract by keyword
python fasta_subset_extractor.py --input db.fasta --keyword "Homo sapiens" --output subset.fasta

# Extract by length range
python fasta_subset_extractor.py --input db.fasta --min-length 50 --max-length 500 --output subset.fasta

# Combine filters
python fasta_subset_extractor.py --input db.fasta --keyword "kinase" --min-length 100 --output subset.fasta
```
