# FASTA Cleaner

Clean a FASTA database by removing duplicates, fixing headers, filtering by length, and removing stop codons.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Remove duplicates and filter by length
python fasta_cleaner.py --input messy.fasta --remove-duplicates --min-length 6 --output clean.fasta

# Remove stop codons and fix headers
python fasta_cleaner.py --input messy.fasta --remove-stop-codons --fix-headers --output clean.fasta

# All cleaning operations
python fasta_cleaner.py --input messy.fasta --remove-duplicates --min-length 6 --remove-stop-codons --fix-headers --remove-invalid-chars --output clean.fasta
```
