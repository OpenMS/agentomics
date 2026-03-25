# FASTA In-Silico Digest Stats

Digest a FASTA database in silico and report peptide statistics.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic trypsin digestion
python fasta_in_silico_digest_stats.py --input db.fasta --enzyme Trypsin --output stats.tsv

# With missed cleavages and length filter
python fasta_in_silico_digest_stats.py --input db.fasta --enzyme Trypsin --missed-cleavages 2 --min-length 7 --max-length 30 --output stats.tsv
```
