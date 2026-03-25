# FASTA Statistics Reporter

Report comprehensive statistics from a FASTA protein database.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic statistics
python fasta_statistics_reporter.py --input db.fasta

# Include tryptic peptide count
python fasta_statistics_reporter.py --input db.fasta --enzyme Trypsin --output stats.json

# With missed cleavages
python fasta_statistics_reporter.py --input db.fasta --enzyme Trypsin --missed-cleavages 2 --output stats.json
```
