# Contaminant Database Merger

Append cRAP contaminant proteins to a target FASTA database with a configurable prefix, removing duplicates.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Add built-in cRAP contaminants
python contaminant_database_merger.py --input target.fasta --add-crap --prefix CONT_ --output merged.fasta

# Add custom contaminant file
python contaminant_database_merger.py --input target.fasta --contaminants custom.fasta --output merged.fasta

# Both
python contaminant_database_merger.py --input target.fasta --add-crap --contaminants extra.fasta --output merged.fasta
```
