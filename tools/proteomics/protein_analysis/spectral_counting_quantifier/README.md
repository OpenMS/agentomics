# Spectral Counting Quantifier

Calculate protein abundances from spectral counts using emPAI and NSAF methods.

## Usage

```bash
python spectral_counting_quantifier.py --input peptide_counts.tsv --fasta db.fasta --method nsaf --output abundances.tsv
python spectral_counting_quantifier.py --input peptide_counts.tsv --fasta db.fasta --method empai --output abundances.tsv
```

## Input Format

TSV file with columns: `protein`, `peptide`, `spectral_count`
