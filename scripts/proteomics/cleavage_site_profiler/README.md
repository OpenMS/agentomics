# Cleavage Site Profiler

Extract P4-P4' windows around cleavage sites from neo-N-terminal peptides and compute position-specific frequencies.

## Usage

```bash
python cleavage_site_profiler.py --input neo_nterm.tsv --fasta reference.fasta --window 4 --output profile.tsv
```

## Input Format

- `neo_nterm.tsv`: columns `peptide`, `protein`
- `reference.fasta`: Reference proteome FASTA file

## Output

- `profile.tsv` - Peptides with cleavage windows
- `profile_frequencies.tsv` - Position-specific amino acid frequencies (P4..P1, P1'..P4')
