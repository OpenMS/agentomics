# Phospho Motif Analyzer

Extract amino acid windows around phosphosites and compute position-specific frequencies.

## Usage

```bash
python phospho_motif_analyzer.py --input phosphosites.tsv --fasta proteome.fasta --window 7 --output motifs.tsv
```

## Input Format

- `phosphosites.tsv`: Tab-separated with columns `peptide`, `protein`, `site` (1-based position)
- `proteome.fasta`: FASTA file with protein sequences

## Output

- `motifs.tsv` - Input rows with added `motif_window` column
- `motifs_frequencies.tsv` - Position-specific amino acid frequencies
