# N-terminal Modification Annotator

Classify N-terminal peptides as protein N-terminus, signal peptide, neo-N-terminus, etc.

## Usage

```bash
python nterm_modification_annotator.py --input nterm_peptides.tsv --fasta reference.fasta --output annotated.tsv
```

## Input Format

- `nterm_peptides.tsv`: columns `peptide`, `protein`
- `reference.fasta`: Reference proteome FASTA file

## Output

- `annotated.tsv` - Peptides with `nterm_type`, `nterm_modification`, `start_position` columns
