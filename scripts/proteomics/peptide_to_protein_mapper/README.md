# Peptide to Protein Mapper

Map peptides to proteins by searching a FASTA database.

## Usage

```bash
python peptide_to_protein_mapper.py --peptides peptides.tsv --fasta db.fasta --output mapped.tsv
```

## Input

- **peptides.tsv** - TSV with a `peptide` column
- **db.fasta** - Protein FASTA database

## Output

TSV with columns: `peptide`, `protein`, `protein_description`, `start`, `end`, `is_unique`
