# Crosslink Classifier

Classify crosslinks as intra-protein, inter-protein, or monolink.

## Usage

```bash
python xl_link_classifier.py --crosslinks links.tsv --fasta proteome.fasta --output classified.tsv
```

## Input Format

- `links.tsv`: columns `peptide1`, `peptide2` (empty peptide2 for monolinks)
- `proteome.fasta`: FASTA file with protein sequences

## Output

- `classified.tsv` - Crosslinks with `link_type`, `proteins1`, `proteins2`, `shared_proteins`
