# Peptide and Protein Quantifier

Roll up peptide-level quantification to protein-level abundances from an annotated consensusXML file using PeptideAndProteinQuant.

## Usage

```bash
# Default (top 3 peptides, include all proteins)
python peptide_protein_quantifier.py --input annotated.consensusXML --output protein_quant.csv

# Top 5 peptides, exclude proteins with fewer
python peptide_protein_quantifier.py --input annotated.consensusXML --output protein_quant.csv \
    --top-n 5 --no-include-all
```

## Dependencies

```
pyopenms
click
```
