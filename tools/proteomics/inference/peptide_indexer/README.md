# Peptide Indexer

Map peptide identifications to proteins in a FASTA database.

## Usage

```bash
python peptide_indexer.py --ids peptides.idXML --fasta database.fasta --output indexed.idXML
```

## Requirements

```
pip install pyopenms click
```

## How it works

Wraps `pyopenms.PeptideIndexing` to match peptide sequences against a FASTA protein database, adding protein accession references to each peptide hit. This is a prerequisite step for protein inference.
