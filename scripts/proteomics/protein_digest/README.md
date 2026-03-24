# Protein In-Silico Digest

Perform in-silico enzymatic digestion of a protein sequence and report
the resulting peptides with their masses.

## Usage

```bash
python protein_digest.py --sequence MKVLWAALLVTFLAGCQAK... --enzyme Trypsin
python protein_digest.py --sequence MKVLWAALLVTFLAGCQAK... --enzyme Lys-C --missed-cleavages 2
python protein_digest.py --list-enzymes
```
