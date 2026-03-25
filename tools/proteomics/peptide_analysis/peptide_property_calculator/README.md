# Peptide Property Calculator

Calculate physicochemical properties of peptide sequences: pI, GRAVY, charge at pH, instability index, amino acid composition.

## Usage

```bash
python peptide_property_calculator.py --sequence PEPTIDEK --ph 7.0
python peptide_property_calculator.py --sequence PEPTIDEK --output properties.json
python peptide_property_calculator.py --input peptides.tsv --output properties.tsv
```

## Input

- `--sequence`: Single peptide sequence
- `--input`: TSV file with a `sequence` column
- `--ph`: pH for charge calculation (default: 7.0)
- `--output`: Output file (.json or .tsv)
