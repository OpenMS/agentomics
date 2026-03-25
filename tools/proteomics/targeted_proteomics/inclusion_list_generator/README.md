# Inclusion List Generator

Generate instrument inclusion lists from peptide data for targeted MS experiments.

## Usage

```bash
python inclusion_list_generator.py --input peptides.tsv --format thermo --charge 2,3 --output inclusion.csv
python inclusion_list_generator.py --input peptides.tsv --format generic --charge 2 --output inclusion.csv
```

## Formats

- **thermo** - Thermo Scientific instrument format
- **generic** - Generic CSV format
