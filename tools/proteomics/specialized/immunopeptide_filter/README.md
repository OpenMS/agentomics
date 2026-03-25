# Immunopeptide Filter

Filter peptides for MHC class I or II binding by length and optional motif.

## Usage

```bash
python immunopeptide_filter.py --input peptides.tsv --class-i --length-range 8-11 --output immunopeptides.tsv
python immunopeptide_filter.py --input peptides.tsv --class-ii --output immunopeptides.tsv
python immunopeptide_filter.py --input peptides.tsv --class-i --motif "^.{1}[LIV]" --output immunopeptides.tsv
```
