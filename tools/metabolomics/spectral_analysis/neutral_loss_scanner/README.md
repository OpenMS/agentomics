# Neutral Loss Scanner

Scan MS2 spectra for characteristic neutral losses from precursor ions in mzML files.

## Usage

```bash
python neutral_loss_scanner.py --input file.mzML --losses 97.977,162.053 --tolerance 0.02
python neutral_loss_scanner.py --input file.mzML --losses 97.977 --tolerance 0.05 --output matches.tsv
```
