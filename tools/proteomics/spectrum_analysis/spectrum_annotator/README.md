# Spectrum Annotator

Annotate observed MS2 spectrum peaks with theoretical fragment ion matches for a given peptide sequence.

## Usage

```bash
python spectrum_annotator.py --mz-list "100.5,200.3,300.1" --intensities "1000,500,200" --sequence PEPTIDEK --charge 2 --tolerance 0.02
python spectrum_annotator.py --mz-list "100.5,200.3" --intensities "1000,500" --sequence PEPTIDEK --charge 1 --output annotation.tsv
```
