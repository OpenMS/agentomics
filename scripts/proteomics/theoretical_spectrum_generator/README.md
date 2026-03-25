# Theoretical Spectrum Generator

Generate theoretical b/y/a/c/x/z fragment ion spectra for a peptide sequence with annotated TSV output.

## Usage

```bash
python theoretical_spectrum_generator.py --sequence PEPTIDEK --charge 2 --ion-types b,y
python theoretical_spectrum_generator.py --sequence PEPTIDEK --charge 1 --ion-types b,y,a --add-losses --output fragments.tsv
```
