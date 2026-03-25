# Peptide Spectral Match Validator

Validate PSMs by recomputing theoretical fragment ions and measuring coverage.

## Usage

```bash
python peptide_spectral_match_validator.py --mzml run.mzML --peptides psms.tsv --output validation.tsv
python peptide_spectral_match_validator.py --mzml run.mzML --peptides psms.tsv --tolerance 0.05 --output validation.tsv
```
