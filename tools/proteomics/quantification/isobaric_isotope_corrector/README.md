# Isobaric Isotope Corrector

Correct isotopic impurities in isobaric labeling quantification data using the manufacturer-provided correction matrix.

## Usage

```bash
# Correct TMT 6-plex data
python isobaric_isotope_corrector.py --input quant.consensusXML --method tmt6plex --output corrected.consensusXML

# Correct iTRAQ 4-plex data
python isobaric_isotope_corrector.py --input quant.consensusXML --method itraq4plex --output corrected.consensusXML
```

## Supported Methods

- `tmt6plex` — TMT 6-plex
- `tmt10plex` — TMT 10-plex
- `itraq4plex` — iTRAQ 4-plex
- `itraq8plex` — iTRAQ 8-plex

## Python API

```python
from isobaric_isotope_corrector import correct_isotope_impurities

n_features = correct_isotope_impurities("quant.consensusXML", "tmt6plex", "corrected.consensusXML")
```

## Pipeline

This tool is step 2 in the isobaric quantification workflow:

1. isobaric_channel_extractor (mzML -> consensusXML)
2. **isobaric_isotope_corrector** (consensusXML -> corrected consensusXML)
3. isobaric_quantifier (consensusXML -> quantified consensusXML)
