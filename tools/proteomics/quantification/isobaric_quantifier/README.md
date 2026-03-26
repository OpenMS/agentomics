# Isobaric Quantifier

Quantify isobaric labeling experiments (TMT, iTRAQ) with optional isotope correction and normalization.

## Usage

```bash
# Quantify TMT 6-plex data
python isobaric_quantifier.py --input quant.consensusXML --method tmt6plex --output quantified.consensusXML

# Quantify iTRAQ 4-plex data
python isobaric_quantifier.py --input quant.consensusXML --method itraq4plex --output quantified.consensusXML
```

## Supported Methods

- `tmt6plex` — TMT 6-plex
- `tmt10plex` — TMT 10-plex
- `tmt16plex` — TMT 16-plex
- `itraq4plex` — iTRAQ 4-plex
- `itraq8plex` — iTRAQ 8-plex

## Python API

```python
from isobaric_quantifier import quantify_isobaric

n_features = quantify_isobaric("quant.consensusXML", "tmt6plex", "quantified.consensusXML")
```

## Pipeline

This tool is step 3 in the isobaric quantification workflow:

1. isobaric_channel_extractor (mzML -> consensusXML)
2. isobaric_isotope_corrector (consensusXML -> corrected consensusXML)
3. **isobaric_quantifier** (consensusXML -> quantified consensusXML)
