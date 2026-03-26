# Isobaric Channel Extractor

Extract reporter-ion channels from isobaric labeling experiments (TMT, iTRAQ).

## Usage

```bash
# TMT 6-plex extraction
python isobaric_channel_extractor.py --input run.mzML --method tmt6plex --output quant.consensusXML

# iTRAQ 4-plex extraction
python isobaric_channel_extractor.py --input run.mzML --method itraq4plex --output quant.consensusXML
```

## Supported Methods

- `tmt6plex` — TMT 6-plex
- `tmt10plex` — TMT 10-plex
- `tmt16plex` — TMT 16-plex
- `itraq4plex` — iTRAQ 4-plex
- `itraq8plex` — iTRAQ 8-plex

## Python API

```python
from isobaric_channel_extractor import extract_channels

n_features = extract_channels("run.mzML", "tmt6plex", "quant.consensusXML")
```

## Pipeline

This tool is step 1 in the isobaric quantification workflow:

1. **isobaric_channel_extractor** (mzML -> consensusXML)
2. isobaric_isotope_corrector (consensusXML -> corrected consensusXML)
3. isobaric_quantifier (consensusXML -> quantified consensusXML)
