# PSM Feature Extractor

Extract rescoring features from PSMs by comparing experimental spectra to theoretical spectra.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python psm_feature_extractor.py --mzml run.mzML --peptides psms.tsv --output features.tsv
```

## PSM TSV Format

The input PSMs file should be a TSV with columns:
- `sequence` (required): peptide sequence
- `charge` (required): charge state
- `rt` (required): retention time in seconds
- `scan_index` (optional): spectrum index in mzML
- `mz` (optional): precursor m/z (calculated from sequence if not provided)
