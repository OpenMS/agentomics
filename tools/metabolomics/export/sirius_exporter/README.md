# SIRIUS Exporter

Export features and MS2 spectra to SIRIUS .ms format for molecular formula identification.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic export
python sirius_exporter.py --features features.tsv --mzml data.mzML --output sirius_input.ms

# Custom tolerances
python sirius_exporter.py --features features.tsv --mzml data.mzML --mz-tolerance 0.02 --rt-tolerance 60 --output sirius_input.ms
```

## Feature TSV Format

The input features file should be a TSV with columns:
- `mz` (required): precursor m/z
- `rt` (required): retention time in seconds
- `charge` (optional): charge state
- `name` (optional): compound name
