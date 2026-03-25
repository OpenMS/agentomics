# Spectral Library Builder

Build a spectral library in MSP format from mzML + peptide identification list.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python spectral_library_builder.py --input run.mzML --peptides identified.tsv --output library.msp
```

## Peptide TSV Format

The input peptide identifications file should be a TSV with columns:
- `sequence` (required): peptide sequence
- `charge` (required): charge state
- `rt` (required): retention time in seconds
- `mz` (optional): precursor m/z (calculated from sequence if not provided)
- `score` (optional): identification score
