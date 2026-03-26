# FLASHDeconv Wrapper

Deconvolve intact-protein mass spectra using the FLASHDeconv algorithm. Takes an mzML file with multiply-charged spectra and produces a TSV of deconvolved monoisotopic masses.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic usage
python flash_deconv.py --input intact.mzML --output masses.tsv

# Custom mass range
python flash_deconv.py --input intact.mzML --output masses.tsv \
    --min-mass 5000 --max-mass 50000
```

## Options

| Option | Default | Description |
|---|---|---|
| `--input` | (required) | Input mzML file |
| `--output` | (required) | Output TSV file |
| `--min-mass` | 5000 | Minimum mass in Da |
| `--max-mass` | 100000 | Maximum mass in Da |

## Python API

```python
from flash_deconv import deconvolve_intact

n_masses = deconvolve_intact("intact.mzML", "masses.tsv", min_mass=5000, max_mass=50000)
```
