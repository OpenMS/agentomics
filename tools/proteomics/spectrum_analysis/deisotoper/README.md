# Deisotoper

Remove isotope peaks from MS2 spectra, keeping only the monoisotopic peak for each isotope envelope. Optionally converts peaks to singly-charged.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic usage
python deisotoper.py --input run.mzML --output deisotoped.mzML

# Custom tolerance and charge range
python deisotoper.py --input run.mzML --output deisotoped.mzML \
    --fragment-tolerance 0.05 --min-charge 1 --max-charge 3
```

## Options

| Option | Default | Description |
|---|---|---|
| `--input` | (required) | Input mzML file |
| `--output` | (required) | Output mzML file |
| `--fragment-tolerance` | 0.1 | Fragment mass tolerance in Da |
| `--min-charge` | 1 | Minimum charge state |
| `--max-charge` | 5 | Maximum charge state |

## Python API

```python
from deisotoper import deisotope

n_spectra = deisotope("run.mzML", "deisotoped.mzML", fragment_tolerance=0.1)
```
