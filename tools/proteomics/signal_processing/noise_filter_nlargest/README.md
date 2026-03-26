# Noise Filter - N Largest

Keep only the N most intense peaks in each spectrum.

Wraps `pyopenms.NLargest`.

## Usage

```bash
# Basic usage
python noise_filter_nlargest.py --input run.mzML --output filtered.mzML

# Keep top 50 peaks
python noise_filter_nlargest.py --input run.mzML --output filtered.mzML --n 50
```

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--input` | str | required | Input mzML file |
| `--output` | str | required | Output mzML file |
| `--n` | int | 100 | Number of most intense peaks to keep per spectrum |

## Python API

```python
from noise_filter_nlargest import filter_nlargest

n_spectra = filter_nlargest("input.mzML", "output.mzML", n=50)
```
