# Noise Filter - Threshold Mower

Remove peaks below an intensity threshold from MS spectra.

Wraps `pyopenms.ThresholdMower`.

## Usage

```bash
# Basic usage
python noise_filter_threshold.py --input run.mzML --output filtered.mzML

# Custom threshold
python noise_filter_threshold.py --input run.mzML --output filtered.mzML --threshold 100.0
```

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--input` | str | required | Input mzML file |
| `--output` | str | required | Output mzML file |
| `--threshold` | float | 100.0 | Minimum intensity threshold |

## Python API

```python
from noise_filter_threshold import filter_threshold

n_spectra = filter_threshold("input.mzML", "output.mzML", threshold=100.0)
```
