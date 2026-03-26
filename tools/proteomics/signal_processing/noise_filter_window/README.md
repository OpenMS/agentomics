# Noise Filter - Window Mower

Remove low-intensity peaks from MS spectra using a sliding window approach. Within each window of a given m/z width, only the N most intense peaks are retained.

Wraps `pyopenms.WindowMower`.

## Usage

```bash
# Basic usage
python noise_filter_window.py --input run.mzML --output filtered.mzML

# Custom window size and peak count
python noise_filter_window.py --input run.mzML --output filtered.mzML --window-size 50.0 --peak-count 3
```

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--input` | str | required | Input mzML file |
| `--output` | str | required | Output mzML file |
| `--window-size` | float | 50.0 | Window size in Da |
| `--peak-count` | int | 3 | Number of peaks to keep per window |

## Python API

```python
from noise_filter_window import filter_window

n_spectra = filter_window("input.mzML", "output.mzML", window_size=50.0, peak_count=3)
```
