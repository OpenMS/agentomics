# Linear Resampler

Resample spectra in an mzML file to uniformly spaced m/z values using linear interpolation.

## Usage

```bash
python linear_resampler.py --input run.mzML --output resampled.mzML --spacing 0.01
```

## API

```python
from linear_resampler import resample_experiment

count = resample_experiment("run.mzML", "resampled.mzML", spacing=0.01)
```
