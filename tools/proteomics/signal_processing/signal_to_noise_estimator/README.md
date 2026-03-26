# Signal-to-Noise Estimator

Estimate signal-to-noise ratios for peaks in MS spectra using the median method.

## Usage

```bash
python signal_to_noise_estimator.py --input run.mzML --output sn_report.tsv
```

## Output

TSV file with columns: `spectrum_index`, `mz`, `intensity`, `sn_ratio`

## API

```python
from signal_to_noise_estimator import estimate_sn

count = estimate_sn("run.mzML", "sn_report.tsv")
```
