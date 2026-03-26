# Peak Type Estimator

Estimate whether spectra in an mzML file contain profile or centroided data.

## Usage

```bash
python peak_type_estimator.py --input run.mzML --output report.tsv
```

## Output

TSV file with columns: `spectrum_index`, `ms_level`, `peak_type`

## API

```python
from peak_type_estimator import estimate_peak_type

result = estimate_peak_type("run.mzML", "report.tsv")
# result = {"type_counts": {"profile": 5, "centroid": 3}, "total_spectra": 8}
```
