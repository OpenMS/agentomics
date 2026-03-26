# Elution Peak Detector

Detect elution peaks from LC-MS data by first finding mass traces and then splitting them into individual chromatographic peaks.

Wraps `pyopenms.MassTraceDetection` + `pyopenms.ElutionPeakDetection` in a two-step workflow.

## Usage

```bash
python elution_peak_detector.py --input run.mzML --output peaks.mzML \
    --width-filtering auto
```

## Options

| Option | Default | Description |
|---|---|---|
| `--input` | required | Input mzML file with MS1 spectra |
| `--output` | required | Output mzML file (peaks as chromatograms) |
| `--width-filtering` | auto | Peak width filtering mode (auto/fixed/off) |
| `--mass-error-ppm` | 10.0 | Mass error tolerance in ppm |
| `--noise-threshold` | 1000.0 | Noise intensity threshold |

## Python API

```python
from elution_peak_detector import detect_elution_peaks

peaks = detect_elution_peaks(
    "run.mzML",
    "peaks.mzML",
    width_filtering="auto",
)
print(f"Detected {peaks} elution peaks")
```
