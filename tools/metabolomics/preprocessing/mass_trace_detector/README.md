# Mass Trace Detector

Detect mass traces in LC-MS data using MassTraceDetection.

A mass trace is a series of connected peaks across consecutive spectra sharing a similar m/z value. This tool wraps `pyopenms.MassTraceDetection` and stores detected traces as chromatograms in an output mzML file.

## Usage

```bash
python mass_trace_detector.py --input run.mzML --output traces.mzML \
    --mass-error-ppm 10.0 --noise-threshold 1000
```

## Options

| Option | Default | Description |
|---|---|---|
| `--input` | required | Input mzML file with MS1 spectra |
| `--output` | required | Output mzML file (traces as chromatograms) |
| `--mass-error-ppm` | 10.0 | Mass error tolerance in ppm |
| `--noise-threshold` | 1000.0 | Noise intensity threshold |

## Python API

```python
from mass_trace_detector import detect_mass_traces

traces = detect_mass_traces(
    "run.mzML",
    "traces.mzML",
    mass_error_ppm=10.0,
    noise_threshold=1000.0,
)
print(f"Detected {traces} mass traces")
```
