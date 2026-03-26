# Internal Calibrator

Recalibrate m/z values in an mzML file using internal calibrants from peptide identifications (idXML). Wraps pyopenms `InternalCalibration` with linear and quadratic transformation models.

## Features

- Load calibrant references from idXML peptide identifications
- Build global linear or quadratic calibration models
- Report before/after ppm error statistics (median and MAD)
- Supports weighted model variants

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### CLI

```bash
# Linear calibration
python internal_calibrator.py --input run.mzML --ids peptides.idXML --output calibrated.mzML --model linear

# Quadratic calibration
python internal_calibrator.py --input run.mzML --ids peptides.idXML --output calibrated.mzML --model quadratic

# Custom ppm tolerance for calibrant matching
python internal_calibrator.py --input run.mzML --ids peptides.idXML --output calibrated.mzML --tol-ppm 10.0
```

### Python API

```python
from internal_calibrator import calibrate_mz

result = calibrate_mz("run.mzML", "peptides.idXML", "calibrated.mzML", model="linear")
print(f"Before: {result['before_ppm_median']:.2f} ppm")
print(f"After:  {result['after_ppm_median']:.2f} ppm")
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input mzML file |
| `--ids` | required | Peptide identifications idXML file |
| `--output` | required | Output calibrated mzML file |
| `--model` | `linear` | Model type: `linear`, `linear_weighted`, `quadratic`, `quadratic_weighted` |
| `--tol-ppm` | `20.0` | Tolerance in ppm for calibrant matching |

## Output

The function returns a dictionary with:

- `success`: whether calibration succeeded
- `n_calibrants`: number of calibrant points found
- `before_ppm_median`: median ppm error before calibration
- `before_ppm_mad`: MAD of ppm error before calibration
- `after_ppm_median`: median ppm error after calibration
- `after_ppm_mad`: MAD of ppm error after calibration
