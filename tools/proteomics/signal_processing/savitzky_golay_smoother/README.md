# Savitzky-Golay Smoother

Apply Savitzky-Golay smoothing to mass spectrometry data in mzML format using pyopenms's SavitzkyGolayFilter. Preserves peak shape better than simple moving average filters by fitting local polynomials.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic usage with defaults (frame_length=11, polynomial_order=3)
python savitzky_golay_smoother.py --input run.mzML --output smoothed.mzML

# Custom parameters
python savitzky_golay_smoother.py --input run.mzML --output smoothed.mzML --frame-length 15 --polynomial-order 4
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input mzML file path |
| `--output` | required | Output mzML file path |
| `--frame-length` | 11 | Number of data points in smoothing frame (must be odd) |
| `--polynomial-order` | 3 | Order of the fitting polynomial |
