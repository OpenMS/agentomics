# Gaussian Smoother

Apply Gaussian smoothing to mass spectrometry data in mzML format using pyopenms's GaussFilter. Reduces high-frequency noise in profile spectra by convolving with a Gaussian kernel.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic usage with default width (0.2 Th)
python gaussian_smoother.py --input run.mzML --output smoothed.mzML

# Custom Gaussian width
python gaussian_smoother.py --input run.mzML --output smoothed.mzML --gaussian-width 0.5
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input mzML file path |
| `--output` | required | Output mzML file path |
| `--gaussian-width` | 0.2 | Width of the Gaussian kernel in Thomson |
