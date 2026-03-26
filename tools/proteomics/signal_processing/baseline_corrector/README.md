# Baseline Corrector

Remove baseline drift from mass spectrometry data in mzML format using pyopenms's MorphologicalFilter. Uses top-hat filtering to estimate and subtract the baseline from profile spectra.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic usage with default structuring element length (3.0 Th)
python baseline_corrector.py --input run.mzML --output corrected.mzML

# Custom structuring element length
python baseline_corrector.py --input run.mzML --output corrected.mzML --struct-element-length 5.0
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input mzML file path |
| `--output` | required | Output mzML file path |
| `--struct-element-length` | 3.0 | Length of the structuring element in Thomson |
