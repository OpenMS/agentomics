# Spectrum Normalizer

Normalize mass spectra in mzML format using pyopenms's Normalizer. Supports normalization to maximum intensity (`to_one`) or to total ion current (`to_TIC`).

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Normalize to max intensity = 1.0 (default)
python spectrum_normalizer.py --input run.mzML --output normalized.mzML

# Normalize to TIC (sum of intensities = 1.0)
python spectrum_normalizer.py --input run.mzML --output normalized.mzML --method to_TIC
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input mzML file path |
| `--output` | required | Output mzML file path |
| `--method` | to_one | Normalization method: `to_one` or `to_TIC` |
