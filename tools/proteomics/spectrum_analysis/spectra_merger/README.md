# Spectra Merger

Merge MS spectra from an mzML file using pyopenms `SpectraMerger`. Supports block-wise merging of consecutive spectra and precursor-based merging of MS2 spectra sharing the same precursor m/z.

## Features

- **Block merge**: combine consecutive spectra in fixed-size blocks (e.g. averaging 3 sequential MS1 scans)
- **Precursor merge**: combine MS2 spectra with matching precursor m/z into single consensus spectra
- Configurable m/z binning width and unit (ppm or Da)

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### CLI

```bash
# Block-wise merge: combine every 3 consecutive spectra
python spectra_merger.py --input run.mzML --output merged.mzML --mode block --block-size 3

# Precursor-based merge: merge MS2 spectra by precursor m/z
python spectra_merger.py --input run.mzML --output merged.mzML --mode precursor

# Custom binning width
python spectra_merger.py --input run.mzML --output merged.mzML --mode block --block-size 5 --mz-binning-width 0.01 --mz-binning-unit Da
```

### Python API

```python
from spectra_merger import merge_spectra

# Block merge
n_spectra = merge_spectra("run.mzML", "merged.mzML", mode="block", block_size=3)

# Precursor merge
n_spectra = merge_spectra("run.mzML", "merged.mzML", mode="precursor")
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input mzML file |
| `--output` | required | Output merged mzML file |
| `--mode` | `block` | Merging mode: `block` or `precursor` |
| `--block-size` | `3` | Spectra per block (block mode only) |
| `--mz-binning-width` | `10.0` | Bin width for peak alignment |
| `--mz-binning-unit` | `ppm` | Unit for binning width: `ppm` or `Da` |
