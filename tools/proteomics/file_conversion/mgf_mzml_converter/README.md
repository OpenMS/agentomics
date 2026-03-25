# MGF ↔ mzML Converter

Bidirectional conversion between MGF (Mascot Generic Format) and mzML.
Consolidates `mgf_to_mzml_converter` and `mzml_to_mgf_converter` into a single tool.

## Features

- **MGF → mzML**: preserves title, precursor m/z, charge, retention time, fragment peaks
- **mzML → MGF**: full filter suite (charge, RT range, m/z range, min peaks,
  min base-peak intensity, MS level)
- **Auto-detection** of conversion direction from file extensions
- Conversion statistics printed to stdout

## Installation

```bash
pip install pyopenms
```

## CLI Usage

```bash
# Auto-detect direction from file extensions
python mgf_mzml_converter.py --input spectra.mgf --output spectra.mzML
python mgf_mzml_converter.py --input run.mzML --output spectra.mgf

# Explicit direction
python mgf_mzml_converter.py --direction mgf2mzml --input spectra.mgf --output out.mzML
python mgf_mzml_converter.py --direction mzml2mgf --input run.mzML --output out.mgf

# mzML → MGF with filters
python mgf_mzml_converter.py --input run.mzML --output out.mgf \
    --charge 2 3 \
    --rt-min 600 --rt-max 1800 \
    --mz-min 400 --mz-max 1200 \
    --min-peaks 5 \
    --min-intensity 1000
```

## Options

| Option | Direction | Description |
|---|---|---|
| `--direction` | both | `mgf2mzml` or `mzml2mgf` (auto if omitted) |
| `--ms-level` | mzml2mgf | MS level to export (default: 2) |
| `--min-peaks` | mzml2mgf | Minimum fragment peaks (default: 1) |
| `--charge` | mzml2mgf | Keep only these precursor charges (repeatable) |
| `--rt-min` / `--rt-max` | mzml2mgf | Retention time window in seconds |
| `--mz-min` / `--mz-max` | mzml2mgf | Precursor m/z window |
| `--min-intensity` | mzml2mgf | Minimum base-peak intensity |
