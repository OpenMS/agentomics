# mzML to MGF Converter

Convert MS2 spectra from mzML format to MGF (Mascot Generic Format).

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Convert MS2 spectra
python mzml_to_mgf_converter.py --input run.mzML --ms-level 2 --output spectra.mgf

# Filter by minimum number of peaks
python mzml_to_mgf_converter.py --input run.mzML --min-peaks 10 --output spectra.mgf
```
