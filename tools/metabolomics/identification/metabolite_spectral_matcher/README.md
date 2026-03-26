# Metabolite Spectral Matcher

Match experimental MS/MS spectra against a spectral library for metabolite identification.

Wraps `pyopenms.MetaboliteSpectralMatching` to compare query spectra against library spectra and return match scores.

## Usage

```bash
python metabolite_spectral_matcher.py --input spectra.mzML \
    --library lib.mzML --output matches.mzTab --precursor-tol 0.1
```

## Options

| Option | Default | Description |
|---|---|---|
| `--input` | required | Input mzML file with query MS/MS spectra |
| `--library` | required | Spectral library mzML file |
| `--output` | required | Output mzTab file |
| `--precursor-tol` | 0.1 | Precursor mass tolerance in Da |

## Python API

```python
from metabolite_spectral_matcher import match_spectra

matches = match_spectra(
    "spectra.mzML",
    "library.mzML",
    "matches.mzTab",
    precursor_tol=0.1,
)
print(f"Matched {matches} spectra")
```
