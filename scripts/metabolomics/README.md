# Metabolomics Scripts

A collection of agentic-created Python scripts for common metabolomics tasks
using [pyopenms](https://pyopenms.readthedocs.io/).

## Scripts

### `metabolite_feature_detection.py`
Detect small-molecule features (isotope envelopes) in centroided LC-MS data
using the `FeatureFinderMetabo` pipeline (mass tracing → elution peak
detection → feature finding).  Output is written as a featureXML file.

```bash
python metabolite_feature_detection.py --input sample.mzML
python metabolite_feature_detection.py --input sample.mzML --output features.featureXML --noise 1e5
```

### `mass_accuracy_calculator.py`
Calculate the mass accuracy (ppm error) between a theoretical value derived
from a peptide sequence or molecular formula and one or more observed m/z
values.

```bash
# Peptide sequence
python mass_accuracy_calculator.py --sequence PEPTIDEK --observed 803.4560

# Molecular formula
python mass_accuracy_calculator.py --formula C6H12O6 --observed 181.0709

# Multiple observed values at charge 2
python mass_accuracy_calculator.py --sequence ACDEFGHIK --charge 2 \
    --observed 554.2478 554.2480 554.2482
```

### `isotope_pattern_matcher.py`
Generate a theoretical isotope distribution for any molecular formula and
optionally compute a cosine similarity score against observed peaks.

```bash
# Generate pattern for glucose
python isotope_pattern_matcher.py --formula C6H12O6

# Compare against observed peaks
python isotope_pattern_matcher.py --formula C6H12O6 \
    --peaks 181.0709,100.0 182.0742,6.7 183.0775,0.4
```
