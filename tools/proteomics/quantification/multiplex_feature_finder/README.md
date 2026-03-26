# Multiplex Feature Finder

Detect multiplex-labeled peptide features (e.g. SILAC light/heavy pairs) from an mzML file using the FeatureFinderMultiplexAlgorithm.

## Usage

```bash
# Basic SILAC light/heavy (Lys8)
python multiplex_feature_finder.py --input run.mzML --output features.featureXML --labels "[][Lys8]"

# SILAC with Lys8 and Arg10
python multiplex_feature_finder.py --input run.mzML --output features.featureXML --labels "[][Lys8,Arg10]"

# Triple SILAC
python multiplex_feature_finder.py --input run.mzML --output features.featureXML --labels "[][Lys4][Lys8]"

# Custom parameters
python multiplex_feature_finder.py --input run.mzML --output features.featureXML \
    --labels "[][Lys8]" --charge-low 2 --charge-high 5 --mz-tolerance 6.0 --intensity-cutoff 500.0
```

## Supported Labels

Arg6, Arg10, Lys4, Lys6, Lys8, Leu3, Dimethyl0, Dimethyl4, Dimethyl6, Dimethyl8, ICPL0, ICPL4, ICPL6, ICPL10

## Dependencies

```
pyopenms
click
```
