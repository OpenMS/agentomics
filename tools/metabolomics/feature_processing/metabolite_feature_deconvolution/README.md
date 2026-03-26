# Metabolite Feature Deconvolution

Group metabolite features by adduct relationships (e.g. [M+H]+, [M+Na]+, [M+K]+) into consensus feature groups. This helps identify features that originate from the same underlying metabolite.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic usage
python metabolite_feature_deconvolution.py --input features.featureXML --output grouped.consensusXML
```

## Options

| Option | Default | Description |
|---|---|---|
| `--input` | (required) | Input featureXML file |
| `--output` | (required) | Output consensusXML file |
| `--top` | 20 | Number of top groups to print |

## Python API

```python
from metabolite_feature_deconvolution import deconvolve_adducts

n_groups = deconvolve_adducts("features.featureXML", "grouped.consensusXML")
```
