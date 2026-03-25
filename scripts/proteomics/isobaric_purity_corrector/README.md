# Isobaric Purity Corrector

Correct TMT/iTRAQ reporter ion quantification for isotopic impurity using a purity correction matrix.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python isobaric_purity_corrector.py --input quant.tsv --label TMT16plex \
    --purity-matrix purity.csv --output corrected.tsv
```

### Input format

Quantification TSV with channel columns matching the labeling scheme:

```
spectrum_id	126	127N	127C
spec1	1000.0	50.0	30.0
```

Purity matrix CSV (N x N, no headers):

```
0.95,0.03,0.02
0.02,0.94,0.04
0.01,0.03,0.96
```

### Parameters

| Flag | Description |
|------|-------------|
| `--input` | Input quantification TSV |
| `--label` | Labeling scheme (TMT6plex, TMT10plex, TMT16plex, iTRAQ4plex, etc.) |
| `--purity-matrix` | Purity correction matrix CSV |
| `--output` | Output corrected TSV |
