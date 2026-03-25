# SILAC Half-Life Calculator

Fit exponential decay to SILAC H/L ratios for protein turnover analysis.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python silac_halflife_calculator.py --input hl_ratios.tsv \
    --timepoints 0,6,12,24,48 --output halflives.tsv
```

### Input format

Tab-separated file with a `protein_id` column and one column per timepoint:

```
protein_id	t0	t6	t12	t24	t48
P12345	10.0	7.5	4.2	1.8	0.5
```

### Parameters

| Flag | Description |
|------|-------------|
| `--input` | Input TSV with protein IDs and H/L ratios |
| `--timepoints` | Comma-separated timepoint values |
| `--output` | Output half-lives TSV |
