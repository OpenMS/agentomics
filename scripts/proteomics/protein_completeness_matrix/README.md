# Protein Completeness Matrix

Compute data completeness per protein and per sample from a quantification matrix.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python protein_completeness_matrix.py --input quant_matrix.tsv \
    --min-completeness 0.5 --output completeness.tsv
```

### Input format

Tab-separated quantification matrix with `protein_id` as the first column:

```
protein_id	sample1	sample2	sample3
P12345	100.5	NA	95.2
P67890	200.1	180.3	190.7
```

### Parameters

| Flag | Description |
|------|-------------|
| `--input` | Input quantification matrix TSV |
| `--min-completeness` | Minimum completeness fraction to retain a protein (default: 0.0) |
| `--output` | Output completeness TSV |
