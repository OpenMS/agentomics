# Top-Down Coverage Calculator

Compute per-residue bond cleavage coverage from fragment ions in top-down proteomics.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python topdown_coverage_calculator.py --sequence PROTEINSEQ \
    --fragments observed.tsv --tolerance 10 --output coverage.tsv
```

### Input format

The fragments file is a TSV with a `mass` column containing observed fragment ion masses:

```
mass
300.1589
401.2066
```

### Parameters

| Flag | Description |
|------|-------------|
| `--sequence` | Protein amino acid sequence |
| `--fragments` | TSV with observed fragment ion masses |
| `--tolerance` | Tolerance in ppm (default: 10) |
| `--output` | Output coverage TSV |
