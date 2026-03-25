# Coefficient of Variation Calculator

Calculate CV% (coefficient of variation) across replicates for each feature.

## Usage

```bash
python coefficient_of_variation_calculator.py --input matrix.tsv --groups groups.tsv --output cv_report.tsv
```

## Input Files

- **matrix.tsv** - Quantification matrix (rows=features, columns=samples)
- **groups.tsv** - Group assignments with columns: `sample`, `group`
