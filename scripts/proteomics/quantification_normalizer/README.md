# Quantification Normalizer

Normalize quantification matrices using median, quantile, or total intensity normalization.

## Usage

```bash
python quantification_normalizer.py --input matrix.tsv --method median --output normalized.tsv
python quantification_normalizer.py --input matrix.tsv --method quantile --output normalized.tsv
python quantification_normalizer.py --input matrix.tsv --method total_intensity --output normalized.tsv
```

## Methods

- **median** - Shift columns so all have the same median
- **quantile** - Force all columns to have the same distribution
- **total_intensity** - Scale columns to the same total intensity
