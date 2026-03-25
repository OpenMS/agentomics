# Intensity Distribution Reporter

Report per-sample intensity statistics from a quantification matrix.

## Usage

```bash
python intensity_distribution_reporter.py --input matrix.tsv --output intensity_stats.tsv
```

## Output Columns

`sample`, `n_values`, `n_missing`, `mean`, `median`, `sd`, `min`, `max`, `q1`, `q3`
