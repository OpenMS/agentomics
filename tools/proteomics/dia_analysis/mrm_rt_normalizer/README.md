# MRM RT Normalizer

Normalize retention times using MRMRTNormalizer outlier removal and linear fitting.

## Usage

```bash
python mrm_rt_normalizer.py --input rt_pairs.tsv --output model.tsv --method iterative
python mrm_rt_normalizer.py --input rt_pairs.tsv --output model.tsv --method ransac
```

## Input Format

Tab-separated file with columns:
- `irt` (or `reference_rt`) - Reference/iRT retention time
- `measured_rt` - Measured retention time

## Output

TSV file with comment lines containing model parameters and cleaned RT pairs:
```
# slope=1.100000
# intercept=5.000000
# r_squared=0.999000
# outliers_removed=2
reference_rt    measured_rt
0.0    5.1
1.0    6.2
...
```

## Dependencies

- pyopenms
- click
