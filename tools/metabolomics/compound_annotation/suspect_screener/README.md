# Suspect Screener

Match detected features against a suspect screening list by exact mass within a given ppm tolerance. Results are ranked by mass error.

## Usage

```bash
python suspect_screener.py --input features.tsv --suspects suspect_list.csv --ppm 5 --output matches.tsv
```

### Input formats

**features.tsv** (tab-separated):
```
feature_id	mz	rt	intensity
F1	180.0634	120.5	5000
```

**suspect_list.csv** (comma-separated):
```
name,formula,exact_mass
Glucose,C6H12O6,180.0634
```

If `exact_mass` is empty, it will be computed from the formula using pyopenms.
