# Kovats Retention Index Calculator

Calculate Kovats retention indices from n-alkane standard retention times for GC-MS data.

Formula: RI = 100n + 100 * (log(RT_x) - log(RT_n)) / (log(RT_{n+1}) - log(RT_n))

## Usage

```bash
python kovats_ri_calculator.py --input features.tsv --standards alkane_rts.tsv --output ri_values.tsv
```

### Input formats

**alkane_rts.tsv** (tab-separated):
```
carbon_number	rt
8	2.0
10	5.0
12	10.0
```

**features.tsv** (tab-separated):
```
feature_id	rt
F1	3.5
F2	7.0
```
