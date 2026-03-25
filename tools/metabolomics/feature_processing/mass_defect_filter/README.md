# Mass Defect Filter

Compute mass defect and Kendrick mass defect for features, then filter by mass defect range.

## Usage

```bash
python mass_defect_filter.py --input features.tsv --mdf-min 0.1 --mdf-max 0.3
python mass_defect_filter.py --input features.tsv --mdf-min 0.1 --mdf-max 0.3 --kendrick-base CH2 --output filtered.tsv
```
