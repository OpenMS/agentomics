# Formula Validator - Seven Golden Rules

Apply the Seven Golden Rules (Kind & Fiehn, 2007) to filter molecular formulas for plausibility.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Apply all rules
python formula_validator_golden_rules.py --input formulas.tsv --rules all --output validated.tsv

# Apply specific rules
python formula_validator_golden_rules.py --input formulas.tsv --rules rdbe,hc,oc --output validated.tsv
```

## Rules

- **rdbe**: RDBE must be non-negative
- **hc**: H/C ratio 0.2-3.1
- **nc**: N/C ratio <= 1.3
- **oc**: O/C ratio <= 1.2
- **sc**: S/C ratio <= 0.8
- **pc**: P/C ratio <= 0.3

## Input format

Tab-separated file with a `formula` column.

## Output format

Tab-separated file with per-rule pass/fail flags and overall valid column.
