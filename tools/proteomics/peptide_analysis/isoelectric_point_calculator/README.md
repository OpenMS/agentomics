# Isoelectric Point Calculator

Calculate pI for peptides and proteins using Henderson-Hasselbalch with multiple pKa sets.

## Usage

```bash
python isoelectric_point_calculator.py --sequence ACDEFGHIK --pk-set lehninger --output pi.json
python isoelectric_point_calculator.py --fasta proteins.fasta --output pi.tsv
python isoelectric_point_calculator.py --sequence PEPTIDEK --charge-curve --output curve.json
```

## Available pKa Sets

- `lehninger` (default)
- `emboss`
- `stryer`
- `solomon`
