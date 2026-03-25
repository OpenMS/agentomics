# Metabolite Class Predictor

Predict compound class from molecular formula using mass defect, element ratios (H:C, O:C), and RDBE (Ring and Double Bond Equivalents).

## Usage

```bash
python metabolite_class_predictor.py --input formulas.tsv --output predictions.tsv
```

### Input format

**formulas.tsv** (tab-separated):
```
formula
C6H12O6
C16H32O2
C2H5NO2
```

### Classification rules

Uses van Krevelen-style heuristic analysis:
- Carbohydrates: O:C 0.6-1.2, H:C 1.5-2.5
- Lipids: H:C > 1.5, O:C < 0.3
- Amino acids/Peptides: contains N, moderate ratios
- Alkaloids: contains N, high RDBE
- Terpenoids: moderate H:C, low O:C
- Phenolics/Polyketides: higher RDBE, moderate ratios
