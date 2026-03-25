# Van Krevelen Data Generator

Compute H:C and O:C ratios from molecular formulas and classify compounds into biochemical classes (lipids, carbohydrates, amino acids, nucleotides) for Van Krevelen diagram analysis.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic ratio computation
python van_krevelen_data_generator.py --input formulas.tsv --output van_krevelen.tsv

# With biochemical class assignment
python van_krevelen_data_generator.py --input formulas.tsv --classify --output van_krevelen.tsv
```

## Input format

Tab-separated file with a `formula` column:

```
formula
C6H12O6
C16H32O2
C3H7NO2
```

## Output format

Tab-separated file with columns: formula, C, H, O, hc_ratio, oc_ratio, and optionally class.
