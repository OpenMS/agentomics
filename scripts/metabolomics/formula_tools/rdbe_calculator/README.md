# RDBE Calculator

Calculate Ring and Double Bond Equivalents (RDBE) for molecular formulas using the standard formula: RDBE = (2C + 2 - H + N + P) / 2.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python rdbe_calculator.py --input formulas.tsv --output rdbe.tsv
```

## Input format

Tab-separated file with a `formula` column:

```
formula
C6H6
C6H12O6
C10H8
```

## Output format

Tab-separated file with columns: formula, C, H, N, P, rdbe.
