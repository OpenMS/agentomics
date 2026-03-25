# Kendrick Mass Defect Analyzer

Compute Kendrick Mass Defect (KMD) for configurable base units (CH2, CF2, C2H4O) and group features into homologous series.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# From formula column
python kendrick_mass_defect_analyzer.py --input features.tsv --base CH2 --output kmd.tsv

# From m/z column with custom tolerance
python kendrick_mass_defect_analyzer.py --input features.tsv --base CF2 --kmd-tolerance 0.01 --output kmd.tsv
```

## Input format

Tab-separated file with either a `formula` or `mz` column:

```
formula
C16H32O2
C18H36O2
C20H40O2
```

## Output format

Tab-separated file with columns: formula/exact_mass, kendrick_mass, nominal_kendrick_mass, kmd, series.
