# MID Natural Abundance Corrector

Correct mass isotopomer distributions (MIDs) for natural 13C abundance using a correction matrix approach.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python mid_natural_abundance_corrector.py --input isotopologues.tsv --formula C6H12O6 --tracer 13C --output corrected.tsv
```

## Input format

Tab-separated file with columns: sample, M0, M1, M2, ... (fractional abundances):

```
sample	M0	M1	M2	M3	M4	M5	M6
ctrl_1	0.92	0.06	0.02	0.0	0.0	0.0	0.0
labeled_1	0.10	0.05	0.15	0.30	0.25	0.10	0.05
```

## Output format

Tab-separated file with corrected MID values (M0_corrected, M1_corrected, ...).
