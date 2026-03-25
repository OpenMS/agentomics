# RNA Mass Calculator

Calculate mass, molecular formula, and isotope patterns for RNA sequences.

## Usage

```bash
python rna_mass_calculator.py --sequence AAUGC --charge 2
python rna_mass_calculator.py --sequence AAUGCAAUGG --charge 3 --output mass.json
python rna_mass_calculator.py --sequence AAUGC --isotopes 5
```

## Options

- `--sequence` - RNA sequence (A, C, G, U characters)
- `--charge` - Charge state for m/z calculation (default: 1)
- `--isotopes` - Number of isotope peaks to calculate (default: 0 = off)
- `--output` - Output JSON file path (optional; prints to stdout if omitted)
