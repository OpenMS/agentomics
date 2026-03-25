# Isotope Pattern Fit Scorer

Score observed vs theoretical isotope patterns using cosine similarity. Detect Cl/Br halogenation from enhanced M+2 peaks.

## Usage

```bash
python isotope_pattern_fit_scorer.py --observed "180.063:100,181.067:6.5,182.070:0.5" --formula C6H12O6 --output fit.json
```

### Output

JSON file with:
- `cosine_similarity`: Score between 0 and 1
- `observed_peaks` and `theoretical_peaks`: Peak lists
- `halogen_detection`: M+2 excess analysis with Cl/Br flagging
