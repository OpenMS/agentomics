# Isotope Pattern Analyzer

Generate theoretical isotope distributions for molecular formulas, score observed
isotope patterns using cosine similarity, and detect halogenation (Cl/Br).

This tool consolidates `isotope_pattern_matcher`, `isotope_pattern_scorer`, and
`isotope_pattern_fit_scorer` into a single, improved utility.

## Features

- Theoretical isotope pattern generation via pyopenms `CoarseIsotopePatternGenerator`
- Cosine similarity scoring between observed and theoretical patterns
- **Da or ppm m/z tolerance** — choose your preferred unit
- Halogen (Cl/Br) detection from M+2 peak enhancement
- JSON output with per-peak detail
- Terminal bar-chart preview of the theoretical distribution
- Optional numpy acceleration for cosine computation

## Installation

```bash
pip install pyopenms
```

## CLI Usage

```bash
# Generate and display the isotope pattern for glucose
python isotope_pattern_analyzer.py --formula C6H12O6

# Score observed peaks against the formula (colon-separated format)
python isotope_pattern_analyzer.py --formula C6H12O6 \
    --observed "180.063:100,181.067:6.5,182.070:0.5" \
    --output result.json

# Use legacy comma-separated format (one --peaks flag per peak)
python isotope_pattern_analyzer.py --formula C6H12O6 \
    --peaks 180.063,100.0 --peaks 181.067,6.5 \
    --output result.json

# Use ppm tolerance
python isotope_pattern_analyzer.py --formula C6H12O6 \
    --observed "180.063:100,181.067:6.5" \
    --tolerance 10 --tolerance-unit ppm

# Detect halogenation (chlorinated compound example)
python isotope_pattern_analyzer.py --formula C6H5Cl \
    --observed "112.007:100,113.011:5.5,114.004:33.0" \
    --output halogen_result.json
```

## Output JSON Structure

```json
{
  "formula": "C6H12O6",
  "cosine_similarity": 0.9987,
  "n_peaks_compared": 3,
  "tolerance": 0.05,
  "tolerance_unit": "da",
  "peaks": [
    {"peak_index": 0, "obs_mz": 180.063, "theo_mz": 180.0634, "obs_intensity": 100.0, "theo_intensity": 100.0},
    ...
  ],
  "theoretical_pattern": [...],
  "halogen_detection": {
    "m2_ratio_observed": 0.5,
    "m2_ratio_theoretical": 0.42,
    "m2_excess": 0.08,
    "halogen_flag": false,
    "possible_halogen": "none"
  }
}
```

## Halogen Detection Thresholds

| M+2 excess above theoretical | Interpretation                  |
|------------------------------|---------------------------------|
| < 10 %                       | No halogenation detected        |
| 10–20 %                      | Cl (weak signal)                |
| 20–70 %                      | Cl                              |
| > 70 %                       | Br                              |
