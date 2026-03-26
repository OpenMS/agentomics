# Posterior Error Probability Estimator

Estimate posterior error probabilities (PEP) for peptide-spectrum matches using a mixture model approach.

## Usage

```bash
python posterior_error_probability.py --input peptides.idXML --output pep.idXML
```

## Requirements

```
pip install pyopenms click
```

## How it works

Wraps `pyopenms.PosteriorErrorProbabilityModel` to fit a two-component mixture model (correct + incorrect) to the PSM score distribution, then assigns PEP values between 0 and 1 to each hit.
