# Spectral Entropy Scorer

Compute spectral entropy and entropy-based similarity between mass spectra, implementing the method from Li & Fiehn (2021).

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python spectral_entropy_scorer.py --query query_peaks.tsv --library lib_peaks.tsv \
    --tolerance 0.02 --output scores.tsv
```

## Input format

Tab-separated files with columns: spectrum_id, mz, intensity.

```
spectrum_id	mz	intensity
Q1	100.05	1000
Q1	150.08	500
Q1	200.12	250
```

## Output format

Tab-separated file with columns: query_id, library_id, query_entropy, entropy_similarity.
