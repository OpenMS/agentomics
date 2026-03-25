# Spectrum Similarity Scorer

Compute cosine similarity between MS2 spectra from MGF files using custom MGF reader and pyopenms SpectrumAlignment.

## Usage

```bash
python spectrum_similarity_scorer.py --query query.mgf --library reference.mgf --tolerance 0.02
python spectrum_similarity_scorer.py --query query.mgf --library ref.mgf --tolerance 0.02 --output scores.tsv
```
