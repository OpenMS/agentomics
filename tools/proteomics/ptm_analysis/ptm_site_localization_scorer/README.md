# PTM Site Localization Scorer

Score PTM site localization confidence using fragment ion coverage comparison.

## Usage

```bash
python ptm_site_localization_scorer.py --mz-list "200.1,300.2,400.3" --intensities "100,200,150" \
    --peptide "PEPS(Phospho)TIDEK" --tolerance 0.02 --output scores.tsv
```
