# Phosphosite Class Filter

Classify phosphosites into Class I/II/III by localization probability and report enrichment efficiency.

## Usage

```bash
python phosphosite_class_filter.py --input phosphosites.tsv --class1-threshold 0.75 --output classified.tsv
```

## Input Format

Tab-separated file with columns: `peptide`, `protein`, `site`, `localization_prob`, `modification`

## Output

- `classified.tsv` - Input rows with added `site_class` and `valid_peptide` columns
- `classified_summary.tsv` - Summary counts and enrichment efficiency
