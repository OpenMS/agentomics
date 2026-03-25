# Phospho Enrichment QC

Compute phospho-enrichment efficiency and pSer/pThr/pTyr ratios from search results.

## Usage

```bash
python phospho_enrichment_qc.py --input search_results.tsv --output enrichment.tsv
```

## Input Format

Tab-separated file with a `sequence` column containing peptide sequences with modification annotations.

## Output

- `enrichment.tsv` - Enrichment metrics including efficiency and residue ratios
