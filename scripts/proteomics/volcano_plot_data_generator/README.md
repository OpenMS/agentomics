# Volcano Plot Data Generator

Generate volcano plot data from differential expression results.

## Usage

```bash
python volcano_plot_data_generator.py --input de_results.tsv --fc-threshold 1.0 --pvalue 0.05 --output volcano.tsv
```

## Output Columns

- `feature` - Feature identifier
- `log2fc` - Log2 fold change
- `pvalue` - P-value
- `neg_log10_pvalue` - -log10(p-value) for plotting
- `regulation` - Classification: `up`, `down`, or `ns`
