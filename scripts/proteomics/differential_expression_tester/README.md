# Differential Expression Tester

Perform t-tests with Benjamini-Hochberg FDR correction on quantification matrices.

## Usage

```bash
python differential_expression_tester.py --input matrix.tsv --design design.tsv --test ttest --output de_results.tsv
python differential_expression_tester.py --input matrix.tsv --design design.tsv --test welch --output de_results.tsv
```

## Input Files

- **matrix.tsv** - Quantification matrix (rows=features, columns=samples)
- **design.tsv** - Experimental design with columns: `sample`, `condition`

## Output

TSV with columns: `feature`, `log2fc`, `pvalue`, `adj_pvalue`
