# DIA-NN Result Converter

Convert DIA-NN report.tsv to a standardized TSV format.

## Usage

```bash
python diann_result_converter.py --input report.tsv --output standardized.tsv
```

## Column Mapping

| DIA-NN Column | Standard Column |
|---|---|
| Stripped.Sequence | peptide |
| Modified.Sequence | modified_peptide |
| Precursor.Charge | charge |
| Precursor.Mz | mz |
| RT | rt |
| Protein.Group | protein |
| Q.Value | qvalue |
| Precursor.Quantity | intensity |
