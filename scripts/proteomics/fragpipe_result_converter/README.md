# FragPipe Result Converter

Convert FragPipe psm.tsv to a standardized TSV format.

## Usage

```bash
python fragpipe_result_converter.py --input psm.tsv --output standardized.tsv
```

## Column Mapping

| FragPipe Column | Standard Column |
|---|---|
| Peptide | peptide |
| Modified Peptide | modified_peptide |
| Charge | charge |
| Observed M/Z | mz |
| Retention | rt |
| Protein | protein |
| Hyperscore | score |
| Intensity | intensity |
| Spectrum File | raw_file |
