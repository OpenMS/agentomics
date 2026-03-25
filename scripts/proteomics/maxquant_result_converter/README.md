# MaxQuant Result Converter

Convert MaxQuant evidence.txt to a standardized TSV format.

## Usage

```bash
python maxquant_result_converter.py --input evidence.txt --output standardized.tsv
```

## Column Mapping

| MaxQuant Column | Standard Column |
|---|---|
| Sequence | peptide |
| Modified sequence | modified_peptide |
| Charge | charge |
| m/z | mz |
| Retention time | rt |
| Proteins | protein |
| Score | score |
| PEP | pep |
| Intensity | intensity |
