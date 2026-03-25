# iRT Calculator

Convert observed retention times to indexed retention times (iRT) using reference peptides and linear regression.

## Usage

```bash
python irt_calculator.py --input identifications.tsv --reference irt_standards.tsv --output irt_converted.tsv
```

## Input Format

- `--reference`: TSV with columns `sequence`, `observed_rt`, `irt`
- `--input`: TSV with columns `sequence`, `rt`
