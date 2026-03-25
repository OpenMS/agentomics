# HDX Back-Exchange Estimator

Estimate per-peptide back-exchange from fully deuterated controls.

## Usage

```bash
python hdx_back_exchange_estimator.py --peptides peptides.tsv --fully-deuterated fd.tsv --max-backexchange 40 --output report.tsv
```

## Input Format

- `peptides.tsv`: columns `sequence`, `centroid_mass` (undeuterated reference)
- `fd.tsv`: columns `sequence`, `centroid_mass` (fully deuterated controls)

## Output

- `report.tsv` - Per-peptide back-exchange estimates with threshold flags
