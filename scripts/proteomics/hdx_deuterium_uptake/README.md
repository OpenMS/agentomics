# HDX Deuterium Uptake Calculator

Calculate deuterium uptake from HDX-MS data: mass shift, fractional uptake, and back-exchange correction.

## Usage

```bash
python hdx_deuterium_uptake.py --peptides peptides.tsv --undeuterated ref.tsv --timepoints 0,10,60 --output uptake.tsv
```

## Input Format

- `peptides.tsv`: columns `sequence`, `timepoint`, `centroid_mass`
- `ref.tsv`: columns `sequence`, `centroid_mass` (undeuterated reference)

## Output

- `uptake.tsv` - Per-peptide mass shift and fractional uptake at each timepoint
