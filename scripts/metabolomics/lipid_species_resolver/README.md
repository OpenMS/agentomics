# Lipid Species Resolver

From sum composition (e.g. PC 36:2), enumerate all possible acyl chain combinations and compute exact masses using pyopenms.

## Usage

```bash
python lipid_species_resolver.py --input lipids.tsv --output resolved.tsv
python lipid_species_resolver.py --input lipids.tsv --lipid-class PC --output resolved.tsv
```

### Input format

**lipids.tsv** (tab-separated):
```
lipid
PC 36:2
PE 34:1
TG 54:3
```

### Example output

For PC 36:2, the tool enumerates combinations like (16:0/20:2), (18:1/18:1), (18:0/18:2), etc.
