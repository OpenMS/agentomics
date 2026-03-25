# Crosslink Distance Validator

Validate crosslinks against PDB structure distances by computing CA-CA distances.

## Usage

```bash
python xl_distance_validator.py --crosslinks links.tsv --pdb structure.pdb --max-distance 30 --output distances.tsv
```

## Input Format

- `links.tsv`: columns `peptide1`, `peptide2`, `chain1`, `residue1`, `chain2`, `residue2`
- `structure.pdb`: Standard PDB format file

## Output

- `distances.tsv` - Crosslinks with computed distances and satisfaction flags
