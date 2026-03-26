# MRM Decoy Generator

Generate decoy transitions from target transition lists using DecoyGenerator.

## Usage

```bash
python mrm_decoy_generator.py --input transitions.tsv --output decoy_transitions.tsv --method shuffle
python mrm_decoy_generator.py --input transitions.tsv --output decoy_transitions.tsv --method reverse
```

## Methods

- `shuffle` - Shuffle amino acid residues (preserving terminal residues)
- `reverse` - Reverse amino acid sequence (preserving enzymatic cutting sites)

## Input/Output Format

Tab-separated file with columns:
- `PrecursorMz` - Precursor m/z value
- `ProductMz` - Product (fragment) m/z value
- `LibraryIntensity` - Library intensity
- `PeptideSequence` - Peptide sequence
- `ProteinName` - Protein name
- `transition_name` - Unique transition identifier
- `transition_group_id` - Transition group identifier
- `PrecursorCharge` - Precursor charge state

## Dependencies

- pyopenms
- click
