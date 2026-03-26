# MRM Feature Finder

Find and score features in MRM/SRM chromatograms using MRMFeatureFinderScoring.

## Usage

```bash
python mrm_feature_finder.py --input chromatograms.mzML --transitions transitions.tsv --output features.featureXML
```

## Transitions TSV Format

Tab-separated file with columns:
- `PrecursorMz` - Precursor m/z value
- `ProductMz` - Product (fragment) m/z value
- `LibraryIntensity` - Library intensity
- `PeptideSequence` - Peptide sequence
- `ProteinName` - Protein name
- `transition_name` - Unique transition identifier (must match chromatogram native IDs)
- `transition_group_id` - Transition group identifier

## Dependencies

- pyopenms
- click
