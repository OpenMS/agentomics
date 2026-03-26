# Chromatogram Extractor

Extract ion chromatograms (XICs) from DIA/SWATH-MS mzML files using a transition list.

## Usage

```bash
python chromatogram_extractor.py --input dia.mzML --transitions transitions.tsv --output chromatograms.mzML
python chromatogram_extractor.py --input dia.mzML --transitions transitions.tsv --output chromatograms.mzML --mz-tol 10.0
```

## Transitions TSV Format

Tab-separated file with columns:
- `PrecursorMz` - Precursor m/z value
- `ProductMz` - Product (fragment) m/z value
- `LibraryIntensity` - Library intensity
- `PeptideSequence` - Peptide sequence
- `ProteinName` - Protein name
- `transition_name` - Unique transition identifier
- `transition_group_id` - Transition group identifier

## Dependencies

- pyopenms
- click
