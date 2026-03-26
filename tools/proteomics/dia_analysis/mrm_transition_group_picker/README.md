# MRM Transition Group Picker

Pick peaks in MRM transition group chromatograms using MRMTransitionGroupPicker.

## Usage

```bash
python mrm_transition_group_picker.py --input chromatograms.mzML --output picked.mzML
python mrm_transition_group_picker.py --input chromatograms.mzML --output picked.mzML --min-peak-width 5.0
```

## Input Format

mzML file with chromatograms. Native IDs should follow the convention `groupid_index`
(e.g., `peptide_0`, `peptide_1`) to enable automatic grouping.

## Dependencies

- pyopenms
- click
