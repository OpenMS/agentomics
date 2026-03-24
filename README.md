# agentomics

A repository of agentic-created scripts using [pyopenms](https://pyopenms.readthedocs.io/)
for proteomics and metabolomics users.

## Requirements

```bash
pip install pyopenms
```

Or install all dependencies at once:

```bash
pip install -r requirements.txt
```

## Repository Structure

```
agentomics/
├── scripts/
│   ├── proteomics/          # Scripts for proteomics workflows
│   │   ├── peptide_mass_calculator.py
│   │   ├── protein_digest.py
│   │   ├── spectrum_file_info.py
│   │   └── feature_detection_proteomics.py
│   └── metabolomics/        # Scripts for metabolomics workflows
│       ├── metabolite_feature_detection.py
│       ├── mass_accuracy_calculator.py
│       └── isotope_pattern_matcher.py
└── tests/
    ├── test_proteomics.py
    └── test_metabolomics.py
```

## Proteomics Scripts

| Script | Description |
|--------|-------------|
| [`peptide_mass_calculator.py`](scripts/proteomics/peptide_mass_calculator.py) | Monoisotopic/average masses and b/y fragment ion series for peptide sequences |
| [`protein_digest.py`](scripts/proteomics/protein_digest.py) | In-silico enzymatic protein digestion (Trypsin, Lys-C, …) |
| [`spectrum_file_info.py`](scripts/proteomics/spectrum_file_info.py) | Summary statistics for mzML files (spectra, RT range, TIC) |
| [`feature_detection_proteomics.py`](scripts/proteomics/feature_detection_proteomics.py) | FeatureFinderCentroided-based peptide feature detection from LC-MS/MS data |

See [`scripts/proteomics/README.md`](scripts/proteomics/README.md) for usage examples.

## Metabolomics Scripts

| Script | Description |
|--------|-------------|
| [`metabolite_feature_detection.py`](scripts/metabolomics/metabolite_feature_detection.py) | FeatureFinderMetabo-based metabolite feature detection from LC-MS data |
| [`mass_accuracy_calculator.py`](scripts/metabolomics/mass_accuracy_calculator.py) | Compute m/z mass accuracy (ppm error) for sequences or molecular formulas |
| [`isotope_pattern_matcher.py`](scripts/metabolomics/isotope_pattern_matcher.py) | Generate theoretical isotope distributions and score against observed peaks |

See [`scripts/metabolomics/README.md`](scripts/metabolomics/README.md) for usage examples.

## Running Tests

```bash
pip install pytest pyopenms
python -m pytest tests/ -v
```
