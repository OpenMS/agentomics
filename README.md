# agentomics

A repository of agentic-created tools using [pyopenms](https://pyopenms.readthedocs.io/) for proteomics and metabolomics.

All code in this repo is written by AI agents. See [AGENTS.md](AGENTS.md) for the contributor guide.

## Requirements

```bash
pip install pyopenms
```

## Scripts

### Proteomics

| Script | Description |
|--------|-------------|
| [`peptide_mass_calculator`](scripts/proteomics/peptide_mass_calculator/) | Monoisotopic/average masses and b/y fragment ions for peptide sequences |
| [`protein_digest`](scripts/proteomics/protein_digest/) | In-silico enzymatic protein digestion |
| [`spectrum_file_info`](scripts/proteomics/spectrum_file_info/) | Summary statistics for mzML files |
| [`feature_detection_proteomics`](scripts/proteomics/feature_detection_proteomics/) | Peptide feature detection from LC-MS/MS data |

### Metabolomics

| Script | Description |
|--------|-------------|
| [`mass_accuracy_calculator`](scripts/metabolomics/mass_accuracy_calculator/) | m/z mass accuracy (ppm error) for sequences or formulas |
| [`isotope_pattern_matcher`](scripts/metabolomics/isotope_pattern_matcher/) | Theoretical isotope distributions and cosine similarity scoring |
| [`metabolite_feature_detection`](scripts/metabolomics/metabolite_feature_detection/) | Metabolite feature detection from LC-MS data |

## Validation

Each script is validated in an isolated venv. See [AGENTS.md](AGENTS.md) for validation commands.

## License

BSD 3-Clause — see [LICENSE](LICENSE).
