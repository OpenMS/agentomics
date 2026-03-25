# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Agentomics is a collection of standalone CLI tools built with [pyopenms](https://pyopenms.readthedocs.io/) for proteomics and metabolomics workflows. These tools fill gaps not yet covered by OpenMS/pyopenms. All code in this repo is agentic-only development — written entirely by AI agents.

## Commands

```bash
# Install dependencies for a specific tool
pip install -r tools/proteomics/peptide_analysis/peptide_mass_calculator/requirements.txt

# Lint a specific tool
ruff check tools/proteomics/peptide_analysis/peptide_mass_calculator/

# Run tests for a specific tool
PYTHONPATH=tools/proteomics/peptide_analysis/peptide_mass_calculator python -m pytest tools/proteomics/peptide_analysis/peptide_mass_calculator/tests/ -v

# Lint all tools
ruff check tools/

# Run all tests across all tools
for d in tools/*/*/*/; do PYTHONPATH="$d" python -m pytest "$d/tests/" -v; done

# Run a script directly
python tools/proteomics/peptide_analysis/peptide_mass_calculator/peptide_mass_calculator.py --sequence PEPTIDEK --charge 2
python tools/metabolomics/formula_tools/isotope_pattern_matcher/isotope_pattern_matcher.py --formula C6H12O6
```

## Architecture

### Per-Tool Directory Structure

Each tool is a self-contained directory under `tools/<domain>/<topic>/<tool_name>/`:

```
tools/<domain>/<topic>/<tool_name>/
├── <tool_name>.py        # The tool (importable functions + click CLI)
├── requirements.txt      # pyopenms + script-specific deps
├── README.md             # Usage examples
└── tests/
    ├── conftest.py       # requires_pyopenms marker + sys.path setup
    └── test_<tool_name>.py
```

Domains: `proteomics/`, `metabolomics/`

Proteomics topics: `spectrum_analysis/`, `peptide_analysis/`, `protein_analysis/`, `fasta_utils/`, `file_conversion/`, `quality_control/`, `targeted_proteomics/`, `identification/`, `ptm_analysis/`, `structural_proteomics/`, `specialized/`, `rna/`

Metabolomics topics: `formula_tools/`, `feature_processing/`, `spectral_analysis/`, `compound_annotation/`, `drug_metabolism/`, `isotope_labeling/`, `lipidomics/`, `export/`

### Key Patterns

- pyopenms import wrapped in try/except with user-friendly error message
- Mass-to-charge: `(mass + charge * PROTON) / charge` with `PROTON = 1.007276`
- Every script has dual interface: importable functions + click CLI + `__main__` guard
- Tests use `@requires_pyopenms` skip marker from conftest.py
- File-I/O scripts use synthetic test data generated with pyopenms objects

## Contributing

See `AGENTS.md` for the full AI contributor guide. Two Claude Code skills are available:

- **`contribute-script`** — guided workflow for adding a new script
- **`validate-script`** — validate any script in an isolated venv (ruff + pytest)
