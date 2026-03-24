# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Agentomics is a collection of standalone CLI tools built with [pyopenms](https://pyopenms.readthedocs.io/) for proteomics and metabolomics workflows. These tools fill gaps not yet covered by OpenMS/pyopenms. All code in this repo is agentic-only development — written entirely by AI agents.

## Commands

```bash
# Install dependencies for a specific script
pip install -r scripts/proteomics/peptide_mass_calculator/requirements.txt

# Lint a specific script
ruff check scripts/proteomics/peptide_mass_calculator/

# Run tests for a specific script
PYTHONPATH=scripts/proteomics/peptide_mass_calculator python -m pytest scripts/proteomics/peptide_mass_calculator/tests/ -v

# Lint all scripts
ruff check scripts/

# Run all tests across all scripts
for d in scripts/*/*/; do PYTHONPATH="$d" python -m pytest "$d/tests/" -v; done

# Run a script directly
python scripts/proteomics/peptide_mass_calculator/peptide_mass_calculator.py --sequence PEPTIDEK --charge 2
python scripts/metabolomics/isotope_pattern_matcher/isotope_pattern_matcher.py --formula C6H12O6
```

## Architecture

### Per-Script Directory Structure

Each script is a self-contained directory under `scripts/<domain>/<tool_name>/`:

```
scripts/<domain>/<tool_name>/
├── <tool_name>.py        # The tool (importable functions + argparse CLI)
├── requirements.txt      # pyopenms + script-specific deps
├── README.md             # Usage examples
└── tests/
    ├── conftest.py       # requires_pyopenms marker + sys.path setup
    └── test_<tool_name>.py
```

Domains: `proteomics/`, `metabolomics/`

### Key Patterns

- pyopenms import wrapped in try/except with user-friendly error message
- Mass-to-charge: `(mass + charge * PROTON) / charge` with `PROTON = 1.007276`
- Every script has dual interface: importable functions + argparse CLI + `__main__` guard
- Tests use `@requires_pyopenms` skip marker from conftest.py
- File-I/O scripts use synthetic test data generated with pyopenms objects

## Contributing

See `AGENTS.md` for the full AI contributor guide. Two Claude Code skills are available:

- **`contribute-script`** — guided workflow for adding a new script
- **`validate-script`** — validate any script in an isolated venv (ruff + pytest)
