# AGENTS.md — AI Contributor Guide

This file instructs AI agents (Claude Code, GitHub Copilot, Cursor, Gemini, etc.) how to contribute tools to the agentomics repository.

## Project Purpose

Agentomics is a collection of standalone CLI tools built with [pyopenms](https://pyopenms.readthedocs.io/) for proteomics and metabolomics workflows. These tools fill gaps not yet covered by OpenMS/pyopenms. All code in this repo is written by AI agents.

## Contribution Requirements

Every tool must be a **self-contained directory** under `scripts/<domain>/<topic>/<tool_name>/`:

```
scripts/<domain>/<topic>/<tool_name>/
├── <tool_name>.py        # The tool itself
├── requirements.txt      # pyopenms + any tool-specific deps (no version pins)
├── README.md             # Brief description + CLI usage examples
└── tests/
    ├── conftest.py       # Shared test config (see below)
    └── test_<tool_name>.py
```

### Topics

**Proteomics topics:** `spectrum_analysis/`, `peptide_analysis/`, `protein_analysis/`, `fasta_utils/`, `file_conversion/`, `quality_control/`, `targeted_proteomics/`, `identification/`, `ptm_analysis/`, `structural_proteomics/`, `specialized/`, `rna/`

**Metabolomics topics:** `formula_tools/`, `feature_processing/`, `spectral_analysis/`, `compound_annotation/`, `drug_metabolism/`, `isotope_labeling/`, `lipidomics/`, `export/`

### Rules

- `<domain>` is `proteomics` or `metabolomics`
- `<topic>` is one of the topic directories listed above
- `requirements.txt` always includes `pyopenms` with no version pin — builds against latest
- No cross-tool imports — each tool is fully independent
- No `__init__.py` files — these are NOT Python packages
- No tools that duplicate functionality already in OpenMS/pyopenms TOPP tools

## Code Patterns

### Tool structure

Every tool must have:

1. **Module docstring** with description, features, and usage examples
2. **pyopenms import guard:**
   ```python
   import sys
   try:
       import pyopenms as oms
   except ImportError:
       sys.exit("pyopenms is required. Install it with:  pip install pyopenms")
   ```
3. **Importable functions** as the primary interface (with type hints and numpy-style docstrings)
4. **`main()` function** with argparse CLI
5. **`if __name__ == "__main__": main()`** guard
6. **`PROTON = 1.007276`** constant where mass-to-charge calculations are needed

### Test structure

Every `tests/conftest.py` must contain:

```python
import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    import pyopenms  # noqa: F401
    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False

requires_pyopenms = pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not installed")
```

Test files:
- Decorate test classes with `@requires_pyopenms` from conftest
- Import tool functions inside test methods: `from <tool_name> import <function>`
- For file-I/O tools: generate synthetic data using pyopenms objects, write to `tempfile.TemporaryDirectory()`

## Validation

Every tool must pass validation in an **isolated venv** before it can be merged. Run these commands from the repo root:

```bash
TOOL_DIR=scripts/<domain>/<topic>/<tool_name>
VENV_DIR=$(mktemp -d)
python -m venv "$VENV_DIR"
"$VENV_DIR/bin/python" -m pip install -r "$TOOL_DIR/requirements.txt"
"$VENV_DIR/bin/python" -m pip install pytest ruff
"$VENV_DIR/bin/python" -m ruff check "$TOOL_DIR/"
PYTHONPATH="$TOOL_DIR" "$VENV_DIR/bin/python" -m pytest "$TOOL_DIR/tests/" -v
rm -rf "$VENV_DIR"
```

Both ruff and pytest must pass with zero errors.

## Linting

Ruff is configured in `ruff.toml` at the repo root:
- Line length: 120
- Rules: E (pycodestyle errors), F (pyflakes), W (pycodestyle warnings), I (isort)

## What NOT to Do

- Do not add cross-tool imports
- Do not add dependencies to a shared/root requirements file
- Do not create tools that duplicate existing pyopenms CLI tools or OpenMS TOPP tools
- Do not pin pyopenms to a specific version
- Do not add `__init__.py` files
- Do not add tools that don't actually use pyopenms (pure stats/math tools belong elsewhere)
