---
name: contribute-script
description: Guide creation of a new pyopenms script contribution — scaffolding through validation
---

# Contribute Script

Guide an AI agent through creating a new pyopenms CLI tool for the agentomics repo. Follow every step — this is a rigid skill.

## Prerequisites

Read `AGENTS.md` in the repo root for the full contributor guide and code patterns.

## Steps

### 1. Understand the tool

Ask the user:
- What does this tool do? What pyopenms functionality does it use?
- What gap in OpenMS/pyopenms does it fill?

### 2. Determine the domain

Ask: Is this a **proteomics** or **metabolomics** tool? If neither fits, discuss whether a new domain directory is needed.

### 3. Pick a name

Choose a descriptive snake_case name for the tool (e.g. `peptide_mass_calculator`, `isotope_pattern_matcher`). Confirm with the user.

### 4. Create a feature branch

```bash
git checkout -b add/<tool_name>
```

### 5. Scaffold the directory

```bash
mkdir -p scripts/<domain>/<tool_name>/tests
```

Create these files:

**`requirements.txt`:**
```
pyopenms
```
Add any additional dependencies the script needs (one per line, no version pins).

**`tests/conftest.py`:**
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

### 6. Write the script

Create `scripts/<domain>/<tool_name>/<tool_name>.py` following these patterns:

- Module-level docstring with description, supported features, and CLI usage examples
- pyopenms import guard:
  ```python
  try:
      import pyopenms as oms
  except ImportError:
      sys.exit("pyopenms is required. Install it with:  pip install pyopenms")
  ```
- `PROTON = 1.007276` constant where mass-to-charge calculations are needed
- Importable functions as the primary interface (with type hints and numpy-style docstrings)
- `main()` function with argparse CLI
- `if __name__ == "__main__": main()` guard

### 7. Write tests

Create `scripts/<domain>/<tool_name>/tests/test_<tool_name>.py`:

- Import `requires_pyopenms` from conftest
- Decorate test classes with `@requires_pyopenms`
- Use `from <tool_name> import <function>` inside test methods
- For file-I/O scripts: generate synthetic data using pyopenms objects in test fixtures, write to `tempfile.TemporaryDirectory()`
- Cover: basic functionality, edge cases, key parameters

### 8. Write README

Create `scripts/<domain>/<tool_name>/README.md` with a brief description and CLI usage examples.

### 9. Validate

Invoke the `validate-script` skill on the new script directory. Both ruff and pytest must pass.

### 10. Commit

```bash
git add scripts/<domain>/<tool_name>/
git commit -m "Add <tool_name>: <brief description>"
```
