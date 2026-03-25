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

### 3. Pick a topic

Choose the topic directory under the selected domain from the options documented in `AGENTS.md`. Confirm the topic with the user before scaffolding files.

### 4. Pick a name

Choose a descriptive snake_case name for the tool (e.g. `peptide_mass_calculator`, `isotope_pattern_matcher`). Confirm with the user.

### 5. Create a feature branch

```bash
git checkout -b add/<tool_name>
```

### 6. Scaffold the directory

```bash
mkdir -p tools/<domain>/<topic>/<tool_name>/tests
```

Create these files:

**`requirements.txt`:**
```
pyopenms
click
```
Add any additional dependencies the script needs (one per line, no version pins).

**`tests/conftest.py`:**
```python
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
```

### 7. Write the script

Create `tools/<domain>/<topic>/<tool_name>/<tool_name>.py` following these patterns:

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
- `main()` function with click CLI
- `if __name__ == "__main__": main()` guard

### 8. Write tests

Create `tools/<domain>/<topic>/<tool_name>/tests/test_<tool_name>.py`:

- Add `pytest.importorskip("pyopenms")` at module level (after `import pytest`)
- Use `from <tool_name> import <function>` inside test methods
- For file-I/O scripts: generate synthetic data using pyopenms objects in test fixtures, write to `tempfile.TemporaryDirectory()`
- Cover: basic functionality, edge cases, key parameters

### 9. Write README

Create `tools/<domain>/<topic>/<tool_name>/README.md` with a brief description and CLI usage examples.

### 10. Validate

Invoke the `validate-script` skill on the new script directory. Both ruff and pytest must pass.

### 11. Commit

```bash
git add tools/<domain>/<topic>/<tool_name>/
git commit -m "Add <tool_name>: <brief description>"
```
