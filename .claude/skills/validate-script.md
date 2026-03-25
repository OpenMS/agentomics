---
name: validate-script
description: Validate a pyopenms script in an isolated venv — runs ruff lint and pytest
---

# Validate Script

Validate any script in the agentomics repo by running ruff and pytest in a fresh isolated venv.

## Steps (follow exactly — rigid skill)

1. **Identify the script directory.** If the user provided a path, use it. Otherwise, ask which script to validate. The path should be `tools/<domain>/<topic>/<tool_name>/`.

2. **Verify the directory structure.** Confirm it contains:
   - `<tool_name>.py`
   - `requirements.txt`
   - `tests/` directory with at least one `test_*.py` file

3. **Create a temporary venv and run validation.** Execute these commands:

   ```bash
   SCRIPT_DIR=<path-to-script-directory>
   VENV_DIR=$(mktemp -d)
   python -m venv "$VENV_DIR"
   "$VENV_DIR/bin/python" -m pip install -r "$SCRIPT_DIR/requirements.txt"
   "$VENV_DIR/bin/python" -m pip install pytest ruff
   "$VENV_DIR/bin/python" -m ruff check "$SCRIPT_DIR/"
   PYTHONPATH="$SCRIPT_DIR" "$VENV_DIR/bin/python" -m pytest "$SCRIPT_DIR/tests/" -v
   rm -rf "$VENV_DIR"
   ```

4. **Report results.** Summarize pass/fail for both ruff and pytest. If either fails, show the relevant error output so the user can fix it.

5. **Clean up.** Ensure the temporary venv is removed even if validation fails.
