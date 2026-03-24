# AI Contributor Skills & Plans — Design Spec

## Purpose

Define the skills, contributor docs, and CI pipeline that enable AI agents to contribute validated, self-contained pyopenms scripts to the agentomics repo. Every contribution must build against the latest pyopenms, resolve its own dependencies, and pass linting + tests in isolation.

## Per-Script Directory Structure

Every script is a self-contained package under `scripts/<domain>/<tool_name>/`:

```
scripts/proteomics/peptide_mass_calculator/
├── peptide_mass_calculator.py
├── requirements.txt
├── README.md
└── tests/
    ├── conftest.py
    └── test_peptide_mass_calculator.py
```

Rules:

- `requirements.txt` always includes `pyopenms` with no version pin (builds against latest) plus any script-specific dependencies
- Tests live inside each script's own `tests/` directory (not a top-level `tests/`)
- Each script is self-contained — no cross-script imports
- These directories are NOT Python packages — no `__init__.py` files
- Each `tests/` directory includes a `conftest.py` that defines the `requires_pyopenms` marker. The `conftest.py` also adds the parent script directory to `sys.path` as a fallback for cases where `PYTHONPATH` is not set (e.g. running `pytest` directly without the validation wrapper):
  ```python
  import sys, os, pytest
  sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
  try:
      import pyopenms
      HAS_PYOPENMS = True
  except ImportError:
      HAS_PYOPENMS = False
  requires_pyopenms = pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not installed")
  ```
- Existing 7 scripts live on the `origin/copilot/add-agentic-scripts-for-proteomics` branch (not yet on `main`). Migration involves copying the files from that branch and restructuring them into per-script directories. The existing monolithic test files (`tests/test_proteomics.py` covering 4 scripts, `tests/test_metabolomics.py` covering 3 scripts) must be decomposed so each script's test class moves into its own `tests/test_<tool_name>.py`. The single root `requirements.txt` is replaced by per-script `requirements.txt` files

### Test Data Strategy

Scripts that process files (mzML, featureXML) need test data. Strategy:

- **Prefer synthetic data** — generate minimal test data in fixtures using pyopenms objects (e.g. `MSExperiment`, `FeatureMap`) directly in test setup
- If synthetic generation is impractical, include small test fixture files in the script's `tests/` directory
- File-I/O tests that cannot run without external data use an additional `@pytest.mark.skipif` marker (e.g. `requires_test_data`)

## Validation Pipeline

For any given script directory, validation runs:

```bash
python -m venv /tmp/validate_<tool>
/tmp/validate_<tool>/bin/python -m pip install -r <script>/requirements.txt
/tmp/validate_<tool>/bin/python -m pip install pytest ruff
/tmp/validate_<tool>/bin/python -m ruff check <script>/
PYTHONPATH=<script> /tmp/validate_<tool>/bin/python -m pytest <script>/tests/ -v
rm -rf /tmp/validate_<tool>
```

Note: Commands use direct venv binary paths (`/tmp/.../bin/python`) instead of `source activate` to avoid platform-specific activation scripts. On Windows, substitute `bin/` with `Scripts/`. CI runs on Ubuntu so this is not a concern there.

Rules:

- Each script is validated in isolation — its own fresh venv, no shared state
- `PYTHONPATH=<script>` is the primary mechanism for test imports; `conftest.py` provides a `sys.path` fallback for direct `pytest` invocation
- ruff configuration lives in repo root as `ruff.toml`:
  ```toml
  line-length = 120
  [lint]
  select = ["E", "F", "W", "I"]
  ```
  (pycodestyle errors/warnings, pyflakes, isort)
- The same validation logic is used by: the Claude Code `validate-script` skill, the GitHub Actions CI pipeline, and documented in `AGENTS.md` for other AI agents to replicate
- Validation must pass before a contribution is considered complete

## Claude Code Skills

Two skills in `.claude/skills/`:

### `contribute-script`

Guides an AI agent through creating a new script end-to-end. Rigid — follow exactly, no skipping steps.

1. **Ask what the tool does** — what pyopenms functionality does it wrap, what gap does it fill
2. **Determine domain** — proteomics or metabolomics (or prompt if a new domain is needed)
3. **Scaffold directory** — create `scripts/<domain>/<tool_name>/` with `requirements.txt`, empty `README.md`, empty test file
4. **Write the script** — following established patterns:
   - pyopenms try/except import with user-friendly error message
   - `PROTON = 1.007276` constant where mass-to-charge calculations are needed
   - Importable functions as primary interface
   - `main()` function with argparse CLI
   - `if __name__ == "__main__"` guard
   - Type hints in function signatures
   - Numpy-style docstrings
5. **Write tests** — pytest with `@requires_pyopenms` marker, covering core functionality
6. **Write README** — CLI usage examples
7. **Validate** — invoke the `validate-script` logic (fresh venv, ruff, pytest)
8. **Commit** — on a feature branch named `add/<tool_name>`

### `validate-script`

Standalone validation — can be invoked on any script directory. Rigid.

1. Detect target script directory (from argument or prompt user)
2. Create temporary venv
3. Install from the script's `requirements.txt`
4. Install `pytest` and `ruff`
5. Run `ruff check` on the script directory
6. Run `python -m pytest` on the script's `tests/` directory
7. Report results — pass/fail with details
8. Clean up temporary venv

## AGENTS.md

Platform-agnostic contributor guide at repo root for any AI agent (Copilot, Cursor, Gemini, etc.). Contents:

1. **Project purpose** — agentic-only pyopenms tools for proteomics/metabolomics that don't yet exist in OpenMS
2. **Contribution requirements:**
   - Self-contained directory under `scripts/<domain>/<tool_name>/`
   - Must include: script `.py`, `requirements.txt`, `README.md`, `tests/` with pytest tests
   - Must use latest pyopenms (no version pinning)
   - Must pass ruff + pytest in an isolated venv
3. **Code patterns to follow:**
   - pyopenms import guard (try/except with install message)
   - Dual interface: importable functions + argparse CLI + `__main__` guard
   - Numpy-style docstrings, type hints
   - `@requires_pyopenms` test skip marker
4. **Validation steps** — exact shell commands to create a venv, install, lint, test
5. **What not to do:**
   - No cross-script imports
   - No adding deps to a shared requirements file
   - No scripts that duplicate existing pyopenms functionality

## GitHub Actions CI

`.github/workflows/validate.yml`:

- **Trigger:** Pull requests that touch anything under `scripts/`
- **Detection job:** Diffs against base branch to identify changed script directories, outputs them as a JSON matrix. Outputs a `has_changes` flag — the validation matrix is conditional on this flag so PRs that only touch non-script files don't produce an empty matrix error.
- **Validation matrix:** For each changed script directory, a parallel job that:
  1. Checks out the repo
  2. Sets up Python 3.11 (pinned for pyopenms wheel availability — update when pyopenms supports newer versions)
  3. Creates a venv
  4. `pip install -r <script>/requirements.txt`
  5. `pip install pytest ruff`
  6. `ruff check <script>/`
  7. `PYTHONPATH=<script> python -m pytest <script>/tests/ -v`
- **Designed for branch protection** — can be set as a required status check to block merges on failure
- Uses the same `ruff.toml` at repo root as local validation

## Updated CLAUDE.md

After implementation, `CLAUDE.md` must reflect:

- New per-script directory structure and how to navigate it
- Per-script test commands: `PYTHONPATH=scripts/<domain>/<tool> python -m pytest scripts/<domain>/<tool>/tests/ -v`
- Reference to the two Claude Code skills (`contribute-script`, `validate-script`)
- Reference to `AGENTS.md` for the full contributor guide
- Ruff lint command: `ruff check scripts/`

## Deliverables

1. Merge content from `origin/copilot/add-agentic-scripts-for-proteomics` branch, then restructure into per-script directories
2. Create `ruff.toml` at repo root
3. Create per-script `conftest.py` with `requires_pyopenms` marker and `PYTHONPATH`/`sys.path` setup
4. Create `.claude/skills/contribute-script.md`
5. Create `.claude/skills/validate-script.md`
6. Create `AGENTS.md` at repo root
7. Create `.github/workflows/validate.yml`
8. Update `CLAUDE.md` to reference new structure and skills
9. Update root `README.md` to reflect new structure
