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
    └── test_peptide_mass_calculator.py
```

Rules:

- `requirements.txt` always includes `pyopenms` with no version pin (builds against latest) plus any script-specific dependencies
- Tests live inside each script's own `tests/` directory (not a top-level `tests/`)
- Each script is self-contained — no cross-script imports
- Tests use the `@requires_pyopenms` skip marker pattern (`pytest.mark.skipif`)
- Existing 7 scripts are migrated to this structure

## Validation Pipeline

For any given script directory, validation runs:

```bash
python -m venv /tmp/validate_<tool>
source /tmp/validate_<tool>/bin/activate
pip install -r <script>/requirements.txt
pip install pytest ruff
ruff check <script>/
python -m pytest <script>/tests/ -v
deactivate
rm -rf /tmp/validate_<tool>
```

Rules:

- Each script is validated in isolation — its own fresh venv, no shared state
- ruff configuration lives in repo root as `ruff.toml` (line length 120, standard Python rules)
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
8. **Commit** — on a feature branch

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
- **Detection job:** Diffs against base branch to identify changed script directories, outputs them as a JSON matrix
- **Validation matrix:** For each changed script directory, a parallel job that:
  1. Checks out the repo
  2. Sets up Python (latest stable 3.x)
  3. Creates a venv
  4. `pip install -r <script>/requirements.txt`
  5. `pip install pytest ruff`
  6. `ruff check <script>/`
  7. `python -m pytest <script>/tests/ -v`
- **Designed for branch protection** — can be set as a required status check to block merges on failure
- Uses the same `ruff.toml` at repo root as local validation

## Deliverables

1. Migrate existing 7 scripts to per-script directory structure
2. Create `ruff.toml` at repo root
3. Create `.claude/skills/contribute-script.md`
4. Create `.claude/skills/validate-script.md`
5. Create `AGENTS.md` at repo root
6. Create `.github/workflows/validate.yml`
7. Update `CLAUDE.md` to reference new structure and skills
8. Update root `README.md` to reflect new structure
