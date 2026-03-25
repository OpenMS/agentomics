# Design: Rename scripts/ to tools/ and Migrate argparse to click

**Date:** 2026-03-25
**Status:** Approved

## Summary

Rename the `scripts/` directory to `tools/` and convert all 123 CLI tools from `argparse` to `click` for a cleaner, more declarative CLI interface.

## Phase 1 — Directory Rename

Atomic rename via `git mv scripts/ tools/`.

Update all references to `scripts/` in:
- `CLAUDE.md` — directory structure, commands, architecture
- `AGENTS.md` — contribution requirements, validation, code patterns
- `README.md` — tool structure, catalog links, examples
- `.github/workflows/validate.yml` — CI paths and grep pattern
- `.claude/skills/contribute-script.md` — scaffolding paths
- `.claude/skills/validate-script.md` — validation paths
- `docs/superpowers/specs/2026-03-24-ai-contributor-skills-design.md`
- `docs/superpowers/plans/2026-03-24-ai-contributor-skills.md`

Single commit: "Rename scripts/ to tools/"

## Phase 2 — Click Migration

### Conversion Rules

Each tool's `main()` function is converted from argparse to click decorators:

| argparse pattern | click equivalent |
|---|---|
| `import argparse` | `import click` |
| `argparse.ArgumentParser(description="...")` | `@click.command()` |
| `parser.add_argument("--foo", required=True, help="...")` | `@click.option("--foo", required=True, help="...")` |
| `parser.add_argument("--foo", type=int, default=1)` | `@click.option("--foo", type=int, default=1)` |
| `parser.add_argument("--foo", action="store_true")` | `@click.option("--foo", is_flag=True)` |
| `args = parser.parse_args()` | removed — click injects as function params |
| `args.foo` | `foo` (function parameter) |
| `if __name__ == "__main__": main()` | unchanged — click commands are callable |

### What stays the same

- Module docstrings
- pyopenms import guard
- Importable functions (the library interface)
- Test structure and conftest.py
- `PROTON` constant where used
- `if __name__ == "__main__": main()` guard

### Requirements update

Each tool's `requirements.txt` gets `click` added.

### Execution

Two parallel agents, one per domain:
- **Agent 1:** `tools/proteomics/` (~89 tools)
- **Agent 2:** `tools/metabolomics/` (~34 tools)

No file conflicts since domains don't overlap.

## Phase 3 — CI and Docs Update

- Add `click` to the shared venv install in `.github/workflows/validate.yml`
- Update AGENTS.md code pattern: "`main()` function with click CLI"
- Update README.md: "click interface" instead of "argparse interface"
- Update CLAUDE.md: click references in architecture section

## Directory Structure (after)

```
tools/<domain>/<topic>/<tool_name>/
├── <tool_name>.py        # The tool (importable functions + click CLI)
├── requirements.txt      # pyopenms + click + tool-specific deps
├── README.md             # Usage examples
└── tests/
    ├── conftest.py       # requires_pyopenms marker + sys.path setup
    └── test_<tool_name>.py
```
