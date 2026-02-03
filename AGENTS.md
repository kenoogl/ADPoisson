# Repository Guidelines

## Project Structure & Module Organization
This repository is currently a Kiro-style, spec-driven workspace. Most content is documentation rather than application code.

- `.kiro/specs/<feature>/`: Feature specs (e.g., `.kiro/specs/adpoisson/requirements.md`, `spec.json`).
- `.kiro/settings/`: Kiro templates and rules used to generate or validate specs.
- `.gemini/commands/kiro/`: CLI command definitions for spec workflows.
- `README.md`, `GEMINI.md`: High-level project context and process notes.

There is no application source tree or tests directory yet. If you add code, create a clear top-level folder (for example, `src/` and `tests/`) and document it here.

## Build, Test, and Development Commands
No build/test scripts or package manifests are present (no `package.json`, `pyproject.toml`, `go.mod`, etc.). If you introduce a build or runtime, document the exact commands and expected outputs here and in `README.md`.

## Coding Style & Naming Conventions
- Documentation is Markdown. Keep headings consistent and concise.
- Spec feature folders are lowercase and short (e.g., `adpoisson`).
- For spec documents, follow the language defined in the feature’s `spec.json` (currently `"language": "ja"`).

## Testing Guidelines
No testing framework is configured. If you add executable code, add a corresponding test setup and document:
- test directory location
- test naming conventions
- the command to run tests

## Commit & Pull Request Guidelines
Git history contains only `"first commit"`, so no established commit convention exists. Use short, imperative messages and include scope when helpful (e.g., `spec: add adpoisson requirements`).

For pull requests, include:
- a brief summary of changes
- the spec path being updated (for example, `.kiro/specs/adpoisson/requirements.md`)
- any required phase approvals (requirements → design → tasks)

## Process Notes
Follow the Kiro workflow described in `GEMINI.md`. Specs move through requirements, design, and tasks before implementation.
