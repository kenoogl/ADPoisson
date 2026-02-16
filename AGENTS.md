# AGENT.md

## Purpose
This repository is designed for AI-assisted research analysis.

## Global Rules
- config.yaml is ground truth for experiments.
- run_summary.json contains results only.
- Do not modify results files.
- Separate facts from evaluation.

---

## Experimental Analysis Mode

### Priority Files
1. intent.md
2. experiments/<exp>/config.yaml
3. results/<exp>/run_summary.json
4. logs/<exp>.json
5. experiments/<exp>/notes.md (human thoughts, optional)
6. results/<exp>/ (generated results)

### AI Operating Rules
1. config.yaml is the single source of truth for experiment settings.
2. run_summary.json contains result summary only.
3. Do not infer solver settings from run_summary.json.
4. When comparing experiments, use run_summary.json fields.
5. Prioritize objective metrics over human notes.
6. Distinguish facts from evaluation.

### Analysis Rules
- Compare run_summary across experiments.
- Check convergence flag first.
- Avoid using human notes as ground truth.

---

## Code Development Mode

### Project Structure & Module Organization
This repository is currently a Kiro-style, spec-driven workspace. Most content is documentation rather than application code.

- `.kiro/specs/<feature>/`: Feature specs (e.g., `.kiro/specs/adpoisson/requirements.md`, `spec.json`).
- `.kiro/settings/`: Kiro templates and rules used to generate or validate specs.
- `.gemini/commands/kiro/`: CLI command definitions for spec workflows.
- `README.md`, `GEMINI.md`: High-level project context and process notes.

There is no application source tree or tests directory yet. If you add code, create a clear top-level folder (for example, `src/` and `tests/`) and document it here.

### Build, Test, and Development Commands
No build/test scripts or package manifests are present (no `package.json`, `pyproject.toml`, `go.mod`, etc.). If you introduce a build or runtime, document the exact commands and expected outputs here and in `README.md`.

### Coding Style & Naming Conventions
- Documentation is Markdown. Keep headings consistent and concise.
- Spec feature folders are lowercase and short (e.g., `adpoisson`).
- For spec documents, follow the language defined in the feature’s `spec.json` (currently `"language": "ja"`).

### Testing Guidelines
No testing framework is configured. If you add executable code, add a corresponding test setup and document:
- test directory location
- test naming conventions
- the command to run tests

### Commit & Pull Request Guidelines
Git history contains only `"first commit"`, so no established commit convention exists. Use short, imperative messages and include scope when helpful (e.g., `spec: add adpoisson requirements`).
For pull requests, include:
- a brief summary of changes
- the spec path being updated (for example, `.kiro/specs/adpoisson/requirements.md`)
- any required phase approvals (requirements → design → tasks)

### Source of Truth
- src/ contains implementation.
- tests/ define correctness.
- config.yaml defines runtime parameters.

### Rules
- Do not modify config structure unless requested.
- Preserve reproducibility.
- Avoid breaking run_exp interface.
- Suggest refactoring only when complexity increases.
- Perform write/execute actions only after explicit user intent such as "実行して", "作成して", "適用して", or equivalent.


### Process Notes
Follow the Kiro workflow described in `GEMINI.md`. Specs move through requirements, design, and tasks before implementation.

### Execution Policy
- Do not create files, edit files, or run commands when the user asks only to "show", "review", or "confirm".

---

## When Both Apply
- Never alter experimental results to fix code.
- Code fixes require new experiment runs.

