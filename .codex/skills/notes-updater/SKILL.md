# notes-updater

## Purpose
Update `experiments/<exp>/notes.md` for experiment analysis while preserving human decisions.

## Inputs
- Required: `<exp_name>`
- Files:
  - `experiments/<exp>/config.yaml`
  - `results/<exp>/run_summary.json`
  - `ai_context/intent.md` (optional)
  - `experiments/<exp>/notes.md`

## Rules
1. Read only the files above.
2. Ensure notes header is exactly:
   - `# Experiment Notes: <exp>`
   - Never leave `<exp>` placeholder.
3. Update only:
   - `## 1. AI Summary (Facts)`
   - `## 2. AI Analysis (Evaluation)`
4. Never modify:
   - `## 3. Human Thoughts (Decision)`
5. Output language: Japanese.
6. Separate facts and evaluation.
7. If required section markers are missing, stop and report error (do not rewrite entire file).

## Output Structure
- `## 1. AI Summary (Facts)`
  - Run settings (from config)
  - runtime / iterations / converged / residual / error (from run_summary)
  - anomalies (fact-only)
- `## 2. AI Analysis (Evaluation)`
  - Alignment with Intent
  - Hypotheses
  - Comparison (with known baseline if available)
  - Minimal Next Experiment (exactly one)

## Procedure
1. Validate file existence (`config.yaml`, `run_summary.json`, `notes.md`).
2. Parse config and summary.
3. Verify notes structure and header.
4. Replace section 1/2 block only.
5. Keep section 3 unchanged.
6. Report updated file path.

## Failure Conditions
- Missing config/summary/notes
- Broken notes structure (missing section markers)
- Header mismatch that cannot be safely fixed