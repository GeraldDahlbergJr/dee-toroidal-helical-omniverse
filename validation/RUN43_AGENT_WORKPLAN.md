# Run 43 independent validation — work allocation

Status: VALIDATION PENDING. The saved candidate's COMPLETE marker is **not** proof of inner-boundary convergence or physical validity.

## Agent 1 — Inner-boundary convergence
- Identify the inner-boundary observable, mesh spacings and consistent boundary discretization.
- Calculate observed order p = log(|Q_h-Q_h/2|/|Q_h/2-Q_h/4|)/log(2) only for genuinely comparable grids.
- Require p >= 1.5 if this is the predeclared acceptance threshold; report uncertainty and grid independence.
- If the required three-grid data are missing, report INDETERMINATE, not PASS.

## Agent 2 — Residual and stability
- Inspect actual residual history, tolerance, boundary constraint violations, solver exit conditions and numerical conditioning.
- Compare norm definitions and scaling across meshes. Report failures and missing data.

## Agent 3 — Reproducibility
- Preserve run ID, commit SHA, environment, parameters, checkpoint identifiers, hashes and artifact provenance.
- Reproduce diagnostics independently without overwriting original artifacts.
- Compare outputs within documented numerical tolerances.

## Agent 4 — PRD readiness
- Compile methods, figures, numerical error budget, limitations and citation trail.
- Do not claim peer review, novelty, physical validation or PRD acceptance prematurely.

## Decision gate
Only label Run 43 NUMERICALLY VALIDATED when Agents 1–3 have traceable evidence and explicit PASS results. PRD submission readiness is a separate review. Keep all failures and inconclusive findings visible.

## Source
Completed candidate continuation workflow run 38463009735, preserved artifact run43-auto-stage-38463009735-1. Inspect the actual artifact before any numerical conclusion.
