# Verification record

- Recovered local HEAD and remote checkpoint branch both equal
  `4f8bcbe8e2ad2e33ce4ac100eb70b6d2c0b4ac73`.
- Prior product-rule report/script/result SHA256 manifest: all three PASS.
- Exact four-node certificate: PASS under SymPy 1.14.0. Rational matrix
  equalities tested on every STF basis column; symbolic diagonal-norm
  obstruction verified; all six normal Dirichlet Gram minors positive.
- New report/script/result SHA256 manifest: PASS.
- `git diff --exit-code 4f8bcbe8e2ad2e33ce4ac100eb70b6d2c0b4ac73 -- .`:
  PASS before staging additions; no tracked checkpoint file changed.
- Frozen audit Git blob unchanged:
  `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`.
- `git diff --check`: PASS. No solver imports/executions, nonlinear solves,
  GMRES, full grids, agents, workflow dispatches, or PRD submissions.
- New branch: `research/run43-adjoint-flux-necessary-20261105`.
  Existing workflow push branch filters do not select this branch. Commit
  includes `[skip ci]` as an additional explicit instruction not to run CI.

Remote publication is verified after committing by matching `git ls-remote`
for the new branch to the local commit SHA and fetching the exact commit through
GitHub. The resulting SHA is returned in the session; this record does not
claim a self-referential final commit hash.
