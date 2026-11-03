# Run 43 coordinated original-source replay — November 3, 2026

Four coordinated agents separately reviewed operator consistency, boundary mathematics, solver implementation, and independent physical validation. Each used its own write directory. The source repository, prior baseline, rejected candidates, checkpoints and physical source definitions were preserved. No merge, candidate promotion, acceptance change or Physical Review D submission was performed.

The delivery branch descends from `6aefba3e96ae911cd0895e32abe7b172c4af3914`, retaining the entire prior Run 43 evidence. The original remote refs and complete Git history were also saved to a verified local Git bundle (SHA256 in session-provenance.json). The frozen audit remains Git blob `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`.

## Equations and correction basis

The original versioned DEE action gives `E=(Pi^T K Pi+psi^-4 D^T K D)/2+U` and `S_i=-Pi^T K D_i`. The maximal conformally flat initial-data reduction of the Einstein constraints is

    Delta psi + |LW|^2 psi^-7/8 + 2 pi G E(psi) psi^5 = 0
    div(LW) = Delta W + grad(div W)/3 = 8 pi G psi^6 S

with the original fixed sources, G=.001, periodic angles, and exact Dirichlet psi=1,W=0 on r=.05,.25. These runs solve those initial constraints; they do not establish full Einstein–DEE spacetime or matter evolution. Available equation provenance is the repository's recovered-formulation transcription; no independently inspectable author manuscript archive was found in the selected tree.

The independently checked solver–audit identity is

    (T-E)_i = (sum_j G_j G_j-L) W_i + sum_j [G_j,G_i] W_j.

Smooth Cartesian derivative equivalence does not make these finite-grid operators equal. The frozen endpoint composition gives `(DrDr-Drr)f_1=(-f0+3f1-3f2+f3)/(4h^2)=h f'''/4+O(h^2)`; the outer odd term has the opposite sign. The radial vector symbol is `I+er er^T/3`. The existing boundary repair cancels this radial defect without an adjustable coefficient. No additional third-normal-derivative boundary condition follows from the Dirichlet data. Angular/product-rule discrepancies remain.

The operator agent independently reproduced the full finite-grid coefficient identity at relative errors 5.16e-16 and 5.94e-16. Its exact continuum analytic-field consistency check gives frozen inner-M2 orders .309947903,.795480313,1.104801455. These are analytic operator diagnostics, never substituted for original-source acceptance.

## Fresh numerical results

All three variants were genuinely solved afresh at n=8,16,32,64, from the original sources with no code or parameter changes. Each grid completed two nonlinear iterations and four successful GMRES solves. Separate-process validation reproduces the unchanged gate. All twelve new field NPZ files are byte-identical to their corresponding preserved historical originals.

| Variant | Band | Metric | 8→16 | 16→32 | 32→64 |
|---|---|---|---:|---:|---:|
| baseline | inner | H | 1.740847585 | 1.929825463 | 1.996175084 |
| baseline | inner | M2 | 1.875425682 | 1.970278594 | 1.994731881 |
| baseline | outer | H | 1.745663903 | 1.963989315 | 2.018376693 |
| baseline | outer | M2 | 2.229142986 | 2.132462281 | 2.071256761 |
| fourth | inner | H | 1.740849102 | 1.929826071 | 1.996174863 |
| fourth | inner | M2 | 0.330580315 FAIL | 0.855566995 FAIL | 1.112652071 FAIL |
| fourth | outer | H | 1.745667434 | 1.963989971 | 2.018377225 |
| fourth | outer | M2 | 1.852947017 | 1.970457203 | 1.992557652 |
| local | inner | H | 1.740847590 | 1.929825391 | 1.996175110 |
| local | inner | M2 | 0.234124077 FAIL | 0.727095846 FAIL | 1.116755252 FAIL |
| local | outer | H | 1.745663928 | 1.963989215 | 2.018376479 |
| local | outer | M2 | 1.998783011 | 1.532308823 | 1.394913922 FAIL |

Baseline: **12/12 PASS**, preserved reference. Fourth-order: **9/12 FAIL**, rejected. Local boundary correction: **8/12 FAIL**, rejected. Every order is required to be finite and at least 1.5 in the fixed inner [.05,.10] and outer [.20,.25] bands. No threshold, band, metric, derivative, volume weight, or frozen row selection changed.

All fields are finite with positive psi and exact original radial boundary values. Independent direct four-dimensional stress reconstruction agrees with the original energy and momentum projections to roundoff (maximum energy difference 1.07e-14, momentum difference 1.12e-16). These checks use actual solved metric samples, not manufactured physical forcing.

## Scientific decision

The requested new passing correction is **not established**. Both failed candidates remain rejected. The original-source fourth-order fields independently satisfy `M_frozen=psi^-10(HW-b)+psi^-10(TW-HW)` within 5.15e-18 in two implementations. First-row algebraic residuals are about 1e-13 while operator-mismatch residuals are about 1e-3, so reducing linear tolerances cannot resolve this discrepancy. The boundary correction is derived from the original constraints, but fails four unchanged orders. No source-derived next correction with evidence of improved accuracy and a 12/12 result is established. This is the genuine scientific blocker; there was no execution blocker. A finite-grid audit failure cannot be relabeled as a pass because a separate diagnostic improves. No fitted coefficient, altered source, compensating audit defect or unsupported boundary condition is used. Failure to establish a justified correction is a scientific blocker, not an impossibility proof.

## Evidence and reproduction

Fresh regression: 45 tests passed on untouched sources. All 629 pre-existing source/evidence files were rehashed: zero changed or missing. Original stress sources, Run38 source checkpoint, action transcription, and frozen audit also match the November 2 branch byte for byte. Complete field arrays, all 24 nonlinear checkpoints, GMRES histories, executed source snapshots, logs, settings, source hashes and exit statuses accompany the report. Independent validation checks original-source integrity, actual 4D stress projections, positivity and exact boundary values in addition to the unchanged frozen audit. Reproduction instructions and SHA256 manifests accompany the completed evidence.


The new evidence is under `research/results/coordinated-replay-20261103-1831/`:

- `baseline/`, `fourth/`, `local/`: all genuinely solved arrays, checkpoints, settings, reports and executed solver snapshots.
- `solver/`: exact execution commands, environment, source bundle, timings, residuals, solve/audit exit statuses, preservation checks and replay driver.
- `operator/`: two independent physical row decompositions, complete sparse coefficient matrices, source snapshots, analytic operator diagnostics, logs and hash manifests.
- `boundary/`: original-action-to-constraint derivation and exact rational signed-stencil certificate.
- `validation/`: independently reloaded twelve-check reports, actual 4D stress checks, original-source integrity and bitwise reproduction evidence.
- `coordination/`: original refs, session provenance, unchanged frozen audit snapshot and 45-test regression log.

From a checkout of the delivery commit, use a fresh tag:

```bash
bash research/results/coordinated-replay-20261103-1831/replay-all.sh unique-run43-replay
```

The replay refuses existing solve/output paths, verifies preserved evidence, performs all twelve genuine solves and independent audits, and retains scientific failures. Its final exit 1 is expected when either rejected candidate still fails; solve completion is not scientific acceptance. The unchanged frozen auditor itself returns 0 for baseline and 1 for each failed candidate. Independent-validation process exit 0 records successful diagnostic execution; its JSON PASS/FAIL is the scientific result. No failing JSON is relabeled.

`SHA256SUMS` covers the complete new evidence and this report. The verified original-history Git bundle remains at `/workspace/run43-investigation/coordination/preserved-repository.bundle`, with its digest recorded in session provenance; it is not duplicated into Git because complete original history is already retained in branch ancestry.
