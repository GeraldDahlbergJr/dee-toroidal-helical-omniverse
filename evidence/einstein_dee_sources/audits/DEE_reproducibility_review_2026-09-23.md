# DEE reproducibility review — 23 September 2026

Prepared for Gerald Dahlberg Jr. from `DEE_LaTeX_v2(5).zip` and the public repository [`dee-toroidal-helical-omniverse`](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse), `main` at commit `7e988c0` (22 September 2026). This is a provenance check, not scientific validation or a full referee report.

## What is independently checkable now

| Claim or calculation | Result | Limit |
| --- | --- | --- |
| Manuscript Eq. (18), illustrative parameters in the “Illustrative normalized point” section | Independently recalculated: ρ² = 0.945; ρ = 0.972111; effective radial mass squared = 7.56; kinetic matrix `[[1.2835, 0.0945], [0.0945, 1.189]]`; determinant = 1.51715125. This agrees with the rounded text. | Algebraic homogeneous branch only; no coupled stability result. |
| `Weak-field_scaling_check.csv` versus `Linearized_Einstein_solve_diagnostics.csv` | All six rows agree with multiplication by `10⁻⁶`, subject to rounding in the scaling CSV. Largest listed scaled metric component is 0.001157492, consistent with the manuscript’s rounded `1.16 × 10⁻³`. | Checks arithmetic on archived outputs, not how the unscaled fields were computed. |
| GitHub tests | `python -m unittest discover -s tests -v`: 10 tests passed in the inspected checkout. These cover formulas, kinetic matrix and NEC diagnostic, radial solver boundary behavior, and a smooth-source grid refinement example. | They do not load or regenerate the 15 archived CSVs. |

## The provenance gap

The ZIP contains `main.tex`, `main.pdf`, an earlier conversion source, and 15 archived exploratory CSVs. The GitHub repository contains `dee_field_model.py`, `nec_scan.py`, `linearized_einstein.py`, `dee_evolution.py`, `dee_geometry.py`, and tests. No checked-in program writes the archived CSVs, and there is no recorded run configuration tying a manuscript value to a command and Git commit. The manuscript itself calls the CSVs archival working outputs and says their numerical values require independent reproduction before journal submission.

| Archived results | What is needed to rerun them |
| --- | --- |
| `Linearized_Einstein_solve_diagnostics.csv`, `Weak-field_scaling_check.csv` | The radial grid and outer radius; explicit tensor source profiles and Gaussian envelope; field, phase, and coupling parameters; normalization of coordinates and `G`; component basis and NEC null vectors; outer boundary conditions; saved raw profiles, residuals, and refinement runs. The current radial integrator takes a supplied source but does not construct this source from the manuscript fields. |
| `Direct_null-vector_NEC_test_with_DEE.csv`, `DEE_phase_NEC_summary.csv`, both compensated/worst-null searches, reduced-metric and flat-cylinder condition tables | Exact metric functions and switching profile, sampled parameter domains and seeds, tetrad/null-vector normalization, curvature and stress-tensor computation, numerical differentiation method, sampling meshes, and condition thresholds. No corresponding generator is in the repository. |
| `Dynamic_collapse_test.csv`, `Ten_closest_approaches.csv`, `Broad-pass_summary.csv` | Initial conditions, ranges or seed, integrator and tolerances, collapse threshold, detailed geometry or force model, and classification rules. `dee_evolution.py` is a different driven relaxation demonstration, not an implementation of manuscript Eq. (31) or these archived scans. |

## Review-ready status

The analytic example and the **arithmetic of the published weak-field rescaling** are checked. The archived unscaled gravitational calculation and the reduced-model scans remain **unreproduced from available source**. The phrase “Linearized Einstein–DEE consistency test: PASS” describes an exploratory checkpoint under the specified scaling; it should not be presented as an independently reproduced numerical result or a nonlinear result. The manuscript already says that nonlinear non-collapse remains unproved.

Recommended sequence:

1. Add a generator for one clearly specified weak-field case, with a versioned configuration and raw radial profiles. Derive the source from the stated matter model, or explicitly label a chosen source as illustrative.
2. Rerun at several radial resolutions and report residuals of Eq. (34), boundary behavior, and a convergence table. Compare newly generated diagnostics with the archival CSVs; preserve any discrepancies instead of fitting undocumented parameters to the old numbers.
3. Rebuild the remaining numerical experiments as separate, named models, each with a command, configuration, output, and scope statement. If a run cannot be reconstructed, remove its precise numbers from the review manuscript or label them as historical exploratory observations.
4. Update the manuscript’s data/code statement only after each reported value points to its reproducing run. Check derivations, dimensions, references, and the journal’s disclosure requirements before submission.

Nothing in this audit proves the model physically valid. It establishes which parts of the current package can be checked and what work remains to make the numerical claims reviewable.
