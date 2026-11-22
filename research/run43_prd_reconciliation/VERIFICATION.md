# Run 43 verification matrix

Equation and table identifiers are those of the recovered referee report. Each row checks the cited evidence or the underlying calculation. The manuscript side of every comparison remains open until the exact draft can be examined. The referee's reported 106 matching entries describe its earlier review.

Use these Git revisions for all paths in the matrix.

| Source revision | Commit |
| --- | --- |
| Frozen replay | [ff07102de88a4b63686a4f71f9fa851602ef7e0b](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/tree/ff07102de88a4b63686a4f71f9fa851602ef7e0b) |
| Referee source cutoff and v2 replay | [19a0e7f94a16d759761a5fe55470263b837c9adf](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/tree/19a0e7f94a16d759761a5fe55470263b837c9adf) |
| Coupled recovery | [a86b8f046908f852722d1d034f4e0711d5fc9466](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/tree/a86b8f046908f852722d1d034f4e0711d5fc9466) |
| First-row study | [7368b51fb9ad4026826294837e173ec7e22f8c41](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/tree/7368b51fb9ad4026826294837e173ec7e22f8c41) |

## Equations

| Eq. | Source-side check at referee cutoff | Manuscript comparison |
| --- | --- | --- |
| 1 | `reproduce/dee_stress_energy_3p1.py` sets the radial kinetic entry to 1, phase diagonals to 1 + 0.30ρ² and 1 + 0.20ρ², off-diagonal to 0.10ρ², and potential to (ρ² − 1)². | Open. Check grouping and square in the printed potential. |
| 2 | Same file implements half the sum of temporal and metric-contracted spatial kinetic terms plus the potential, and the negative momentum covector contraction. The conformally flat spatial inverse introduces ψ⁻⁴. | Open. Check index placement, metric and signs. |
| 3 | `reproduce/solve_toroidal_helical_coupled_constraints.py` implements the maximal conformally flat Hamiltonian with +ψ⁻⁷ times squared strain divided by 8 and +2πGEψ⁵. G is 10⁻³. | Open. Check conventions and identification of the benchmark source. |
| 4 | Same file has momentum source 8πGψ⁶S. Physical residual scales the divergence of strain by ψ⁻¹⁰, giving a source subtraction of 8πGψ⁻⁴S in background Cartesian components. | Open. Verify definition of S and each residual label. |
| 5 | `reproduce/diagnose_toroidal_boundary_convergence.py` and `research/audit_run43_arrays.py` give normalized RMS with weight rqψ⁶, q = R + r cos θ, radial bands r ≤ 0.10 and r ≥ 0.20 excluding both walls, and finite log₂(R_N/R_2N) ≥ 1.5 on three doublings. The benchmark fixes ψ = 1 and W = 0 at r = 0.05 and 0.25, with periodic angles. | Open. Distinguish the 12 frozen checks from prospective v2 interpretability. |
| 6 | `research/run43_frozen_boundary/REPORT.md` derives the first-interior composed radial row (3, −5, 1, 1)/(4h²) at offsets (−1, 0, 1, 2). Its normalized Taylor moments give 0, 0, 1, 1/4, 5/24 through degree four. Reflection changes the sign of the cubic term. The radial momentum multiplier is ψ⁻¹⁰(I + eᵣeᵣᵀ/3). | Open. Verify first-row indexing and wall exclusion. |
| 7 | Same report, weighted-band limit. For exact samples of fixed smooth periodic fields, ψ positive and bounded away from zero, nonsingular geometry, identical source, simultaneous angular refinement with steps O(h), and quadrature and denominator errors O(h), a first-row defect hB + O(h²) with O(h) mass yields squared RMS C²h³ + O(h⁴). Remaining rows contribute O(h⁴). C > 0 gives leading h^(3/2). | Open. Require stated hypotheses and distinguish computed-field residuals. |
| 8 | `research/run43_continuum_auditor_v1/auditor.py` constructs Cartesian contractions using direct derivative stencils and analytic metric derivatives. Mixed terms enter both ordered contractions. | Open. Check ordered indices and derivative definitions. |
| 9 | `research/run43_continuum_auditor_v2/diagnostics.py` computes pointwise operator gain times the global maximum nodal Euclidean norm of W with conformal scaling, adds the source norm, and uses weighted band RMS of this scale. The floor is 64 times binary64 epsilon times that RMS. | Open. Define the norm and state empirical scope. |
| 10 | V2 diagnostics and measurements. Composing a centered first difference gives curvature response relative to exact of [sin(kΔ)/(kΔ)]². It vanishes at the Nyquist cosine. The FFT recovers the sampled cosine curvature, while a Nyquist sine can vanish on every original node. | Open. Verify equation and accompanying limits. |

## Tables

| Table | Entries reported by referee | Source data and checks | Manuscript comparison |
| --- | ---: | --- | --- |
| 1 | 36 | Baseline and coupled `independent-audit.json` replayed from NPZ arrays. Hybrid `research/run43_hybrid/measurements.json` RMS orders independently recomputed against `research/run43_frozen_momentum/results/frozen-acceptance.json`. All 36 source-side comparisons verified. | Open. Verify all rounded cells and restrict caption to these three field families. |
| 2 | 6 | Referee-cutoff `research/run43_frozen_boundary/results/deficits.csv`, `rr_first` signed shares on hybrid and coupled. The three percentages are 6.97, 69.69, 108.25 for hybrid and 10.61, 81.77, 112.81 for coupled. The first-row study uses its own allocation and gives finest shares 98.3%, 102.8% and 97.4% for hybrid, coupled and reference. Cite the precise ledger and definition. | Open. Verify row labels and definitions. Signed shares above 100% arise through cancellation. |
| 3 | 32 | Referee-cutoff `research/run43_continuum_auditor_v1/results/comparison.csv` and coupled `research/run43_coupled_recovery/results/comparison.csv`. Identify field, region, evaluator, grid and reconstruction. The coupled CSV has identical Git blob `70fce7fe773881078d3bd4b52e50526caa175cfb` at both cited revisions. | Open. Verify cells and add the alternate reconstruction interpretation below the table. |
| 4 | 12 | V2 `measurements.json`, `archives.regions`, N=64 inner-band triples of historical v1 RMS, independent RMS and vector gap. Replayed archive assertions against pinned inputs. Independent reconstruction values and hybrid/4N summaries are inherited. | Open. Verify the four triples and which data support which checks. |
| 5 | 20 | V1 `results/manufactured.csv`, 64→128 regional orders. V2 `measurements.json`, `angular.cosine_modes`, N=64 full operator ratios across both axes, with zero at mode 32. | Open. Match the twelve selected regional rows and eight cosine ratios to the two panels and distinguish manufactured evidence from physical residuals. |

## Results that govern the manuscript claims

| Field family and evaluator | Evidence | Scope |
| --- | --- | --- |
| Exact November 3 frozen baseline | 12/12, minimum **1.7408475852068372**. Full replay JSON matches the archived JSON byte for byte. | Historical finite-grid acceptance. |
| Strong coupled fields, frozen gate | 9/12. Inner momentum orders 0.855454968, 0.971826621, 1.127653118. Full raw-array replay JSON matches byte for byte. | Three inner momentum failures. |
| Quadratic hybrid fields, frozen gate | 9/12. Inner momentum orders 0.474713503, 0.881282789, 1.123133441. All orders recomputed from archived RMS. | Three inner momentum failures. Archived RMS recalculation. |
| Same hybrid fields, substituted operator | 12/12 diagnostic combining six frozen Hamiltonian and six v1 momentum comparisons. V1 inner momentum orders 1.868352, 1.568025, 1.629676. | Compact radial second derivative and analytic metric rule. Frozen-gate result remains 9/12. |
| Coupled fields, v1 mechanical ratios | 6/12 from the six frozen Hamiltonian checks. | V2 withholds interpretability of the momentum ratios based on roundoff-limited residuals. |

The first-interior radial term explains a measured discrepancy under angular refinement. At N=64 with 8N angles and 8N radial reference intervals, the inner-band frozen, independently reconstructed and vector-gap RMS are approximately 1.996e−4, 5.187e−8 and 1.996e−4, respectively. The 4N-angle term decomposition attributes 96.24% of the finest reference gap by signed projection to the composed radial vector Laplacian. Coarse metric terms can cancel other contributions. The coupled frozen first-row RMS is about 7.03e−4 at N=64. Independently evaluated quartic traces and Fourier/angular derivatives diagnose finite-grid sensitivity, with radial stencil choices and join traces contributing unresolved uncertainty. [Angular sweep](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/5996f9f5de797aba7d67fc0025b7007ec40399fb/research/run43_asymptotic/REPORT.md) · [Term decomposition](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/a86b8f046908f852722d1d034f4e0711d5fc9466/research/run43_frozen_momentum/REPORT.md).

V2 classifies residuals at or below its empirical floor as roundoff-limited, those above the floor through ten times the floor as near-floor, and those above that margin as above-floor. The 27 manufactured arithmetic cases split 9/9/9. All 24 coupled interior outcomes are roundoff-limited and all 24 native-N reference outcomes are above-floor. The N=64 coupled inner floor is about 5.52e−13, against a v1 residual RMS of 1.87e−15. The angular screen reports underresolved-risk or screen-clear-on-samples, with the latter limited to sampled content. Hybrid and 4N v2 floors and screens lack raw-input verification. The v2 replay differed at six extended-precision-derived values on two N=64 coupled rows. The largest relative change was approximately 1.21e−7. Classification, hashes, scores and the manuscript-level rounded values reported by the audit were consistent. The replay did not regenerate the independent quartic reconstruction.

| N=64 coupled region | Quantity | Archived value | Replay value |
| --- | --- | ---: | ---: |
| Inner band | Extended residual RMS | 1.8719640710190286e−15 | 1.8719640720680373e−15 |
| Inner band | Float64 versus extended RMS | 3.652626830956602e−18 | 3.652626450921317e−18 |
| Inner band | Evaluation error / epsilon scale | 0.00042369967857899235 | 0.00042369963449542304 |
| Inner band excluding first row | Extended residual RMS | 1.913248285981953e−15 | 1.9132482870548276e−15 |
| Inner band excluding first row | Float64 versus extended RMS | 3.4656777067828628e−18 | 3.4656772881004694e−18 |
| Inner band excluding first row | Evaluation error / epsilon scale | 0.00040430446711160544 | 0.00040430441826829795 |

The final draft must identify the evidence revision for each claim, avoid conflating these scopes, check every in-text citation and reproduce the manuscript's stated commands. [Reproduction](REPRODUCTION.md) covers the audit commands actually exercised. [Data availability and verified source citations](DATA_AVAILABILITY.md) provides entries for manuscript reconciliation.
