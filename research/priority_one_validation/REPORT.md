# Run 43 priority-one independent physical validation

Starting checkpoint: `5ed79712e8170e1d794b358b16918094f79431f8`.
Branch: `run43-independent-physical-validation`.

**Overall verdict: NOT YET ESTABLISHED.** The available independent physical validation requirements do not pass. No claim of independently validated Einstein–DEE initial data or PRD submission readiness follows. Three bounded investigations completed; **zero corrective experiments and zero new DEE solves**. All 952 preexisting tracked files, including accepted archives and historical evidence, remain byte-identical to the checkpoint. No manuscript work or submission was undertaken.

## 1. Quantitative root-cause report

The most defensible explanation of the frozen-versus-independent discrepancy is **shared radial boundary discretization error and cancellation in the frozen diagnostic**. This mechanism is verified by an exact discrete identity and fresh saved-field measurements. It explains why the frozen gate can pass without independently demonstrating continuum accuracy. It does **not** uniquely explain every continuum residual or establish whether the large radial reconstruction spread is actual unresolved field structure, interpolation amplification, or both.

The accepted solver uses the same composed Cartesian tensor divergence as the frozen momentum audit, with a prescribed angular stabilization Q. Consequently the frozen physical residual F is psi^-10 QW plus solve error. In the unchanged inner band (.05,.10], fresh measurements give:

| n | frozen vector RMS F | RMS(F−psi^-10 QW) | supplementary fourth-derivative RMS |
|---:|---:|---:|---:|
| 8 | 3.26003804e-4 | 3.6966e-13 | 1.75433877e-3 |
| 16 | 8.88511970e-5 | 7.0131e-13 | 1.78650363e-3 |
| 32 | 2.26751582e-5 | 4.7213e-13 | 1.00113689e-3 |
| 64 | 5.68952745e-6 | 1.2127e-13 | 4.05160645e-4 |

Frozen vector orders are 1.875426, 1.970279, 1.994732. Supplementary fourth-derivative orders are −.026211, .835500, 1.305073; every Cartesian component fails p>=1.5 on all three pairs. These supplementary operators are diagnostic comparisons, not replacement acceptance criteria or certified continuum oracles.

At first interior radial row, the original composed closure obeys

    (Dr Dr−Drr)f = (−f0+3f1−3f2+f3)/(4h²)
                 = h f'''(r1)/4 + O(h²).

The vector principal factor is I+er er^T/3. Replacing only endpoint D1 weights by the historically defined matched closure gives residual B and an exact difference E=F−B:

    E1 = psi1^-10 (h/4)(I+er er^T/3) delta3 W.

Fresh relative identity errors range from 2.25e-15 to 2.03e-14. At n=64, ||E||=4.19321230e-4 and ||B||=4.18900911e-4, their weighted cosine is −.99990836, and cancellation fraction −2<E,B>/(||E||²+||B||²)=.99990786. Each term is about 74 times ||F||. These signed contributions cancel; their norms cannot be summed into an error bound. Matched-closure residual orders also remain below threshold, so the endpoint term is not a verified complete explanation of continuum error.

The current checkpoint's nonpolynomial evaluator qualification passes its stated measured-grid tests: all H,Mx,My,Mz conservative fixture orders on 16→32 and 32→64 exceed 1.5 (range 5.288–6.117), with smallest error/floating-bound ratio 418.1. The earlier polynomial and coefficient-construction failures remain preserved historical evidence; they are not the current qualification outcome. Qualification on this smooth fixture does not certify arbitrary saved-field reconstruction error.

Fresh decomposition of the saved-field sensitivity identifies its measured source:

| n=64 region | radial 9−7 RMS / R9, Mx,My,Mz | full sensitivity / R9, Mx,My,Mz | largest angular-cut contribution / R9 |
|---|---|---|---:|
| inner band | .745610, .745610, .744184 | .750743, .750743, .749390 | 7.45e-14 |
| outer band | .778982, .778982, .783075 | .779455, .779455, .783450 | 2.65e-13 |
| inner first row | .735584, .735584, .734143 | .735584, .735584, .734143 | 3.14e-14 |

The envelope takes a pointwise maximum before RMS. Thus its norm need not equal the largest separate difference norm. The first four interior radial planes account for 98.8320%, 98.8320%, 99.1389% of the squared weighted 9−7 momentum difference over the whole interior. Saved W angular power above the n//3 cutoff is 3.45722e-31; that of psi−1 is 1.05045e-25. Resolved high angular modes do not explain the reported spread. Unsampled aliases remain unbounded. Boundary localization does not by itself prove a radial checkerboard or distinguish true field structure from differentiation conditioning. No particular reconstruction width has been shown closest to the unknown continuum field.

The independent action/source review finds no demonstrated sign, normalization, conformal-power or tensor-component error. With Pi=+n·dPhi,

    E=.0127638+.068614375 (r^-2+9q^-2) psi^-4,
    S_cart=−.030316 (e_theta/r−3 e_phi/q), G=.001.

The phase kinetic eigenvalues are (1.1305957881,1.3419042119). A new deterministic 512-sample direct 4D stress projection comparison, including nonzero lapse/shift and general positive spatial metrics, agrees with production at maximum absolute differences (E,S_i,S_ij,trace)=(7.37e-14,3.29e-14,1.24e-12,2.93e-13). These are observed algebraic discrepancies, not a new physical tolerance or error certification.

The physical continuum equations reviewed are

    H=−8 psi^-5 Delta psi−psi^-12 |LW|²−16piG E,
    M^i=psi^-10 partial_j(LW)^ij−8piG psi^-4 delta^ij S_j.

Momentum entries called Mx,My,Mz in the diagnostics are contravariant physical Cartesian components. The conformal source factor in the solve is psi^6 S_i. Cartesian W components differentiate as Cartesian scalar components on the chart; a second insertion of toroidal vector connection terms would double-count the coordinate rotation. No coordinate/component mistake was found. This is consistency with the specified conformally flat, maximal initial-data model, not validation of general Einstein–DEE evolution or of the proposed action experimentally.

## 2. Evidence-backed discrepancy matrix

| Classification | Finding | Evidence | Scope / consequence |
|---|---|---|---|
| Confirmed implementation error, historical evaluator | Floating monomial weight construction violated constant moments | `run43_final_validation/REPORT.md`, `REPAIR_REPORT.md` | Already repaired before starting checkpoint; cannot explain the current spread by itself |
| Confirmed discretization consistency defect | Composed radial closure is locally O(h) on first interior row | Agent 1 exact identity, relative agreement <=2.03e-14 | Solver and frozen audit share the defect; not a missing continuum equation term |
| Confirmed evaluator limitation | Frozen residual essentially measures stabilization on its own solved discrete system | Agent 1 F−psi^-10 QW RMS <=7.02e-13 across four grids | Frozen 12/12 is insufficient for independent continuum accuracy |
| Confirmed evaluator limitation | Baseline envelope is empirical, not a certified error enclosure | Agent 2 decomposition and current qualification scope | Neither 9-point nor 7-point reconstruction is certified closer to continuum |
| Confirmed numerical residual requirement failure | Current independent nominal inner momentum orders and sensitivity fail unchanged criteria | Current saved audit, recomputed envelope, matrix below | Failure of available numerical evidence; does not prove continuum model has no solution |
| Confirmed numerical sensitivity | Spread is radial and concentrated near inner surface | Agent 2 9−7/cutoff decomposition and radial localization | Resolved angular-cutoff artifact is excluded as its measured cause |
| Unresolved possibility | True unresolved radial field structure vs polynomial reconstruction amplification | No baseline-field derivative/error enclosure | Do not label either a verified complete physical root cause |
| Unresolved possibility | Unsampled angular aliases, one-sided continuum compatibility and derivative/seam accuracy | Unique angular samples, no certified continuum trace/derivative bounds | Low saved high-mode power alone cannot bound these |
| Not supported as explanation | Wrong source sign, normalization, conformal powers or tensor basis | Agent 3 derivation and 512 direct stress checks | No correction justified on this basis |
| Unresolved accuracy specification | Absolute physical accuracy and solution-error bounds | No predeclared physical absolute residual tolerance/stability estimate recovered | Algebraic solver tolerances cannot substitute |

## 3. Controlled corrective experiment decision

**NOT INITIATED: no new, specifically justified correction.** The exact boundary defect is known, but replacing composed derivatives, matched endpoints, local edge corrections, direct Hessians, fourth-vector, compatible closures and radial-filter variants already have preserved experiments. The coordinated replay records fourth-vector 9/12 and local 8/12 versus baseline 12/12; other failures are documented in the agent-1 report's historical references. Repeating those experiments without a new mechanism is not justified. Changing source data, adding a third-derivative boundary condition, excluding a row, relaxing p>=1.5 or discarding the sensitivity budget would change the problem or criteria.

No experimental solver branch or candidate was created because the conditional experiment prerequisite was not met. The new branch contains diagnosis and evidence only. No fallback candidate, proof search or repair loop was initiated.

Evidence missing before a correction can be selected: a quantified link from the measured saved-field radial behavior to a stable consistent correction distinct from rejected variants, together with an independently justified baseline radial reconstruction/error estimate. Finer or independently sampled derivative evidence could help, but no such solve is authorized by this report.

## 4. Independent validation evidence matrix

Component order in every vector cell is **H, Mx, My, Mz**. Current independent nine-point/Fourier residuals use normalized physical-volume weights J psi^6. These are point-sampled RMS diagnostics, not certified continuum integrals. At n=64:

| Region | RMS | 32→64 nominal p | sensitivity/RMS | Verdict under existing necessary criteria |
|---|---|---|---|---|
| inner band | 8.535156e-4,1.860788e-3,1.860788e-3,2.486626e-3 | 1.803981,1.027781,1.027781,.990610 | .268370,.750743,.750743,.749390 | FAIL |
| inner first row | 2.335634e-3,8.839260e-3,8.839260e-3,1.181274e-2 | 1.606925,.498995,.498995,.461755 | .470951,.735584,.735584,.734143 | FAIL |
| outer band | 5.272748e-6,1.105690e-4,1.105690e-4,1.245985e-4 | 2.035908,2.312865,2.312865,2.299152 | .016570,.779455,.779455,.783450 | FAIL: momentum sensitivity |

Exact values and all resolutions/regions/pairs are retained in the checkpoint's `run43_nonpolynomial_validation/results.json`, the new decomposition and `checkpoint_verification.json`. Rounded report cells do not define acceptance. No requirement was weakened: p>=1.5 and sensitivity <=10% remain necessary; sensitivity-adjusted rates still lack certification because the envelope is empirical.

| Claim | Verdict | Evidence / limit |
|---|---|---|
| Checkpoint/accepted archive preservation | PASS | 952/952 tracked blobs unchanged, 31/31 input manifest checks, 4/4 accepted archive/result SHA256 checks |
| Continuum equations and source algebra | PASS, scoped | Agent 3 independent derivation and new direct 4D projection comparison |
| Existing evaluator manufactured qualification | PASS, scoped archived result | Nonpolynomial fixture, all four components; new arithmetic check reproduces lower orders |
| Frozen numerical acceptance | PASS, archived only | Historical 12/12 preserved; no solve rerun or promotion to continuum evidence |
| Radial nodal Dirichlet traces | PASS | Exactly psi=1 and all W=0 at both endpoints, four grids |
| Angular periodic model prescription | PASS, scoped implementation | Periodic chart/operators; independent continuous seam/derivative accuracy NOT YET ESTABLISHED |
| H continuum convergence/uncertainty everywhere required | FAIL for available necessary gate | Inner-band H sensitivity .268> .10 at n=64; nominal p alone insufficient |
| Mx, My, Mz continuum convergence/uncertainty | FAIL for available necessary gate | Every inner momentum component below p=1.5 and far above .10 sensitivity |
| Boundary continuum constraint/normal-derivative accuracy | NOT YET ESTABLISHED | Nodal Dirichlet pass is insufficient; no extra endpoint PDE condition imposed |
| Certified arbitrary-field reconstruction uncertainty | NOT YET ESTABLISHED | Empirical spread is not a bound; fixture bound not transferred to baseline |
| Absolute residual and solution accuracy | NOT YET ESTABLISHED | Physical absolute tolerance and stability/error control missing |
| Complete independent physical validation | NOT YET ESTABLISHED | Relevant requirements do not all pass |

The finite annulus r∈[.05,.25], R=1, psi=1,W=0 at radial surfaces and periodic angles is an artificial benchmark. Its boundary data do not prescribe a physical exterior, asymptotic flatness, traction or an extra normal-derivative condition. Global phase winding admissibility needs a target-space convention; the prescribed local phase gradients suffice for the scoped constraint calculation. No full evolution, global nonlinear stability, generic nonaxisymmetric solution or observed physical confirmation is validated here.

## 5. Evidence and reproduction

Read the three independent investigation reports: [momentum/boundary](agent1/REPORT.md), [evaluator/uncertainty](agent2/REPORT.md), [equations/sources](agent3/REPORT.md). Their scripts, raw JSON and process logs are in those directories. New purpose: isolate cancellation and stabilization in the accepted baseline, decompose the existing envelope, and extend source checks to arbitrary lapse/shift. These are not repeated historical DEE solves.

Use Python 3.12 with NumPy 2.3.5, SciPy 1.17.0 and SymPy 1.14.0 (recorded environment); SciPy is needed by agent 1 imports, SymPy by agent 2's evaluator. From a checkout of this branch, run the bounded diagnostics as described in [REPRODUCE.md](REPRODUCE.md). Every script writes only within this new evidence directory. **Do not execute old renderers/audits that overwrite historical evidence.** No network or solver execution is needed to reproduce the diagnosis.

`checkpoint_verification.json` records every starting tracked file's Git blob and SHA256, accepted NPZ hashes/shapes/traces, input manifest verification, independently recomputed fixture and saved-audit rate arithmetic, and software versions. Its PASS is a preservation/arithmetic result, not a physical-validation pass. `SHA256SUMS` seals the new report/scripts/raw evidence, excluding itself.

Commit and remote identity are recorded in the final delivery message after push and `git ls-remote` comparison; the commit cannot contain its own SHA. This branch preserves the supplied checkpoint as its parent, with no baseline archive edits.

**STOP.** No additional solve, correction candidate, theoretical extension or submission is authorized by these deliverables. Further work awaits separate authorization.
