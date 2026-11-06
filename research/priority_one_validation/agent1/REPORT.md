# Agent 1: momentum and radial boundary diagnosis

Checkpoint: `5ed79712e8170e1d794b358b16918094f79431f8`. No baseline edits or DEE solves. The new calculation evaluates only the accepted replay fields under `research/results/coordinated-replay-20261103-1831/baseline`; SHA-256 values are checked against their saved result manifests. Its purpose differs from historical fourth-vector budgets: quantify cancellation on the accepted baseline itself, and test whether its frozen residual is simply the prescribed stabilizer.

## Quantitative findings

The accepted operator in `research/run43_consistent_momentum.py`, `TensorGrid`, is the same composed Cartesian tensor divergence used by the frozen audit, minus the fixed angular fourth-difference term Q. Consequently its frozen momentum residual F should equal psi^-10 QW plus linear-solve error. Fresh read-only evaluation confirms this on the inner band:

| n | frozen RMS F | RMS(F−psi^-10 QW) | fourth-derivative audit RMS |
|---:|---:|---:|---:|
| 8 | 3.26003804e-4 | 3.6966e-13 | 1.75433877e-3 |
| 16 | 8.88511970e-5 | 7.0131e-13 | 1.78650363e-3 |
| 32 | 2.26751582e-5 | 4.7213e-13 | 1.00113689e-3 |
| 64 | 5.68952745e-6 | 1.2127e-13 | 4.05160645e-4 |

Thus the frozen order measures the vanishing angular stabilization defect on these solved fields. It does not independently bound the continuum momentum error. Frozen orders are 1.87542568, 1.97027859, 1.99473188; the alternative fourth-derivative audit gives −0.02621148, 0.83549959, 1.30507331, all below 1.5. Each Cartesian component independently falls below 1.5 for all three pairs in that supplementary audit (component orders in diagnostics.json). This failure is not an omitted vector component in the frozen vector RMS.

The radial first derivative uses centered interior stencils and a three-point one-sided boundary closure. At the first interior row,

    (Dr Dr − Drr)f = (−f0+3f1−3f2+f3)/(4h²)
                   = h f'''(r1)/4 + O(h²).

The composed operator therefore has a confirmed locally first-order boundary consistency defect. Fixed Dirichlet values do not force f''' to vanish. For vector Cartesian components its principal factor is I+er er^T/3. The basis er is independent of r, so this factor is analytic and is not a numerical component conversion.

Change only the endpoint D1 weights to the historically defined matched closure; call that supplementary residual B and E=F−B. On the same saved fields, exact algebra predicts E only on row 1 in the inner band:

    E1 = psi1^-10 (h/4)(I+er er^T/3) delta3 W.

The measured relative identity error is 2.25e−15, 3.04e−15, 6.86e−15, 2.03e−14. This verifies the closure contribution independently of a Taylor estimate of W'''.

| n | RMS E | RMS B | cos(E,B) | cancellation fraction |
|---:|---:|---:|---:|---:|
| 8 | 1.17645319e-3 | 1.00779065e-3 | −.96717687 | .95571147 |
| 16 | 1.34286672e-3 | 1.30996677e-3 | −.99806376 | .99775679 |
| 32 | 8.80092806e-4 | 8.75989961e-4 | −.99967746 | .99966654 |
| 64 | 4.19321230e-4 | 4.18900911e-4 | −.99990836 | .99990786 |

Cancellation fraction is −2<E,B>/(||E||²+||B||²), with the frozen physical-volume weights and unchanged band (.05,.10]. These are signed vector contributions; their norms must not be added as errors. At n64 each contribution is about 74 times F. Strong cancellation explains why the frozen norm is a poor independent accuracy measure. The calculation does not establish that the endpoint term alone explains the full continuum discrepancy: matched-closure residual orders are −.37833424,.58054398,1.06430532, also below threshold.

## Coordinates, tensor handling, and physical boundary data

The audited unknown is explicitly W_cartesian. `Audit.grad` differentiates each Cartesian scalar component and projects derivatives with er, et/r, ep/q; `Audit.tensor` symmetrizes the two Cartesian indices and subtracts 2/3 of the Cartesian trace. It correctly differentiates the varying basis through subsequent tensor-component derivatives. No missing connection term is demonstrated by this inspection. Finite-difference product-rule and commutator defects remain real truncation effects, as documented by the historical exact identity T−expanded=(sum GjGj−L)Wi+sum [Gj,Gi]Wj. The present new exact endpoint identity establishes a closure effect without diagnosing a coordinate error.

Prescribed radial data psi=1,W=0 at r=.05,.25 and periodic angles are finite-annulus benchmark data. They supply no extra third-normal derivative condition. A toroidal annulus excludes the r=0 coordinate axis, so no axis regularity condition may be invented to remove the defect. This report does not infer a matched exterior or asymptotic flatness.

## Corrective experiment decision

No new candidate is justified by this investigation. The specific obvious consistency correction—replace composed radial second derivatives near the boundary—is already tested, with failures preserved in `RUN43_LOCAL_BOUNDARY_20261103.md` (8/12), `RUN43_COORDINATED_REPLAY_20261103.md` (fourth-vector 9/12, local 8/12), and `RUN43_COMPATIBLE_CLOSURE_20261103.md` (direct Hessian and compatible closures). Matched endpoint, global radial, second Hessian, radial filter and coordinate-KO variants also already exist. Repeating one without a newly justified mechanism would violate the requested bounded experiment rule.

Confirmed: known local O(h) consistency defect; frozen diagnostic shares it with solver; frozen residual equals angular stabilizer to ~1e−13; extremely strong endpoint/remainder cancellation. Confirmed numerical failure: supplementary fourth-derivative inner-band momentum convergence is below 1.5 on all pairs and all components. Unresolved: absolute continuum error, degree of unresolved radial field structure versus differentiation conditioning, and a stable correction satisfying both frozen and independent continuum criteria. No coordinate, source, tensor-component or physical boundary prescription error is established here. Those possibilities cannot be called verified root causes.

Reproduce from repository root:

```sh
OPENBLAS_NUM_THREADS=1 python research/priority_one_validation/agent1/diagnose_baseline.py
```

Output: `diagnostics.json`, plus stdout (preserved as `diagnostics.log`). Historical archives are read only. Fourth-derivative and closure-subtraction diagnostics are supplementary and never substitute for the frozen gate.
