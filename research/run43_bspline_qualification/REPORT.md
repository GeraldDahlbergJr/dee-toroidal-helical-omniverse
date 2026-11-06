# Run 43 — single bounded manufactured qualification

**Overall verdict: FAILED qualification (resource gate). Method accuracy at n=32 is undetermined.**

The single authorized execution completed n=8 and n=16, then stopped on a memory allocation failure at n=32 before its first solve. No retry, correction candidate, original-source solve, or manuscript work followed. The 4 GiB address-space ceiling denied another 54 MiB array during independently derived analytic forcing evaluation on the 192³ assembly quadrature points. This is a resource failure of this implementation, not evidence that degree-four Galerkin approximation cannot converge.

Base: `5a83752a0e8df98e9db5c25041aebcd3e87ed6f9`. Branch: `run43-independent-physical-validation`. Run computation elapsed 54.7154 s (including symbolic construction); session began 14:57:04 UTC. One CPU affinity and single-thread BLAS/OpenMP were enforced. Address space, a stricter ceiling than resident memory alone, was capped at 4 GiB. Exact failure-time RSS was not captured; last completed-grid peak RSS was 675620 KiB. No memory-peak value is inferred from the denied allocation.

## Criteria

FAIL below means qualification evidence is absent where explicitly indicated; it does not assert an observed accuracy defect.

| Predeclared criterion | Result | Evidence / limitation |
|---|---|---|
| Every H/M component, every region, both pairs: independent order ≥1.5 | FAIL — incomplete | 8→16 passes provisionally; 16→32 unavailable |
| Relative linearized algebraic residual ≤1e-10 | PASS for completed solves; FAIL for full three-grid qualification | Maximum 1.87628e-12; n=32 never solved |
| No unexpected discrete nullspace | FAIL — not fully established | Scalar/vector principal forms coercive; full coupled Jacobian not certified |
| Radial/angular traces ≤100 eps times fixed field scales | PASS at n=8,16; FAIL for full three-grid qualification | All measured defects zero; no n=32 field |
| Explicit quadrature and evaluator uncertainty, ≤10% per metric | FAIL — incomplete | Residual-integration comparisons small; assembly-quadrature and total evaluator/solver bounds absent |
| n=8,16,32 within 4 GiB | FAIL | n=32 allocation denied under cap |
| One thread; ≤30 minutes; no additional resolution/candidate | PASS | One execution; only authorized n attempted; stopped at first resource failure |

## Independent strong-residual table

RMS uses J psi^6; components are physical-scaled manufactured constraints as defined in PROTOCOL.md. p_low here includes only twice the observed residual-integration quadrature difference. It is **not** a total uncertainty-certified order. All n=32 values and 16→32 orders are unavailable.

| Region | Component | n=8 | n=16 | n=32 | p_low 8→16 |
|---|---|---:|---:|---|---:|
| inner | H | 8.529809627e-03 | 9.230528358e-04 | unavailable | 3.208028 |
| inner | Mx | 1.120614081e-03 | 9.227598785e-05 | unavailable | 3.602189 |
| inner | My | 9.740362265e-04 | 7.763592829e-05 | unavailable | 3.649176 |
| inner | Mz | 1.107028521e-02 | 7.944824337e-04 | unavailable | 3.800533 |
| bulk | H | 5.738399271e-03 | 5.282241434e-04 | unavailable | 3.441426 |
| bulk | Mx | 8.078425909e-04 | 5.987986960e-05 | unavailable | 3.753917 |
| bulk | My | 1.188383746e-03 | 5.180181437e-05 | unavailable | 4.519850 |
| bulk | Mz | 8.763856058e-03 | 5.421313791e-04 | unavailable | 4.014851 |
| outer | H | 5.852719666e-03 | 7.503065336e-04 | unavailable | 2.963555 |
| outer | Mx | 3.806839179e-04 | 2.049145512e-05 | unavailable | 4.215467 |
| outer | My | 7.408226594e-04 | 1.884070319e-05 | unavailable | 5.297199 |
| outer | Mz | 2.687195612e-03 | 1.259751258e-04 | unavailable | 4.414891 |
| first_inner | H | 6.544311936e-03 | 4.478391122e-04 | unavailable | 3.869189 |
| first_inner | Mx | 1.099033623e-03 | 7.740829018e-05 | unavailable | 3.827599 |
| first_inner | My | 7.802027515e-04 | 6.700553490e-05 | unavailable | 3.541486 |
| first_inner | Mz | 1.167417347e-02 | 6.332199339e-04 | unavailable | 4.204470 |
| first_outer | H | 1.743064237e-03 | 2.054724929e-04 | unavailable | 3.084608 |
| first_outer | Mx | 3.601539529e-04 | 1.134453441e-05 | unavailable | 4.988506 |
| first_outer | My | 7.325910309e-04 | 1.168487456e-05 | unavailable | 5.970288 |
| first_outer | Mz | 2.448231812e-03 | 5.386034588e-05 | unavailable | 5.506373 |

## Field errors and traces

| Field | RMS n=8 | RMS n=16 | Sample max n=8 | Sample max n=16 | Trace tolerance |
|---|---:|---:|---:|---:|---:|
| psi | 2.912058550e-07 | 4.221890108e-09 | 1.464590848e-06 | 2.140922217e-08 | 2.220446049e-14 |
| Wx | 4.094534731e-07 | 4.585745849e-09 | 2.067303345e-06 | 2.745756256e-08 | 5.551115123e-17 |
| Wy | 2.598128024e-06 | 2.856980574e-08 | 1.004953201e-05 | 1.316875233e-07 | 3.885780586e-17 |
| Wz | 4.386608334e-06 | 4.308036344e-08 | 1.715109984e-05 | 1.940178753e-07 | 4.996003611e-17 |

All radial and angular maximum absolute trace defects are 0 in both completed grids. Fixed characteristic field scales [1,.0025,.00175,.00225] were encoded before execution. The trace representation is compatible with the required tolerance; it was not relaxed. Interior psi is positive (by fixture design and small measured field errors); no certified whole-function lower bound for the numerical field was computed. Field maxima are sampled, not uniform error bounds.

## Quadrature, evaluator and conditioning assessment

Assembly uses 6-point Gauss per span and exact geometry factors. Independent residual integration uses 6 and 8 points per span. The integration uncertainty proxy is twice their RMS difference.
- n=8: largest integration proxy / scored RMS = 2.514627424e-05; largest raw basis derivative disagreement (orders 0–2 combined) = 1.000444172e-11; maximum linear relative residual = 1.603545925e-12.
- n=16: largest integration proxy / scored RMS = 1.510764304e-06; largest raw basis derivative disagreement (orders 0–2 combined) = 2.182787284e-11; maximum linear relative residual = 1.876279924e-12.

These differences are assessments, not rigorous enclosures. No assembly 6-versus-higher-order operator comparison was completed before the stop. The basis discrepancy combines different derivative units and is not a physical residual error bound. The evaluator uses independent local-power differentiation and Cartesian chain rules, but shares the analytic forcing generator, chart constants and tensor-contraction helper with assembly; it is not a separately reproduced experiment. No certified propagation bound from coefficient/evaluator/linear-solve errors into each physical metric was completed. Consequently the required total 10% uncertainty budget is NOT established despite small integration-rule differences.

The open/periodic tensor basis is linearly independent. For radial-zero periodic vector fields, integral |LW|²/2=integral |grad W|²+|div W|²/3; Poincare excludes a vector principal kernel, and likewise the scalar stiffness kernel. Positive quadrature weights alone do not replace a proof of sampled-gradient injectivity or local coupled Jacobian invertibility. CG convergence without breakdown for completed right-hand sides is useful evidence but does not certify absence of all discrete null modes. No full coupled spectrum/local inverse estimate was obtained.

## Reproduction and artifacts

See PROTOCOL.md for analytic solution, independent forcing derivation, geometry, spline construction, weak equations, residual definitions and fixed tolerances. `experiment.py` is the assembly/fixture implementation; `evaluator.py` is the independent strong evaluator; `run.py` is the bounded driver. Saved coefficient arrays represent psi-1 and Cartesian W, in that order. `results.json` and `execution.log` preserve the actual failed execution unchanged.

Reproduction commands (in a disposable checkout; not executed again in this authorization):

```sh
python -m pip install numpy==2.3.5 scipy==1.17.0 sympy==1.14.0
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  timeout 1500s python research/run43_bspline_qualification/run.py \
  > research/run43_bspline_qualification/execution.log 2>&1
```

Driver selects the lowest available CPU and enforces the address-space ceiling. The historical invocation used CPU 0 and timeout 1550s, with ample session time remaining; the observed allocation stop occurred after 55 s. Python 3.12.14. Results and timings can depend on platform allocation overhead. `results.json` retains the driver status INCONCLUSIVE because its generic exception handler records an incomplete run; this report explicitly adjudicates the resource criterion as FAIL and overall qualification as FAILED. No data field has been rewritten to hide the exception.

The implementation was left exactly as executed after failure; no memory correction or further experiment is authorized here. Only report/handoff/provenance work followed. The protocol prose was written while the single run was executing; all inherited numerical requirements and code trace scales were fixed before the run. No retrospective threshold changes were made.

A fixture pass would qualify a method for further testing only. This incomplete fixture does not establish independent physical validation of the original Einstein–DEE solution. **STOP; await further authorization.**
