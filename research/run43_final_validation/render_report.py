"""Render recorded diagnostics without reevaluating any fields."""
import json
from pathlib import Path
p=Path(__file__).resolve().parent
r=json.loads((p/'results.json').read_text());q=json.loads((p/'qualification.json').read_text())
lines=['# Run 43 final bounded validation — NOT YET ESTABLISHED','',
'Reference: `b8022efa971711bba4421fc3468cffd2b2608fb4`. One agent; zero new DEE solves. No equations, solver, candidates, or saved inputs changed. No repair after qualification failure.','',
'## Decision','',
'The independent evaluator FAILED qualification. Its analytic polynomial fixture has nonzero radial, both angular and mixed derivatives in the toroidal chart, and nonzero exact momentum in all Cartesian components. The numerical Hamiltonian errors grow with refinement. A necessary constant-annihilation check exposes radial coefficient construction error exceeding the roundoff bound for application of the stored matrix. These results prevent independent continuum validation. This is an evaluator failure, not a demonstrated failure of the accepted baseline. Numerical work stopped without repair.','',
'| Evidence | Status | Quantitative basis |','|---|---|---|',
'| Accepted baseline identity | VALIDATED | Four archive SHA256 values equal saved result hashes; full input manifest retained |',
'| Equations and parameters recovered | VALIDATED | Original source definitions and specification hashed; G=.001, R=1, a=.05, b=.25, rho²=.945, Pi=(0,.12,-.04), phase gradients .25 and .20 times (dtheta-3dphi) |',
'| Independent evaluator qualification | FAILED | All four constant-annihilation tests fail; analytic fixture errors below |',
'| Radial Dirichlet traces | VALIDATED | psi=1 and W=0 exactly at both saved radial endpoints, all four grids |',
'| Angular periodic trace equality | NOT YET ESTABLISHED | Endpoint-exclusive samples contain no independent duplicate 2pi traces |',
'| Frozen numerical gate | VALIDATED (archived evidence only) | Preserved baseline independent-audit.json records 12/12; not rerun or promoted to continuum evidence |',
'| Continuum residual convergence p>=1.5 | NOT YET ESTABLISHED | Evaluator unqualified; provisional conservative orders do not supply a valid pass |',
'| Evaluator uncertainty <=10% | FAILED for this evaluator | Radial-width and angular-cutoff sensitivity exceeds budget; no certified error bound |',
'| Absolute physical residual accuracy | NOT YET ESTABLISHED | No predeclared physical absolute tolerance recovered; algebraic rtol=1e-10, atol=1e-13 are not physical accuracy requirements |',
'| Absolute solution accuracy | NOT YET ESTABLISHED | Residual convergence alone would not bound solution error |','',
'## Qualification','',
'| n | max H fixture error | max Mx error | max My error | max Mz error | constant Drr defect | dot-product roundoff bound |','|---|---|---|---|---|---|---|']
for a,b in zip(r['fixture']['rows'],q['necessary_consistency_check']):
 lines.append('| '+str(a['n'])+' | '+' | '.join(f'{v:.8e}' for v in a['max_abs_error']+[b['constant_second_derivative_defect'],b['stored_matrix_dot_product_roundoff_bound']])+' |')
lines+=['','The bound is gamma_(n+1) times the maximum absolute row sum, with gamma_m=m*eps/(1-m*eps). It bounds rounding in the stored-matrix dot product, not coefficient construction. Failure means the constructed matrix cannot be treated as annihilating constants up to application roundoff. No physical tolerance was invented.','',
'## Boundary audit','',
'All saved arrays have shape psi=(n+1,n,n), W=(3,n+1,n,n), finite entries and positive psi. Maximum radial trace errors are exactly zero. Periodicity of Fourier reconstruction follows by construction and cannot independently verify missing endpoint samples. Adjacent last/first angular samples are different spatial points; their differences are recorded in results.json and are not seam jumps.','',
'## Provisional residual and uncertainty tables','',
'**UNQUALIFIED DIAGNOSTICS — do not use these values as physical pass/fail evidence.** H=-8 psi^-5 Delta psi - psi^-12 |LW|² -16 pi G E. M_i=psi^-10 div(LW)_i -8 pi G psi^-4 S_i. Components are Cartesian. Geometry and kinetic projections are reconstructed inside evaluate.py without importing project code.','',
'Radial derivatives differentiate local degree-eight interpolants (nine nodes); angular derivatives differentiate real periodic trigonometric interpolants. The scalar Laplacian uses the analytic toroidal metric formula; momentum differentiates the reconstructed Cartesian strain.','',
'RMS weights are J psi^6, normalized by the sum of those weights over each sampled region. Normalized magnitudes divide each residual RMS by the RMS of the sum of absolute equation terms for that component. This is a dimensionless cancellation measure, not an absolute physical accuracy requirement.','',
'Uncertainty is the RMS of the pointwise maximum of |R9-R7|, |R7-R5| and |R9-R9_cut|, where cut removes angular modes above n//3. It is an empirical sensitivity envelope, not a guaranteed bound, and does not include unknown aliases or fully bound coefficient error.','',
'Each cell lists H, Mx, My, Mz in that order. First rows mean radial indices 1 and n-1; fixed bands exclude the radial boundary itself.','',
'| n | region | RMS | normalized RMS | uncertainty RMS | uncertainty/residual |','|---|---|---|---|---|---|']
for a in r['rows']:
 for name,v in a['regions'].items():
  fmt=lambda key:', '.join(f'{x:.6e}' for x in v[key])
  lines.append(f"| {a['n']} | {name} | {fmt('rms')} | {fmt('normalized_rms')} | {fmt('uncertainty_rms')} | {fmt('uncertainty_fraction')} |")
lines+=['','## Provisional orders','',
'Nominal p=log2(R_coarse/R_fine). Conservative sensitivity estimate p_low=log2((R_coarse-U_coarse)/(R_fine+U_fine)); null means the numerator interval contains zero. These are not certified lower orders because U is unqualified. A numerical convergence pass would require every relevant pair and component to meet p_low>=1.5 and both uncertainty fractions<=.10 with a qualified evaluator. No such pass is claimed.','',
'| region | pair | nominal H,Mx,My,Mz | conservative H,Mx,My,Mz |','|---|---|---|---|']
for o in r['orders']:
 fmt=lambda key:', '.join('null' if v is None else f'{v:.6f}' for v in o[key])
 lines.append(f"| {o['region']} | {o['pair']} | {fmt('nominal')} | {fmt('conservative')} |")
lines+=['','## Provenance and reproduction','',
'Accepted input directory: `research/results/coordinated-replay-20261103-1831/baseline/`. Selection is supported by `research/RUN43_COORDINATED_REPLAY_20261103.md`, which preserves the baseline and rejects fourth/local candidates. No candidate fields were loaded. `input_sha256.json` hashes baseline archives, result metadata, checkpoints, settings, executed solver snapshot, audits, source equations and original physical parameters. Frozen diagnostic Git blob: `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`.','',
'Environment: Python 3.12, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0. SciPy is not imported by the evaluator. From repository root:','',
'```sh','OPENBLAS_NUM_THREADS=1 python research/run43_final_validation/evaluate.py','OPENBLAS_NUM_THREADS=1 python research/run43_final_validation/test_manufactured.py','python research/run43_final_validation/render_report.py','```','',
'The qualification command intentionally exits 1 and writes qualification.json/log evidence. This reproduces the failed evaluator, not a new solve or repair. `results.json` contains component maxima as well as all tables; no nonfinite JSON values are used. Historical solves and audits are not rerun.','',
'Final outcome: independent continuum validation NOT YET ESTABLISHED. No repair or further investigation initiated. STOP.','']
(p/'REPORT.md').write_text('\n'.join(lines))
