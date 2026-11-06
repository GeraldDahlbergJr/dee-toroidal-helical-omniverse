# Run 43 final bounded validation — NOT YET ESTABLISHED

Reference: `b8022efa971711bba4421fc3468cffd2b2608fb4`. One agent; zero new DEE solves. No equations, solver, candidates, or saved inputs changed. No repair after qualification failure.

## Decision

The independent evaluator FAILED qualification. Its analytic polynomial fixture has nonzero radial, both angular and mixed derivatives in the toroidal chart, and nonzero exact momentum in all Cartesian components. The numerical Hamiltonian errors grow with refinement. A necessary constant-annihilation check exposes radial coefficient construction error exceeding the roundoff bound for application of the stored matrix. These results prevent independent continuum validation. This is an evaluator failure, not a demonstrated failure of the accepted baseline. Numerical work stopped without repair.

| Evidence | Status | Quantitative basis |
|---|---|---|
| Accepted baseline identity | VALIDATED | Four archive SHA256 values equal saved result hashes; full input manifest retained |
| Equations and parameters recovered | VALIDATED | Original source definitions and specification hashed; G=.001, R=1, a=.05, b=.25, rho²=.945, Pi=(0,.12,-.04), phase gradients .25 and .20 times (dtheta-3dphi) |
| Independent evaluator qualification | FAILED | All four constant-annihilation tests fail; analytic fixture errors below |
| Radial Dirichlet traces | VALIDATED | psi=1 and W=0 exactly at both saved radial endpoints, all four grids |
| Angular periodic trace equality | NOT YET ESTABLISHED | Endpoint-exclusive samples contain no independent duplicate 2pi traces |
| Frozen numerical gate | VALIDATED (archived evidence only) | Preserved baseline independent-audit.json records 12/12; not rerun or promoted to continuum evidence |
| Continuum residual convergence p>=1.5 | NOT YET ESTABLISHED | Evaluator unqualified; provisional conservative orders do not supply a valid pass |
| Evaluator uncertainty <=10% | FAILED for this evaluator | Radial-width and angular-cutoff sensitivity exceeds budget; no certified error bound |
| Absolute physical residual accuracy | NOT YET ESTABLISHED | No predeclared physical absolute tolerance recovered; algebraic rtol=1e-10, atol=1e-13 are not physical accuracy requirements |
| Absolute solution accuracy | NOT YET ESTABLISHED | Residual convergence alone would not bound solution error |

## Qualification

| n | max H fixture error | max Mx error | max My error | max Mz error | constant Drr defect | dot-product roundoff bound |
|---|---|---|---|---|---|---|
| 8 | 2.48175044e-06 | 1.20097844e-10 | 1.20091585e-10 | 3.05405675e-11 | 3.02050466e-07 | 1.14336477e-09 |
| 16 | 6.16187618e-06 | 2.72929918e-10 | 2.72907792e-10 | 7.08411749e-11 | 8.07456672e-07 | 8.63875602e-09 |
| 32 | 2.72248513e-05 | 1.20200143e-09 | 1.19997593e-09 | 2.61176367e-10 | 3.48128378e-06 | 6.70773997e-08 |
| 64 | 1.94955139e-04 | 8.15440683e-09 | 8.14436065e-09 | 1.68537111e-09 | 2.46278942e-05 | 5.28488604e-07 |

The bound is gamma_(n+1) times the maximum absolute row sum, with gamma_m=m*eps/(1-m*eps). It bounds rounding in the stored-matrix dot product, not coefficient construction. Failure means the constructed matrix cannot be treated as annihilating constants up to application roundoff. No physical tolerance was invented.

## Boundary audit

All saved arrays have shape psi=(n+1,n,n), W=(3,n+1,n,n), finite entries and positive psi. Maximum radial trace errors are exactly zero. Periodicity of Fourier reconstruction follows by construction and cannot independently verify missing endpoint samples. Adjacent last/first angular samples are different spatial points; their differences are recorded in results.json and are not seam jumps.

## Provisional residual and uncertainty tables

**UNQUALIFIED DIAGNOSTICS — do not use these values as physical pass/fail evidence.** H=-8 psi^-5 Delta psi - psi^-12 |LW|² -16 pi G E. M_i=psi^-10 div(LW)_i -8 pi G psi^-4 S_i. Components are Cartesian. Geometry and kinetic projections are reconstructed inside evaluate.py without importing project code.

Radial derivatives differentiate local degree-eight interpolants (nine nodes); angular derivatives differentiate real periodic trigonometric interpolants. The scalar Laplacian uses the analytic toroidal metric formula; momentum differentiates the reconstructed Cartesian strain.

RMS weights are J psi^6, normalized by the sum of those weights over each sampled region. Normalized magnitudes divide each residual RMS by the RMS of the sum of absolute equation terms for that component. This is a dimensionless cancellation measure, not an absolute physical accuracy requirement.

Uncertainty is the RMS of the pointwise maximum of |R9-R7|, |R7-R5| and |R9-R9_cut|, where cut removes angular modes above n//3. It is an empirical sensitivity envelope, not a guaranteed bound, and does not include unknown aliases or fully bound coefficient error.

Each cell lists H, Mx, My, Mz in that order. First rows mean radial indices 1 and n-1; fixed bands exclude the radial boundary itself.

| n | region | RMS | normalized RMS | uncertainty RMS | uncertainty/residual |
|---|---|---|---|---|---|
| 8 | interior | 6.881107e-03, 2.532324e-03, 2.532324e-03, 2.085058e-03 | 1.376235e-02, 4.081880e-01, 4.081880e-01, 2.873708e-01 | 6.333420e-03, 1.598470e-03, 1.598470e-03, 1.279079e-03 | 9.204071e-01, 6.312265e-01, 6.312265e-01, 6.134502e-01 |
| 8 | inner_band | 1.601937e-02, 3.325541e-03, 3.325541e-03, 3.326488e-03 | 1.598663e-02, 4.587109e-01, 4.587109e-01, 3.501586e-01 | 1.517157e-02, 1.684682e-03, 1.684682e-03, 1.342485e-03 | 9.470766e-01, 5.065887e-01, 5.065887e-01, 4.035742e-01 |
| 8 | outer_band | 2.616504e-03, 3.151276e-03, 3.151276e-03, 2.222232e-03 | 1.199969e-02, 6.196179e-01, 6.196179e-01, 4.521858e-01 | 2.030749e-03, 2.141776e-03, 2.141776e-03, 1.700454e-03 | 7.761306e-01, 6.796534e-01, 6.796534e-01, 7.652007e-01 |
| 8 | inner_first_row | 2.120014e-02, 4.893480e-03, 4.893480e-03, 4.859915e-03 | 1.671845e-02, 6.640189e-01, 6.640189e-01, 5.106266e-01 | 2.310388e-02, 2.566333e-03, 2.566333e-03, 2.043031e-03 | 1.089798e+00, 5.244393e-01, 5.244393e-01, 4.203842e-01 |
| 8 | outer_first_row | 3.587764e-03, 4.257201e-03, 4.257201e-03, 2.877523e-03 | 1.805237e-02, 9.853881e-01, 9.853881e-01, 9.932047e-01 | 2.750179e-03, 2.929472e-03, 2.929472e-03, 2.326323e-03 | 7.665440e-01, 6.881215e-01, 6.881215e-01, 8.084466e-01 |
| 16 | interior | 3.280370e-03, 2.471043e-03, 2.471043e-03, 2.871958e-03 | 6.019761e-03, 3.842990e-01, 3.842990e-01, 3.718371e-01 | 2.387314e-03, 1.815542e-03, 1.815542e-03, 2.075333e-03 | 7.277576e-01, 7.347270e-01, 7.347270e-01, 7.226195e-01 |
| 16 | inner_band | 8.551176e-03, 5.544722e-03, 5.544722e-03, 6.689482e-03 | 7.146690e-03, 5.580451e-01, 5.580451e-01, 5.242190e-01 | 6.280714e-03, 3.973508e-03, 3.973508e-03, 4.698226e-03 | 7.344854e-01, 7.166289e-01, 7.166289e-01, 7.023303e-01 |
| 16 | outer_band | 9.000830e-05, 2.063098e-03, 2.063098e-03, 2.130998e-03 | 4.258778e-04, 4.268331e-01, 4.268331e-01, 4.183958e-01 | 1.460257e-05, 1.614118e-03, 1.614118e-03, 1.693502e-03 | 1.622358e-01, 7.823760e-01, 7.823760e-01, 7.946988e-01 |
| 16 | inner_first_row | 1.606566e-02, 1.243157e-02, 1.243157e-02, 1.498991e-02 | 8.862585e-03, 9.907056e-01, 9.907056e-01, 1.000000e+00 | 1.430033e-02, 8.691277e-03, 8.691277e-03, 1.026069e-02 | 8.901178e-01, 6.991292e-01, 6.991292e-01, 6.845064e-01 |
| 16 | outer_first_row | 5.779625e-05, 3.875794e-03, 3.875794e-03, 3.984676e-03 | 3.062333e-04, 9.998008e-01, 9.998008e-01, 1.000000e+00 | 2.794224e-05, 3.024328e-03, 3.024328e-03, 3.173207e-03 | 4.834612e-01, 7.803118e-01, 7.803118e-01, 7.963525e-01 |
| 32 | interior | 1.099594e-03, 1.431985e-03, 1.431985e-03, 1.851006e-03 | 1.927731e-03, 2.250754e-01, 2.250754e-01, 2.393455e-01 | 5.437386e-04, 1.061074e-03, 1.061074e-03, 1.364248e-03 | 4.944903e-01, 7.409817e-01, 7.409817e-01, 7.370303e-01 |
| 32 | inner_band | 2.980297e-03, 3.793936e-03, 3.793936e-03, 4.940989e-03 | 2.270886e-03, 3.702256e-01, 3.702256e-01, 3.619358e-01 | 1.483013e-03, 2.802524e-03, 2.802524e-03, 3.631365e-03 | 4.976056e-01, 7.386851e-01, 7.386851e-01, 7.349470e-01 |
| 32 | outer_band | 2.163674e-05, 5.493830e-04, 5.493830e-04, 6.132344e-04 | 1.042916e-04, 1.184309e-01, 1.184309e-01, 1.285510e-01 | 1.101941e-06, 4.284988e-04, 4.284988e-04, 4.815812e-04 | 5.092915e-02, 7.799636e-01, 7.799636e-01, 7.853135e-01 |
| 32 | inner_first_row | 7.114252e-03, 1.249190e-02, 1.249190e-02, 1.626870e-02 | 3.181615e-03, 9.974918e-01, 9.974918e-01, 1.000000e+00 | 4.941481e-03, 9.022059e-03, 9.022059e-03, 1.168622e-02 | 6.945889e-01, 7.222329e-01, 7.222329e-01, 7.183253e-01 |
| 32 | outer_first_row | 1.267342e-05, 1.451236e-03, 1.451236e-03, 1.615517e-03 | 6.932102e-05, 4.672105e-01, 4.672105e-01, 5.736114e-01 | 2.970055e-06, 1.127811e-03, 1.127811e-03, 1.267571e-03 | 2.343530e-01, 7.771386e-01, 7.771386e-01, 7.846223e-01 |
| 64 | interior | 3.090123e-04, 6.734663e-04, 6.734663e-04, 8.986489e-04 | 5.295507e-04, 1.055927e-01, 1.055927e-01, 1.158187e-01 | 8.247690e-05, 5.057843e-04, 5.057843e-04, 6.736426e-04 | 2.669049e-01, 7.510166e-01, 7.510166e-01, 7.496171e-01 |
| 64 | inner_band | 8.535313e-04, 1.860788e-03, 1.860788e-03, 2.486626e-03 | 6.211073e-04, 1.798654e-01, 1.798654e-01, 1.789175e-01 | 2.290595e-04, 1.396973e-03, 1.396973e-03, 1.863454e-03 | 2.683669e-01, 7.507428e-01, 7.507428e-01, 7.493904e-01 |
| 64 | outer_band | 5.221961e-06, 1.105690e-04, 1.105690e-04, 1.245985e-04 | 2.540456e-05, 2.344510e-02, 2.344510e-02, 2.585474e-02 | 5.860497e-07, 8.618362e-05, 8.618362e-05, 9.761676e-05 | 1.122279e-01, 7.794554e-01, 7.794554e-01, 7.834503e-01 |
| 64 | inner_first_row | 2.335757e-03, 8.839260e-03, 8.839260e-03, 1.181274e-02 | 9.324416e-04, 9.926482e-01, 9.926482e-01, 9.999981e-01 | 1.099978e-03, 6.502015e-03, 6.502015e-03, 8.672242e-03 | 4.709302e-01, 7.355836e-01, 7.355836e-01, 7.341432e-01 |
| 64 | outer_first_row | 9.760793e-07, 4.111812e-04, 4.111812e-04, 4.618640e-04 | 5.421575e-06, 1.002818e-01, 1.002818e-01, 1.182820e-01 | 2.230354e-06, 3.197660e-04, 3.197660e-04, 3.621977e-04 | 2.285013e+00, 7.776767e-01, 7.776767e-01, 7.842084e-01 |

## Provisional orders

Nominal p=log2(R_coarse/R_fine). Conservative sensitivity estimate p_low=log2((R_coarse-U_coarse)/(R_fine+U_fine)); null means the numerator interval contains zero. These are not certified lower orders because U is unqualified. A numerical convergence pass would require every relevant pair and component to meet p_low>=1.5 and both uncertainty fractions<=.10 with a qualified evaluator. No such pass is claimed.

| region | pair | nominal H,Mx,My,Mz | conservative H,Mx,My,Mz |
|---|---|---|---|
| interior | [8, 16] | 1.068782, 0.035342, 0.035342, -0.461947 | -3.371335, -2.198560, -2.198560, -2.617825 |
| interior | [16, 32] | 1.576888, 0.787104, 0.787104, 0.633725 | -0.879802, -1.927248, -1.927248, -2.012960 |
| interior | [32, 64] | 1.831235, 1.088338, 1.088338, 1.042480 | 0.505737, -1.668728, -1.668728, -1.691590 |
| inner_band | [8, 16] | 0.905622, -0.737526, -0.737526, -1.007895 | -4.128836, -2.536242, -2.536242, -2.520991 |
| inner_band | [16, 32] | 1.520667, 0.547420, 0.547420, 0.437095 | -0.975128, -2.069812, -2.069812, -2.106013 |
| inner_band | [32, 64] | 1.803940, 1.027781, 1.027781, 0.990610 | 0.467860, -1.716324, -1.716324, -1.731889 |
| outer_band | [8, 16] | 4.861439, 0.611124, 0.611124, 0.060480 | 2.485265, -1.864972, -1.864972, -2.873762 |
| outer_band | [16, 32] | 2.056575, 1.908928, 1.908928, 1.797019 | 1.729526, -1.123010, -1.123010, -1.323344 |
| outer_band | [32, 64] | 2.050819, 2.312865, 2.312865, 2.299152 | 1.821955, -0.702757, -0.702757, -0.755216 |
| inner_first_row | [8, 16] | 0.400094, -1.345076, -1.345076, -1.624989 | null, -3.182170, -3.182170, -3.164146 |
| inner_first_row | [16, 32] | 1.175196, -0.006984, -0.006984, -0.118107 | -2.771710, -2.524048, -2.524048, -2.563428 |
| inner_first_row | [32, 64] | 1.606822, 0.498995, 0.498995, 0.461755 | -0.661083, -2.144478, -2.144478, -2.160366 |
| outer_first_row | [8, 16] | 5.955965, 0.135414, 0.135414, -0.469635 | 3.288220, -2.377660, -2.377660, -3.698887 |
| outer_first_row | [16, 32] | 2.189170, 1.417210, 1.417210, 1.302466 | 0.932363, -1.598817, -1.598817, -1.829006 |
| outer_first_row | [32, 64] | 3.698664, 1.819436, 1.819436, 1.806456 | 1.597516, -1.176339, -1.176339, -1.243887 |

## Provenance and reproduction

Accepted input directory: `research/results/coordinated-replay-20261103-1831/baseline/`. Selection is supported by `research/RUN43_COORDINATED_REPLAY_20261103.md`, which preserves the baseline and rejects fourth/local candidates. No candidate fields were loaded. `input_sha256.json` hashes baseline archives, result metadata, checkpoints, settings, executed solver snapshot, audits, source equations and original physical parameters. Frozen diagnostic Git blob: `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`.

Environment: Python 3.12, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0. SciPy is not imported by the evaluator. From repository root:

```sh
OPENBLAS_NUM_THREADS=1 python research/run43_final_validation/evaluate.py
OPENBLAS_NUM_THREADS=1 python research/run43_final_validation/test_manufactured.py
python research/run43_final_validation/render_report.py
```

The qualification command intentionally exits 1 and writes qualification.json/log evidence. This reproduces the failed evaluator, not a new solve or repair. `results.json` contains component maxima as well as all tables; no nonfinite JSON values are used. Historical solves and audits are not rerun.

Final outcome: independent continuum validation NOT YET ESTABLISHED. No repair or further investigation initiated. STOP.
