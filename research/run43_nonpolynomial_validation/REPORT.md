# Run 43 non-polynomial evaluator qualification and accepted-baseline audit

Publication-readiness verdict: **NOT READY for PRD submission on this evidence. STOP.**

Continued from `5799e440b5f27bfe81c3ce5c44814a7d2d652157` on `validation/run43-final-bounded`. One agent; zero new DEE solves. All 931 preexisting tracked files are byte-identical to that commit, including the accepted solver, equations, archives and all historical validation evidence. All 31 original input hashes and all four accepted archive/result SHA256 comparisons pass. New work is confined to this directory.

## Decision and final evidence matrix

| Evidence | Status | Basis |
|---|---|---|
| Baseline and historical evidence preservation | PASS | 931/931 original Git blobs unchanged; 31/31 manifest hashes; 4/4 accepted archive hashes |
| Non-polynomial analytic fixture | PASS | Exponential radial factor; both angles and mixed derivatives; Cartesian W and exact momentum nonzero in all components |
| Independent exact continuum residuals | PASS | Analytic gradient/Hessian and div LW identity; 12/12 arbitrary-function symbolic checks |
| Floating-point separation | PASS | All required errors exceed 10U; smallest E/U is 418.1 at n=64 H |
| Unchanged radial primitive checks | PASS | 8/8 constant and 72/72 degree 0–8 derivative checks under the original arithmetic bounds |
| Evaluator p>=1.5 qualification | PASS | 8/8 required component/pair tests; lower orders 5.288–6.117 |
| Saved-baseline continuum convergence | REQUIREMENT NOT MET | Inner-band and inner-first-row nominal momentum orders below 1.5 on 32->64; sensitivity-adjusted component/pair gate fails |
| Baseline uncertainty <=10% | REQUIREMENT NOT MET | All n=64 momentum sensitivity fractions approximately 0.735–0.785; envelope is empirical, not certified |
| Radial Dirichlet traces | PASS | psi=1 and Cartesian W=0 exactly at both radial endpoints on all four saved grids |
| Independent angular seam equality | NOT ESTABLISHED | Endpoint-exclusive saved arrays contain no duplicate 2pi traces |
| Archived frozen numerical gate | PRESERVED ONLY | Prior 12/12 result not rerun and not promoted to continuum evidence |
| Absolute physical residual accuracy | NOT ESTABLISHED | No predeclared physical absolute tolerance; algebraic solver tolerances are not physical criteria |
| Absolute solution accuracy | NOT ESTABLISHED | Residual convergence would still require a stability/error estimate |
| Publication readiness | NOT READY | Evaluator qualification passes; physical baseline convergence/uncertainty and absolute accuracy evidence remain insufficient |

## Fixture qualification

Exact definitions, independent symbolic derivation, and the complete floating-point bound are in [CONTINUUM_AND_ROUNDOFF.md](CONTINUUM_AND_ROUNDOFF.md). The fixture is not a physical solution and its exact continuum residual is not assumed to vanish. All values below are maximum absolute errors at all grid nodes, including radial endpoints. Component order throughout: **H, Mx, My, Mz**. U bounds the discrepancy between the measured error and the exact discrete truncation norm, including sampling, coefficient, metric, matter, derivative, nonlinear and residual-assembly errors.

| n | Error E | Floating bound U | E/U | Exact continuum residual maxima |
|---|---|---|---|---|
| 8 | 7.712181e+01, 5.207673e+00, 5.992201e+00, 4.979135e+00 | 1.182234e-08, 6.036773e-12, 6.579754e-12, 6.689527e-12 | 6.523395e+09, 8.626585e+11, 9.107029e+11, 7.443179e+11 | 7.756492e+02, 5.121838e+01, 5.990402e+01, 4.938409e+01 |
| 16 | 4.942634e+00, 3.171744e-01, 3.704203e-01, 3.160322e-01 | 7.973260e-08, 9.667616e-11, 1.078263e-10, 9.688761e-11 | 6.199013e+07, 3.280792e+09, 3.435343e+09, 3.261843e+09 | 7.756492e+02, 5.217031e+01, 6.221824e+01, 5.215089e+01 |
| 32 | 1.265102e-01, 7.990644e-03, 9.226643e-03, 7.883793e-03 | 5.825562e-07, 1.519403e-09, 1.691558e-09, 1.564701e-09 | 2.171640e+05, 5.259068e+06, 5.454522e+06, 5.038531e+06 | 7.756492e+02, 5.279989e+01, 6.221824e+01, 5.215089e+01 |
| 64 | 1.853805e-03, 1.153810e-04, 1.329166e-04, 1.136647e-04 | 4.434219e-06, 1.938780e-08, 2.145528e-08, 1.988422e-08 | 4.180680e+02, 5.951217e+03, 6.195054e+03, 5.716326e+03 | 7.760612e+02, 5.285717e+01, 6.221824e+01, 5.215089e+01 |

| Pair | Nominal p | Conservative lower p | Pass every component |
|---|---|---|---|
| [16, 32] | 5.287954e+00, 5.310821e+00, 5.327213e+00, 5.325038e+00 | 5.287948e+00, 5.310820e+00, 5.327213e+00, 5.325038e+00 | True |
| [32, 64] | 6.092620e+00, 6.113834e+00, 6.117212e+00, 6.116034e+00 | 6.089167e+00, 6.113591e+00, 6.116979e+00, 6.115782e+00 | True |

The requirement remains p>=1.5 on 16->32 and 32->64 for every component. Both errors on each pair must additionally exceed 10U. No rate is inferred from roundoff-dominated data. The smallest separation on the required grids is H at n=64: E=1.853805e-3 versus U=4.434219e-6. These are measured-grid rates, not a proof of an asymptotic rate at arbitrary resolution.

## Independent saved-baseline continuum audit

Only `research/results/coordinated-replay-20261103-1831/baseline/fields_{8,16,32,64}.npz` were evaluated. No candidate, checkpoint field, new solve or equation modification was used. Checkpoint files were hashed only for preservation. The unchanged standalone evaluator reconstructs the toroidal metric, physical coefficients, Cartesian strain and continuum residual independently of the accepted solver.

RMS weights are J psi^6, normalized by their sum in each sampled region. Normalized residuals divide by the RMS of the sum of absolute equation terms. This dimensionless cancellation ratio is not a physical absolute tolerance. First rows are radial indices 1 and n-1; bands are fixed physical radial widths excluding the endpoints. Maxima and term-scale RMS are also retained in results.json.

Uncertainty U here is the RMS of the pointwise maximum of |R9-R7|, |R7-R5| and |R9-R9_cut|, with angular modes above n//3 removed in R9_cut. This is an empirical reconstruction-sensitivity envelope, **not** a certified upper error bound. The fixture arithmetic bound is not transferred to arbitrary baseline fields or unknown aliases. The required sensitivity budget is U/R<=0.10; failing it cannot be cured by calling the envelope conservative.

| n | Region | Residual RMS | Normalized RMS | Sensitivity RMS U | U/R |
|---|---|---|---|---|---|
| 8 | interior | 6.881112e-03, 2.532324e-03, 2.532324e-03, 2.085058e-03 | 1.376236e-02, 4.081880e-01, 4.081880e-01, 2.873708e-01 | 6.333420e-03, 1.598470e-03, 1.598470e-03, 1.279079e-03 | 9.204065e-01, 6.312265e-01, 6.312265e-01, 6.134502e-01 |
| 8 | inner_band | 1.601938e-02, 3.325541e-03, 3.325541e-03, 3.326488e-03 | 1.598664e-02, 4.587109e-01, 4.587109e-01, 3.501586e-01 | 1.517157e-02, 1.684682e-03, 1.684682e-03, 1.342485e-03 | 9.470760e-01, 5.065887e-01, 5.065887e-01, 4.035742e-01 |
| 8 | outer_band | 2.616509e-03, 3.151276e-03, 3.151276e-03, 2.222232e-03 | 1.199971e-02, 6.196179e-01, 6.196179e-01, 4.521858e-01 | 2.030754e-03, 2.141776e-03, 2.141776e-03, 1.700454e-03 | 7.761309e-01, 6.796534e-01, 6.796534e-01, 7.652007e-01 |
| 8 | inner_first_row | 2.120016e-02, 4.893480e-03, 4.893480e-03, 4.859915e-03 | 1.671847e-02, 6.640189e-01, 6.640189e-01, 5.106266e-01 | 2.310388e-02, 2.566333e-03, 2.566333e-03, 2.043031e-03 | 1.089798e+00, 5.244393e-01, 5.244393e-01, 4.203842e-01 |
| 8 | outer_first_row | 3.587771e-03, 4.257201e-03, 4.257201e-03, 2.877523e-03 | 1.805241e-02, 9.853881e-01, 9.853881e-01, 9.932047e-01 | 2.750185e-03, 2.929472e-03, 2.929472e-03, 2.326323e-03 | 7.665444e-01, 6.881215e-01, 6.881215e-01, 8.084466e-01 |
| 16 | interior | 3.280362e-03, 2.471043e-03, 2.471043e-03, 2.871958e-03 | 6.019745e-03, 3.842990e-01, 3.842990e-01, 3.718371e-01 | 2.387314e-03, 1.815542e-03, 1.815542e-03, 2.075333e-03 | 7.277594e-01, 7.347270e-01, 7.347270e-01, 7.226195e-01 |
| 16 | inner_band | 8.551153e-03, 5.544722e-03, 5.544722e-03, 6.689482e-03 | 7.146672e-03, 5.580451e-01, 5.580451e-01, 5.242190e-01 | 6.280714e-03, 3.973508e-03, 3.973508e-03, 4.698226e-03 | 7.344873e-01, 7.166289e-01, 7.166289e-01, 7.023303e-01 |
| 16 | outer_band | 9.000638e-05, 2.063098e-03, 2.063098e-03, 2.130998e-03 | 4.258687e-04, 4.268331e-01, 4.268331e-01, 4.183958e-01 | 1.460268e-05, 1.614118e-03, 1.614118e-03, 1.693502e-03 | 1.622405e-01, 7.823760e-01, 7.823760e-01, 7.946988e-01 |
| 16 | inner_first_row | 1.606559e-02, 1.243157e-02, 1.243157e-02, 1.498991e-02 | 8.862550e-03, 9.907056e-01, 9.907056e-01, 1.000000e+00 | 1.430033e-02, 8.691277e-03, 8.691277e-03, 1.026069e-02 | 8.901212e-01, 6.991292e-01, 6.991292e-01, 6.845064e-01 |
| 16 | outer_first_row | 5.778459e-05, 3.875794e-03, 3.875794e-03, 3.984676e-03 | 3.061715e-04, 9.998008e-01, 9.998008e-01, 1.000000e+00 | 2.794246e-05, 3.024328e-03, 3.024328e-03, 3.173207e-03 | 4.835624e-01, 7.803118e-01, 7.803118e-01, 7.963525e-01 |
| 32 | interior | 1.099605e-03, 1.431985e-03, 1.431985e-03, 1.851006e-03 | 1.927750e-03, 2.250754e-01, 2.250754e-01, 2.393455e-01 | 5.437385e-04, 1.061074e-03, 1.061074e-03, 1.364248e-03 | 4.944855e-01, 7.409817e-01, 7.409817e-01, 7.370303e-01 |
| 32 | inner_band | 2.980326e-03, 3.793936e-03, 3.793936e-03, 4.940989e-03 | 2.270909e-03, 3.702256e-01, 3.702256e-01, 3.619358e-01 | 1.483012e-03, 2.802524e-03, 2.802524e-03, 3.631365e-03 | 4.976007e-01, 7.386851e-01, 7.386851e-01, 7.349470e-01 |
| 32 | outer_band | 2.162252e-05, 5.493830e-04, 5.493830e-04, 6.132344e-04 | 1.042230e-04, 1.184309e-01, 1.184309e-01, 1.285510e-01 | 1.102523e-06, 4.284988e-04, 4.284988e-04, 4.815812e-04 | 5.098960e-02, 7.799636e-01, 7.799636e-01, 7.853135e-01 |
| 32 | inner_first_row | 7.114388e-03, 1.249190e-02, 1.249190e-02, 1.626870e-02 | 3.181676e-03, 9.974918e-01, 9.974918e-01, 1.000000e+00 | 4.941480e-03, 9.022059e-03, 9.022059e-03, 1.168622e-02 | 6.945755e-01, 7.222329e-01, 7.222329e-01, 7.183253e-01 |
| 32 | outer_first_row | 1.249651e-05, 1.451236e-03, 1.451236e-03, 1.615517e-03 | 6.835327e-05, 4.672105e-01, 4.672105e-01, 5.736114e-01 | 2.971632e-06, 1.127811e-03, 1.127811e-03, 1.267571e-03 | 2.377970e-01, 7.771386e-01, 7.771386e-01, 7.846223e-01 |
| 64 | interior | 3.090070e-04, 6.734663e-04, 6.734663e-04, 8.986489e-04 | 5.295416e-04, 1.055927e-01, 1.055927e-01, 1.158187e-01 | 8.247540e-05, 5.057843e-04, 5.057843e-04, 6.736426e-04 | 2.669046e-01, 7.510166e-01, 7.510166e-01, 7.496171e-01 |
| 64 | inner_band | 8.535156e-04, 1.860788e-03, 1.860788e-03, 2.486626e-03 | 6.210958e-04, 1.798654e-01, 1.798654e-01, 1.789175e-01 | 2.290575e-04, 1.396973e-03, 1.396973e-03, 1.863454e-03 | 2.683695e-01, 7.507428e-01, 7.507428e-01, 7.493904e-01 |
| 64 | outer_band | 5.272748e-06, 1.105690e-04, 1.105690e-04, 1.245985e-04 | 2.565165e-05, 2.344510e-02, 2.344510e-02, 2.585474e-02 | 8.736692e-08, 8.618362e-05, 8.618362e-05, 9.761676e-05 | 1.656952e-02, 7.794554e-01, 7.794554e-01, 7.834503e-01 |
| 64 | inner_first_row | 2.335634e-03, 8.839260e-03, 8.839260e-03, 1.181274e-02 | 9.323925e-04, 9.926482e-01, 9.926482e-01, 9.999981e-01 | 1.099969e-03, 6.502015e-03, 6.502015e-03, 8.672242e-03 | 4.709508e-01, 7.355836e-01, 7.355836e-01, 7.341432e-01 |
| 64 | outer_first_row | 2.930433e-06, 4.111812e-04, 4.111812e-04, 4.618640e-04 | 1.627711e-05, 1.002818e-01, 1.002818e-01, 1.182820e-01 | 3.322370e-07, 3.197660e-04, 3.197660e-04, 3.621977e-04 | 1.133747e-01, 7.776767e-01, 7.776767e-01, 7.842084e-01 |

## Baseline convergence and uncertainty limits

Nominal p=log2(Rcoarse/Rfine). Sensitivity-adjusted p=log2((Rcoarse-Ucoarse)/(Rfine+Ufine)); null means the coarse interval contains zero. These adjusted values are not certified lower orders because baseline U is empirical. A necessary numerical pass requires adjusted p>=1.5 and U/R<=.10 at both ends for every relevant component/pair/region. The audit does not meet it. In particular the 32->64 inner-band momentum nominal orders are 1.028, 1.028, 0.991; inner-first-row orders are 0.499, 0.499, 0.462.

| Region | Pair | Nominal p | Sensitivity-adjusted p | All components meet p and budget |
|---|---|---|---|---|
| interior | [8, 16] | 1.068787e+00, 3.534195e-02, 3.534195e-02, -4.619470e-01 | -3.371321e+00, -2.198560e+00, -2.198560e+00, -2.617825e+00 | False |
| interior | [16, 32] | 1.576870e+00, 7.871038e-01, 7.871038e-01, 6.337248e-01 | -8.798248e-01, -1.927248e+00, -1.927248e+00, -2.012960e+00 | False |
| interior | [32, 64] | 1.831273e+00, 1.088338e+00, 1.088338e+00, 1.042480e+00 | 5.057900e-01, -1.668728e+00, -1.668728e+00, -1.691590e+00 | False |
| inner_band | [8, 16] | 9.056270e-01, -7.375260e-01, -7.375260e-01, -1.007895e+00 | -4.128815e+00, -2.536242e+00, -2.536242e+00, -2.520991e+00 | False |
| inner_band | [16, 32] | 1.520649e+00, 5.474199e-01, 5.474199e-01, 4.370946e-01 | -9.751516e-01, -2.069812e+00, -2.069812e+00, -2.106013e+00 | False |
| inner_band | [32, 64] | 1.803981e+00, 1.027781e+00, 1.027781e+00, 9.906104e-01 | 4.679124e-01, -1.716324e+00, -1.716324e+00, -1.731889e+00 | False |
| outer_band | [8, 16] | 4.861472e+00, 6.111242e-01, 6.111242e-01, 6.048022e-02 | 2.485291e+00, -1.864972e+00, -1.864972e+00, -2.873762e+00 | False |
| outer_band | [16, 32] | 2.057493e+00, 1.908928e+00, 1.908928e+00, 1.797019e+00 | 1.730352e+00, -1.123010e+00, -1.123010e+00, -1.323344e+00 | False |
| outer_band | [32, 64] | 2.035908e+00, 2.312865e+00, 2.312865e+00, 2.299152e+00 | 1.936695e+00, -7.027567e-01, -7.027567e-01, -7.552160e-01 | False |
| inner_first_row | [8, 16] | 4.001005e-01, -1.345076e+00, -1.345076e+00, -1.624989e+00 | null, -3.182170e+00, -3.182170e+00, -3.164146e+00 | False |
| inner_first_row | [16, 32] | 1.175163e+00, -6.983595e-03, -6.983595e-03, -1.181071e-01 | -2.771776e+00, -2.524048e+00, -2.524048e+00, -2.563428e+00 | False |
| inner_first_row | [32, 64] | 1.606925e+00, 4.989951e-01, 4.989951e-01, 4.617552e-01 | -6.609359e-01, -2.144478e+00, -2.144478e+00, -2.160366e+00 | False |
| outer_first_row | [8, 16] | 5.956259e+00, 1.354136e-01, 1.354136e-01, -4.696353e-01 | 3.288413e+00, -2.377660e+00, -2.377660e+00, -3.698887e+00 | False |
| outer_first_row | [16, 32] | 2.209160e+00, 1.417210e+00, 1.417210e+00, 1.302466e+00 | 9.480508e-01, -1.598817e+00, -1.598817e+00, -1.829006e+00 | False |
| outer_first_row | [32, 64] | 2.092340e+00, 1.819436e+00, 1.819436e+00, 1.806456e+00 | 1.545647e+00, -1.176339e+00, -1.176339e+00, -1.243887e+00 | False |

## Boundary audit and absolute-tolerance limitations

All four arrays have the expected shapes, finite values, positive psi, and exactly zero radial trace errors. Saved angular grids exclude the duplicate 2pi endpoint. Adjacent last/first samples are distinct physical points; their differences in results.json are not seam jumps. Fourier reconstruction is periodic by construction and cannot independently establish equality of unsaved physical seam traces. No radial boundary equation has been substituted for a continuum residual.

| n | max |psi-1| at radial endpoints | max |W| at radial endpoints |
|---|---|---|
| 8 | 0.0 | 0.0 |
| 16 | 0.0 | 0.0 |
| 32 | 0.0 | 0.0 |
| 64 | 0.0 | 0.0 |

No predeclared absolute physical residual tolerance was recovered. Solver rtol=1e-10 and atol=1e-13 control algebraic solves and do not authorize an absolute continuum-accuracy claim. Even a successful residual-convergence test would not bound absolute solution error without stability/coercivity and source/interpolation-error estimates. The current audit additionally fails its convergence/uncertainty requirements. This does not prove nonexistence or failure of the continuum DEE model; it shows that publication readiness of the accepted baseline has not been established.

## Reproduction and process evidence

From repository root, use the following commands in order. The qualifier and symbolic checks must pass before the audit. No command runs a DEE solve. The historical failed qualification and its evaluator are preserved; do not run historical renderers to describe this attempt.

```sh
OPENBLAS_NUM_THREADS=1 python research/run43_nonpolynomial_validation/qualify.py
python research/run43_nonpolynomial_validation/check_identities.py
OPENBLAS_NUM_THREADS=1 python research/run43_nonpolynomial_validation/audit.py
python research/run43_nonpolynomial_validation/preserve.py
python research/run43_nonpolynomial_validation/render_report.py
OPENBLAS_NUM_THREADS=1 python research/run43_nonpolynomial_validation/verify_evidence.py
```

Environment: python 3.12.14, numpy 2.3.5, sympy 1.14.0, mpmath 1.3.0. Qualification and audit process logs, machine-readable component results, identity checks and preservation evidence are retained. SHA256SUMS covers this evidence directory. The existing evaluator source hash is recorded in preservation.json.

Implementation note: the first high-precision point-seed prototype already exhibited the reported resolved truncation sequence. Strengthening its seeds to interval arithmetic initially hit the tan(pi/2) interval pole in the Fourier-matrix bound, producing an incomplete log with NaNs and a JSON serialization exception. `interval_seed_debug.log` preserves that development failure; it is not a completed qualification result. Evaluating cot as cos/sin removes that representational pole. The fixture, existing evaluator, equations and p threshold were unchanged. Final qualification uses interval seeds throughout and contains only finite results. No physical/evaluator repair cycle was initiated after a failed rate test.

**Final verdict: NOT READY for publication on the available independent continuum evidence. All authorized evidence work is complete. STOP; no additional repair or solve is initiated.**
