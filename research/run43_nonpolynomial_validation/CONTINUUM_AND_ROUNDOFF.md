# Analytic fixture and floating-point enclosure

One supplemental fixture; the polynomial fixture and its failed evidence remain unchanged. No DEE equation, accepted solver, grid, archive, or original evaluator is modified. The fixture parameters below were fixed before its first evaluation; no amplitude, frequency, stencil, or convergence criterion was tuned afterward.

## Exact continuum definition

The toroidal chart is

\[
x=(1+r\cos\theta)\cos\phi,\quad y=(1+r\cos\theta)\sin\phi,\quad z=r\sin\theta,
\quad 1/20\le r\le1/4.
\]

Put \(q=1+r\cos\theta\), \(s=5(r-1/20)\), and

\[
F=e^{12(s-1)}[1+\sin\theta/5+3\cos\phi/20+\sin(\theta+\phi)/10],
\quad \psi=1+F/50,\quad W=(10,13,-9)F/1000.
\]

W components are Cartesian. F is smooth, periodic in both angles, non-polynomial in r, and has nonzero r, theta, phi, r-theta, r-phi, and theta-phi derivatives. All W components, and all exact momentum residual components, are nonzero. This fixture tests the differential residual evaluator; it is not a solution of the physical equations or a homogeneous-boundary manufactured solve. In particular it need not obey the accepted solution's radial Dirichlet conditions.

The orthonormal Cartesian vectors are

\[
e_r=(\cos\theta\cos\phi,\cos\theta\sin\phi,\sin\theta),\quad
 e_\theta=(-\sin\theta\cos\phi,-\sin\theta\sin\phi,\cos\theta),\quad
 e_\phi=(-\sin\phi,\cos\phi,0).
\]

For an arbitrary scalar f, its exact gradient is
\(g=e_rf_r+e_\theta f_\theta/r+e_\phi f_\phi/q\).
Its symmetric Hessian in this orthonormal frame has entries

\[
h_{rr}=f_{rr},\quad h_{\theta\theta}=f_{\theta\theta}/r^2+f_r/r,
\quad h_{\phi\phi}=f_{\phi\phi}/q^2+\cos\theta f_r/q-\sin\theta f_\theta/(rq),
\]
\[
h_{r\theta}=f_{r\theta}/r-f_\theta/r^2,\quad
h_{r\phi}=f_{r\phi}/q-\cos\theta f_\phi/q^2,\quad
h_{\theta\phi}=f_{\theta\phi}/(rq)+\sin\theta f_\phi/q^2.
\]

The Cartesian Hessian is \(H_f=\sum_{a,b}e_a h_{ab}e_b^T\); \(\Delta f=\operatorname{tr}h\).
All derivatives of F are elementary analytic derivatives: \(F_r=60F\), \(F_{rr}=3600F\), \(F_{r\theta}=60F_\theta\), and \(F_{r\phi}=60F_\phi\); angular derivatives differentiate the displayed three trigonometric terms exactly. `check_identities.py` verifies all nine Cartesian Hessian entries by symbolic chain-rule differentiation of the Cartesian gradient of an arbitrary smooth f, including the frame derivatives. It also independently verifies the Cartesian identity
\(\operatorname{div}LW=\Delta W+\nabla\operatorname{div}W/3\).
No numerical derivative enters those checks.

With \(c=(10,13,-9)/1000\), the exact quantities used in the continuum residual are

\[
(LW)_{ij}=c_i g_j+c_j g_i-\tfrac23\delta_{ij}(c\cdot g),\quad
(\operatorname{div}LW)_i=c_i\Delta F+\tfrac13\sum_j c_j(H_F)_{ij},\quad
\Delta\psi=\Delta F/50.
\]

The physical coefficients are independently reconstructed using
\(G=1/1000\), \(\rho^2=189/200\),
\(K=\left(\begin{smallmatrix}1+.30\rho^2&.10\rho^2\\.10\rho^2&1+.20\rho^2\end{smallmatrix}\right)\),
\(v=(.12,-.04)\), \(a=(.25,.20)\), \(d=e_\theta/r-3e_\phi/q\),
\(A=v^TKv/2+(\rho^2-1)^2\), \(B=(a^TKa)|d|^2/2\), and \(S_i=-(v^TKa)d_i\).
Thus the following formulas establish the exact, generally nonzero continuum residual, rather than assuming it is zero:

\[
H=-8\psi^{-5}\Delta\psi-\psi^{-12}\sum_{ij}(LW)_{ij}^2-16\pi G(A+B\psi^{-4}),
\quad M_i=\psi^{-10}(\operatorname{div}LW)_i-8\pi G\psi^{-4}S_i.
\]

`qualification.json` records the sampled exact residual maxima for each component and grid. Analytic formulas define the exact continuum functions; finite-precision values of those functions are enclosed, not treated as exact floating-point oracles.

## Floating-point error bound

The existing evaluator uses FFT angular differentiation, exact-rational-derived nine-point radial weights, the toroidal scalar Laplacian, and two successive Cartesian gradients for momentum. Qualification compares its residual with the analytic continuum residual above at every node, including the radial endpoints.

The bound does **not** assume a vendor-specific FFT error constant. A separate direct-matrix reference encloses the same ideal real trigonometric/radial collocation operator. For even n, the Fourier derivative matrices have, for j unequal k,

\[
D^{(1)}_{jk}=\tfrac12(-1)^{j-k}\cot[\pi(j-k)/n],\quad
D^{(2)}_{jk}=-\tfrac12(-1)^{j-k}\csc^2[\pi(j-k)/n],
\]

and diagonal entries 0 and \(-(n^2+2)/12\), respectively. D1 uses the same zero Nyquist derivative as the evaluator. The cotangent is evaluated as cos/sin to avoid an interval pole in tan at quarter turns. Radial reference entries come directly from exact rational uniform-grid Lagrange weights, without importing the evaluator's radial constructor.

All coordinates, trigonometric/exponential seeds, pi, rational physical constants and matrix coefficients start as `mpmath.iv` intervals at 80 decimal digits, then become outward binary64 balls. Ball arithmetic propagates errors through sampling, metric factors, all derivative applications, strain, divergence, powers of psi, matter terms and residual assembly. Product bounds include both first-order input errors and their product; reciprocal bounds use the strictly positive lower absolute endpoint. A matrix application adds coefficient/input uncertainty and \(\gamma_m |D||u|\), with \(\gamma_m=m\epsilon/(1-m\epsilon)\). Here epsilon is NumPy's machine epsilon (twice the usual unit roundoff), so this is conservative for ordinary IEEE binary64 operations and BLAS reductions. Positive bound arithmetic is inflated by 1024 epsilon and rounded outward with `nextafter`; contractions have at most 65 terms. The tiny positive floor prevents subnormal bound underflow. These are standard floating-point-model bounds, not a formal proof of hardware or library correctness; `mpmath.iv` supplies the transcendental interval seeds.

Let Rstar be the exact discrete residual of the analytically sampled fixture, Cstar its exact continuum residual, Rd the direct-reference center with radius bd, Rc the continuum center with radius bc, and Rfft the actually evaluated existing FFT residual. The triangle inequality gives

\[
|R_{fft}-R^*|\le |R_{fft}-R_d|+b_d,
\quad |C^*-R_c|\le b_c.
\]

Accordingly, if E is max|Rfft-Rc|, then the true discrete truncation norm T=max|Rstar-Cstar| lies in [max(0,E-U),E+U], where U includes max|Rfft-Rd|, max bd, max bc, and subtraction/reduction rounding. The observed FFT/direct discrepancy includes all differences in the existing evaluator's sampled geometry, physical constants and FFT operations; it is not assumed to vanish. Reference radii and FFT discrepancies are reported separately. The same bound covers the reference's input sampling error, even though the evaluator receives the fixture's binary64 samples.

For each required pair 16->32 and 32->64, each of H, Mx, My, Mz must have E>10U at both ends and

\[
p_{low}=\log_2[(E_{coarse}-U_{coarse})/(E_{fine}+U_{fine})]\ge1.5.
\]

The extra 10U condition prevents extracting a rate from roundoff-dominated data. No absolute manufactured tolerance or rate exemption is introduced. A resolved lower order is reported only on the required pairs; no asymptotic claim beyond the measured grids follows.

## Scope of the saved-baseline uncertainty

The fixture bound qualifies this evaluator on this analytic test. It does not certify unresolved angular modes or a uniform error bound for arbitrary saved fields. The baseline audit retains the previously specified empirical 5/7/9-point radial and n//3 angular-cutoff sensitivity envelope and the 10% budget. Those sensitivities are explicitly **not** guaranteed upper bounds. The conservative baseline order computed from them is a sensitivity-adjusted estimate, not a certified mathematical lower order. Failure of this necessary budget or of the unchanged p>=1.5 criterion prevents publication readiness. No amount of passing a manufactured evaluator test establishes absolute accuracy of the physical baseline.
