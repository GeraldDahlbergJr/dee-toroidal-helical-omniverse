# New physical test protocol: weak-field Einstein–DEE consistency

Date: 2026-11-08 UTC. Status: **DESIGN ONLY — EXECUTION BLOCKED**.
Scientific verdict: **NOT YET ESTABLISHED**.

Preserved reference:
[21683875ed93cd413d3b4472ba360c6e21a70a9d](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/commit/21683875ed93cd413d3b4472ba360c6e21a70a9d).
See [FORMAL_CORRECTION_EQ34.md](FORMAL_CORRECTION_EQ34.md).
All existing repository content, especially the original manuscript and frozen
Run 43 evidence, remains unchanged. No scientific code or PDE solve is executed
in this phase. Stop after committing/pushing these two documents; explicit
author approval is required before any physical test execution.

## 1. Question, claim and evidence separation

Test whether a specified nontrivial DEE matter configuration and complete metric
satisfy the weak-field expansion of the Einstein–matter system, including gauge,
all components, matter closure and conservation, on a declared domain with
compatible boundaries. A successful result is limited to that configuration,
domain and perturbative order. It establishes neither nonlinear survival nor
global non-collapse. Vacuum or zero-source solutions are controls, not DEE
physical evidence. Manufactured sources are implementation controls only.

This is a new experiment. It does not reproduce the historical Gaussian-source
checkpoint. The missing historical Gaussian form/width/parameters, generator,
grid/domain, component basis and full profiles must remain UNKNOWN. The
historical epsilon_s=10^-6 and illustrative polynomial parameters are not
default inputs for this experiment.

## 2. Established definitions versus proposed choices

“Established” below means specified by the preserved model or a mathematical
definition, not experimentally established physics.

| Item | Status |
|---|---|
| Three-scalar action, stress tensor and Euler–Lagrange equations | Preserved model definitions: recovered main.tex Eqs. (3)–(9) |
| Polynomial potential/kinetic functions | Preserved optional realization, Eqs. (14)–(15); coefficients for this test unresolved |
| Restricted twist operator | Derived under the restrictions in the correction record |
| Complete tensor equations, trace reversal, covariant differentiation | Mathematical definitions, with conventions fixed below |
| Stationary cylindrical annular experiment, finite boundaries, thresholds and refinement rules below | New proposed test-design assumptions, requiring approval |
| Particular source profiles, dimensional scales, couplings and boundary values | Unresolved; must be independently specified and frozen before execution |
| Historical Gaussian run | Incomplete provenance; no guessed reconstruction permitted |

## 3. Complete metric, basis and gauge

Use $c=1$, $G>0$, signature $(-,+,+,+)$ and
$R^\rho{}_{\sigma\mu\nu}=\partial_\mu\Gamma^\rho_{\nu\sigma}
-\partial_\nu\Gamma^\rho_{\mu\sigma}
+\Gamma^\rho_{\mu\lambda}\Gamma^\lambda_{\nu\sigma}
-\Gamma^\rho_{\nu\lambda}\Gamma^\lambda_{\mu\sigma}$.
Coordinates are $(t,r,\phi,z)$, with background
$g^{(0)}=\operatorname{diag}(-1,1,r^2,1)$.

The proposed stationary, axisymmetric, z-independent metric retains all ten
covariant symmetric functions:

$$
\begin{split}
ds^2={}&(-1+h_{tt})dt^2+(1+h_{rr})dr^2+
(r^2+h_{\phi\phi})d\phi^2+(1+h_{zz})dz^2\\
&+2h_{tr}dt\,dr+2h_{t\phi}dt\,d\phi+2h_{tz}dt\,dz\\
&+2h_{r\phi}dr\,d\phi+2h_{rz}dr\,dz+2h_{\phi z}d\phi\,dz .
\end{split}
$$

No component may be set to zero solely because a solver omits it. Any symmetry
or gauge reduction must be derived and the eliminated equations evaluated.
Define $h=g^{(0)\mu\nu}h_{\mu\nu}$ and
$\bar h_{\mu\nu}=h_{\mu\nu}-g^{(0)}_{\mu\nu}h/2$.
Impose background harmonic gauge
$C_\nu=\nabla^{(0)\mu}\bar h_{\mu\nu}=0$; this is a linearized condition,
not automatically the exact nonlinear harmonic condition.

Compute geometric equations in the coordinate basis using background
connections, or transform the entire tensor to Cartesian components.
For reporting use the background orthonormal coframe
$(dt,dr,r\,d\phi,dz)$, dual frame
$(\partial_t,\partial_r,r^{-1}\partial_\phi,\partial_z)$.
Use a separately constructed physical orthonormal tetrad of g for physical
energy conditions. Label coordinate, background-frame and physical-frame
quantities explicitly. Dimensionless weak-field size is a frame norm, not
the maximum of dimensionally different coordinate entries.

## 4. Matter closure and perturbative ordering

Fields are $\rho$, $\Theta$, $\Psi$ with

$$
\mathcal L_m=-\tfrac12(\nabla\rho)^2-V
-\tfrac12K_\Theta(\nabla\Theta)^2-\tfrac12K_\Psi(\nabla\Psi)^2
-\Lambda\nabla_\alpha\Theta\nabla^\alpha\Psi ,
$$

$$
\begin{split}
T_{\mu\nu}={}&\partial_\mu\rho\partial_\nu\rho
+K_\Theta\partial_\mu\Theta\partial_\nu\Theta
+K_\Psi\partial_\mu\Psi\partial_\nu\Psi\\
&+\Lambda(\partial_\mu\Theta\partial_\nu\Psi+
\partial_\mu\Psi\partial_\nu\Theta)+g_{\mu\nu}\mathcal L_m .
\end{split}
$$

Contractions use the stated metric at the stated expansion order. Candidate
stationary matter has $\rho=\rho(r)$ and
$\Theta=\omega_\Theta t+m_\Theta\phi+k_\Theta z+\theta(r)$,
$\Psi=\omega_\Psi t+m_\Psi\phi+k_\Psi z+\psi(r)$.
The radial profiles are unknown fields, not arbitrary prescribed sources.
If using the preserved polynomial realization, specify
$V=\lambda_V(\rho^2-v^2)^2/4$,
$K_\Theta=1+a\rho^2$, $K_\Psi=1+b\rho^2$, $\Lambda=c_{\rm mix}\rho^2$.
Here $c_{\rm mix}$ denotes the manuscript mixing coefficient, not light speed.

Evaluate all three matter residuals:

$$
E_\rho=\Box_g\rho-V'-\tfrac12K_\Theta'(\nabla\Theta)^2
-\tfrac12K_\Psi'(\nabla\Psi)^2-\Lambda'\nabla\Theta\cdot\nabla\Psi ,
$$

$$
E_\Theta=\nabla_\mu J_\Theta^\mu,\quad
E_\Psi=\nabla_\mu J_\Psi^\mu,\quad
J_\Theta^\mu=K_\Theta\nabla^\mu\Theta+\Lambda\nabla^\mu\Psi,\quad
J_\Psi^\mu=K_\Psi\nabla^\mu\Psi+\Lambda\nabla^\mu\Theta .
$$

A weak-gravity family must have a declared dimensionless parameter epsilon,
$h=\epsilon h_1+O(\epsilon^2)$, and a consistent expansion of the fields and
source. One permissible NEW bookkeeping choice is to multiply the entire
matter action by epsilon, giving $T^{\rm total}=\epsilon T$ and unchanged
matter Euler–Lagrange equations for epsilon>0. This changes the gravitational
coupling of matter and is a new modeling assumption requiring explicit approval;
it is not a recovered historical procedure. Alternatively specify a physically
defined field/coupling family and derive its ordering. Freeze exactly one
choice before execution. Leading matter must solve its flat-background
equations, and source corrections must be included at the claimed order.
Never multiply a stress tensor by a radial Gaussian or other envelope without
deriving the additional action, equations, forces and stress required to
preserve conservation. An externally supported source tests that augmented
model only, and all support stress must be included.

## 5. Domain, topology and boundary specification

Proposed domain: $t\in\mathbb R$, $r\in[r_{\rm in},R]$ with
$0<r_{\rm in}<R<\infty$, $\phi\in S^1$ of period $2\pi$, $z\in\mathbb R$
with translation invariance; report extensive quantities per unit z length.
This annular choice avoids silently assuming a regular axis in the presence
of winding. It tests a finite cylindrical region, not an isolated globally
asymptotically flat object.

Transport phases are proposed to be angular-valued of period $2\pi$, with
integer $m_\Theta,m_\Psi$. A different field topology needs explicit approval.
Time and z gradients need not be periodic.

At both radial boundaries prescribe values of $\rho,\theta,\psi$ and a
gauge-compatible gravitational boundary operator/data set. Proposed metric
data are traces of a declared interior-core/exterior continuation satisfying
the same perturbative equations and gauge, with matching induced metric and
extrinsic curvature. If a jump is intended, give its surface stress and
junction conditions. Do not impose ten arbitrary zero metric values.
Matter flux matching must include $J_\Theta^r,J_\Psi^r$ and radial stress;
a reflecting alternative must explicitly impose and verify zero normal
current/energy flux and derive its compatible gravitational data.

The actual boundary functions, continuation, flux choice, length scales and
values are **UNRESOLVED**. Freeze the exact independent boundary conditions,
count them against the constrained radial system, show compatibility, and
demonstrate that the chosen problem is not overdetermined before approval.
An annular result excludes claims about the omitted core. Any later axis
extension must enforce smooth Cartesian tensors, scalar regularity and allowed
winding; $H=r^2 f(r^2)$ is the twist condition. No asymptotic-flatness condition
may be borrowed from a localized 3D source for an infinite cylinder.

## 6. Componentwise residuals and conservation

Evaluate the full linearized Einstein tensor directly, independently of the
solver's reduced equations:

$$
\begin{split}
G^{(1)}_{\mu\nu}=\tfrac12[&
\nabla^{(0)}_\alpha\nabla^{(0)}_\mu h^\alpha{}_\nu+
\nabla^{(0)}_\alpha\nabla^{(0)}_\nu h^\alpha{}_\mu
-\Box_0 h_{\mu\nu}-\nabla^{(0)}_\mu\nabla^{(0)}_\nu h\\
&-g^{(0)}_{\mu\nu}
(\nabla^{(0)}_\alpha\nabla^{(0)}_\beta h^{\alpha\beta}-\Box_0 h)] .
\end{split}
$$

Report separate residuals
$R^{(1)}_{\mu\nu}=G^{(1)}_{\mu\nu}-8\pi G T^{[1]}_{\mu\nu}$
for tt, tr, tphi, tz, rr, rphi, rz, phiphi, phiz, zz.
Here $T^{[1]}$ is the complete source retained at first gravitational order,
including any declared amplitude. In harmonic gauge cross-check against
$-\Box_0\bar h_{\mu\nu}/2$; never assume Eq. (34).

Also evaluate $R^{\rm full}_{\mu\nu}=G_{\mu\nu}[g]-8\pi G T_{\mu\nu}[g,q]$.
For a first-order construction this full residual is expected to be
$O(\epsilon^2)$, not identically zero. Freeze the expansion and verify
that scaling; smallness alone cannot turn this into a nonlinear certificate.

Check all four $C_\nu$, all four
$\nabla^{(0)\mu}T^{[1]}_{\mu\nu}$, and all four
$\nabla_g^\mu T_{\mu\nu}$ at the appropriate retained order. Verify the
Noether identity
$\nabla^\mu T_{\mu\nu}=E_\rho\partial_\nu\rho+
E_\Theta\partial_\nu\Theta+E_\Psi\partial_\nu\Psi$
with the adopted residual sign conventions. Check local current divergences
and integrated radial flux balances, including boundary flux and any surface
sources. Check the linearized Bianchi identity independently. On t=constant
slices report Hamiltonian and all three momentum constraints using
$K_{ij}=-\mathcal L_n\gamma_{ij}/2$,
$E=T_{\mu\nu}n^\mu n^\nu$, $S_i=-\gamma_i{}^\mu T_{\mu\nu}n^\nu$,
and the stated metric/source order.

For every residual report a maximum norm, volume-weighted RMS norm, location
of the maximum, and boundary/interface norms. Do not pool tensor components
or discard a zero-source equation. Include points between solver nodes and
one-sided boundary evaluations.

## 7. Proposed numerical acceptance contract

All thresholds below are NEW protocol choices, not historical measurements.
Before execution freeze a manifest containing units, length L, gravitational
scale, field scales, epsilon values, source scales and boundary scales.
Dimensionless residuals use these fixed scales; never choose a denominator
after viewing the result. For each Einstein component use
$S_{ab}=\max(\epsilon/L^2,\|8\pi G T^{[1]}_{ab}\|_\infty)$ in the background
frame; zero-source components retain the first floor. Matter residual scales
follow the fixed action/field units, currents use their manifest current scale
divided by L, and gauge uses epsilon/L. Conservation uses a fixed stress
scale divided by L. Boundary errors use the prescribed field/flux scales.

Proposed requirements, for both maximum and RMS norms separately:

- First-order Einstein, matter, gauge, conservation, constraints and boundary
  residuals at their claimed retained order: normalized error <=10^-6.
- For the full Einstein residual of a first-order family:
  $L^2\|R^{\rm full}_{ab}\|/\epsilon^2<=10$ and observed epsilon order >=1.8
  over epsilon, epsilon/2, epsilon/4, after discretization errors are bounded.
  This coefficient bound is a preregistered test criterion, not a theorem.
  Apply analogous declared second-order budgets to full matter/conservation
  residuals when only first-order fields are constructed.
- Weak field: $\max_{ab}|h_{\hat a\hat b}|<=10^{-3}$ and spatial metric positive
  definite with one timelike metric eigen-direction. Check inverse and lapse.
- Use at least three radial resolutions N, 2N, 4N (N fixed in manifest).
  Require order >=1.8 for a declared second-order method, or a preregistered
  method-specific target. If errors reach a precision floor, demonstrate
  a higher-precision bound instead of reporting a spurious order.
- Independently bound discretization, quadrature, boundary/matching,
  differentiation, algebraic and floating-point errors. Their sum must be
  <=10% of each applicable acceptance threshold. A refinement difference or
  two quadrature rules agreeing is diagnostic, not a rigorous total-error bound.
- Independent fields must agree within their combined justified error bounds
  and <=10^-6 of their fixed field scale. Independent residual decisions must
  agree component by component.

The epsilon experiment must use the same approved physical family and boundary
scaling, not retune the model. A failure to meet error bounds is inconclusive
(**NOT YET ESTABLISHED**), not a numerical PASS. A demonstrated incompatible
equation or inadmissible field rejects that candidate configuration. Changing
a threshold, source or boundary after inspection requires a new protocol version
and a new experiment identifier.

## 8. Physical admissibility

Require finite smooth fields in the stated domain, rho>0 with a positive lower
bound exceeding its uncertainty, positive-definite kinetic matrix
$K_\Theta>0$, $K_\Psi>0$, $K_\Theta K_\Psi-\Lambda^2>0$ with uncertainty margins,
and $\lambda_V>0$ if using the polynomial stable-potential realization.
The homogeneous $m_{\rm eff}^2$ formula alone is not a coupled stability test.

Check the complete null cone in a physical tetrad:
$k^{\hat a}=(1,\mathbf n)$, $|\mathbf n|=1$.
Establish the minimum of $T_{\hat a\hat b}k^{\hat a}k^{\hat b}$ by an analytic
quadratic-form argument with certified coefficient bounds or a globally bounded
optimization; two azimuthal directions or finite angular sampling is insufficient.
Require NEC>=0 within a justified enclosure; an enclosure straddling zero is
unresolved unless the action-based proof closes the sign. Report observer
energy density and fluxes. Any additional WEC/DEC claim needs its own proof.

Require finite curvature invariants, no signature change, spacelike t slices,
and spacelike closed phi orbits ($g_{\phi\phi}>0$). These are local admissibility
checks, not a proof of global causality or absence of horizons. Report all
boundary support and fluxes; do not hide unphysical stresses outside the domain.
A finite annular weak-field test cannot establish global regularity or non-collapse.

## 9. Independent reproducibility and verification

Before running, publish a versioned manifest with every parameter, exact
profile/boundary expression or input-file hash, equations and perturbative order,
mesh, differentiation, quadrature, precision, solver stopping rule, tolerances,
software versions, commands and expected artifact paths. No null field is zero.
Record whether inputs are recovered definitions or new assumptions.

After approval, implementation controls may include flat/vacuum cases and the
exact H=r² twist identity, followed by manufactured convergence cases; they do
not count as the physical test. The physical source must come from the approved
matter equations and stress tensor rather than from the desired metric residual.

Require a separate verifier implementation that does not import the solver's
reduced operator or reuse its derivative/residual routine. Derive tensors in
Cartesian components or by an independent connection/curvature calculation;
verify coordinate transformations and gauge explicitly. Independently construct
the approved solution (different discretization or analytic solution where
available) and cross-check the same input manifest. If independent personnel
are unavailable, disclose that limitation; a second routine alone does not
constitute independent scientific review.

Archive full component arrays, fields, boundary data, error budgets, convergence
tables, environment lock information and SHA256 manifests under a NEW experiment
directory, never inside frozen evidence. A clean-environment verifier must
reproduce the inputs and decisions. Independent review must inspect matter
closure, boundary compatibility, error accounting and every omitted component.
No Run 43 score, manufactured forcing or historical PASS substitutes for this.

## 10. Unresolved definitions and execution gate

The following must be settled explicitly, versioned and approved before execution:

1. Whether the proposed stationary annular symmetry/topology is the intended
   physical question, and whether a nontrivial admissible solution exists.
2. Physical units and L; G; all potential/kinetic coefficients; phase frequencies,
   winding numbers and axial gradients; radial field data and positive-rho branch.
3. Weak-gravity parameter family and scaling of the complete action, fields,
   sources and boundaries; retained matter and metric perturbative orders.
4. Numerical r_in and R; core/exterior continuation or declared boundary support;
   exact constraint/gauge-compatible gravitational boundary data, matter values,
   current/energy flux conditions and any surface stress.
5. Numerical method, N, epsilon sequence, precision, dimensional residual scales,
   algebraic tolerance and independent error-enclosure/verification method.
6. Independent verifier/reviewer and final acceptance-contract approval.

Historical missing Gaussian data, the missing Eq. (34)-specific audit and missing
later-manuscript TeX remain provenance gaps. They are not filled by these new
choices and are not prerequisites to a separately approved new experiment.

No physical result is claimed by this protocol. This document is deliberately
not an executable configuration while the above definitions remain unresolved.
**Scientific verdict: NOT YET ESTABLISHED. Await explicit approval before
executing the physical test.**
