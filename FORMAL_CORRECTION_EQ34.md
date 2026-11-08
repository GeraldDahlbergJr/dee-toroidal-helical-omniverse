# Formal correction record: manuscript Equation (34)

Date: 2026-11-08 UTC. Scientific verdict: **NOT YET ESTABLISHED**.

## Evidence and scope

Preserved evidence reference: [21683875ed93cd413d3b4472ba360c6e21a70a9d](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/commit/21683875ed93cd413d3b4472ba360c6e21a70a9d).
All source citations below refer to that commit. This is a new correction record,
not an alteration of the historical manuscript or an assertion that a missing
historical audit has been recovered.

Authoritative recovered draft:
[evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex](evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex),
Eqs. (32)–(38), Appendix C. The exact TeX corresponding to the later July 21
PDF remains unrecovered.

## Incorrect universal reduction

Equation (34) writes, for every indexed component,

$$
\bar h_{\mu\nu}''+\frac1r\bar h_{\mu\nu}'
=-16\pi G T_{\mu\nu}(r).
$$

Stationarity, axisymmetry and translation invariance along z do not justify
applying this scalar radial Laplacian universally to cylindrical tensor
components. Cylindrical coordinate basis vectors vary with position; a flat
background still has nonzero cylindrical connection coefficients. In Cartesian
components a rotationally invariant tensor need not have radial scalar
components. Orthonormal cylindrical components also require connection terms.

The covariant background equation in linearized harmonic gauge is
$\nabla^{(0)\alpha}\nabla^{(0)}_\alpha\bar h_{\mu\nu}
=-16\pi G T_{\mu\nu}$, with
$\nabla^{(0)\mu}\bar h_{\mu\nu}=0$. Its radial reduction must be derived for the
declared tensor basis and complete ansatz. Equations (35) and (C1), obtained
from the scalar reduction, consequently cannot be used as universal tensor
solutions. This does not invalidate their use for an appropriately derived
scalar sector or the implementation of that documented scalar equation.

## Correct restricted phi-z twist operator

Conventions: signature $(-,+,+,+)$, units $c=1$, flat background
$ds_0^2=-dt^2+dr^2+r^2d\phi^2+dz^2$, trace reversal
$\bar h_{\mu\nu}=h_{\mu\nu}-g^{(0)}_{\mu\nu}h/2$.
For an isolated stationary, axisymmetric, z-independent symmetric perturbation
$h_{\phi z}=\bar h_{\phi z}=H(r)$, its line element is
$2H\,d\phi\,dz$. The trace vanishes. With $x=r\cos\phi$, $y=r\sin\phi$,

$$
\bar h_{xz}=-yH/r^2,\qquad \bar h_{yz}=xH/r^2.
$$

These components have zero Cartesian divergence in this restricted ansatz.
Direct Cartesian differentiation yields

$$
\Delta\bar h_{xz}=-\frac{y}{r^2}(H''-H'/r),\qquad
\Delta\bar h_{yz}=\frac{x}{r^2}(H''-H'/r).
$$

Transforming the source by the same covariant tensor rule gives

$$
\boxed{H''-H'/r=-16\pi G T_{\phi z}}.
$$

For the orthonormal amplitude $F=H/r$, the corresponding operator is
$F''+F'/r-F/r^2$, with source $T_{\hat\phi\hat z}=T_{\phi z}/r$.
These are different representations of the same restricted equation.
A smooth axis branch requires $H=r^2 f(r^2)$ with smooth f; $H(0)=0$ alone
does not establish smoothness. Other components, time or z dependence, and
other gauges require a fresh coupled derivation.

The preserved derivations are
[docs/TENSOR_TWIST_INDEPENDENT_DERIVATION.md](docs/TENSOR_TWIST_INDEPENDENT_DERIVATION.md)
and [docs/TENSOR_TWIST_HARMONIC_GAUGE_CHECK.md](docs/TENSOR_TWIST_HARMONIC_GAUGE_CHECK.md).
They establish this restricted operator, not a complete Einstein–DEE solution.

## Exact counterexample H(r)=r²

For $H=r^2$, $H'=2r$ and $H''=2$, hence

$$
\mathcal L_{\rm twist}[r^2]=2-2=0,\qquad
\mathcal L_{\rm scalar}[r^2]=2+2=4.
$$

Its Cartesian components are $h_{xz}=-y$, $h_{yz}=x$, which are smooth and
harmonic, with zero divergence. Thus the restricted linearized vacuum equation
is satisfied, whereas the universal scalar reduction would incorrectly infer
$T_{\phi z}=-1/(4\pi G)$ for nonzero G. One may multiply H by an arbitrarily
small constant on a finite domain to remain weak field. This local operator
counterexample asserts neither asymptotic flatness nor satisfaction of an
unrelated outer homogeneous boundary condition, and is not an exact nonlinear
vacuum metric claim. No numerical experiment is needed for this identity.

## Correction to interpretation of the historical PASS

The historical Equation (38) “Linearized Einstein--DEE consistency test: PASS”
is preserved as a reported result. Small perturbations after source scaling and
positive sampled directional NEC contractions do **not** certify all ten
linearized Einstein equations, harmonic gauge, source conservation, matter
equations, boundary compatibility, or the complete nonlinear Einstein equations.
The universal Eq. (34) reduction is incorrect. The available historical summary
therefore cannot establish complete weak-field Einstein–DEE consistency.

The historical Gaussian envelope was an acknowledged modeling choice, not
derived from the matter action. Its exact generator/configuration, width/form,
domain/grid, basis and complete source/metric profiles remain unrecovered.
Do not replace those missing inputs with illustrative parameters or later
manufactured checkpoints. The missing prior Eq. (34)-specific audit is still
unrecovered; this document makes its own explicit correction.

## Preservation and next step

Only this record and NEW_PHYSICAL_TEST_PROTOCOL.md are added on a separate
research branch descended directly from the evidence reference. No original
manuscript, evidence seal, archived output, solver, or frozen Run 43 evidence
is edited or reinterpreted as a physical certification. Run 43 manufactured
initial-data/certification records are distinct from the historical Gaussian
checkpoint. Their verdict remains **NOT YET ESTABLISHED**.

This phase is documentation and test design only. No PDE solve or physical
consistency test is authorized or performed. Execute the new protocol only
after the unresolved definitions are fixed in a versioned manifest and the
author explicitly approves execution.
