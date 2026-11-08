# Nonlinear DEE Validation Protocol

## Status

This document defines a falsifiable validation program for extending the repository beyond its verified stationary, axisymmetric, z-independent, flat-background linearized harmonic-gauge sector.

It does **not** assume that the nonlinear DEE construction is a solution of the Einstein equations. The purpose is to determine whether a specified nonlinear metric/matter ansatz survives the equations and consistency checks.

## Governing equations

Use the full Einstein equations (c=1):

$$
G_{\mu\nu}[g] = 8\pi G\,T_{\mu\nu}.
$$

No terms quadratic or higher order in the metric perturbation or its derivatives may be dropped in the nonlinear test.

The geometric identity

$$
\nabla_\mu G^{\mu\nu}=0
$$

implies that an accepted matter model must satisfy

$$
\nabla_\mu T^{\mu\nu}=0.
$$

## Stage A — Specify the nonlinear ansatz

Before numerical claims are made, freeze an explicit coordinate system and metric ansatz. Record every nonzero component and every assumed symmetry. The ansatz must be non-degenerate and Lorentzian over the tested domain.

For a stationary axisymmetric benchmark, a Weyl–Lewis–Papapetrou representation is a useful independent reference because the exact Einstein equations reduce to coupled nonlinear equations for metric/twist potentials. It is a benchmark, not an assumption that DEE must use that gauge.

## Stage B — Independent geometry construction

From the metric itself compute, rather than prescribe:

1. inverse metric $g^{\mu\nu}$;
2. Christoffel symbols $\Gamma^\rho{}_{\mu\nu}$;
3. Riemann tensor $R^\rho{}_{\sigma\mu\nu}$;
4. Ricci tensor $R_{\mu\nu}$;
5. Ricci scalar $R$;
6. Einstein tensor $G_{\mu\nu}$.

Cross-check the implementation on exact metrics with known answers (at minimum Minkowski and one nontrivial exact GR solution appropriate to the implementation).

## Stage C — Nonlinear residual

For a proposed DEE metric and matter configuration define

$$
E_{\mu\nu}=G_{\mu\nu}-8\pi G T_{\mu\nu}.
$$

Report componentwise and normed residuals. A small residual at one grid resolution is insufficient; residual convergence under refinement is required.

## Stage D — Bianchi / matter consistency

Compute the covariant divergence of both sides independently:

$$
B^\nu=\nabla_\mu G^{\mu\nu},\qquad
C^\nu=\nabla_\mu T^{\mu\nu}.
$$

The discrete implementation must demonstrate convergence of these residuals toward zero at the expected numerical order. The geometric Bianchi identity is also a stringent code check.

## Stage E — Weak-field bridge

Introduce an amplitude parameter $\epsilon$:

$$
g_{\mu\nu}(\epsilon)=\eta_{\mu\nu}+\epsilon h_{\mu\nu}+O(\epsilon^2).
$$

Evaluate the full nonlinear residual at a sequence of decreasing $\epsilon$. The first-order coefficient must reproduce the independently verified linearized equations. In particular, for the restricted stationary axisymmetric $\phi z$ sector, the appropriate first-order limit must recover

$$
H''-\frac{1}{r}H'=-16\pi G T_{\phi z}
$$

under the same conventions and assumptions used in the linear verification.

A practical numerical test is to verify that the difference between the nonlinear tensor and its first-order approximation scales as $O(\epsilon^2)$.

## Stage F — Nonlinear contribution test

Demonstrate that the nonlinear calculation is genuinely nonlinear by isolating the difference

$$
N_{\mu\nu}(\epsilon)=G_{\mu\nu}[\eta+\epsilon h]-\epsilon G^{(1)}_{\mu\nu}[h].
$$

For a smooth weak-field family this remainder should begin at second order unless symmetry removes that coefficient. Fit its amplitude scaling rather than assuming it.

## Stage G — Boundary, axis, and signature checks

For every numerical configuration verify:

- regularity at the symmetry axis where applicable;
- finite curvature invariants in the claimed regular domain;
- determinant $\det g\neq0$;
- Lorentzian signature throughout the domain;
- stated boundary/exterior conditions;
- no hidden coordinate singularity is being interpreted as physical curvature.

## Stage H — Independent benchmark against exact stationary-axisymmetric GR

Where the chosen ansatz overlaps a standard stationary-axisymmetric vacuum sector, compare its equations/residuals with the Weyl–Lewis–Papapetrou/Ernst formulation. Exact stationary-axisymmetric vacuum GR contains nonlinear coupling between the gravitational and twist potentials; this provides a strong check that a purported nonlinear implementation has not accidentally retained only its linear part.

## Acceptance criteria

A nonlinear DEE configuration survives this validation stage only if all applicable criteria pass:

1. exact/benchmark metrics reproduce their known Einstein tensors to numerical tolerance;
2. $E_{\mu\nu}$ converges toward zero with grid refinement for the proposed solution;
3. Bianchi and matter-conservation residuals converge toward zero;
4. the $\epsilon\to0$ limit recovers the already verified linearized sector;
5. nonlinear-minus-linear scaling is consistent with the measured perturbative order;
6. signature, determinant, axis regularity, curvature, and boundary checks pass;
7. results can be regenerated from versioned inputs without manual tuning.

## Failure criteria

The ansatz does not survive in its tested form if a persistent nonzero Einstein residual remains under refinement, conservation cannot be satisfied by the stated matter model, the weak-field limit disagrees with the verified linear sector under matching conventions, or the metric loses its required signature/regularity in the claimed domain.

A failure of a particular ansatz is not by itself a proof that every possible DEE formulation is impossible; it falsifies the tested formulation and identifies what must be changed.

## Current checkpoint

The repository has already established a numerical and analytic checkpoint for the restricted linear twist sector. The next implementation task is therefore **not** to claim nonlinear validation, but to encode an explicit nonlinear DEE metric and matter closure and subject them to Stages B–H above.
