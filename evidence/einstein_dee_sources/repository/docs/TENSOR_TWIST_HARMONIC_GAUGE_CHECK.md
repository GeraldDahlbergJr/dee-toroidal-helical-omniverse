# Harmonic-gauge field-equation check for the tensor twist sector

## Scope

This note independently checks the previously derived Cartesian transformation of the stationary, axisymmetric cylindrical perturbation

$$
\bar h_{\phi z}=H(r)
$$

against the linearized Einstein equation on a flat background. It does not assume the reduced radial operator in advance.

## 1. Linearized equation and conventions

Use signature $(-,+,+,+)$, trace reversal

$$
\bar h_{\mu\nu}=h_{\mu\nu}-\frac12\eta_{\mu\nu}h,
$$

and harmonic gauge

$$
\partial^\mu\bar h_{\mu\nu}=0.
$$

The linearized Einstein equation is then

$$
\Box\bar h_{\mu\nu}=-16\pi G T_{\mu\nu}.
$$

For the off-diagonal $\phi z$ component, the background metric has $\eta_{\phi z}=0$, so trace reversal does not change that component: $\bar h_{\phi z}=h_{\phi z}$.

We restrict this check to stationary, axisymmetric, $z$-independent fields. The d'Alembertian therefore reduces to the spatial transverse Laplacian when evaluated in Cartesian components.

## 2. Cartesian representation

With $x=r\cos\phi$, $y=r\sin\phi$ and

$$
d\phi=\frac{-y\,dx+x\,dy}{r^2},
$$

the cylindrical covariant component gives

$$
\bar h_{xz}=-\frac{y}{r^2}H(r),\qquad
\bar h_{yz}= \frac{x}{r^2}H(r).
$$

Define $F(r)=H(r)/r$. Then

$$
\bar h_{xz}=-\sin\phi\,F(r),\qquad
\bar h_{yz}=\cos\phi\,F(r).
$$

These are azimuthal $m=1$ Cartesian-component patterns, not radial scalars.

## 3. Harmonic-gauge check

For $\nu=z$,

$$
\partial^\mu\bar h_{\mu z}
=\partial_x\bar h_{xz}+\partial_y\bar h_{yz}=0
$$

for an axisymmetric azimuthal field. For $\nu=x$ or $y$, the only added component in this restricted ansatz carries a $z$ index and its $z$ derivative vanishes. Thus this isolated stationary axisymmetric twist perturbation is compatible with the harmonic condition within the stated restricted ansatz.

This is a compatibility check, not a claim that an arbitrary larger DEE perturbation automatically satisfies harmonic gauge; additional components must be checked together when present.

## 4. Apply the Cartesian Laplacian

For either $m=1$ angular factor, the scalar Cartesian Laplacian acting on the Cartesian component amplitude is

$$
\nabla^2[F(r)e^{\pm i\phi}]
=\left(F''+\frac1rF'-\frac{F}{r^2}\right)e^{\pm i\phi}.
$$

Substituting $F=H/r$ gives

$$
F''+\frac1rF'-\frac{F}{r^2}
=\frac1r\left(H''-\frac1rH'\right).
$$

The source transforms in exactly the same way:

$$
T_{xz}=-\frac{y}{r^2}T_{\phi z},\qquad
T_{yz}=\frac{x}{r^2}T_{\phi z}.
$$

Therefore the common angular factors cancel from the Cartesian linearized Einstein equations and yield

$$
\boxed{H''-\frac1rH'=-16\pi G\,T_{\phi z}.}
$$

Thus the tensor-aware operator is recovered from the harmonic-gauge linearized field equation independently of the earlier reduced-coordinate guess.

## 5. Axis regularity

Smooth Cartesian metric components require

$$
H(r)=O(r^2)
$$

for the regular axisymmetric branch used here. The homogeneous function $H=r^2$ obeys the derived operator exactly. A constant $H$ would make the Cartesian components singular and is excluded by regularity.

## 6. Analytic identities

For

$$
\mathcal L_T[H]=H''-\frac1rH',
$$

one obtains

$$
\mathcal L_T[r^2]=0,\qquad
\mathcal L_T[r^4]=8r^2,\qquad
\mathcal L_T[r^6]=24r^4.
$$

These provide manufactured-solution tests independent of the scalar radial solver.

## 7. Numerical verification protocol

The accompanying `tensor_twist.py` discretizes the independently derived BVP with centered second-order differences. `tests/test_tensor_twist.py` checks the analytic polynomial identities and a manufactured solution

$$
H_{\rm exact}(r)=a r^2+b r^4,
$$

whose source is fixed analytically by

$$
T_{\phi z}(r)=-\frac{8b r^2}{16\pi G}.
$$

Grid refinement at 41, 81, 161 and 321 points is required to show observed order greater than 1.9 at every refinement.

## 8. Interpretation boundary

This result verifies the operator for this **stationary, axisymmetric, z-independent, flat-background linearized harmonic-gauge sector**. It does not by itself establish nonlinear Einstein dynamics, matter closure, global stability, exterior matching, or experimental validity. If additional perturbation components or different gauge conditions are introduced, the coupled system must be re-derived rather than assuming this isolated equation remains unchanged.
