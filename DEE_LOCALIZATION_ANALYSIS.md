# DEE analytical localization analysis

Starting commit: 14b8a1384d15c5e9fc296c2daf0f8d750b4dbee1.
Date: 2026-11-07 America/Los_Angeles (2026-11-08 UTC).
New localization verdict: **OBSTRUCTION IDENTIFIED**, within the scope below.

Preserved, independent verdicts:
**MATTER EXISTENCE: ESTABLISHED** for the previously verified extended exact solution.
**COMPLETE EINSTEIN–DEE SPACETIME: NOT YET ESTABLISHED**.
Neither verdict is revised by this investigation.

## 1. Foundation, equations and scope

Read at the starting commit: the original recovered manuscript
evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex (blob a20c0ad7b3469487e50fa2be3a193d1aacc92ab6), particularly Eqs. (3)–(9), (14)–(15); DEE_MATTER_FIELD_DERIVATION.md; DEE_MATTER_EXISTENCE_ANALYSIS.md; DEE_MATTER_INDEPENDENT_CHECK.md; FORMAL_CORRECTION_EQ34.md; NEW_PHYSICAL_TEST_PROTOCOL.md. Their contents are preserved unchanged. The previously proved exact solutions and local annular existence result remain established within their original scopes.

Use their normalized units, signature (-,+,+,+), and prescribed flat metric
ds²=-dt²+dr²+r²dphi²+dz². No gravitational field is constructed.
Q=(Theta,Psi)^T, M=[[K_Theta,Lambda],[Lambda,K_Psi]].
The original matter density is

$$\mathcal L=-\tfrac12(\partial\rho)^2-V(\rho)
-\tfrac12\partial_\mu Q^T M(\rho)\partial^\mu Q.$$

With Q=w t+m phi+k z+q(r), u=q'(r), the unchanged equations are

$$\rho''+\rho'/r=V_{,\rho}+
\tfrac12[-w^TM_{,\rho}w+u^TM_{,\rho}u+
m^TM_{,\rho}m/r^2+k^TM_{,\rho}k],$$
$$rM u=C,\qquad J^t=-Mw,\quad J^r=Mu,\quad
J^\phi=Mm/r^2,\quad J^z=Mk.$$

The two-vector C is constant. Stress is the manuscript tensor; in particular

$$\varepsilon=\tfrac12\rho'^2+V+
\tfrac12(w^TMw+u^TMu+m^TMm/r^2+k^TMk),$$
$$p_r=\rho'^2+u^TMu+\mathcal L,\qquad
p_\phi=m^TMm/r^2+\mathcal L,$$
$$p_r'+(p_r-p_\phi)/r=0.$$

The principal no-go theorem assumes:
(a) this stationary ansatz on the complete transverse plane, including a regular axis;
(b) smooth, positive rho and globally defined phases, with no defects/surface sources;
(c) C² couplings and potential, M strictly positive along the solution;
(d) rho tends to a finite ordinary vacuum rho_v with M(rho_v) strictly positive and finite;
(e) U=V-V(rho_v)>=0 along the field range, U(rho_v)=U'(rho_v)=0;
(f) finite nonnegative transverse excitation energy defined below.
The original polynomial realization with lambda_V>0, v>0, rho_v=v and positive M(v) satisfies the potential/vacuum assumptions. A positive-semidefinite A in M=I+rho²A guarantees positivity globally but is not necessary for the theorem.

This is not a theorem for arbitrary unspecified V,M of Eq. (3), singular/punctured cores, degenerating target metrics, external boundaries, nonstationary fields, different phase ansatzes, or curved spacetime.

## 2. Three meanings of finite energy and a legitimate vacuum subtraction

An ordinary vacuum has constant rho=rho_v, constant Q, V'(rho_v)=0, and
T_vac=-V(rho_v) g. Subtracting this fixed vacuum tensor changes neither matter equations nor conservation; define U=V-V(rho_v) in the excitation energy. For the polynomial vacuum V(v)=0, no subtraction is needed.

$$e=\varepsilon-V(\rho_v)=\tfrac12\rho'^2+U+
\tfrac12(w^TMw+u^TMu+m^TMm/r^2+k^TMk).$$

Under the stated positivity assumptions e>=0.

| Notion | Definition | What it does not imply |
|---|---|---|
| Finite local density | e finite at each regular point, including axis if claimed | Integrability at infinity |
| Finite annular energy per unit z | $2pi\int_{r_{\rm in}}^R r e\,dr<infinity$ for fixed 0<r_in<R | A valid core or exterior |
| Finite total transverse energy per unit z | $E_\perp=2pi\int_0^\infty r e\,dr<infinity$ | Finite total 3D energy for an infinite z-invariant nonzero source |

The two old exact transport families have constant positive density on the plane and infinite E_perp. They have finite energy on every fixed bounded annulus. They are still valid matter solutions.

With positive M, a nonzero constant phase gradient has nonzero stress and is not the ordinary vacuum. Subtracting its full tensor defines a relative-to-transport-background problem, analyzed separately in section 7. It cannot be called a zero-stress vacuum or an isolated finite-energy source.

## 3. Positive energy bounds remove constant transport gradients

Continuity and M(rho_v)>0 imply constants 0<delta<=L<infinity and r_* such that
delta I<=M(r)<=L I for all r>=r_*.

If w or k is nonzero, then
$e>=\delta(|w|²+|k|²)/2$ on the tail, and

$$E_\perp(r_*,R)\ge
\frac{\pi\delta}{2}(|w|^2+|k|^2)(R^2-r_*^2).$$

Thus finite E_perp forces w=k=0. Even w=+/-k (common null gradients) cannot cancel this bound: the amplitude equation uses Lorentzian contractions, but energy contains the sum of the positive time and axial terms.

For nonzero winding m,

$$E_\perp(r_*,R)\ge\pi\delta|m|^2\log(R/r_*),$$

so finite E_perp forces m=0 under this ordinary-vacuum asymptotic condition, independently of the axis obstruction. This conclusion applies to each angular-valued phase (integer m), or already follows from single-valuedness for globally real-valued phases.

Using u=M^-1 C/r, the radial phase-gradient contribution is

$$E_u=\pi\int\frac{C^TM^{-1}C}{r}\,dr
\ge\frac{\pi|C|^2}{L}\log(R/r_*).$$

Consequently finite tail energy also requires C=0. A regular axis independently gives C=lim r M u=0. Positivity/invertibility then implies u=0 everywhere in a connected domain. Hence all phase derivatives vanish; additive phase constants are not nontrivial transport.

These bounds are necessary and do not by themselves prove a nontrivial amplitude profile exists.

## 4. Radial stress/Pohozaev obstruction to an amplitude lump

After section 3, the only equation is

$$\rho''+\rho'/r=U_{,\rho},\quad
e=\rho'^2/2+U,\quad p_r=\rho'^2/2-U,\quad
p_\phi=-\rho'^2/2-U.$$

Conservation gives the exact identity

$$\boxed{(r^2p_r)'=r(p_r+p_\phi)=-2rU.}$$

Axis regularity and bounded stress imply r²p_r->0 as r->0.
Finite E_perp implies the existence of radii R_j->infinity with
R_j² e(R_j)->0: otherwise r²e has a strictly positive eventual lower bound and
integral r e dr diverges logarithmically. More formally integrability implies
liminf r²e=0, which suffices. Since |p_r|<=e,
R_j²p_r(R_j)->0. Integrating the identity from 0 to R_j gives

$$R_j^2p_r(R_j)=-2\int_0^{R_j}rU\,dr.$$

The left side tends to zero and U>=0. Thus integral rU=0 and U=0 everywhere by continuity. Now r²p_r is identically zero; with U=0 this forces rho'=0 everywhere away from the axis, hence everywhere by continuity. With the stated asymptotic value, rho=rho_v.

This proof requires neither an assumed pointwise boundary limit at every large radius nor a guessed exponential tail. It avoids unsupported integration-by-parts boundary omissions. It also works if the nonnegative potential has several zeros: the stress identity still removes the radial gradient.

**Theorem:** under assumptions (a)–(f), the only smooth stationary solution in the stated ansatz is the ordinary vacuum with constant phases. Nontrivial localized stationary matter is obstructed in that scope.

For the original polynomial realization this is a rigorous no-go result for any coefficient set with positive M along the solution and positive finite M at the vacuum, not only the previously used A>0 subset.

## 5. Independent variational interpretation

For static transverse fields in two dimensions, the derivative part of the energy is invariant under Q_s(r)=Q(r/s), rho_s(r)=rho(r/s), whereas the potential energy scales as s². Thus

$$E(s)=E_{\rm grad}+s^2 E_U,\qquad E'(1)=2E_U.$$

A stationary finite-energy solution with legitimate vanishing boundary contributions must have E_U=0. This agrees with the stress proof, including its exceptional feature that in two dimensions the virial identity alone does not remove arbitrary sigma-model gradient energy. Here the first integrals and radial stress identity remove the remaining radial gradients within the stated ansatz. No universal ban on every two-dimensional sigma-model soliton is claimed.

A cutoff variation yields the same relation without assuming global dilation is an admissible compactly supported variation; the subsequence boundary argument in section 4 supplies the rigorous tail control. On an annulus the boundary term cannot be dropped:

$$R^2p_r(R)-r_{\rm in}^2p_r(r_{\rm in})
=-2\int_{r_{\rm in}}^R rU\,dr.$$

Finite-annular energy therefore does not trigger the complete-plane no-go theorem.

## 6. Winding: two independent obstructions, and the annular alternative

A smooth circle-valued phase on a disk has zero degree on every contractible circle. Thus nonzero winding cannot extend through a regular axis when Theta/Psi remain globally defined phase fields. In addition, a uniformly positive M near r=0 gives a logarithmically divergent integral of the m/r gradient energy and a 1/r² local-density divergence.

On a genuine annulus, integer winding is mathematically admissible and the preserved locally Lipschitz radial ODE gives local solutions with finite annular energy and positive M,rho on a sufficiently short interval. This is a valid alternative domain within the original action, not a regular disk or a proof of localization. If the exterior tends to an ordinary positive-M vacuum, the logarithmic tail obstruction persists. If a disk core with continuous phases is demanded, the topological obstruction returns. A puncture requires an explicitly specified core/inner boundary and flux/support accounting.

For the original polynomial realization M(0)=I: sending the amplitude to zero alone cannot make angular-gradient energy disappear. Allowing a kinetic eigenvalue to vanish at a limiting core would fail strict positive-definiteness there and would also require a justified target-space completion for phases to become undefined. This is not automatically licensed by calling rho an amplitude.

## 7. Subtracting an extended transport background does not give a polynomial radial lump

For completeness retain fixed finite w,k, but require a smooth positive rho tending to a positive homogeneous stationary value rho_b, on a regular full plane. Axis regularity gives m=C=u=0. Write
A=[[a,c_mix],[c_mix,b]], M=I+rho²A and

$$d=k^TAk-w^TAw,\qquad
\lambda_V(\rho_b^2-v^2)+d=0.$$

The original radial equation becomes exactly

$$\boxed{\rho''+\rho'/r=
\lambda_V\rho(\rho^2-\rho_b^2).}$$

The reduced effective potential is
$U_{\rm eff}=V+\frac12\rho^2d$, up to an irrelevant constant, and

$$U_{\rm eff}(\rho)-U_{\rm eff}(\rho_b)
=\frac{\lambda_V}{4}(\rho^2-\rho_b^2)^2\ge0.$$

This effective potential controls the stationary amplitude equation, not the physical energy density. In fact

$$\varepsilon-\varepsilon_b=\tfrac12\rho'^2+
[U_{\rm eff}(\rho)-U_{\rm eff}(\rho_b)]
+(w^TAw)(\rho^2-\rho_b^2),$$

so relative energy need not be positive.

There is nonetheless a stronger maximum-principle obstruction requiring no energy-sign argument. If rho exceeds rho_b somewhere, continuity and convergence to rho_b ensure a positive global maximum above rho_b is attained at finite r (possibly the axis, interpreted in Cartesian coordinates). At this maximum Delta rho<=0 while the displayed equation gives Delta rho>0. If 0<rho<rho_b anywhere, a positive global minimum below rho_b is attained, with Delta rho>=0 but RHS<0. Both are contradictions. Therefore rho=rho_b everywhere. The phases retain their unchanged extended gradients, so the only solution here is the reference homogeneous background, with zero relative excitation. This includes the previously established common-null-transport background at rho_b=v.

This theorem assumes positivity, regularity on the full plane, lambda_V>0, rho_b>0 and the fixed-gradient polynomial realization. It does not apply to an annulus with nonconstant boundary data, a solution crossing rho=0, an asymptotic rho_b=0 branch, or a different potential. The latter alternatives have not been established as admissible localized configurations here.

## 8. What is and is not resolved

Established: ordinary-vacuum, positive-M, regular-axis stationary localization is obstructed as stated; constant transport gradients have area/logarithmic energy tails; finite-radius vacuum patching without added sources fails the matching/uniqueness conditions detailed in DEE_BOUNDARY_FEASIBILITY.md; homogeneous-background subtraction does not yield a positive regular polynomial radial lump.

Not established: localized solutions for arbitrary unspecified V,M of the original general action, nonstandard vacua with limiting kinetic degeneracy, different global topology or phase dependence, sign-changing relative potentials, or nonstationary fields. None is selected as a remedy or claimed to exist. Changing the original polynomial functions would be a new model choice, even though Eq. (3) leaves general functions open; it cannot be presented as a result for the already studied polynomial realization.

No fabricated Gaussian data, numerical parameters, gravitational construction, PDE integration or numerical certification are used. The formal Eq. (34) correction and old matter-existence/spacetime verdicts stand unchanged.
