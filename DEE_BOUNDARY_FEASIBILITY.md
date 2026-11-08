# DEE radial boundary and matching feasibility

Starting commit: 14b8a1384d15c5e9fc296c2daf0f8d750b4dbee1.
Date: 2026-11-07 America/Los_Angeles (2026-11-08 UTC).
Scope: original three-scalar matter action on prescribed flat cylindrical spacetime.

The matter-existence verdict remains ESTABLISHED for the old extended exact family.
The complete Einstein–DEE spacetime verdict remains NOT YET ESTABLISHED.
This document gives matter matching conditions, not gravitational junction construction.

## 1. Conditions derived from the original action

Let an interface be r=R>0 and orient a common normal n=+partial_r from the interior (-) to the exterior (+). Write [f]=f_+-f_-; use this same normal on both sides when comparing fluxes. The action's variation has boundary integrand
$-n_\mu(\partial^\mu\rho\,\delta\rho+J^{\mu T}\delta Q)$ for an outward normal on each region. Combining the two regions gives interface terms
$[\rho']\delta\rho+[J^r]^T\delta Q$.

For a source-free interface in the original theory:

$$[\rho]=0,\quad [Q]=0\pmod{2pi},\quad [\rho']=0,\quad [J^r]=0.$$

Phase continuity modulo 2pi means continuity of the circle-valued fields and of the appropriate local representatives up to constant shifts. A genuine field-value jump yields distributional derivatives whose quadratic energy is not an admissible finite-action layer when the kinetic metric is nondegenerate. Integrating the Euler–Lagrange equations across a thin layer independently recovers [rho']=0 and [J^r]=0; smooth bulk force terms do not produce delta distributions.

Continuity on the entire timelike cylinder for every t,z,phi implies continuity of the tangential phase derivatives, hence identical w,k,m across the interface in the stated ansatz. They cannot be switched off radially by a piecewise choice of constants. At fixed rho and positive M, continuity of J^r also implies continuity of u.

Distributional stress conservation, with no surface stress, requires

$$[T^r{}_\nu]=0\quad(\nu=t,r,phi,z).$$

In particular [p_r]=0. With the field and derivative conditions above, this traction matching follows from the original stress formula; it should still be checked explicitly. Continuity of arbitrary tangential stresses is not an extra general interface equation, although smooth same-action field matching makes them continuous here.

## 2. Regular cylindrical axis

A smooth positive axisymmetric rho has rho=f(r²), so rho'(0)=0. Smooth radial phase representatives likewise have q=q_0+O(r²). Bounded M and regular derivatives imply

$$C=\lim_{r\to0} r M u=0.$$

With M invertible, u=0 for the entire connected source-free radial domain. A globally defined circle-valued phase on a disk requires m=0; real-valued phases do also. A uniform positive kinetic lower bound independently excludes m!=0 by local and integrated energy divergence. For M=I+rho²A, M(0)=I at field value rho=0; amplitude collapse alone does not regularize this winding.

For general smooth V,M after m=C=0, put

$$F(\rho)=V_{,\rho}+\tfrac12[-w^TM_{,\rho}w+k^TM_{,\rho}k].$$

The axis equation (r rho')'=r F(rho) gives the local expansion

$$\rho(r)=\rho_c+\tfrac14F(\rho_c)r^2+O(r^4)$$

when the functions/solution have the corresponding smoothness. This follows by integrating the equation from zero; it is not a freely assignable derivative. rho_c is the free central amplitude; q_0 contains two irrelevant phase offsets. Positivity of rho and M must persist beyond the axis; axis regularity alone gives no global existence.

## 3. Finite annulus with specified interior and exterior

The annular radial system has the two amplitude data, two phase constants C, and two phase offsets at a regular radius. Fixed w,m,k are symmetry parameters. Equivalent Dirichlet data require solving for compatible C and rho; arbitrary simultaneous values, derivatives and fluxes at both edges generally overdetermine it. Local ODE existence proves only a small-interval solution, not solvability of arbitrary two-end boundary data.

There are two explicitly assessable continuations:

1. **The exact extended solution on all three regions.** Choose the old non-null axial family rho=rho_0, Q=kz+Q_0 in the core r<r_in, annulus and exterior r>R. All field/normal-current/traction data match, with no layer or support. Annular energy is pi(R²-r_in²)epsilon_0. The global transverse energy is infinite. This is a specified, admissible continuation demonstrating local annular feasibility, not localization.

2. **Ordinary vacuum core and ordinary vacuum exterior, same original theory and no shells.** At the outer edge the vacuum gives rho=v, rho'=0, Q constant, J^r=0; tangential continuity forces w=k=m=0 throughout the annulus. Smooth matching then gives C=0, and amplitude uniqueness at R forces rho=v on the connected annulus. The core consequently matches the same trivial solution. No nontrivial compact annular lump results. If the two vacua differ, a source-free profile which is exactly at one equilibrium with zero derivative at finite R still cannot depart from it by smooth ODE uniqueness.

An unspecified core/exterior is not boundary data. An annulus with arbitrary reflecting walls is a boundary-value problem with physical walls or prescribed boundary driving; its finite energy does not establish a free, globally conserved configuration without accounting for the walls and exterior.

For static phases the integrated radial balance is

$$R^2p_r(R)-r_{\rm in}^2p_r(r_{\rm in})
=-2\int_{r_{\rm in}}^R r U\,dr.$$

Boundary traction is a real term, not a disposable integration constant. If zero radial traction is imposed at both boundaries and U>=0, this identity already forces U=0; in the scalar-only reduction it then forces rho'=0. For a winding annulus, transverse phase gradient terms cancel in p_r+p_phi, but topology, phase boundary values and stress matching still need checking.

## 4. Asymptotic ordinary vacuum

The reference vacuum is rho->rho_v finite, Q derivatives->0, V'(rho_v)=0,
M(rho_v)>0, U=V-V(rho_v)>=0 and U(rho_v)=0. If V(rho_v) is nonzero, subtract the constant vacuum tensor only for defining excitation energy; this does not set its gravitational effect to zero or solve a gravitational equation.

In the restricted phase ansatz finite transverse energy forces constant w=k=m=0 and C=0 by the area and logarithmic bounds in DEE_LOCALIZATION_ANALYSIS.md. Remaining necessary amplitude conditions are rho->rho_v,
integral r rho'^2 dr finite and integral r U dr finite. A rigorous virial proof uses a sequence R_j with R_j²e(R_j)->0, rather than imposing an unproved pointwise rate.

For the polynomial potential at v>0, a formal small tail eta=rho-v would obey
eta''+eta'/r=2lambda_V v² eta+O(eta²).
A decaying linearized tail is only a necessary asymptotic possibility. It neither supplies a regular center nor defeats the global no-go theorem. No fitted tail or numerical shooting solution is claimed.

If M loses its positive lower bound at infinity, the tail bounds must be revisited using the actual eigenvalues. Necessary integrability tests include
integral r w^TMw dr, integral r k^TMk dr,
integral m^TMm dr/r and integral C^TM^-1 C dr/r.
They are not existence proofs. No degenerate-vacuum solution is supplied.

## 5. Finite-radius vacuum matching and no compacton result

Suppose the exterior r>=R is exactly the ordinary vacuum, V,M are C² near it, M is invertible, and there is no interface action. Continuity and flux matching give at R
rho=rho_v, rho'=0, C=0, w=k=m=0. On every interval r>0 the amplitude equation is an ordinary second-order equation with locally Lipschitz force. Its unique solution with these initial values is rho=rho_v. Connected inward continuation remains the vacuum down to the regular axis.

Therefore a nontrivial smooth compactly supported departure from vacuum is impossible in this scope. The result is a finite-radius uniqueness theorem; it does not invoke a gravitational equation or numerical solve. Non-Lipschitz potentials, a degenerate kinetic sector, discontinuous fields or explicit shell sources fall outside it.

The old non-null axial family cannot be spliced to vacuum: rho_0 differs from v, k differs from zero on the boundary, and its radial pressure -epsilon_0 differs from zero. Even the old null family, whose radial traction is zero, cannot be cut to constant phases because its nonzero tangential phase derivatives fail continuity. Zero traction alone is not matter matching.

## 6. What an actual surface layer would require

This section identifies missing support; it does not add it to the original model.

Let a static thin cylindrical surface carry a tangential surface stress S^ab, with coordinate components on (t,phi,z), no normal delta stress, and let the bulk stress be piecewise smooth. Write $p_s=S_{\hat\phi\hat\phi}=R^2S^{\phi\phi}$. Direct cylindrical distributional conservation gives

$$[T^{rr}]-R S^{\phi\phi}
=[p_r]-p_s/R=0,$$
$$D_a S^a{}_b+[T^r{}_b]=0,\qquad b=t,phi,z.$$

Equivalently [p_r]=p_s/R. The circumference pressure term is derived from
Gamma^r_phi_phi=-R. Our sign convention [outside minus inside] makes a fluid with positive interior pressure and vacuum outside require negative p_s, i.e. surface tension. This supplies an independent sign check.

For the old non-null axial interior and vacuum exterior, [p_r]=epsilon_0, so a necessary formal stress balance is p_s=R epsilon_0 (compressive surface pressure). This condition alone does not repair its field-value and tangential-gradient discontinuities. The surface energy density sigma and axial stress are not determined by this radial balance; an explicit physical surface action/equation of state and energy checks would be required.

If surface phase charge/current j^a exists, total current conservation requires
$D_a j^a+[J^r]=0$. Normal derivative jumps are scalar surface sources:
a distributional bulk equation E_rho=s_rho delta(r-R) has [rho']=s_rho,
and E_Q=s_Q delta(r-R) has [J^r]=s_Q. A consistent surface action must supply the compensating sources with its own variation/sign conventions, as well as its stress and dynamics. In the unaugmented original action these surface terms are absent, so the jumps vanish. Merely declaring a shell is not a solution.

For a smooth finite-thickness layer governed by the same original action, all original field equations and stress conservation must hold throughout it; the complete-plane obstruction still applies under its hypotheses. A hand-applied envelope cannot replace that layer.

## 7. Boundary feasibility conclusions

Regular axis plus ordinary asymptotic vacuum and finite E_perp permits only the trivial solution within the proven scope. An explicitly continued extended exact solution provides compatible annular data but infinite global transverse energy. Annular winding is locally legitimate on a genuine annulus; a smooth disk core and ordinary vacuum tail are obstructed. Exact finite-radius vacuum interfaces without new sources cannot create nontrivial compact objects. Adding walls/shells would be an augmented physical problem requiring explicit support, not an admissible hidden assumption.

No statement here generalizes these restrictions to every V,M, every topology or every nonstationary configuration allowed by the broad original action. Matter matching conditions do not establish metric/extrinsic-curvature matching, harmonic gauge, or any Einstein equation. Those remain outside this investigation.
