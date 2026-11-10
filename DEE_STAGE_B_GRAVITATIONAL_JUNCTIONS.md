# Stage B: gravitational junctions and admissibility

Date: 2026-11-10 UTC. Base: `6b3c46da4962db7f1b48abfee604e0201b9530e2`.
Assumptions and classifications are defined in `DEE_STAGE_B_EXTERIOR_MATCHING.md`.

## Interface geometry

**VERIFIED:** in cylindrical coordinates the inherited metric is

    ds²=(-1+a_t r²)dt²+(1+a_r r²)dr²
        +r²(1+a_phi r²)dphi²+(1+a_z r²)dz²+O(G²),
    a_t=2pi G(2epsilon-B), a_r=2pi G(B-3epsilon),
    a_phi=2pi G(B-epsilon), a_z=-2pi G(B+2epsilon).

The cylindrical conversion includes the transverse Cartesian tensor, not just
its diagonal entries. At the timelike cylinder r=R use the *same* normal
pointing toward increasing r on both sides, n=(g_rr)^(-1/2)partial_r, and
K_ab=(1/2)(g_rr)^(-1/2)partial_r gamma_ab. Using outward normals separately
would reverse a sign and create a false jump. To first order,

    gamma_tt=-1+a_t R², gamma_zz=1+a_z R²,
    gamma_phiphi=R²+a_phi R⁴,
    K_tt=a_t R, K_zz=a_z R,
    K_phiphi=R+(2a_phi-a_r/2)R³.

K_phiphi includes the background curvature R and the normal normalization;
matching coordinate derivatives alone is not the invariant general condition.
For general radial H_t,H_z,H_phi,H_r in this coordinate form,
K_tt=H_t'/2, K_zz=H_z'/2 and
K_phiphi=r+r H_phi+r² H_phi'/2-r H_r/2 at this order.

**VERIFIED:** Einstein gravity with no surface stress requires [gamma_ab]=0
and [K_ab]=0 (Darmois/Israel conditions for a nonnull interface). The collar
construction uses the same metric function on both sides, hence all six
quantities above agree; the whole metric and derivatives agree as well.
Together with the amplitude/phase flux checks, this is simultaneous matter
and gravitational matching with no distributional layer. The Einstein source
remains T^(0)=diag(epsilon,-epsilon,-epsilon,B-epsilon) in Cartesian coordinates.
There is no jump to vacuum. Smoothly continuing a solution of the leading
metric equation preserves that equation without recalculating Stage A½.

**CONDITIONAL (order limit):** u=O(G) changes T by O(G), so G delta T is
O(G²). Leading metric matching and the O(G) matter response are consistent
in this hierarchy, but neither second-order Einstein equations nor their
junction data have been solved. This is not a claim of an exact coupled
Einstein–DEE solution.

## Why an immediate vacuum exterior fails

**OBSTRUCTION IDENTIFIED:** for an interface with [gamma]=[K]=0, the Gauss
and Codazzi constraints require continuity of T_nn and T_na. The inherited
radial normal stress is T_nn=-epsilon+O(G); T_na=0. An immediate vacuum has
T_nn=0. With k!=0 and positive M_0,V>=0, epsilon>0, so these cannot agree
at leading order. This rules out *any* vacuum geometry directly joined to
these interior data without a surface layer, not just a particular vacuum
metric or an asymptotically flat choice. An O(G) matter response cannot cancel
a nonzero zeroth-order pressure in the stated hierarchy. If epsilon=0 this
argument alone is silent; under the physical assumptions that means k=0 and
V(rho_0)=0, the trivial background, not the nontrivial axial family.

**VERIFIED:** a direct independent constraint calculation gives the same
result. The induced metric on each fixed-r cylinder has constant coefficients
in (t,phi,z), so its intrinsic scalar curvature is zero locally. For a spacelike
normal, G_nn=(K²-K_ab K^ab-{}^(3)R)/2. Raising indices with gamma gives

    G_nn^(1)=a_z-a_t=-8pi G epsilon.

It agrees with 8pi G T_nn^(0). Therefore a vacuum side sharing gamma and K
would inherit this nonzero constraint value, contradicting vacuum. The
verifier computes this expression from the full diagonal cylinder metric,
not by asserting the inherited Einstein equations.

**OBSTRUCTION IDENTIFIED (matter view):** smooth scalar traces enforce the
same tangential phase gradient k on the interface and rho=rho_0+O(G).
For positive M, the null contraction includes (l·partial Q)^T M(l·partial Q),
which is nonzero for some null directions when k!=0. This is also incompatible
with a smooth immediate vacuum. A finite responding transition region could
alter these data and is **NOT YET ESTABLISHED**, rather than excluded by the
interface argument.

## Conservation, regularity and energy

**VERIFIED:** the original-action identity is

    nabla_mu T^mu_nu=E_rho partial_nu rho+(partial_nu Q)^T E_Q.

For the radial first-order response, mixed stress variations are

    delta T^r_r=delta T^phi_phi=delta T^t_t=(B/2)h_zz,
    delta T^z_z=2rho_0 C u-(B/2)h_zz.

The radial divergence including connection terms is
partial_r delta T^r_r-B Gamma^z_zr
=B h_zz'/2-B h_zz'/2=0. The t,phi,z divergences vanish by symmetry.
Thus conservation is compatible both in the bulk and distributionally at the
smooth join. As in Stage A½, first-order conservation alone does not test the
amplitude equation, since partial rho_0=0; the exterior amplitude EOM is
checked separately. The linear stress formulas are in mixed cylindrical
components; using partial^mu delta T_mu nu alone would omit necessary terms.

**CONDITIONAL:** finite-radius collar fields and their derivatives are smooth
and rho stays positive if |u|/rho_0<<1. No K_0 singularity occurs away from the
axis. To fill an annular core within this family requires the additional
zero-K_0 condition in the matter file. The center of a continued disk remains
regular because its original I_0 mode is analytic in r².

**VERIFIED (structural energy theorem):** whenever M(rho) is positive definite
and V>=0, the original action obeys NEC and WEC. For any null l,
T_ll=(l·partial rho)²+(l·partial Q)^T M(l·partial Q)>=0. In an observer's
orthonormal frame the energy density is one half the sum of all temporal
and spatial kinetic squares (after positive kinetic diagonalization) plus V,
so it is nonnegative. The energy flux magnitude is bounded by that density
using 2|ab|<=a²+b², giving DEC as well. This is a pointwise property of the
original action, applicable to any actual regular solution satisfying these
hypotheses. It does not prove existence or turn the truncated fields into an
exact solution. Saturated directions must not be assessed by an arbitrary
sign of discarded higher-order terms. **CONDITIONAL:** M positivity on the
collar follows from the strict bound in the companion file, not from boundedness
alone. No universal strong energy condition follows: in the axial background
rho_energy+p_r+p_phi+p_z=B-2epsilon=-2V, negative when V>0.

## Conclusion ledger and reproducibility

| Conclusion | Classification | Scope |
|---|---|---|
| Radial responding collar satisfies scalar EOM and all no-layer junctions | VERIFIED | Retained perturbative equations, finite radius |
| Collar is weak and has positive kinetic matrix | CONDITIONAL | Strict initial margins; sufficiently small collar |
| Frozen exterior for nonzero forcing | OBSTRUCTION IDENTIFIED | Same metric and phases, F!=0 |
| Immediate vacuum join | OBSTRUCTION IDENTIFIED | Nonzero epsilon, no shell, retained interior data |
| Globally bounded same-family weak response | OBSTRUCTION IDENTIFIED | Nonzero forcing; metric also fails for nonzero background |
| Isolated finite-mass solution with retained infinite axial matter | OBSTRUCTION IDENTIFIED | Nondecaying positive axial background |
| Other geometry, transition region or complete spacetime | NOT YET ESTABLISHED | Requires further equations/boundaries |

**VERIFIED (reproduction):** run `python DEE_STAGE_B_VERIFIER.py` from the
repository root with Python 3 and SymPy 1.14.0. It raises on any failed exact
identity, uses no fitted parameters and no numerical PDE. Positivity,
Wronskian nonvanishing, weak-field margins, global asymptotics and scope
restrictions have the analytical proofs stated here and in the companion file.

**VERIFIED (observed run, 2026-11-10 UTC):** the targeted verifier exited 0
with SymPy 1.14.0. All 26 printed exact identity checks were zero, including
both exterior traces, the amplitude equation, both phase currents, the direct
Gauss constraint and all four conservation components. Final output:
`PASS: Stage B exterior Cauchy matching, currents, junction geometry, constraint and conservation.`
