# DEE matter existence and admissibility analysis

Starting commit: 27ae206a7b25b2ba09fdc89cfe1bdbe7cc60dd6f.
Research date: 2026-11-07 America/Los_Angeles (2026-11-08 UTC).

**Verdict: EXISTENCE ESTABLISHED — an explicitly verified, nontrivial admissible solution exists within the stated scope.**

Scope is the original matter action on prescribed flat spacetime, with the original optional polynomial realization, finite local energy, and nonzero axial phase gradients. It is not an isolated finite-energy cylinder, not a winding solution through an axis, not a toroidal core, not a complete Einstein-equation solution, and not proof of gravitational stability or non-collapse.

## 1. Exact, non-null axial family

Use precisely manuscript Eqs. (14)–(15), not different field equations:

$$A=\begin{pmatrix}a&c_{\rm mix}\\c_{\rm mix}&b\end{pmatrix},
\quad M(\rho)=I+\rho^2A,\quad V=\lambda_V(\rho^2-v^2)^2/4.$$

Choose the symbolic open parameter set
$\lambda_V>0$, $v>0$, $a>0$, $b>0$, $0<c_{\rm mix}^2<ab$.
Thus A is positive definite and mixing is nonzero. Let a constant two-vector k have both entries nonzero and obey

$$0<C=k^TAk<\lambda_V v^2.$$

This parameter set is nonempty: for any fixed vector e with both entries nonzero, set k=s e with
$0<|s|<\sqrt{\lambda_Vv^2/(e^TAe)}$.
No fitted physical values or historical Gaussian inputs are needed for this existence proof.

Set

$$\boxed{\rho=\rho_0=\sqrt{v^2-C/\lambda_V}>0,\qquad
\Theta=k_\Theta z+\Theta_0,\qquad
\Psi=k_\Psi z+\Psi_0.}$$

Here w=m=u=0. Constants Theta_0 and Psi_0 are arbitrary. The amplitude branch is nontrivially displaced from the potential minimum: rho_0<v and V(rho_0)>0. Both transport fields have nonzero gradients and their coupling to rho is retained.

### All field equations

The full amplitude equation, evaluated without deleting any terms, is

$$E_\rho=-\rho_0[\lambda_V(\rho_0^2-v^2)+k^TAk]=0.$$

The currents are $J^t=J^r=J^\phi=0$, $J^z=M(\rho_0)k$, each constant in Cartesian coordinates, so both divergences vanish identically. This checks all three Euler–Lagrange equations, not merely the algebraic branch.

### Complete stress and conservation

Write $B=k^TM(\rho_0)k>0$, $V_0=C^2/(4\lambda_V)>0$ and
$\varepsilon_0=V_0+B/2>0$. Then

$$T_{\mu\nu}=\operatorname{diag}
(\varepsilon_0,-\varepsilon_0,-r^2\varepsilon_0,B-\varepsilon_0)$$

in (t,r,phi,z) coordinates, and all six offdiagonal components are zero.
In the orthonormal frame the third entry is -epsilon_0. In Cartesian coordinates this is simply
diag(epsilon_0,-epsilon_0,-epsilon_0,B-epsilon_0), a constant tensor.
Thus all four stress divergences vanish directly. Radial balance is
$p_r'=0$, $p_r=p_\phi=-\varepsilon_0$, so no radial force is missing.

The negative transverse pressure is a tension of this matter solution, not an external confining force. There is no cutoff or boundary wall: the configuration extends over all Minkowski space. A finite annulus is its restriction, with identical same-field continuation on either side. Matching values, normal derivatives, currents and stress is automatic; no jump, shell or supporting stress is introduced. If the annular material were cut off into vacuum, the nonzero traction would require a new matching/support analysis; this proof does not permit that cut.

### Admissibility

M has eigenvalues $1+\rho_0^2\alpha_i>1$, where alpha_i are the positive eigenvalues of A. The full three-field kinetic Hessian diag(1,M) is strictly positive. All local derivatives, currents and frame stress components are finite and smooth. In Cartesian coordinates there is no axis singularity: rho is constant and the phases depend only on z. Circle-valued phases are smooth maps exp(i k z); there is no periodic z identification. No winding quantization is needed for these axial gradients on z in R.

For every null vector n,
$T(n,n)=B(n^z)^2\ge0$ directly. WEC and DEC also hold for this diagonal example: energy epsilon_0>=0, |p_r|=|p_phi|=epsilon_0 and
$|p_z|=|B/2-V_0|\le B/2+V_0=\varepsilon_0$.
These energy conditions do not establish dynamical stability.

All stated local requirements are satisfied with no envelope or artificial source scaling. Its energy per unit z length integrated over the entire transverse plane is infinite; on any finite annulus it is
$\pi(R^2-r_{\rm in}^2)\varepsilon_0$.
It therefore cannot be used as an isolated finite-transverse-energy source without further work.

## 2. Second exact family: common null axial transport

For contrast, still within the same unchanged action, take rho=v,
Q=q(t+z)+Q_0 with a nonzero constant two-vector q. Set m=u=0.
Then all three contractions X_Theta, X_Psi and X_ThetaPsi vanish, as do V(v) and V'(v). Hence E_rho=0 for any smooth coupling functions evaluated at v with V'(v)=0; the phase currents are constant and conserved. For the polynomial family M(v)>0 as above.

With ell_mu=(1,0,0,1),

$$T_{\mu\nu}=N\ell_\mu\ell_\nu,\qquad
N=q^TM(v)q>0,\qquad \ell_\mu\ell^\mu=0.$$

This is a nonzero null-dust-form stress realized by the scalar fields. Its only nonzero coordinate components are T_tt=T_tz=T_zz=N. The physical z energy flux is -N, consistent with dependence on t+z. It is smooth, conserved and has positive finite local density, but is also spatially extended with infinite total/transverse energy. Its stress is stationary although its phases vary linearly in time. This family is an additional exact check, not needed to rescue the non-null family.

## 3. Winding and the constant-amplitude obstruction

For the polynomial model with zero radial currents and positive M, u=0.
If rho=rho_0>0 is constant, the amplitude equation requires for every radius in a nonempty open interval

$$\lambda_V(\rho_0^2-v^2)-w^TAw+k^TAk+m^TAm/r^2=0.$$

Differentiate with respect to r to obtain $-2m^TAm/r^3=0$. Thus a necessary condition is m^TAm=0. If A is positive definite, this forces m=0.
Consequently a constant positive amplitude carrying nonzero winding is impossible for this chosen positive-definite A on an interval. A radius-by-radius use of manuscript Eq. (18) while still declaring rho constant is inconsistent.

This is an obstruction to that specific constant-amplitude winding ansatz, not to the full original action: positive M alone does not imply positive A, and variable rho introduces derivatives. For A with a nonzero null direction the constant-amplitude winding obstruction can change, but axis topology and finite-energy requirements must still be checked. No coefficients were changed to force a winding PASS here.

## 4. Rigorous local annular existence with winding

This result is separate from the exact axial family. Let r_0>0,
rho(r_0)=rho_*>0, rho'(r_0)=p_*, and M(rho_*)>0. Take smooth couplings, V and fixed finite w,m,k,C. Substitute u=M^-1 C/r into the derived amplitude equation:

$$\rho'=p,$$
$$p'=-p/r+V_{,\rho}
+\tfrac12[-w^TM_{,\rho}w+u^TM_{,\rho}u
+m^TM_{,\rho}m/r^2+k^TM_{,\rho}k].$$

Together with Q_radial'=M^-1 C/r this is a first-order system with locally Lipschitz right-hand side on the open set r>0,rho>0,det M>0 (C² couplings suffice locally). Picard–Lindelof yields a unique local solution for arbitrary finite initial phase offsets. By continuity, a sufficiently small closed interval about r_0 retains strictly positive rho and M and bounded derivatives/energy. It satisfies all field equations and hence conserved stress, including when integer m!=0 for angular-valued phases.

This proves local existence on some annulus with data chosen as above. It does not prove a solution extends to preassigned radii, satisfies two prescribed boundary values, matches an unknown core/exterior, remains bounded globally, or has finite integrated energy. For C=0 the phase radial parts are constant. This local theorem does not bypass the axis obstruction. The continuation criterion fails if rho reaches zero, M loses invertibility, r reaches zero with winding, or fields/derivatives become unbounded.

## 5. Regular-axis and localization obstructions

For strictly positive uniform M and m!=0, the frame angular gradient m/r gives energy at least delta |m|²/(2r²). The axis energy diverges and the circle-valued phase cannot extend continuously over the disk. In this polynomial A>=0 family M>=I even when rho tends to zero, so amplitude suppression cannot regularize winding. This is stronger than merely observing a coordinate singularity.

Any proposed radial multiplication f(r)T of a conserved stress has

$$\nabla_\mu(fT^\mu{}_\nu)=f'(r)T^r{}_\nu.$$

For the exact non-null axial family $T^r{}_r=-\varepsilon_0\ne0$, so a nonconstant Gaussian multiplier violates radial conservation. For the null family T^r_nu=0, multiplying the stress alone can preserve conservation, but conservation is not matter closure: changing the scalar amplitude/phases to generate such an envelope must independently solve their field equations. No Gaussian matter profile is claimed.

No global nonexistence theorem for every localized solution of the original three-scalar theory has been established here. A finite-energy localized winding solution, vacuum matching, nonlinear stability and gravitational backreaction remain open, rather than being inferred from the exact extended solution.

## 6. Verdict and unresolved questions

**EXISTENCE ESTABLISHED**, within the explicit matter-only, finite-local-energy, axial-gradient scope, by the non-null exact family in section 1. It exists for a nonempty open set of coefficients/gradients in the manuscript's original polynomial realization. The general unspecified functions of Eq. (3) are not claimed to admit it for every choice.

Other results have different logical status: the winding theorem is local existence only; the constant-amplitude positive-A winding case and the regular-axis uniformly-positive-M winding case have identified obstructions. No boundary-value solution or localized cylindrical object is established.

Unresolved: physical scale/stiffness restoration and empirical motivation for coefficients; a physically selected rather than existential parameter set; nonzero-winding global matter closure and core/exterior matching; finite transverse energy; coupled matter stability; and every gravitational equation/gauge/boundary question. Original missing historical-source inputs remain missing. The formal correction to Eq. (34) stands unchanged; the historical consistency PASS and frozen Run 43 evidence are not used as existence or Einstein-equation proof.

For G>0 the displayed nonzero T on flat spacetime fails G_mu_nu=8pi G T_mu_nu, since flat G_mu_nu=0. Flat space here is explicitly a prescribed matter-test background. This investigation constructs no gravitational field and stops after committing the analytical evidence.
