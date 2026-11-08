# DEE matter-field derivation

Research date: 2026-11-07 America/Los_Angeles (2026-11-08 UTC).
Starting commit: 27ae206a7b25b2ba09fdc89cfe1bdbe7cc60dd6f.
Scope: analytical matter equations on a prescribed flat background only.

## 1. Authority and assumptions

Authoritative recovered source: evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex, Git blob a20c0ad7b3469487e50fa2be3a193d1aacc92ab6. Equations (3)–(9) define the action, equations, currents and stress; (14)–(15) define an optional polynomial realization; (22)–(24) give kinetic positivity. EQUATIONS.json supplies verbatim equation blocks and source lines. Also read: preserved source README and both September 23 audit records; docs/RECOVERED_NONLINEAR_DEE_FORMULATION.md; FORMAL_CORRECTION_EQ34.md; NEW_PHYSICAL_TEST_PROTOCOL.md; NEW_PHYSICAL_TEST_SPECIFICATION_V0_1.md, all at the starting commit. No later manuscript is silently substituted. The exact later July 21 TeX, historical Gaussian generator, and missing earlier Eq. (34)-specific audit remain unrecovered.

Assumptions: real fields; signature (-,+,+,+); coordinates (t,r,phi,z); c_light=1; ds²=-dt²+dr²+r²dphi²+dz²; sqrt(-g)=r for r>0. The symbol rho is a matter scalar, not the coordinate r. V,K_Theta,K_Psi,Lambda depend only on rho. For derivation take C² fields/functions, and stronger smoothness where axis extension is discussed. Constants omega,m,k are real except that m must be integer for 2pi-periodic angular-valued phases. Nonzero m is not allowed for a globally single-valued real-valued phase. Stationarity and axisymmetry mean invariant gradients and stress; a phase itself can shift under a Killing flow. No source envelope, matter-action scaling, external force, gravitational perturbation or gravitational equation is introduced. The present user instruction authorizes this analytical investigation; it does not execute the earlier numerical protocol.

Write Q=(Theta,Psi)^T and

$$M(\rho)=\begin{pmatrix}K_\Theta&\Lambda\\\Lambda&K_\Psi\end{pmatrix},
\quad {\cal L}=-\tfrac12\partial_\mu\rho\partial^\mu\rho
-V-\tfrac12\partial_\mu Q^T M\partial^\mu Q.$$

This is Eq. (3)'s original matter action, in matrix notation.

## 2. Variation before imposing symmetry

Under compactly supported scalar variations, integration by parts gives

$$\delta S_m=\int\sqrt{-g}\,[E_\rho\delta\rho+
E_\Theta\delta\Theta+E_\Psi\delta\Psi]\,d^4x+\text{boundary terms},$$

$$E_\rho=\Box\rho-V_{,\rho}
-\tfrac12K_{\Theta,\rho}(\partial\Theta)^2
-\tfrac12K_{\Psi,\rho}(\partial\Psi)^2
-\Lambda_{,\rho}\partial\Theta\cdot\partial\Psi,$$

$$E_Q=\nabla_\mu(M\partial^\mu Q),\qquad E_\rho=0,\quad E_Q=0.$$

Indeed delta L from a rho variation is
$-\partial^\mu\rho\,\partial_\mu\delta\rho-
[V_{,\rho}+\frac12\partial_\mu Q^TM_{,\rho}\partial^\mu Q]\delta\rho$.
The phase variation is $-(M\partial^\mu Q)^T\partial_\mu\delta Q$.
The boundary integrand is
$-n_\mu[\partial^\mu\rho\,\delta\rho+J^\mu{}^T\delta Q]$.
Fixed field values, compact support, or compatible flux data make the variational problem well defined. Do not impose both arbitrary values and arbitrary fluxes.

The two shift currents, with the manuscript's sign convention, are

$$J^\mu=M\partial^\mu Q.$$

The canonical Noether current for delta Q=+constant can instead be written -J; this overall sign does not change its conservation.

With $T_{\mu\nu}=-2\,\delta S_m/(\sqrt{-g}\delta g^{\mu\nu})$, use
$\delta\sqrt{-g}=-\sqrt{-g}\,g_{\mu\nu}\delta g^{\mu\nu}/2$.
Symmetric inverse-metric variation gives

$$T_{\mu\nu}=\partial_\mu\rho\partial_\nu\rho+
\partial_\mu Q^TM\partial_\nu Q+g_{\mu\nu}{\cal L}.$$

In particular the mixing term contributes
$\Lambda(\partial_\mu\Theta\partial_\nu\Psi+
\partial_\mu\Psi\partial_\nu\Theta)$; there is no missing factor of two.

## 3. Stationary cylindrical reduction

Let $w=(\omega_\Theta,\omega_\Psi)^T$, $m=(m_\Theta,m_\Psi)^T$,
$k=(k_\Theta,k_\Psi)^T$, and $u=(\theta',\psi')^T$, with primes here meaning d/dr. Derivatives of couplings are labeled ,rho. Define

$$X_\Theta=-\omega_\Theta^2+\theta'^2+m_\Theta^2/r^2+k_\Theta^2,$$
$$X_\Psi=-\omega_\Psi^2+\psi'^2+m_\Psi^2/r^2+k_\Psi^2,$$
$$X_{\Theta\Psi}=-\omega_\Theta\omega_\Psi+\theta'\psi'
+m_\Theta m_\Psi/r^2+k_\Theta k_\Psi.$$

Thus

$$W=-w^TMw+u^TMu+m^TMm/r^2+k^TMk,\qquad
{\cal L}=-\rho'^2/2-V-W/2.$$

Since $\Box\rho=(r\rho')'/r$, the first field equation is

$$\boxed{\rho''+\rho'/r-V_{,\rho}
-\tfrac12K_{\Theta,\rho}X_\Theta
-\tfrac12K_{\Psi,\rho}X_\Psi
-\Lambda_{,\rho}X_{\Theta\Psi}=0.}$$

The complete current components (each entry is a two-vector) are

$$J^t=-Mw,\quad J^r=Mu,\quad J^\phi=Mm/r^2,\quad J^z=Mk.$$

All t,phi,z derivatives in their divergence vanish. Consequently the other two field equations are

$$\boxed{[r(K_\Theta\theta'+\Lambda\psi')]'/r=0,\qquad
[r(K_\Psi\psi'+\Lambda\theta')]'/r=0.}$$

Their first integrals are $rMu=C$, a constant two-vector; where M is invertible,

$$u=M^{-1}C/r=
\frac1{rD}\begin{pmatrix}K_\Psi C_\Theta-\Lambda C_\Psi\\
K_\Theta C_\Psi-\Lambda C_\Theta\end{pmatrix},
\qquad D=K_\Theta K_\Psi-\Lambda^2.$$

The radial fluxes per unit z length through a cylinder are 2pi C. Zero normal current sets C=0 and, for positive M, u=0. Radial profiles cannot otherwise be chosen freely.

For the optional original polynomial realization put
$A=\begin{pmatrix}a&c_{\rm mix}\\c_{\rm mix}&b\end{pmatrix}$,
$M=I+\rho^2 A$, $V=\lambda_V(\rho^2-v^2)^2/4$. Then

$$\boxed{\rho''+\rho'/r-
\rho\{\lambda_V(\rho^2-v^2)-w^TAw+u^TAu+
m^TAm/r^2+k^TAk\}=0.}$$

No homogeneous algebraic branch formula is valid for arbitrary variable rho without the rho''+rho'/r term.

## 4. All ten stress components

All entries below are covariant coordinate components; symmetry supplies the rest:

| Component | Expression |
|---|---|
| tt | $w^TMw-{\cal L}$ |
| tr | $w^TMu$ |
| tphi | $w^TMm$ |
| tz | $w^TMk$ |
| rr | $\rho'^2+u^TMu+{\cal L}$ |
| rphi | $u^TMm$ |
| rz | $u^TMk$ |
| phiphi | $m^TMm+r^2{\cal L}$ |
| phiz | $m^TMk$ |
| zz | $k^TMk+{\cal L}$ |

The orthonormal coframe is (dt,dr,r dphi,dz). Divide each covariant phi index by r to get frame components. In particular

$$\varepsilon=T_{\hat t\hat t}
=\tfrac12\rho'^2+V+\tfrac12(w^TMw+u^TMu+m^TMm/r^2+k^TMk),$$

$$p_r=\rho'^2+u^TMu+{\cal L},\quad
p_\phi=m^TMm/r^2+{\cal L},\quad p_z=k^TMk+{\cal L}.$$

With +omega t convention T_ti has the displayed positive bilinear sign; physical energy flux measured by static observers is -T_ti in their frame.

## 5. Conservation, including radial balance

The nonzero background Christoffels are
$\Gamma^r_{\phi\phi}=-r$, $\Gamma^\phi_{r\phi}=\Gamma^\phi_{\phi r}=1/r$.
For mixed stress,

$$\nabla_\mu T^\mu{}_\nu=
r^{-1}\partial_\mu(rT^\mu{}_\nu)
-\tfrac12(\partial_\nu g_{\alpha\beta})T^{\alpha\beta}.$$

The four equations on shell reduce to

$$[rT_{rt}]'/r=0,\quad [rT_{r\phi}]'/r=0,\quad [rT_{rz}]'/r=0,$$
$$\boxed{p_r'+(p_r-p_\phi)/r=0.}$$

For example the last equation is $T_{rr}'+T_{rr}/r-T_{\phi\phi}/r^3=0$.
The three tangential flux integrals are respectively w^T C, m^T C, k^T C.
To check radial balance directly, differentiate p_r using W and obtain

$$p_r'+(p_r-p_\phi)/r=\rho'E_\rho+
u^T\{(Mu)'+Mu/r\}.$$

The explicit derivative of $m^TMm/r^2$ contributes the term canceled by
$(p_r-p_\phi)/r$; dropping it would produce a false force. More generally,

$$\nabla_\mu T^\mu{}_\nu=
E_\rho\partial_\nu\rho+E_Q^T\partial_\nu Q.$$

This is an off-shell identity, with the residual signs defined above, and on-shell conservation follows for all four components. Conservation by itself need not imply all field equations.

## 6. Kinetic positivity and energy conditions

On a flat t slice the velocity Hessian is diag(1,M). Completing the phase square,

$$x^TMx=K_\Theta(x_\Theta+\Lambda x_\Psi/K_\Theta)^2+
(D/K_\Theta)x_\Psi^2.$$

Strict positivity is equivalent to $K_\Theta>0,D>0$ (and implies $K_\Psi>0$).
The eigenvalues are
$(K_\Theta+K_\Psi\pm\sqrt{(K_\Theta-K_\Psi)^2+4\Lambda^2})/2$.
For the polynomial model,
$D=1+(a+b)\rho^2+(ab-c_{\rm mix}^2)\rho^4$.
A sufficient, not necessary, condition at every rho is A positive semidefinite;
strictly positive a,b and c_mix²<ab is a convenient open subset. Kinetic positivity alone says nothing about the sign of V or full dynamical stability.

The Hamiltonian density is the energy expression above, generalized to all spatial gradients. With M positive and V>=0 it is nonnegative. For any null vector n,

$$T_{\mu\nu}n^\mu n^\nu=(n\cdot\partial\rho)^2+
(n\cdot\partial Q)^TM(n\cdot\partial Q)\ge0.$$

This proves the entire null-cone NEC without directional sampling. WEC follows here from the same positive energy argument in any local orthonormal observer frame when V>=0. No coupled perturbative stability claim follows.

## 7. Regularity, topology, boundaries and dimensions

On a closed annulus 0<r_in<=r<=R, demand smooth fields, bounded derivatives and couplings, rho bounded strictly above zero and M bounded strictly positive. Angular-valued phases permit integer winding; patchwise phase representatives differ by constants. Phase gradients, currents and stress must be globally defined. This does not supply an omitted core.

For an axis extension with smooth axisymmetric real rho and radial phase parts, rho=f(r²), theta=g(r²), psi=h(r²) locally (with appropriate differentiability). In particular their first radial derivatives vanish at the axis. Smooth bounded M and currents imply C=lim rMu=0, hence u=0. A nonzero angular winding around a contractible axis cannot extend as a continuous circle-valued phase on the disk. Independently, if M has a uniform positive lower bound near the axis, its winding energy is at least const*|m|²/(2r²), so finite local energy excludes m!=0. For M=I+rho²A with A positive semidefinite, M>=I; setting rho to zero does not remove this divergence. A different target topology or degenerating kinetic matrix would be an additional model/scope assumption, not a cure within strict axis positivity.

At a finite interface without a shell require field continuity, matching rho normal derivative and phase normal currents, and stress traction matching. A cut to vacuum generally fails these requirements. A finite computational annulus can instead be a restriction of an explicitly given smooth continuation, without a material wall. Natural no-current conditions imply C=0; arbitrary Dirichlet plus arbitrary Neumann data overdetermine the second-order equations. There are six local radial integration constants, or two from rho and two phase fluxes plus two phase offsets after first integration; prescribed w,m,k are symmetry parameters.

Finite local density differs from finite transverse energy per unit length
$2pi\int_0^\infty r\varepsilon\,dr$ and from finite total energy on infinite z.
Nonzero constant density has neither finite transverse nor total energy. No decay or isolated-source claim is implicit.

The preserved manuscript is a normalized action and explicitly leaves restoration of reference mass scales open. No numerical physical constants are selected here. All equations require consistent field/length normalization: each contraction multiplying K or Lambda and V must have the same action-density units. For dimensionless angular phases one can restore a common phase stiffness F², use M_dim=F²(I+rho²A), [A]=[rho]^-2, [k]=[omega]=length^-1, [m]=1. With canonical rho and hbar=c_light=1, [rho]=[F]=mass and [lambda_V]=1. This is a dimensional restoration convention, not a changed equation or fitted parameter: the formulas in this package use the original normalized F=1 units. Quantitative SI predictions and any gravitational coupling require a separately frozen scale map.
