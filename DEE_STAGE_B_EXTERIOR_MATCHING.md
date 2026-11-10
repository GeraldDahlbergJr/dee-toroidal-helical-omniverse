# Stage B: exterior matter matching

Date: 2026-11-10 UTC. Base: `6b3c46da4962db7f1b48abfee604e0201b9530e2`.
Branch: `research/stage-b-exterior-matching`.

## Inputs, equations and meaning of verification

**VERIFIED (scope):** scientific inputs read were the Stage A½ bounded response,
boundary compatibility and backreaction check, the symbolic recovery document
for the metric components, and the original covariant action in
`evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex`, Eqs. (3)–(9),
(14)–(15). No prior evidence, manuscript, Run 43 result or verdict is changed.
VERIFIED below means an analytical result with the indicated exact identity
checks; it does not mean independent peer review or nonlinear existence.

The unchanged action is Einstein gravity plus

    L_m = -(grad rho)^2/2 - V(rho) - (grad Q)^T M(rho)(grad Q)/2,
    M=I+rho² A, A=[[a,c],[c,b]], V=lambda_V(rho²-v²)²/4.

**CONDITIONAL (assumptions throughout):** eta=(-,+,+,+), G>0,
lambda_V,rho_0,v>0, M_0 positive definite, real k, and
rho_0²=v²-C/lambda_V>0. Define C=k^T A k, B=k^T M_0 k,
epsilon=V(rho_0)+B/2, m²=2lambda_V rho_0²>0, and
F=2pi G rho_0 C(B+2epsilon). Use static, radial, z-independent u,
rho=rho_0+u, Q=kz+Q_0, and the verified leading metric. This is a
weak-gravity hierarchy: h,u=O(G) on fixed finite domains; Einstein's equation
is checked through O(G), and matter equations through their first response.
It is not an exact linearization around an Einstein solution with flat metric
and nonzero background stress. All omitted residuals start at O(G²) under
these assumptions; they have not been solved.

The inherited original-action equations are

    E_rho = Box_g rho - V'(rho) -rho (grad Q)^T A(grad Q)=0,
    E_Q = nabla_mu[M(rho) partial^mu Q]=0,
    E_rho^(1) = (Delta_perp-m²)u+rho_0 C h_zz = 0,
    h_zz=-2pi G(B+2epsilon)r²,
    (Delta_perp-m²)u=Fr².

## Responding exterior: a finite collar construction

**VERIFIED:** for any radial disk or annular Stage A½ response, let its
outer radius be R and its traces be g=u(R), q=u'(R). Continue the *same
background matter and polynomial metric* into R<r<R+delta. Set

    p=-F(r²/m²+4/m⁴), I=I_0(mr), K=K_0(mr),
    W_R=I(R)K'(R)-I'(R)K(R)=-1/R,
    A_ext=[(g-p(R))K'(R)-(q-p'(R))K(R)]/W_R,
    D_ext=[I(R)(q-p'(R))-I'(R)(g-p(R))]/W_R,
    u_ext=p+A_ext I+D_ext K.

The Bessel Wronskian is nonzero for every R>0. This solves the amplitude
EOM, u_ext(R)=g and u_ext'(R)=q exactly in the retained equation. ODE
uniqueness makes these coefficients equal to the interior coefficients.
Consequently every derivative also agrees across R: the interface is an
artificial division of one smooth response, not an imposed material wall.
The result applies to zero Dirichlet, zero Neumann or other radial traces;
it does not impose both zero traces. Nonconstant u is compulsory on an open
collar if F!=0, since a constant cannot satisfy Lu=Fr² there.

**VERIFIED:** the two phases have identical traces Q=kz+Q_0. Their radial
normal currents vanish separately: g^{rz}=0, partial_r Q=0, so n_mu J^mu=0.
Their divergence vanishes even for varying rho(r): the only current direction
is z and all coefficients are z-independent. In fact this phase statement is
exact for any diagonal, z-independent cylindrical metric. The scalar normal
flux is continuous because u' and the common metric are continuous. There is
no delta source in either scalar equation. The companion gravitational file
verifies the corresponding source-free junctions.

**CONDITIONAL (physical admissibility of the collar):** if the inherited
solution has strict weak-field, small-amplitude and positive-M margins on its
closed domain, continuity gives some delta>0 preserving these margins. No
arbitrary exterior width is promised. Sufficient checks on the enlarged domain
are max|h_hat a hat b|<<1, |u|/rho_0<<1 and
||A|| |2rho_0 u+u²|<lambda_min(M_0); derivatives must retain the fixed-domain
O(G) scaling. Exponentially growing Bessel coefficients can make the allowed
collar narrow. Finite, smooth mathematical continuation alone is insufficient
to prove weak-field admissibility at a prescribed outer radius.

**VERIFIED:** an annulus admits the same inward continuation on any collar
away from r=0. **CONDITIONAL:** filling its whole central hole smoothly in this
same radial family requires D_ext=0, since K_0 is logarithmically singular.
Generic annular boundary data need not meet that condition. A disk solution
has D_ext=0 already. These statements do not construct arbitrary nonradial
elliptic Cauchy continuations.

## Unchanged exterior and isolation obstructions

**OBSTRUCTION IDENTIFIED:** retaining u=0 throughout an exterior open set
with the same metric and fixed phases gives E_rho^(1)=-Fr², nonzero for
F!=0. Even disregarding that bulk failure, joining unchanged rho and its
normal derivative requires both u=0 and partial_n u=0. Green's identity with
I_0(mr)>0 gives F integral I_0(mr)r² d²x=0 if both traces vanish on the
entire disk/annulus boundary, impossible for F!=0. This inherited obstruction
is cited, not unnecessarily recalculated. For an annulus this proof requires
both boundary circles; specifying Cauchy data only at one circle is different.
F=0 removes this amplitude obstruction, but does not remove nonzero background
gravitation or establish isolation. For nonzero k and positive M_0,V>=0,
B+2epsilon=2(B+V)>0, so the nontrivial forcing exception is C=0.

**OBSTRUCTION IDENTIFIED:** extending this same radial family to r=infinity
cannot yield a uniformly bounded response for F!=0: p grows quadratically,
I_0 exponentially, and K_0 decays. No combination cancels p at infinity.
Independently, the polynomial metric grows as r² whenever the physical
background stress is nonzero, so it loses weak-field control even in the
F=0 exception. These are failures of this continuation, not theorems excluding
all nonlinear exteriors.

**OBSTRUCTION IDENTIFIED (isolation under retained axial symmetry):** for
nonzero k, positive M_0 gives B>0 and epsilon>0. The background energy density
is constant over the original plane. Even if one could localize it transversely,
retaining translation invariance in z and a nonzero positive energy per unit
length gives infinite total energy along the axis. An asymptotically flat,
finite-mass isolated object in the usual four-dimensional sense cannot retain
this background to infinity. Einstein's equations also preclude standard
asymptotically flat falloff while T tends to the same nonzero constant.
No asymptotic-flatness boundary condition is imposed here. A finite axial
extent, different asymptotic matter or different topology changes the problem
and has not been constructed; compactifying z also requires phase topology
and does not by itself remove transverse nondecay.

**NOT YET ESTABLISHED:** a transition to another exterior geometry, z-dependent
phase relaxation, a vacuum reached through a finite transition region, a
nonradial continuation, a complete nonlinear spacetime, or stability. These
possibilities are not ruled out by the restricted radial collar theorem.

## Reproduction and verdict

**VERIFIED:** `python DEE_STAGE_B_VERIFIER.py` checks exterior Cauchy matching,
its EOM by Bessel identities, both currents, cylindrical junction geometry,
normal stress constraint, and radial conservation. The independent routes
are closed-form matching versus ODE uniqueness, and direct extrinsic-curvature
expansion versus the Gauss constraint. Analytical sign/asymptotic arguments
above are not replaced by numeric samples. No numerical PDE or nonlinear
solver, prior verifier rerun, or Run 43 execution is involved.

**CONDITIONAL (Stage B result):** a smooth, source-free, physically admissible
*finite collar at the retained perturbative order* exists for radial responses
with strict admissibility margins. **NOT YET ESTABLISHED:** a globally admissible
exterior or complete Einstein–DEE spacetime. The standing historical matter
existence and localization verdicts and Eq. (34) correction keep their prior
scopes. Stage B stops after publication; Stage C is not authorized.
