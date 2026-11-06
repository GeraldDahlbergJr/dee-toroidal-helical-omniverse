# Run 43 — bounded coupled-estimate investigation

Base: `c80dd4b8ca09f8da402c896b1082c1b4320de40a`.
Branch: `research/run43-coupled-estimates-proof`.
One agent; proof only; continuation of Section 7 of
[the conforming mixed report](RUN43_CONFORMING_MIXED_PROOF_20261105.md).
Its containment, variational adjoint and Dirichlet Korn results are accepted
under their documented assumptions, without repeating their proofs.

**Decision:** conditional sufficient estimates can be stated and derived, but
none establishes the actual coupled solution sequence or the twelve finite-grid
ratios. Stop at the missing discrete positivity, nonlinear stability and sampled
error estimates. This is neither an impossibility theorem for this pair nor an
authorization to implement it.

Only the prior report and the original source/audit definitions were consulted:
`reproduce/solve_toroidal_helical_coupled_constraints.py`,
`reproduce/validate_toroidal_helical_dee_source.py`,
`reproduce/dee_stress_energy_3p1.py`, `research/audit_run43_arrays.py`, and
`reproduce/diagnose_toroidal_boundary_convergence.py`. No historical scan,
experiment, solve, sweep or Actions run was performed. All equations, exact
integration assumptions, geometry, boundary conditions and audit are retained.

## 1. Positivity — CONDITIONAL

Write the unchanged energy as E(u)=A+B(x)u^-4, where
A=(PI^T K PI)/2+U and B=spatial/2. These coefficients are nonnegative:
K has scalar entry 1 and angular block with positive diagonal and determinant
(1+.30 rho^2)(1+.20 rho^2)-(.10 rho^2)^2
=1+.50 rho^2+.05 rho^4>0; U=(rho^2-1)^2>=0.
The sources are smooth and bounded on the fixed shell. Thus the Hamiltonian is

    -Delta u = F(u,sigma)
             = u^-7 |sigma|^2/8 + 2 pi G (A u^5+B u),
    u=1 on both radial boundaries, sigma=LW.

For an already positive weak solution for which F belongs to H^-1 and is a
nonnegative distribution, testing with (1-u)_+ gives u>=1 almost everywhere.
With continuity this holds throughout the closed domain. This is a continuum
comparison argument conditional on existence and admissibility, not an
existence proof for a singular equation. It provides no upper barrier.

The same test is generally not in the enriched conforming scalar space.
Exact integration, positive energy and Dirichlet values do not prove a discrete
maximum principle; no appropriate inverse-positive stiffness/nonlinear map is
established for this pair. Positive values at ordinary Q2 nodes would also not
prove positivity of the entire enriched face/cell function. Restricting the
candidate space to positive functions merely states admissibility; it does not
prove a solution exists there, or remains uniformly separated from zero.

A precise sufficient criterion is a continuous solution u>=m>0 and
||u_h-u||_Linf(closed Omega)<=m/2, giving u_h>=m/2 everywhere.
For example, W^(1,p) convergence with p>3 implies this criterion by embedding
on this smooth three-dimensional domain. Neither estimate is established.
Strict positivity on each fixed compact mesh gives a mesh-dependent minimum;
uniform inverse-power bounds require a uniform positive lower bound or an
alternative weighted-source estimate. A uniform lower bound is sufficient,
not logically necessary if the strain vanishes suitably near small u_h.

## 2. Nonlinear coupling — CONDITIONAL

For a positive fixed-mesh candidate all integrals are finite. Uniform weak
solvability is a different question. H1_0 embeds in L6, so L^(6/5) embeds in
H^-1. If u_h>=m>0 and sigma_h is uniformly bounded in L^(12/5), then

    ||u_h^-7 |sigma_h|^2||_H^-1
      <= C m^-7 ||sigma_h||_(12/5)^2.

For two candidates u,v>=m and strains s,t, Holder gives the stability estimate

    ||u^-7 |s|^2-v^-7 |t|^2||_(6/5)
    <= m^-7 (||s||_(12/5)+||t||_(12/5)) ||s-t||_(12/5)
       +7 m^-8 ||u-v||_Linf ||t||_(12/5)^2.

Thus strong L^(12/5) strain convergence and uniform scalar convergence close
this source pairing. Alternatively, with the same lower bound, almost-everywhere
u_h convergence, strong L^(12/5) strain convergence and dominated convergence
for the fixed limiting strain suffice. These are sufficient conditions, not
claims that this Lebesgue exponent is necessary for every H^-1 source.
The accepted L2 strain estimate alone does not supply them; Section 7 already
exhibits the concentration obstruction, which is not repeated here.

The other terms must also be controlled. Uniform H1 bounds on u_h imply
u_h^5 bounded in L^(6/5); strong L6 convergence implies strong L^(6/5)
convergence of u_h^5, using

    ||u^5-v^5||_(6/5) <= C (||u||_6+||v||_6)^4 ||u-v||_6.

H1 weak convergence alone does not give compactness at this critical exponent.
For momentum, u^6 S is in L1 from H1 alone, again insufficient in general for
H^-1. A sufficient bound is u in L^(36/5), giving u^6 S in L^(6/5);
strong convergence in that norm controls the load difference. A uniform upper
bound and suitable scalar convergence would supply these power estimates.

For a bounded continuum scalar, smooth-domain Dirichlet elliptic regularity
for div L gives W in W^(2,p), for finite p when its load is in Lp; this would
supply much stronger continuum strain integrability. It does not automatically
give mesh-uniform L^(12/5) bounds for the actual enriched Galerkin strains.
Such a discrete estimate or an adequate error/inverse estimate remains missing.

Even the displayed source bounds do not prove coupled existence, uniqueness
or stability. The A u^5 term is critical in three dimensions and the momentum
load grows as u^6. No closed invariant interval, a priori upper bound or
invertibility bound for the coupled linearization has been established.
These are genuine missing hypotheses, not properties inferred from Korn.

## 3. Regularity — UNRESOLVED for the actual solutions

To interpret pointwise second derivatives, H2 in three dimensions is
insufficient: its second derivatives are merely L2. A sufficient continuum
condition is u,W in Hs with s>7/2, which embeds in C2; additional smoothness
is needed for a quantitative Taylor remainder. For the consistency bounds
below a sufficient assumption is C4 regularity up to both radial boundaries,
with bounded derivatives in the periodic parameter chart. Smooth exact
geometry transfers these bounds to Cartesian derivatives. Appropriate higher
Lp elliptic estimates and a scalar bound could bootstrap continuum regularity,
but their starting coupled estimates have not been established here.

The discrete functions remain C0 and piecewise smooth; global C2 regularity
is not inferred. What is needed instead is control of their *sampled* errors
under the actual audit operators, including stencils touching each boundary.
H3 regularity and the accepted O(h^2) energy approximation do not provide it.
No relation between finite-element h and audit n was previously fixed;
a theorem for u_h,W_h must state a choice h=h(n) and prove the estimates for
that sequence. No such choice or rate is asserted here.

## 4. Exact audit requirements and sufficient estimates — UNRESOLVED

Let delta=1/n, n=8,16,32,64. Radial and angular spacings are respectively
.20 delta and 2 pi delta. Let P_n sample continuous fields at the frozen nodes,
e_u=P_n u_h-P_n u and e_W=P_n W_h-P_n W. Use exactly the audit's Cartesian
first derivative G_n (including its three-point one-sided radial endpoints),
STF tensor T_n, composed divergence D_n T_n, and scalar face-flux Laplacian Q_n.
All arrays below include the boundary values needed by these stencils.
The residual arrays are precisely

    H_n(v,z)=-8 v^-5 Q_n v-v^-12 |T_n z|^2
             -16 pi G (A+B v^-4),
    M_n(v,z)=v^-10 D_n T_n z-8 pi G v^-4 S.

The contraction is full Cartesian Frobenius. For each band B use its actual
interior mask: [.05,.10] or [.20,.25], excluding the radial endpoints, and the
normalized RMS with weights J v^6. With m<=v,u<=M these norms are uniformly
equivalent to the normalized J-weighted array norm, denoted ||.||_(B,n).
Constants may depend on m,M, the fixed geometry and smooth solution bounds.
Assume also bounded sampled Q_n P_n u, D_n T_n P_n W, T_n P_n W, and
T_n P_n W_h. Subtraction and the Lipschitz bounds for inverse powers give

    ||H_n(P_n u_h,P_n W_h)-H_n(P_n u,P_n W)||_(B,n)
      <= C [||Q_n e_u||_(B,n)+||T_n e_W||_(B,n)+||e_u||_infinity],
    ||M_n(P_n u_h,P_n W_h)-M_n(P_n u,P_n W)||_(B,n)
      <= C [||D_n T_n e_W||_(B,n)+||e_u||_infinity].

The bounded-strain assumption can instead be replaced by appropriate discrete
Holder bounds on its quadratic product. It cannot simply be omitted.
These inequalities state explicitly the sampled error estimates required;
the accepted energy estimate bounds none of their second-difference terms.
A crude sufficient route, not a claimed rate, uses the stencil operator bounds
||Q_n e||_infinity+||D_n T_n e||_infinity<=C delta^-2||e||_infinity
and ||T_n e||_infinity<=C delta^-1||e||_infinity. Thus nodal errors
O(delta^(p+2)) would bound these residual perturbations by O(delta^p).
Sharper estimates could avoid this stringent sufficient condition.

**Consistency at both boundaries matters.** For C4 exact fields the flux
Hamiltonian has O(delta^2) pointwise consistency at all audited interior nodes,
including the first/last interior layers. The tensor's first derivatives have
O(delta^2) errors, but their composition is not uniformly second order there.
For a scalar radial test g, with radial spacing d=.20/n,

    (Dg)_0 = g'_0-d^2 g'''_0/3+O(d^3),
    (Dg)_i = g'_i+d^2 g'''_i/6+O(d^3)  (interior).

At i=1, D(Dg) uses (Dg)_0; the difference in error coefficients yields
+d g'''/4+O(d^2). The outer first-interior layer has the corresponding
signed term. Smooth geometry multipliers do not remove this generic O(d)
possibility in D_n T_n. Away from those layers the composed consistency is
O(delta^2). Each layer has normalized band volume O(delta), since each fixed
physical band has nonzero thickness. Consequently sufficient consistency
bounds for the exact solution are

    ||H_n(P_n u,P_n W)||_(B,n) <= C_H,B delta^2,
    ||M_n(P_n u,P_n W)||_(B,n) <= C_M,B delta^(3/2).

No nonzero leading coefficient is asserted: solution-dependent cancellation
may improve these orders. Both layers are included, separately in their own
band; dropping them changes the frozen audit. Endpoint and fourth-order
quantities remain additional diagnostics, not acceptance substitutes.

**Finite ratios are separate obligations.** For each of the four band/metric
combinations let R_n be the actual positive finite weighted RMS. The exact
necessary and sufficient numerical inequalities are

    0<R_(2n)<=2^(-3/2) R_n<infinity,  n=8,16,32.

Zeros require separate handling because the frozen logarithmic order must be
finite. Upper bounds tending to zero alone cannot prove these inequalities.
For example, a certified expansion R_n=c n^-p(1+eta_n), c>0,
|eta_n|<=epsilon_n<1, is sufficient if, for every pair,

    2^-p (1+epsilon_(2n))/(1-epsilon_n)<=2^(-3/2).

For p=2 a uniform epsilon <=(sqrt(2)-1)/(sqrt(2)+1) suffices.
At the generic boundary-layer p=3/2 there is no slack for a nonzero symmetric
relative-error bound; one needs signed remainder information, a certified
ratio directly, or a proved improved leading order. Scalar-dependent physical
weights and changing masks must be included in any such expansion. No
expansion, positive lower envelope or finite-n remainder constants for the
actual coupled sequence have been proved.

## 5. Feasibility and hard stop

| Requirement | Classification | Precise limit |
| --- | --- | --- |
| Positivity | CONDITIONAL | Continuum comparison and a uniform-error sufficient criterion proved; existence and whole-domain positivity of discrete coupled solutions unproved. |
| Nonlinear coupling | CONDITIONAL | Source pairing/stability inequalities proved under explicit lower-bound and stronger-norm hypotheses; those hypotheses and coupled solvability unproved. |
| Regularity | UNRESOLVED | No established bounded coupled solution or discrete sampled-derivative error theorem. |
| Sampled audit compatibility | UNRESOLVED | No actual-sequence estimates or twelve certified finite-n ratios; boundary composition has only generic order 3/2 in band RMS. |

No requirement is DISPROVEN for actual coupled solutions. What fails is the
inference from energy convergence to nonlinear H^-1 stability or strong sampled
residual convergence. The exact missing bridge consists of a positive bounded
coupled solution sequence, mesh-uniform source control beyond L2 strain,
sampled second-difference error estimates on the unchanged masks, and signed
finite-resolution ratio control. These are not supplied by the completed
linear mixed results. This bounded investigation stops here.

**One recommended next step:** authorize a proof-only attempt to establish a
mesh-uniform positive bounded branch of solutions for this same coupled
Galerkin system, including the stronger strain integrability needed above;
stop if that first prerequisite fails. No implementation or physical solves.

Validation of this report is limited to mathematical/document review, a
single-file diff, base-parent and frozen-audit blob checks, and remote branch
SHA verification. The commit SHA is reported externally to avoid self-reference.
