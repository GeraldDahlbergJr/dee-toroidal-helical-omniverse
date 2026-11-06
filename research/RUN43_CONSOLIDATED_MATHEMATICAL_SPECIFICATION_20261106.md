# Run 43 — consolidated mathematical specification

Specification checkpoint: `86b6cdf09605a8bc7e14cbfa9cfae4902c8e59ff`.
This document consolidates the existing exact-integration conforming mixed
Einstein–DEE initial-data formulation. It defines a mathematical method, not
an implemented solver or a certified positive coupled branch. All referenced
files are interpreted at this checkpoint. No model, space, quadrature, source,
or acceptance criterion is changed. No historical proofs are repeated.

## 1. Authoritative definitions and scope

The definition sources are:

| Reference | Definition supplied |
| --- | --- |
| [Covariant DEE formulation](../docs/RECOVERED_NONLINEAR_DEE_FORMULATION.md), “Covariant action”, “Euler-Lagrange equations”, “Stress-energy” | Original matter and Einstein equations |
| [Stress-energy source](../reproduce/dee_stress_energy_3p1.py), constants, `coeffs`, `potential`, `projections` | Kinetic coefficients, potential, projection signs |
| [Toroidal seed](../reproduce/validate_toroidal_helical_dee_source.py), constants, `metric`, `source` | Shell, helical gradients, normal velocities |
| [Coupled production formulation](../reproduce/solve_toroidal_helical_coupled_constraints.py), `G`, `Grid.__init__`, `energy`, `longitudinal`, `solve_vector`, `solve` | Cartesian conversion, conformal constraints and boundary data |
| [Conforming mixed specification](RUN43_CONFORMING_MIXED_PROOF_20261105.md), §§1–7 | Exact geometry, mesh class, enriched spaces, pairings, weak momentum and scalar candidate |
| [Invariant-set report](RUN43_DECISIVE_INVARIANT_SET_20261106.md), §§1–2 | Actual scalar map, whole-function interval, discrete Green kernel |
| [Analytical reduction](RUN43_TARGETED_ANALYTICAL_COUNTEREXAMPLE_20261106.md), “a coupled large-amplitude test” | Fixed-mesh limiting integral |

The production code supplies the continuum target and physical data; its
finite-difference matrices are **not** the conforming Galerkin operators
defined below. The latter originate in the conforming mixed specification.
The [coupled-estimates report](RUN43_COUPLED_ESTIMATES_PROOF_20261106.md),
[positive-branch report](RUN43_POSITIVE_COUPLED_BRANCH_PROOF_20261106.md), and
[limiting-integral certification](RUN43_LIMITING_INTEGRAL_CERTIFICATION_20261106.md)
supply the result classifications and limitations recorded in §8.

## 2. Exact geometry, calculus and admissible meshes

Let a=1/20, b=1/4, R=1, and let both angles have period 2π. Set

    X(r,θ,φ)=((1+r cosθ)cosφ, (1+r cosθ)sinφ, r sinθ),
    Ω=X((a,b)×(ℝ/2πℤ)²),   q=1+r cosθ,   J=rq.

The orthonormal frame and Euclidean background metric are

    e_r=(cosθ cosφ,cosθ sinφ,sinθ),
    e_θ=(-sinθ cosφ,-sinθ sinφ,cosθ),
    e_φ=(-sinφ,cosφ,0),
    barγ=dr²+r²dθ²+q²dφ².

For a scalar u,

    ∇u=e_r u_r+e_θ u_θ/r+e_φ u_φ/q,
    Δu=u_rr+(1/r+cosθ/q)u_r+u_θθ/r²
         -sinθ u_θ/(rq)+u_φφ/q²,
    dV=J dr dθ dφ.

The outward normals are n_a=-e_r and n_b=e_r; on r=c,
dS=cq(c,θ)dθdφ. The geometry is exact and smooth, with r≥1/20 and
q≥3/4. Angular chart seams are identified, not physical boundaries.
All vector and tensor components below are Cartesian unless stated otherwise.

An admissible mesh is a conforming tensor-product partition of the parameter
shell into bricks, with matching periodic seam faces and no hanging faces.
Each brick is A_K((0,1)³), with affine A_K, and the physical cell map is
F_K=X∘A_K. Parameter side lengths within each brick are uniformly comparable,
with bounded aspect ratios; the maximal parameter side length h tends to zero
for a refining family. Global quasi-uniformity is not imposed. These are exact
curved cells, without faceting or polynomial geometry replacement. Put
B_K=DF_K. Volume pullbacks use |det B_K|=J|det DA_K|; face pullbacks use
|cof(B_K)n_ref| and the corresponding oriented physical normal.

## 3. Original Einstein–DEE equations and frozen physical data

Use signature (-,+,+,+), c=1, and fields Φ=(ρ,Θ,Ψ). With spacetime contractions,

    L_m=-½(∇ρ)²-V(ρ)-½K_Θ(ρ)(∇Θ)²-½K_Ψ(ρ)(∇Ψ)²
         -Λ(ρ)∇_μΘ∇^μΨ,
    S_action=∫sqrt(-g)[R_scalar/(16πG)+L_m] d⁴x,
    G_μν=8πG T_μν,
    T_μν=∂_μρ∂_νρ+K_Θ∂_μΘ∂_νΘ+K_Ψ∂_μΨ∂_νΨ
         +Λ(∂_μΘ∂_νΨ+∂_μΨ∂_νΘ)+g_μν L_m.

The original matter equations are

    □ρ-V'(ρ)-½K_Θ'(ρ)(∇Θ)²-½K_Ψ'(ρ)(∇Ψ)²
       -Λ'(ρ)∇_μΘ∇^μΨ=0,
    ∇_μ(K_Θ∇^μΘ+Λ∇^μΨ)=0,
    ∇_μ(K_Ψ∇^μΨ+Λ∇^μΘ)=0.

The realization is K_Θ=1+.30ρ², K_Ψ=1+.20ρ², Λ=.10ρ²,
V=(4/4)(ρ²-1)², G=.001. This initial-data benchmark freezes the matter
seed; it does not solve these evolution equations. Its seed is

    ρ=sqrt(.945), Π=(0,.12,-.04), Π_A=n^μ∂_μΦ_A,
    Dρ=0, DΘ=.25 dχ, DΨ=.20 dχ, dχ=dθ-3dφ.

Only these phase gradients are needed; no extra global phase-value convention
is imposed. Define the Cartesian covector direction, identified using barγ,

    D=e_θ/r-3e_φ/q,
    K_kin=diag(1,K_phase),
    K_phase=[[1.2835,.0945],[.0945,1.189]].

The normal velocities and covector gradients remain fixed as the metric
changes. For γ=ψ⁴barγ, the original physical projections are

    E=½ΠᵀK_kinΠ+½K_kin,AB γ^ij D_iΦ_A D_jΦ_B+V,
    S_i=-K_kin,AB Π_A D_iΦ_B,
    S_ij=K_kin,AB D_iΦ_A D_jΦ_B
         +γ_ij[½ΠᵀK_kinΠ-½K_kin,AB γ^kl D_kΦ_A D_lΦ_B-V].

The spatial stress trace is γ^ij S_ij. Distinguish it from the vector source
S below, which is the fixed momentum covector raised with **barγ**, not γ.
With P=(.12,-.04), d=(.25,.20), the actual coefficients are

    A=½PᵀK_phase P+(.945-1)²=.0127638,
    B(x)=½dᵀK_phase d |D|²=.068614375(r^-2+9q^-2),
    S(x)=-(PᵀK_phase d)D=-.030316D,
    E(ψ)=A+Bψ^-4.

The terminating decimals here denote the exact parameter algebra of the
specification; executable floating-point evaluation has its ordinary rounding.
These values agree with the original coefficient, seed and Cartesian source
definitions above. A and B are nonnegative and bounded on this fixed shell.

For maximal conformally flat data, tr K=0 and A_phys^ij=ψ^-10(LW)^ij, where

    (LW)_ij=∂_jW_i+∂_iW_j-(2/3)δ_ij∂_kW_k,
    (div T)_i=∂_jT_ij.

The target constraint equations are exactly

    div LW=ΔW+(1/3)∇div W=8πG ψ⁶ S,
    -Δψ=ψ^-7|LW|²/8+2πG(Aψ⁵+Bψ),
    ψ=1, W=0 at r=a,b; periodic Cartesian traces at angular seams.

These are finite-shell Dirichlet data, without exterior matching or a traction
condition. The contraction |T|²=Σ_ij T_ij² is the full Frobenius contraction;
off-diagonal Cartesian entries count twice.

## 4. Exact conforming scalar, displacement and strain spaces

On the reference cube let Q2 contain polynomials of degree ≤2 in each
variable; P_m denotes Cartesian polynomials of total degree ≤m. Define

    U(K)={qhat∘F_K^-1 : qhat∈Q2}+P_2(R³)|K,
    U_h={u∈C⁰(closure Ω): u|K∈U(K), periodic seam traces},
    U_h,0={u∈U_h: u=0 at r=a,b},
    V_h=U_h³,   V_h,0=U_h,0³.

The local scalar pullback is Q2+span{1,F_i,F_iF_j:1≤i≤j≤3}.
Redundant generators are quotiented as functions. Assembly imposes equality
of entire face trace functions, including periodic faces, by coefficient
constraints in their combined trace spans. Zero radial trace is imposed in
the same way. Ordinary 27-node Q2 agreement alone does not define this enriched
space. The scalar candidate is ψ_h∈1+U_h,0, positive throughout closure Ω.

Let E_1=diag(1,-1,0), E_2=diag(1,1,-2), and let E_3,E_4,E_5 have respectively
xy, xz, yz symmetric pairs equal to one and all other entries zero. These
span the Cartesian symmetric trace-free tensors STF(3). For Cartesian unit e_i,

    g(qhat)=B_K^-T∇_ξ qhat,
    D_i(qhat)=e_i⊗g(qhat)+g(qhat)⊗e_i
              -(2/3)(e_i·g(qhat))I,
    That(K)=span{D_i(ξ^α):i=1,2,3; 0≤α_j≤2}
             +span{E_β p(F_K):β=1,…,5; p∈P_1(R³)},
    T(K)={That∘F_K^-1:That∈That(K)},
    T_h={T∈L²(Ω;STF(3)):T|K∈T(K)}.

This is a componentwise Cartesian composition with geometry-dependent
gradient generators, not a Piola map or a polynomial tensor pullback.
There are no tensor trace constraints between cells: T_h is L² conforming,
not H(div) conforming. The existing containment result is L(V_h)⊂T_h.

## 5. Integration, mixed equations and Galerkin operators

All pairings are exact Lebesgue integrals with Euclidean dV on the exact
geometry. Tensor pairings use the full Frobenius product, or its actual Gram
matrix in an independent STF basis. No finite Gaussian rule is specified as
exact for the nonpolynomial geometry and loads. Any future numerical
quadrature requires separate justification of positivity, norm equivalence
and variational error; it is not part of this specification.

For a prescribed positive scalar input v define W_h(v) and σ_h(v) by

    (σ_h,τ_h)=(LW_h,τ_h)                         ∀τ_h∈T_h,
    -(σ_h,Lz_h)/2=(8πG v⁶S,z_h)                 ∀z_h∈V_h,0,
    (σ_h,W_h)∈T_h×V_h,0.

Thus σ_h=LW_h exactly. Write a(u,z)=(Lu,Lz)/2. Define the linear vector
operator K_h acting on scalar multipliers by

    a(K_h c,z_h)=-(8πG cS,z_h)                  ∀z_h∈V_h,0.

Then W_h(v)=K_h(v⁶). The minus sign and all physical constants and S are
included in K_h. In an independent vector basis z_i its stiffness is
C_ij=a(z_j,z_i), and its load is -8πG∫c S·z_i dV.
Equivalently, with tensor mass H_T and the homogeneous strain coefficient
map L_0, C=L_0ᵀH_T L_0/2. The variational div L operator has matrix
M A_D=-C if represented through the exact displacement mass M.

Define the scalar stiffness and inverse Poisson operator by

    b(p,q)=(∇p,∇q),
    b(R_h f,q_h)=(f,q_h)                       ∀q_h∈U_h,0,
    R_h f∈U_h,0.

For an independent real basis φ_i of U_h,0, H_ij=b(φ_j,φ_i),

    g_h(x,y)=φ(x)ᵀH^-1φ(y),
    R_h f(x)=∫_Ω g_h(x,y)f(y)dV_y.

This finite-dimensional kernel is basis independent and continuous, and each
fixed-mesh slice is bounded. Neither its nonnegativity nor inverse positivity
of R_h is assumed. Loads here are integrable functions on a fixed mesh;
the variational inverse also applies to H^-1 loads by dual pairing.

The weak integration identity fixing the sign is

    (div T,z)=-(T,Lz)/2+⟨Tn,z⟩_∂Ω.

For discontinuous T_h it is interpreted cellwise with internal jump terms,
not as the broken classical divergence alone. Dirichlet reactions remain:
for a trace g and a lift z_h∈V_h,

    reaction_h(g)=(8πG v⁶S,z_h)+(σ_h,Lz_h)/2.

The interior equation makes this independent of the chosen discrete lift.
No raw discontinuous cell traction is declared an equilibrated reaction.

## 6. Actual coupled scalar source and admissible interval

Set a_0(x)=2πG A and b_0(x)=2πG B(x), distinct from the shell radii. Define

    F_h(v)=v^-7|LK_h(v⁶)|²/8+a_0v⁵+b_0v,
    𝒯_h(v)=1+R_h F_h(v),
    I_h(m,M)={v∈1+U_h,0:m≤v(x)≤M for every x∈closure Ω},
    0<m≤1≤M<∞.

The scalar map 𝒯_h is distinct from the strain space T_h.
A coupled solution is a fixed point ψ_h=𝒯_h(ψ_h), together with the momentum
solution above. Equivalently b(ψ_h-1,q_h)=(F_h(ψ_h),q_h) for every scalar test.
No projection of the squared strain, unweighted five-component square, or
source substitution is permitted. Whole-function bounds, rather than nodal
bounds, define admissibility. All integrals are finite on a fixed mesh for
such inputs. Membership does not assert existence of a fixed point.

The exact interval-invariance requirement is

    m-1≤∫g_h(x,y)F_h(v)(y)dV_y≤M-1

for every x and v∈I_h(m,M). For m=1 the lower condition is nonnegativity of
this actual-family integral. Positivity of the load alone, a negative kernel
value alone, or a sign of a coefficient matrix inverse alone does not decide
this restricted-family condition.

## 7. Limiting integral: definition and assumptions

Fix one admissible mesh, one interior cell K=F_K((0,1)³), and the mapped Q2
bubble

    η|K=[∏_(i=1)³4ξ_i(1-ξ_i)]∘F_K^-1,   η=0 off K.

It has zero entire face traces, belongs to U_h,0, satisfies 0≤η≤1 and
max η=1. This is the (0,1)³ version of the bubble in the certification report.
Let v_t=1+tη for finite t≥0; then v_t∈I_h(1,1+t). Define

    Z_j=K_h(η^j), τ_j=LZ_j, j=0,…,6, with η⁰=1,
    W_h(v_t)=Σ_(j=0)⁶ binom(6,j)t^jZ_j,
    σ_h(v_t)=Σ_(j=0)⁶ binom(6,j)t^jτ_j,
    E_0={y∈Ω:η(y)=0}.

The existing fixed-mesh limit is

    t^-12 F_h(v_t) → 1_(E_0)|τ_6|²/8 in L¹(Ω),
    t^-12 R_h F_h(v_t)(x) → C_h(x),
    C_h(x)=(1/8)∫_(E_0)g_h(x,y)|LK_h(η⁶)(y)|²dV_y.

Assumptions are exactly the fixed finite-dimensional spaces and operators,
exact integration, the original bounded a_0,b_0,S, the stated nonnegative
bubble, and a fixed physical evaluation point x∈closure Ω. The limit uses
linearity of K_h, L² convergence of the scaled strain polynomial, L¹
convergence of its square, dominated convergence of (1+tη)^-7, and boundedness
of the fixed-mesh kernel slice. It is not a mesh-uniform assertion and
requires no hypothesized sign of g_h.

No particular cell coordinates, evaluation point, or sign certificate has
been selected in this consolidation. If C_h(x)<0 were separately certified,
the response would be negative for sufficiently large finite t. This would
refute lower-endpoint-1 comparison only for intervals containing that input.
It would not refute every m<1 interval or existence of positive solutions.
C_h(x)≥0 would not settle finite-t behavior. For an upper endpoint M fixed
in advance, t≤M-1 is required; the limit alone supplies no such violation.

## 8. Established results, conditional results and outstanding requirements

| Classification | Existing result and its scope |
| --- | --- |
| Proven under exact integration and the specified mesh class | H¹ conformity, exact strain containment, variational adjoint with reactions, Cartesian quadratic reproduction; conforming report §§3–6. |
| Proven linear stability | On V_h,0, ‖Lz‖²₂=2‖∇z‖²₂+(2/3)‖div z‖²₂; fixed-domain Poincare gives a mesh-independent H¹ Korn bound and zero Dirichlet kernel. Fixed-load momentum has a unique solution and stable elimination. |
| Proven approximation with hypotheses | For H³ zero-trace periodic displacement, best H¹/strain approximation is O(h²); scalar interpolation L² error is O(h³). Fixed-load momentum is quasi-optimal; density gives convergence without this regularity. These are not nonlinear coupled rates. |
| Proven actual-source energy bound | For v∈I_h(m,M), ‖σ_h(v)‖₂≤C M⁶ uniformly in h. This controls the scalar load in L¹, not uniformly in H^-1 or L∞ after Poisson inversion; invariant-set report §3. |
| Proven fixed-mesh statements | Continuity of the frozen scalar map on the interval, the Green representation, and the limit in §7. No kernel-integral sign follows. |
| Conditional continuum positivity | An already positive admissible continuum weak solution with nonnegative H^-1 scalar source satisfies ψ≥1. This assumes existence and valid testing; it is not a discrete maximum principle. |
| Conditional nonlinear passage | A uniform lower scalar bound, uniform scalar convergence and strong L^(12/5) strain convergence suffice for H^-1 convergence of the inverse-power squared-strain term. Scalar powers and the ψ⁶ momentum load require their own stronger bounds/convergence; coupled-estimates report §§1–2. |
| Conditional bounded branch | Suitable mesh-uniform stronger strain and scalar inverse estimates, lower comparison and quantitative upper closure would permit an invariant-set argument. Those assumptions are not established; positive-branch and invariant-set reports. |
| Unresolved | Actual coupled existence, uniqueness/continuation, whole-function discrete positivity, mesh-independent upper bound and stronger strain integrability; source positivity and linear Korn stability do not supply them. |
| Unresolved | Sign of the limiting integral on a specified mesh, bubble and point, and any finite admissible violating input. The certification report is inconclusive; no actual-family counterexample is established. |
| Unresolved | Pointwise boundary traction, sampled second-derivative errors, and compatibility with the unchanged strong-residual audit. Weak energy convergence does not imply these estimates. |

The frozen [array audit](audit_run43_arrays.py) requires n=8,16,32,64,
ψ shape (n+1,n,n), Cartesian W shape (3,n+1,n,n), finite arrays, positive
ψ, and exact radial values ψ=1,W=0. Its
[diagnostic source](../reproduce/diagnose_toroidal_boundary_convergence.py)
has checkpoint Git blob `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`.
It uses the original face-flux Hamiltonian and second-order tensor-divergence
momentum residuals, interior nodes in fixed radial bands [.05,.10] and
[.20,.25], and normalized volume-weighted RMS with weights Jψ⁶.
All twelve logarithmic orders (two bands, two metrics, three refinement pairs)
must be finite and ≥1.5. Supplementary endpoint/fourth-order diagnostics do
not replace these criteria. Sampling continuous conforming fields is defined,
but the audit reconstructs its own strain from samples; it need not equal σ_h.
No coupled sequence, relation h=h(n), sampled-error theorem, or finite-pair
ratio certification is provided by the conforming formulation.

## 9. Reproducibility boundary and stop

To instantiate the method one must supply an explicit admissible brick
partition, construct the quotient/trace-constrained spaces, and use the exact
operators and sources above. A kernel-sign certificate would additionally
need a particular physical bubble cell and evaluation point. These choices
are inputs still required for a particular instance, not missing definitions
of the mesh class or operators and not computations performed here.

This change adds only this specification. Verification is limited to source
comparison, document/link review, diff checks, checkpoint parent and remote
branch SHA verification. No simulations, physical solves, numerical
experiments, GMRES, sweeps or GitHub Actions are run. Publication completes
this task. **STOP; no further positivity investigation is authorized.**
