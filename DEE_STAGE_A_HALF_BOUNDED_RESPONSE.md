# Stage A½: bounded amplitude response

Date: 2026-11-10 UTC. Starting commit: `f2da9c5fa8d3f34275afc8f0651366adbf8e65d5`.
Branch: `research/stage-a-half-bounded-response`.

## Inputs and scope

Only the four requested baseline evidence files and the original covariant matter action were read as scientific inputs. The action source is `evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex`, Eqs. (3)–(9), (14)–(15). It has M(rho)=I+rho²A, A=[[a,c],[c,b]], and V=lambda_V(rho²-v²)²/4. No action, prior evidence or manuscript is changed. Previously verified metric, background relation and linearized equations are inputs, not recalculated Einstein equations.

Use eta=diag(-1,1,1,1), Q=kz+Q_0+f, rho=rho_0+u, with static transverse u. Let m=sqrt(2lambda_V rho_0²)>0, C=k^T A k, B=k^T M_0 k, epsilon=V(rho_0)+B/2, and F=2pi G rho_0 C(B+2epsilon). All constants are real, rho_0,lambda_V>0; the baseline branch rho_0²=v²-C/lambda_V>0 is retained. Physical admissibility additionally requires v²>0 and positive-definite M_0. F may have either sign or vanish. Work is strictly first order in the prescribed weak metric, with f=0 unless explicitly stated.

The transverse domains are the centered disk 0<=r<R or concentric annulus 0<r_i<r<r_o, with finite positive radii. This is a transverse boundary problem, with no axial endcaps. Radial formulas require constant data on each boundary circle; nonradial smooth data have the elliptic existence result but generally nonradial solutions.

## Exact solution

For L=Delta_perp-m², the equation is Lu=Fr². Substitution of p=alpha r²+beta yields

    Lp = -m² alpha r² +4alpha-m² beta,
    alpha=-F/m², beta=-4F/m⁴,
    p(r)=-F(r²/m²+4/m⁴).

For r>0, the complete radial solution is

    u=p+A_I I_0(mr)+A_K K_0(mr).

The modified Bessel equation verifies each homogeneous term. I_0(mr) is analytic in r² at the origin; K_0(mr) has a logarithmic singularity. Disk regularity therefore forces A_K=0, with u'(0)=0. Writing I(r)=I_0(mr), K(r)=K_0(mr), a disk Dirichlet value u(R)=g gives

    u_D(r)=p(r)+(g-p(R)) I(r)/I(R).

I(R)>0, so this is well defined. For disk homogeneous Neumann data u'(R)=0,

    u_N(r)=p(r)+[2FR/(m³ I_1(mR))] I(r).

I_1(mR)>0. These are distinct boundary problems, not simultaneous conditions.

On the annulus both I and K are smooth and bounded, including both boundary circles. Put d_j=g_j-p(r_j) and

    D=I(r_i)K(r_o)-I(r_o)K(r_i),
    A_I=[d_i K(r_o)-d_o K(r_i)]/D,
    A_K=[I(r_i)d_o-I(r_o)d_i]/D.

Then u(r_i)=g_i and u(r_o)=g_o. D<0: I/K is strictly increasing, since I'=mI_1>0, K'=-mK_1<0 and I,K>0. Thus no exceptional positive radii make this system singular.

Homogeneous annular Neumann conditions give

    N (A_I,A_K)^T = (2F r_i/m²,2F r_o/m²)^T,
    N=[[m I_1(mr_i),-m K_1(mr_i)],
       [m I_1(mr_o),-m K_1(mr_o)]].

The outward normal is -partial_r at the inner circle; its zero condition is still u'(r_i)=0. Coefficients are exactly N^{-1} times the stated vector. det N is nonzero by the homogeneous energy identity below; equivalently I_1/K_1 is strictly increasing. General normal data q_i,q_o instead require radial derivatives -q_i,q_o. Mixed Dirichlet/Neumann data on the two circles are also coercive and unique.

## Existence, uniqueness and regularity

Rewrite as (-Delta+m²)u=-Fr². The bilinear form integral(grad u dot grad v +m²uv) is coercive on H_0^1 for Dirichlet and on H^1 for Neumann, since m²>0. Smooth Dirichlet traces admit a lifting. Lax–Milgram gives unique weak solutions in both cases; elliptic boundary regularity makes them smooth on the closed smooth domain for smooth data and forcing. For finite C^{2,alpha} traces one obtains C^{2,alpha} regularity. Uniqueness follows by integrating the homogeneous difference against itself: integral(|grad w|²+m²w²)=0, with the boundary term zero. This also proves the annular boundary matrices invertible. Radial symmetry follows from uniqueness and rotational invariance when the boundary data are radial. The disk Bessel solution is smooth in Cartesian coordinates at r=0.

Unlike the massless Neumann problem, no zero-mean forcing condition applies. Integration instead gives

    -m² integral_domain u = F integral_domain r²

for zero Neumann data. Thus the disk mean is -FR²/(2m²) and the annular mean is -F(r_i²+r_o²)/(2m²). Nonzero F requires a nonzero response mean; it does not prevent solvability.

## Bounds and perturbative meaning

Let R be the outer radius. Comparison with constant sub/supersolutions gives

    ||u_D||_infinity <= max(||g||_infinity, |F|R²/m²),
    ||u_N||_infinity <= |F|R²/m²

for homogeneous Neumann data. The Neumann bound also follows by testing the positive parts above and below these constants. Thus zero Dirichlet or homogeneous Neumann solutions are bounded for every fixed m>0 and finite domain. For F>0 zero-Dirichlet u is strictly negative in the interior, with strictly positive outward derivative by the strong maximum principle and Hopf lemma; signs reverse for F<0.

Bounded does not automatically mean perturbatively small. Sufficient conditions include |F|R²/(m²rho_0)<<1 and, for nonzero Dirichlet data, ||g||/rho_0<<1, together with pi|G|(|B|+|epsilon|)R²<<1 for all supplied metric components. Positive M must remain positive: a sufficient condition is ||A|| |2rho_0 u+u²| < lambda_min(M_0) uniformly. For fixed domain, m>0 and couplings, with data O(G), elliptic estimates give u and its derivatives O(G); constants depend on geometry and m. No uniform statement is made as m tends to zero, radii vary without control, or high-gradient boundary data are introduced. The displayed amplitude bound is sufficient, not necessary. Algebraic 1/m powers in p can cancel in a boundary solution; they alone do not diagnose a singular disk Dirichlet massless limit. For zero Neumann data the mean explicitly exhibits the m->0 obstruction when F stays nonzero.

## Both phase currents

The retained first-order current density and equations are

    delta(sqrt(-g) J^mu)=M_0 partial^mu f
      +2rho_0 u A k delta_z^mu
      +(h/2 eta^{mu z}-h^{mu z})M_0 k,
    M_0 Box f= M_0 k D^z-2rho_0 A k partial_z u,
    D^z=partial_mu h^{mu z}-(1/2)partial_z h.

The supplied harmonic metric has h_zz=-2pi G(B+2epsilon)r², h=4pi G(B-4epsilon)r² and h^{tz}=h^{xz}=h^{yz}=0. Thus D^z=0 directly. For every response above partial_z u=0. Each of the two component equations therefore has zero residual for f=0, including nonzero transverse amplitude gradients. The amplitude-dependent current points in z, so transverse variation of its z component produces no divergence. Normal phase flux through the circular boundaries is zero at this order.

Explicitly, additional phases, if independently required by other data, must satisfy

    (1+a rho_0²)Box f_Theta+c rho_0² Box f_Psi
      =(M_0 k)_Theta D^z-2rho_0(a k_Theta+c k_Psi)partial_z u,
    c rho_0² Box f_Theta+(1+b rho_0²)Box f_Psi
      =(M_0 k)_Psi D^z-2rho_0(c k_Theta+b k_Psi)partial_z u.

For positive M_0 this is Box f=k D^z-2rho_0 M_0^{-1}Ak partial_z u. No additional phases are required by this transverse harmonic problem. Nonzero z-dependent phases would also change the amplitude equation through 2rho_0 k^T A partial_z f and cannot be added while silently retaining the fixed-phase forcing. Phase uniqueness or compatibility with unspecified exterior/axial data is not asserted.

## Classification and reproduction

**VERIFIED:** polynomial coefficients, disk/annular Bessel responses, bounded-domain existence/uniqueness and regularity, admissibility of homogeneous Neumann conditions, conditional perturbative bounds, and both fixed-phase current residuals. These are analytical theorems with exact symbolic identity checks, not numerical simulation or external review.

**CORRECTION REQUIRED:** interpreting boundedness as automatic weak-field validity, or a Dirichlet trace match as complete source-free matter matching. See the companion boundary analysis.

**NOT YET ESTABLISHED:** a compatible nontrivial exterior matter/metric continuation, global weak-field asymptotics, nonlinear stability, or full coupled Einstein–DEE existence. Historical status and Eq. (34) correction remain unchanged.

Run `python DEE_STAGE_A_HALF_BOUNDED_VERIFIER.py` (SymPy 1.14.0). All checked residuals must be exactly zero. Positivity, coercivity, Hopf and matching arguments are analytical proofs above and in the companion document, not claims proved by sampling. No baseline verifier rerun, numerical PDE solver, Run 43, parameter fitting or nonlinear solver is used.

Observed verification on 2026-11-10 UTC: Python verifier exited 0 with SymPy 1.14.0. Particular coefficients were exactly alpha=-F/m² and beta=-4F/m⁴; all printed equality residuals were zero, both phase components separately passed, and K_0's axis limit was +infinity. Final output: `PASS: exact bounded response, boundary and both phase-current identities.`
