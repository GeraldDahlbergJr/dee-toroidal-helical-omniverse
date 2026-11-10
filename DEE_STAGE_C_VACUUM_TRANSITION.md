# Stage C: vacuum-transition feasibility

Date: 2026-11-10 UTC. Base: `0a9c4f84c383089d3f74b49fe2a4d5fb05c75856`.
Branch: `research/stage-c-transition-feasibility`. Single-agent analytical work.

**Result:** a finite, smooth transition to vacuum is impossible while retaining
nonzero constant axial phase gradients and a strictly positive kinetic matrix.
Amplitude relaxation alone cannot achieve vacuum. This is an exact obstruction
to the current ansatz, including any regular finite-width nonlinear radial
collar. It is not a no-go theorem for the entire Einstein–DEE theory.
The companion no-go analysis gives hypotheses and proofs.

## Authority and variation of the unchanged action

Source: `evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex`,
Eqs. (3)–(9), (14)–(15); blob
`a20c0ad7b3469487e50fa2be3a193d1aacc92ab6`. Read alongside the prior matter
derivation/existence analysis and both Stage B reports. Their evidence and
scope are preserved. The recovered source, rather than an unrecovered later
manuscript or historical source envelope, is the authority here.
Signature (-,+,+,+), G>0, real fields, smooth Lorentz metric, and C² or smoother
couplings/fields are assumed. Set Q=(Theta,Psi)^T and

    S = integral sqrt(-g) [R/(16pi G) + L],
    L = -1/2 (d rho)^2 - V(rho) - 1/2 d_mu Q^T M(rho) d^mu Q,
    M = [[K_Theta,Lambda],[Lambda,K_Psi]].

Vary before imposing symmetry. The rho terms in delta L are
`-d^mu rho d_mu delta rho - [V' + (dQ)^T M' (dQ)/2] delta rho`;
the phase terms are `-(M d^mu Q)^T d_mu delta Q`. Integration by parts gives

    E_rho = Box_g rho - V' - 1/2 d_mu Q^T M' d^mu Q = 0,
    E_Q = nabla_mu J^mu = 0,             J^mu = M d^mu Q.

The boundary term is `-n_mu[d^mu rho delta rho + J^mu^T delta Q]`.
Symmetric inverse-metric variation, including delta sqrt(-g), gives

    T_mu nu = d_mu rho d_nu rho + d_mu Q^T M d_nu Q + g_mu nu L,
    G_mu nu = 8pi G T_mu nu.

The cross term in T is Lambda times the sum of both mixed phase products.
These are the full equations, with no weak-field truncation. Off shell,

    nabla_mu T^mu_nu = E_rho d_nu rho + E_Q^T d_nu Q.

The verifier derives the radial equations by a reduced variational calculation
and independently checks conservation by differentiating the stress. The
covariant derivation above fixes the equations before that reduction.

## Vacuum means zero stress, not zero amplitude or zero contraction

On a connected open region with M positive definite, T=0 implies d rho=0
and dQ=0. One proof uses the two null directions e_0 +/- e_i in an orthonormal
frame: each null stress is a sum of positive kinetic squares. Their sum forces
both temporal and each spatial derivative to vanish. Then T=-V g forces V=0;
the rho equation also requires V'=0. Conversely these conditions solve all
matter equations and give zero stress. No V>=0 assumption is needed for this
null-contraction proof. If V>=0, the independent energy-density proof is

    T_hat0hat0 = V + 1/2 sum_{a=0}^3 [(d_hat a rho)^2
                                        + d_hat a Q^T M d_hat a Q].

For the original polynomial realization

    M=I+rho² A, A=[[a,c],[c,b]], V=lambda_V(rho²-v²)²/4,

with lambda_V,v>0, the vacua have rho=+v (also -v if the signed scalar is
allowed), constant phases, V=V'=0. If rho is restricted to positive separation,
only +v is in that field domain. At rho=0, V=lambda_V v^4/4>0 and M=I;
zero amplitude is not vacuum. At rho=v strict kinetic positivity means
`1+a v²>0` and `det M(v)>0`. The established axial family uses A positive
definite, hence M>=I at every rho and is safely within this sector. The source
with unspecified functions has a vacuum only if a permitted rho_v satisfies
V=V'=0 and the desired kinetic conditions; its existence is not automatic.
A common null phase wave has (dQ)^T M dQ=0 but nonzero T and is not vacuum.

If M loses positivity, the above theorem must not be applied. Semidefinite
null directions can hide phase gradients from stress, but M' still enters the
rho equation and the principal symbol becomes degenerate. Indefinite M permits
ghosts/cancellations. Neither is a regular strictly-positive-kinetic exterior.
Nothing here changes the functions to exploit those possibilities.

## Exact radial system and the failed transition

For explicit gravitational and current formulas use the static diagonal family

    ds²=-N(r)²dt²+H(r)²dr²+L(r)²dphi²+Z(r)²dz²,
    rho=rho(r), Q=w t+m phi+k z+q(r), S=N H L Z.

N,H,L,Z are positive, finite and smooth on a closed annular transition;
r=0 is treated by regular-axis limits. Constants are two-vectors; m must
have integer components for circle-valued phases. The full equations become

    (S rho'/H²)'/S - V' - 1/2[-w^T M' w/N²
      +q'^T M' q'/H²+m^T M' m/L²+k^T M' k/Z²] = 0,
    (S M q'/H²)'=0,     S M q'/H²=j (constant two-vector).

An exterior with constant phases has j=0. With nonsingular M this forces
q'=0 on the whole connected radial region. A regular axis also sets j=0
if normal currents are bounded (S/H is proportional to r there).
Tangential current components are -Mw/N², Mm/L², Mk/Z². Their divergence
vanishes by symmetry; conservation does not authorize arbitrary k(r).
The established configuration has w=m=0, k!=0 and
rho_0²=v²-k^T A k/lambda_V. Neither varying rho nor the metric at finite
regular radius makes d_z Q=k vanish. Its axial frame energy alone is
at least lambda_min(M)|k|²/(2Z²)>0. Thus no finite vacuum endpoint exists.
In a nondegenerate asymptotic cylindrical frame with bounded Z and a uniform
positive lower bound on M, the same lower bound prevents even T->0.
If Z diverges or eigenvalues tend to zero at infinity this bound is not
uniform: no universal asymptotic theorem is asserted. Such geometry is not
a finite smooth vacuum join; it requires separate Einstein/asymptotic analysis.

Smooth radial tapering Q=k(r)z gives d_r Q=k'(r)z. Its positive kinetic
energy is z-dependent and unbounded in z wherever k'!=0. More fundamentally,
invariance of dQ and d(dQ)=0 force k'=0. A taper therefore leaves the current
ansatz and requires solving both phase PDEs anew. Stationary/translation-
invariant stress alone is weaker than invariant phase gradients: the latter
is an explicit hypothesis here. Other symmetry realizations are not excluded.

## Necessary constraints, regularity and interpretation

At a timelike interface without a shell, require common induced metric gamma
and extrinsic curvature K, using a common normal. Gauss and Codazzi then
require common T_nn and T_na. For a spacelike normal and convention
K_ab=gamma_a^mu gamma_b^nu nabla_mu n_nu,

    G_nn = (K²-K_ab K^ab - R[gamma])/2 = 8pi G T_nn.

In the radial diagonal family R[gamma]=0 locally, and with
alpha=N'/(HN), beta=L'/(HL), chi=Z'/(HZ),

    G_nn=alpha beta+alpha chi+beta chi,
    T_nn=(rho'^2+q'^T M q')/(2H²)-V
         + w^T M w/(2N²)-m^T M m/(2L²)-k^T M k/(2Z²).

The verifier independently calculates G_nn from the Ricci tensor and from K.
T_na=d_n rho d_a rho+(d_n Q)^T M d_a Q. Vacuum requires zero normal
traction; this necessary condition alone is not sufficient for T=0.
In the static axial zero-current sector with V>=0, any point with rho'=0
has negative T_nn if k!=0. A transition could change rho', but cannot evade
the phase-trace obstruction. The prior perturbative collar is not an exact
Einstein interior, and its O(G) response cannot establish a finite vacuum
transition or cancel its zeroth-order stress.

No scalar delta source requires continuity of rho and its normal derivative,
and of phase traces and normal currents; smooth matching is stronger still.
Nonzero winding with M uniformly positive has divergent axis energy and
cannot extend a circle-valued phase over a contractible disk. Since M>=I
in the established parameter family, sending rho to zero does not cure it.
For a truly isolated object in infinite z, nonzero translation-invariant
positive energy per unit length also gives infinite total energy. Changing
only radial boundaries cannot supply finite axial extent.

The rigorous conclusion is a restricted no-go, not a constructed exterior.
A candidate escaping it must change the nonzero axial phase configuration
and, for a radial taper, relax phase-gradient translation invariance; this is
necessary, not sufficient. It may also require time dependence, different
boundary data or topology. Exact vacuum already exists in the unchanged
theory (constant rho=v, Q, any vacuum Einstein metric), so a theory-wide
claim that vacuum itself is forbidden would be false. Connecting nontrivial
matter to it in a different configuration remains unproved.
