# Stage A½: boundary compatibility and limitations

Date: 2026-11-10 UTC. Baseline: `f2da9c5fa8d3f34275afc8f0651366adbf8e65d5`.
Conventions, exact radial solutions and assumptions are in `DEE_STAGE_A_HALF_BOUNDED_RESPONSE.md`.

## What the boundary problem proves

**VERIFIED:** on a bounded transverse disk or annulus the forced equation has a unique smooth solution for smooth Dirichlet data, or for homogeneous Neumann data. Choosing u=0 on the boundary preserves the original amplitude trace rho=rho_0. Choosing normal derivative zero preserves its original normal derivative. Neither choice alone establishes matching to an unchanged exterior. A finite mathematical boundary supplies data; the original action supplies no material wall that enforces arbitrary data. An annular theorem assumes a genuine inner boundary and does not construct its core.

For the original scalar kinetic term, source-free matching across a regular interface requires continuity of rho and of its normal flux n^mu partial_mu rho. Since the background rho_0 is constant, to first order these conditions against an unchanged amplitude exterior are u=0 and partial_n u=0. Metric corrections multiplying the background amplitude derivative vanish. The two phases likewise require compatible traces and normal currents; their fixed values Q=kz+Q_0 and zero radial flux are compatible locally in this transverse ansatz. This fact does not enforce amplitude matching or gravitational junction conditions.

## Overdetermination: disk

The regular disk response has only one homogeneous coefficient. Requiring u(R)=0 fixes it. Its radial derivative is

    u_D'(R)=(F/m³)[(t²+4) I_1(t)/I_0(t)-2t], t=mR>0.

The bracket is strictly positive. For example the series coefficient of t^{2j+1} in (t²+4)I_1(t)-2t I_0(t) is zero at j=0 and, for j>=1, is

    [j/(j+1)]/[2^{2j-1}(j-1)!j!] > 0.

Equivalently the strong maximum principle and Hopf lemma give u_D'(R) with the sign of F. Hence for F!=0 it cannot vanish. Regularity at the center plus both u(R)=0 and u'(R)=0 is inconsistent, not merely a count of conditions. The zero Neumann solution instead has a nonzero boundary trace and nonzero mean for nonzero F.

On an annulus the two homogeneous coefficients allow two independent boundary conditions total. One can prescribe u and u' at one boundary as an ODE initial-value problem, giving a smooth solution across any finite annulus; this does not generally satisfy prescribed data at the other boundary or provide a regular disk core. Mixed data at separate boundaries are mathematically valid. Requiring both original amplitude and derivative at both boundaries is inconsistent for F!=0.

## Integral obstruction, including the annulus

Let v(r)=I_0(mr)>0, so Lv=0. Green's identity gives

    F integral_domain v r² d²x
       = integral_boundary (v partial_n u-u partial_n v) dl.

If u=partial_n u=0 on the entire boundary, the right side is zero while the left side has the nonzero sign of F. This proves impossibility on either domain without relying on equation counting. The radial verifier checks the associated identity

    [r(vu'-uv')]' = F v r³

on every proposed solution, and explicitly verifies the positive integral primitives. For zero Dirichlet data on an annulus, Hopf also forces nonzero outward derivative on each boundary. Coordinate derivatives at inner and outer boundaries have opposite signs, since the inner outward normal is -partial_r.

**Exception:** F=0. Under zero Dirichlet or zero Neumann data uniqueness gives u=0, and both unchanged traces then match at this order. This is a coefficient exception (including C=0 or B+2epsilon=0 where allowed), not a construction of an exterior or full spacetime. With G>0, rho_0>0, positive M_0, V>=0 and nonzero k, B+2epsilon=2(B+V)>0, so C=0 is the remaining forcing exception. m=0, R=0 or a collapsed annulus are outside the assumed positive-mass nondegenerate domain.

## Exterior and global limitations

**CORRECTION REQUIRED:** retaining the original frozen amplitude throughout an open exterior in the same supplied metric fails its amplitude equation there whenever Fr²!=0. Even an attempted join of traces cannot cure that bulk residual. An exterior responding to gravity could have matching nonzero derivative and trace data; its actual matter equations, metric continuation, gravitational constraints and junction conditions must be solved. Such a continuation is **NOT YET ESTABLISHED** here. Matching to flat original matter while changing the exterior metric is also unproved. No surface action or artificial scalar, phase or stress source is introduced to repair a mismatch.

**VERIFIED:** the forced equation, if artificially extended with the same Fr² to the whole plane, has no bounded radial solution for F!=0. A regular-center radial solution is p+A_I I_0(mr). For A_I!=0 it grows exponentially; for A_I=0 it has quadratic growth. On an exterior annulus, removing I_0 avoids exponential growth but K_0 decays and cannot cancel p's quadratic growth. These conclusions concern this radial equation; extending the weak metric itself to infinity is already invalid since its components grow as r². The bounded-domain solution supplies no asymptotically flat isolated body or finite global energy.

The statements remain distinct:

| Question | Classification | Result |
|---|---|---|
| Local bounded transverse response | VERIFIED | Smooth unique Dirichlet/Neumann solutions under the stated hypotheses |
| Small amplitude and weak metric | VERIFIED, conditional | Explicit sufficient bounds in the response document; finite domain alone is insufficient |
| Both first-order phase currents with fixed phases | VERIFIED | Each divergence residual is zero in the given metric |
| Match both unchanged amplitude and derivative on the entire boundary | CORRECTION REQUIRED for F!=0 | Green identity and disk derivative give an obstruction; F=0 is the exception |
| Actual responding exterior and gravitational matching | NOT YET ESTABLISHED | No exterior solution or junction construction supplied |
| Bounded global radial response for unchanged nonzero forcing | CORRECTION REQUIRED | Polynomial/exponential asymptotics forbid it |
| Full coupled Einstein–DEE solution and global stability | NOT YET ESTABLISHED | Only leading weak-gravity matter response has been checked |

With fixed domain and data O(G), u=O(G), delta T=O(G), and its gravitational source G delta T is O(G²). Therefore the baseline leading metric can consistently be used for this order's local matter response. This order counting does not solve higher-order Einstein equations. Standing MATTER EXISTENCE: ESTABLISHED and LOCALIZATION: OBSTRUCTION IDENTIFIED retain their previous restricted scopes; COMPLETE EINSTEIN–DEE SPACETIME remains NOT YET ESTABLISHED. No historical result is promoted to a complete solution.

## Evidence limits

The new verifier uses exact SymPy expressions for p, both Bessel modes, disk and annular Dirichlet/Neumann conditions, radial residuals, the matching identity and each of the two phase-current equations. Bessel recurrences and derivative identities are explicit. Nonvanishing determinants, existence and matching obstructions rely on the analytical proofs, not floating-point evaluation. The checks are independent algebraic routes by one agent, not independent peer review. Only the three new evidence files are committed; no prior calculations are rerun except checking the inherited phase formula on the newly proposed response. No other research stage follows remote verification.
