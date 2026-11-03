# Independent Run 43 boundary mathematics review

This review changes no repository source, baseline, physical source, checkpoint,
frozen audit, acceptance band, row selection, or threshold. Reference branch:
`origin/research/run43-boundary-continuation-20261103`, commit
`6aefba3` (full immutable commit recorded by the coordinator). This is an
independent algebra check, not a new physical solution or acceptance result.

## Reduction from the original action

The original model is documented in
`docs/RECOVERED_NONLINEAR_DEE_FORMULATION.md`, sections Covariant action,
Stress-energy, and Euler-Lagrange equations. With scalar fields phi^A and kinetic
matrix K_AB, L_m=-K_AB grad(phi^A).grad(phi^B)/2-U. Consequently
T_mn=K_AB partial_m phi^A partial_n phi^B+g_mn L_m. For normal derivatives
Pi^A and spatial covectors D_i phi^A this gives

    E = (K_AB Pi^A Pi^B + K_AB gamma^ij D_i phi^A D_j phi^B)/2 + U,
    S_i = -K_AB Pi^A D_i phi^B.

Taking the normal-normal and normal-spatial projections of G=8 pi G T gives
R+K^2-K_ij K^ij=16 pi G E and D_j(K^ij-gamma^ij K)=8 pi G S^i.
For the explicitly limited maximal, conformally flat initial-data ansatz
gamma_ij=psi^4 delta_ij, K=0, K^ij=psi^-10 A^ij, and A=LW, these become

    Delta psi + (A:A) psi^-7/8 + 2 pi G E(psi) psi^5 = 0,
    partial_j A^ij = Delta W_i + partial_i div W/3 = 8 pi G psi^6 S_i.

In this benchmark the original matter fields, their spatial covectors, and their
normal velocities remain fixed; E is recomputed with gamma^-1=psi^-4 delta.
These projections do not solve the matter evolution equations or all spacetime
Einstein equations. This scope is explicit in the original solver's docstring,
`reproduce/solve_toroidal_helical_coupled_constraints.py` on the reference branch.

## Boundary conditions and consistency

The finite-shell model prescribes psi=1 and W=0 at r=.05 and .25, with periodic
angles. These are model inputs, not exterior matching conditions derived from
the action. The second-order elliptic constraints do not grant an additional
independent W'''=0 boundary condition. A concrete counterexample is
f(r)=a+b-r-ab/r=(r-a)(b-r)/r, a=.05,b=.25: f(a)=f(b)=0 while
f'''=6ab/r^4 is nonzero. Tangential derivatives of W vanish at the Dirichlet
boundary but the normal third derivative need not.

Let D_r be centered inside and the three-point one-sided derivative at endpoints.
At the first interior row, direct matrix multiplication gives

    D_r D_r f = (3f0-5f1+f2+f3)/(4h^2),
    (D_r D_r-D_rr)f = (-f0+3f1-3f2+f3)/(4h^2)
                       = h f'''/4 + h^2 f''''/8 + O(h^3).

The direct centered D_rr contributes h^2 f''''/12, yielding the documented
DrDr expansion f''+h f'''/4+5h^2 f''''/24. At the outer adjacent plane the
odd term reverses sign. The accompanying executable uses exact rational matrix
arithmetic to verify the signed coefficients and derivative moments; its JSON
includes the script hash. No production forcing is constructed by this test.

The unit radial Cartesian vector e_r is independent of r. Thus the second-normal
derivative coefficient of Delta W+grad div W/3 is I+e_r e_r^T/3, with positive
eigenvalues 1,1,4/3. This uniquely justifies the existing local correction

    P_edge (I+e_r e_r^T/3)(D_rr-D_r D_r)W.

It restores the direct second-order normal stencil at the adjacent planes. It
does not fix the angular or finite-product-rule differences between expanded
and divergence-of-tensor discretizations. Coefficient tuning or adding an audit
defect to a higher-order physical solve would have no accuracy derivation here.

## Finite-grid gate and decision

One smooth O(h) plane has O(h^1.5) fixed-volume-band RMS asymptotically, because
its radial volume fraction is O(h). This statement alone does not imply the
specified n=8,16,32,64 pairwise orders reach 1.5. The physical first-row radius
changes at each refinement, and the leading coefficient can vary rapidly; the
counterexample has f''' proportional to r^-4. No row may be excluded to hide it.

The reference branch records original-source local-boundary M2 inner orders
0.234124077,0.727095846,1.116755252 and outer orders
1.998783011,1.532308823,1.394913922. With the six passing H checks this is 8/12,
so reject. Fourth-order inner M2 orders are
0.330580315,0.855566995,1.112652071; its H and outer M2 checks pass, giving 9/12,
so reject. These numbers describe preserved historical executions only; any
fresh execution must be reported separately by the solver and validation agents.

The frozen audit blob remains
`5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`. The Nov 2 12/12 baseline remains
the successful reference. No mathematically justified new boundary condition or
free correction parameter is established by this independent review. A fresh
replay of the two rejected candidates can establish reproducibility, not change
their acceptance status absent genuinely passing unchanged checks.

Reproduce this independent exact algebra certificate with:

    python /workspace/run43-investigation/boundary/verify_boundary_algebra.py
