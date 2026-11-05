# Run 43 — controlled six-row screen

Base: 80c48400b80c8fd391f10cdd06a4f1d9c53ec490.
One agent, one wider family, exact symbolic arithmetic, 16 radial nodes.

Boundary moments: **INFEASIBLE**. Six free boundary parameters remain after
constant/affine and SBP constraints, but all affected and transition quadratic
moments obey sum i H_ii R_i=-1/2. They cannot simultaneously vanish.
Dirichlet strain kernel: **UNRESOLVED** for a moment-feasible construction,
because no such member exists. A checked illustrative member detects the old
mode but has other exact kernel modes; this is no simultaneous remedy.
Adjoint identity: **PASS**, by the full weighted STF Green identity.
Coercivity: **UNRESOLVED** for a consistent construction.

See derivation.md, check_exact.py and certificate.json. Reproduce only this new
certificate with `python research/run43_six_row_controlled/check_exact.py`.
Historical evidence was read without rerunning it. Production, physical
sources, frozen audit and historical files are unchanged. No solver, GMRES,
physical convergence sweep or GitHub Actions was run.

Next necessary step: separately authorize a discretization family that relaxes
the trapezoidal-norm/centered-transition restrictions; widening this tested
block alone failed. Physical solves remain unauthorized.
