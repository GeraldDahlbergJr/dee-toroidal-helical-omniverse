# Independent physical review — Run 43, November 3

This investigation preserves the frozen numerical gate. Passing that gate is a
scoped initial-constraint result, not a proof of full Einstein–DEE evolution.

## Equation provenance and source interpretation

The repository's `docs/RECOVERED_NONLINEAR_DEE_FORMULATION.md` transcribes an
archived author manuscript. The source action has three fields with kinetic
matrix `diag(1,KTheta,KPsi)` plus symmetric Theta/Psi coupling Lambda, and a
potential. Its variation gives `G_mn=8 pi G T_mn`, the rho scalar equation, and
two conserved transport currents. The original manuscript itself is not present
as an independently inspectable archive in the current main tree; the versioned
transcription is the available equation provenance. This does not authorize
inventing additional physical conditions.

For normal field velocities Pi and spatial covectors D, direct stress projection
gives `E=(Pi^T K Pi + gamma^ij D_i^T K D_j)/2+V` and
`S_i=-Pi^T K D_i`. With `gamma=psi^4 delta`, fixed Pi,D therefore require
`E=temporal/2+psi^-4 spatial/2+V`, while the covariant S_i stays fixed. Holding
E fixed as psi changes would change the original physical problem. This review
uses an independent direct four-dimensional stress construction at actual
solved metric samples in `verify_fresh.py`, retaining the mixed kinetic term.

For maximal conformally flat initial data, `A^ij=psi^-10 (L W)^ij` gives
`Delta psi+|L W|^2 psi^-7/8+2 pi G E psi^5=0` and
`partial_j(L W)^ij=8 pi G psi^6 S_i`. Consequently the independently measured
physical momentum is `psi^-10 div(L W)-8 pi G psi^-4 S_i`.

The shell has r in [.05,.25], major radius 1, periodic angles and psi=1,W=0
Dirichlet boundaries. These are artificial benchmark boundary data, not a
derived physical exterior, asymptotic flatness, an axis condition, or an extra
condition on the third normal derivative. The stress seed is axisymmetric;
Cartesian vector components' phi dependence does not establish generic
nonaxisymmetric physics. Solving the initial constraints is legitimate without
also solving matter evolution, but cannot demonstrate its later survival.

## Independent accuracy limitation

The November 3 preserved documents report 12/12 for the baseline,
9/12 for the fourth-vector candidate and 8/12 for both the older matched-endpoint
candidate and the later local-boundary candidate. They are distinct sources;
the fresh selected boundary execution uses `run43_local_boundary.py`, whose
correction only changes first/last interior radial planes. These are historical
claims pending the fresh executions. The original continuation also reports baseline
supplementary inner momentum orders -0.16187,.69641,1.13020 using different
derivatives. The baseline's frozen pass therefore cannot alone establish robust
boundary accuracy across independent operators. The endpoint's improved
supplementary orders do not replace its failed frozen gate.

Discrete composed derivatives and direct second derivatives need not coincide.
For the original centered first derivative and one-sided endpoint closure,
`(Dr Dr-Drr)f` at row 1 equals
`(-f0+3f1-3f2+f3)/(4 h^2)`, hence `h f'''/4+O(h^2)`. Dirichlet values do not
force this term to vanish. Including the original first interior row in the
frozen band is essential; no row or threshold is changed in this review.
An accurate solution can retain an O(h) audit truncation term at this row.
This observation does not prove a 12/12 consistent method cannot exist, and
does not justify adding an audit defect to the physical equations.

## Fresh validation protocol

`verify_fresh.py` reads immutable completed NPZ fields, verifies file hashes
against execution result manifests, grid shapes, finiteness, positive psi,
exact Dirichlet boundary values, and the exact frozen audit Git blob. It
recalculates all twelve H/M2 pairs using the unchanged audit and writes only
inside this validation directory. It independently checks direct action-derived
4D sources at four actual samples at every n. It uses no manufactured field in
physical acceptance and never fits a threshold or removes a row.

Full Einstein–DEE evolution, current conservation over time, curvature and
trapped-surface survival, global phase regularity, and an absolute continuum
error bound are not established by these initial-constraint runs. The available
fixed boundary data and original source must remain intact rather than being
adjusted after inspection to make a numerical gate pass.
