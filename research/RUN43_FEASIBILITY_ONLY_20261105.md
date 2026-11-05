# Run 43 — feasibility after two exact obstructions

Base: `194ed894cc35e0671b172569703f5dfd164b8fd5`.
One agent; existing evidence and symbolic reasoning only. No candidate is
constructed, implemented or tested; no certificate is rerun.

## Necessary changes within the analyzed family

The four-row family fixes uniform nodal locations, trapezoidal radial norm,
scalar SBP with endpoint traces, constant/affine exactness, centered rows beyond
the block, and the collocated Cartesian STF strain with its weighted adjoint.
Its existing [derivation](run43_four_row_boundary/derivation.md) and
[certificate](run43_four_row_boundary/certificate.json) give quadratic residuals
2(3t-2)(3t-1) and 4(2t-1)(3t-1). Their common zero is t=1/3, where the next
affected rows have residuals 1/18 and -1/6. Thus no parameter in this family
has even vanishing pointwise boundary truncation error. Lowering the requested
positive consistency order does not help: the quadratic defect is O(1).
All rows affected by the composition, including transition rows, are necessary
moment constraints; satisfying just the first two is insufficient.

Consequently at least one defining discretization restriction must change.
If the nodes, norm, scalar SBP identity, polynomial requirements, centered
exterior rows and same strain/adjoint composition are all retained, the
permitted boundary coupling must extend beyond the analyzed four-row support.
This is a conditional necessity, not a proof that any particular larger width
works. The evidence determines no smallest feasible width or number of new
unknowns. Merely adding formal parameters that eliminate back to this same
one-parameter family changes nothing.

Independently, at t=1/3 the first two interior derivative columns coincide.
For u=(0,1,-1,0,...,0), D u=0; for angular-constant Cartesian W=u c,
L W=0 exactly, although W is nonzero and has zero radial trace. With
L0=L I and positive definite tensor weight Ht, the existing identity is

    -M A = L0^T Ht L0 / 2,
    ker(L0^T Ht L0) = ker(L0).

Strict coercivity in a norm positive on every admissible nonzero vector
therefore requires changing the strain map or admissible discrete space so
this mode is detected or excluded by a consistent approximation. Reweighting
the same zero strain, changing a positive norm, or changing boundary traction
cannot remove it: E I=0 kills boundary work. A quotient norm hiding the mode
does not establish the requested coercivity on the original Dirichlet space.
Injectivity is necessary for fixed-grid strict positivity; uniform coercivity
additionally requires a mesh-independent discrete Korn estimate in a stated
norm. Removing this one mode alone proves neither that estimate nor absence
of other modes.

These two requirements must hold simultaneously: consistent boundary moments
and a strain space without spurious Dirichlet kernel. They are not logically
equivalent. The earlier [endpoint certificate](run43_boundary_strain/certificate.json)
already illustrates that a positive normal Gram matrix can coexist with failed
boundary consistency.

## Possible changes, not established remedies

| Option | What the evidence warrants |
| --- | --- |
| Wider boundary stencil/block | Necessary only under the retained restrictions stated above; may permit different moments and distinguish the coincident columns. Neither sufficiency nor minimum width is known. |
| Additional boundary or auxiliary degrees of freedom | Not proved necessary. They could enlarge the approximation space, but must be coupled consistently without adding physical boundary data or leaving an uncontrolled kernel. |
| Different discrete strain/space | A change in its action or domain is necessary to eliminate the exhibited kernel. A different placement, space or derivative realization is possible; changing the continuum STF formula is not required by the evidence. |
| Compatible quadrature and adjoint structure | Required jointly if the exact Green/energy identity is retained. A different quadrature is not itself proved necessary and cannot alone repair the fixed strain kernel. |

The [quadrature no-go report](RUN43_QUADRATURE_NO_GO_20261105.md) and its
[certificate](run43_quadrature_no_go_certificate.json) further exclude repairing
the original expanded momentum operator and original collocated strain by
*any* positive nodal quadrature alone, including geometry-dependent weights.
The paired angular modes have equal pointwise strain energy but unequal
expanded work. Thus retaining that expanded matrix, strain and nodal quadrature
class together is impossible for the requested energy identity. At least one
must change; a dense matrix inner product is outside that certificate's scope
and remains unproved. An adjoint definition can supply an algebraic identity
without supplying boundary consistency or coercivity, as the failed closures
already demonstrate. No exact Green identity is asserted necessary for every
convergent discretization; it is a requirement of the energy-compatible family
being assessed here.

## Physical restrictions and decision

Preserve the original continuum momentum equation div(LW)=8 pi G psi^6 S,
Cartesian STF strain, Hamiltonian coupling, original matter sources and
G=.001, toroidal shell a=1/20 to b=1/4, radial psi=1 and W=0, and angular
periodicity. Modifying these would change the original physical problem.
Neither zero traction, an extra normal-derivative condition, W'''=0, nor a
zero-net-source condition may be added as a closure device. Dirichlet data
allow nonzero traction and boundary reaction; see the
[necessary-adjoint report](RUN43_ADJOINT_FLUX_NECESSARY_20261105.md).
Numerical stencil widths, approximation spaces and quadrature rules are
discretization choices, provided their consistency with that same problem is
proved. The frozen audit remains a separate unchanged acceptance requirement;
redesigning a discrete operator cannot silently redefine it.

Decision: the analyzed four-row family is infeasible. Escaping both
obstructions requires a genuine departure from its moment restrictions and a
strain/domain change that removes the spurious kernel, with compatible positive
pairings if the energy identity is required. Existing certificates do not select
a successful redesign or establish a numerical minimum. No alternative is
authorized or advanced here. Only this report is added; physical sources,
production solver, frozen audit and historical evidence remain unchanged.
Stop pending separate authorization.
