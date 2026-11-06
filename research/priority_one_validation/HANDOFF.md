# Run 43 continuity handoff

Development branch: `run43-bspline-validation-completion`.
Evidence parent: `5a18cb1c29e7511ca2a984fa4aafccf04757715b`.
The current verified commit is the tip of the development branch. Verify it
with `git rev-parse HEAD` and
`git ls-remote origin refs/heads/run43-bspline-validation-completion`.

## Current decision

Independent physical validation remains NOT YET ESTABLISHED. PRD readiness
is on hold. See the [streamed qualification report](../run43_bspline_completion/REPORT.md)
for the quantitative gate results, including the resource outcome and the
remaining total uncertainty criterion. The single authorized memory correction
and qualification rerun are spent. This handoff authorizes no further run.

## Preserved evidence

- The accepted frozen baseline and its independent numerical replay passed
  12/12. The independent nonpolynomial evaluator fixture passed.
- The baseline independent continuum audit did not establish validation.
  Finest inner-band momentum orders were approximately 1.03, 1.03 and 0.99.
  Momentum radial reconstruction sensitivity was approximately 74–79%.
- The composed radial closure has a first-order first-interior-row defect.
  Its contribution and the remaining momentum residual cancel strongly in
  the frozen diagnostic. No source or stress-energy algebra correction was
  established.
- Rejected approaches include matched endpoint and local/global radial
  substitutions, direct Hessians, fourth-vector, compatible closures,
  coordinate KO, and radial filters. Preserve their results.
- The archived degree-four B-spline fixture hit 4 GiB in symbolic forcing
  evaluation at n=32. The new implementation streams the original forcing
  and weak forms with one mathematical discretization. Small-grid equivalence
  and the n=32 preflight are preserved in the new evidence directory.

## Next decision and commands

The missing certification is a per-metric enclosure of numerical evaluator
and quadrature uncertainty on the archived spline fields. Assess propagated
assembly and solver effects according to the frozen protocol. Start from
the saved coefficients and the exact q6 Jacobian certificate. No additional
DEE solve is needed for this audit. Original-source execution requires a
complete qualification pass and predeclared physical acceptance tolerances.

From the repository root:

    (cd research/run43_bspline_completion && sha256sum -c SHA256SUMS)

Use the [new reproduction commands](../run43_bspline_completion/REPRODUCE.md)
in a disposable directory. The previous [diagnosis](REPORT.md) and its
[reproduction commands](REPRODUCE.md) remain available. The historical
[fixture protocol](../run43_bspline_qualification/PROTOCOL.md) and
[resource failure report](../run43_bspline_qualification/REPORT.md) are
unchanged.
