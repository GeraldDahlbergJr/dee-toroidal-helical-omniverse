# Agent 2 — reconstruction sensitivity diagnosis

Scope: read-only saved-field diagnostics at starting checkpoint
`5ed79712e8170e1d794b358b16918094f79431f8`. No DEE solve or accepted
archive/evaluator modification. The current primary evidence is
`research/run43_nonpolynomial_validation`, whose audit imports the unchanged
`Evaluator` from `research/run43_final_validation/evaluate.py`. The earlier
polynomial qualification failed; the later nonpolynomial qualification passed
its stated measured-grid test. Those outcomes must not be conflated.

## Quantitative finding

The reported 74–79% momentum envelope is a **radial reconstruction sensitivity**,
not the angular cutoff contribution. Fresh `decompose.py` evaluates the same
four archived fields and separates each envelope ingredient. Components in
this table are Mx, My, Mz; ratios divide a difference RMS by nine-point residual
RMS, using the original physical volume weights.

| n=64 region | 9-minus-7 radial ratio | 7-minus-5 radial ratio | angular-cut ratio (largest component) | full envelope ratio |
|---|---|---|---|---|
| inner band | .745610, .745610, .744184 | .275652, .275652, .276427 | 7.45e-14 | .750743, .750743, .749390 |
| outer band | .778982, .778982, .783075 | .195259, .195259, .193355 | 2.65e-13 | .779455, .779455, .783450 |
| inner first row | .735584, .735584, .734143 | .263735, .263735, .264362 | 3.14e-14 | .735584, .735584, .734143 |

The full envelope takes the pointwise maximum before RMS, so its value need
not equal the maximum of the three separately computed RMS values. This is
not double-counting by summation. At n=32, the inner-band radial 9-minus-7
ratios are .732332, .732332, .728390; the angular-cut ratios are below 4e-12.
The radial sensitivity remains a large fraction of the residual on refinement.

At n=64 the fraction of saved Cartesian W Fourier power above the n//3
angular cutoff is 3.45722e-31; for psi-1 it is 1.05045e-25. The first ratio is
unweighted total sample power over all radii and Cartesian components.
Both the power and derivative-sensitive cutoff experiment demonstrate that
resolved high angular modes do not account for this envelope. They do not
bound frequencies already aliased into the samples.

A final bounded localization diagnostic, recorded in the same JSON, finds
98.8320%, 98.8320%, 99.1389% of the interior weighted squared 9-minus-7
momentum difference within the first four interior radial planes at n=64.
The last four planes contain 1.11146%, 1.11146%, .80425%. This demonstrates
boundary localization of the measured sensitivity; it does not identify a
unique smooth continuum interpolant or prove a gridscale parity mode.

## What this does and does not establish

Confirmed evaluator limitation: the saved-field uncertainty is an empirical
spread among reconstructions, not a certified error enclosure. The smooth
nonpolynomial fixture qualifies the evaluator on that fixture and those grids;
it does not bound radial truncation on these solved fields. No available result
establishes whether the 9-point or 7-point continuum residual is closer to the
unknown exact field residual. Therefore changing the reference width or
removing the spread would not be a justified physical correction.

Confirmed numerical diagnostic failure: under the already declared evaluator,
n=64 inner-band residual RMS is (1.860788e-3, 1.860788e-3, 2.486626e-3), and
32-to-64 nominal orders are (1.027781, 1.027781, .990610). Inner-first-row
orders are (.498995, .498995, .461755). These do not meet p>=1.5 even before
sensitivity adjustment. The <=10% sensitivity requirement also fails. This is
a failure of the present independent validation evidence; it is not proof that
the continuum DEE model has no solution.

Not established: a unique physical origin of the radial sensitivity, an absolute
physical residual bound, a bound on unsampled aliases, or a validated corrected
field. A radial closure defect or grid-dependent boundary error is plausible,
but assigning causal responsibility requires the separate operator investigation
and a demonstrated identity or independent error estimate. This diagnosis
alone justifies no new corrective solve.

Evidence needed to resolve the remaining ambiguity: an independent radial
reconstruction with a justified error estimate on these saved fields (or
independently sampled field/derivative data); a quantitative link between the
accepted discrete equations and the first-row continuum defect; and a declared
physical absolute tolerance plus stability/error evidence for any absolute
solution-accuracy claim. An angular cutoff adjustment is not supported as a
correction by these measurements.

## Reproduction

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python research/priority_one_validation/agent2/decompose.py > research/priority_one_validation/agent2/decompose.log 2>&1
```

The script verifies each input NPZ against its archived result SHA256, performs
no solves, writes only `decomposition.json` in this directory, and checks that
the freshly reconstructed envelope matches the current saved audit. Exact
ratios, all four resolutions, all five regions and input hashes are retained
in that JSON. Baseline arrays, preexisting evaluators and reports are preserved.
