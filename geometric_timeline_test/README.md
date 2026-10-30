# Geometric timeline and closed-input baseline

Prepared for Gerald Ted Dahlberg Jr., October 29, 2026.

This implements a geometric index around the existing local helical–toroidal
trajectory benchmark. It does not use GitHub or alter Run 43. It is not an
autonomous solution of the DEE field equations: the original smooth transition
schedule remains in `test_transition.py`. Instrumentation alone does not remove
that ansatz. The timeline observes the computed trajectory and never changes it
to produce a desired stage or hide a failure.

## Timeline

| Index | Observed stage | Entry condition |
|---|---|---|
| 0 | Counterclockwise helical motion | Positive major-angle rate; full winding and axial rates |
| 1 | Transition toward toroidal circulation | Measured winding rate decreases after at least one full CCW pass |
| 2 | Toroidal circulation | Winding and axial rates are near zero; CCW circulation continues |
| 3 | Restore axial propagation | Winding and axial rates increase after the toroidal hold |
| 4 | Axial propagation restored | Winding and axial rates return near their baseline values |

CCW means viewed from positive z toward the origin. Major-angle and poloidal
angles are distinct and remain lifted (unwrapped). The domain is a translating
toroidal tube throughout; this does not change the domain topology.

Stages are inferred from numerically differentiated trajectory rates, not from
clock-time cutoffs. Fractions are normalized to the existing benchmark's initial
rates. A fraction tolerance of 1e-4 labels a rate as fully active or paused.
The three-turn toroidal guard allows 0.1 turn for stage-detection timing; it is
an instrumentation tolerance, not a proof of exactly three turns. The original
analytic benchmark independently checks its exact three-turn schedule.
The index is currently offline: neighboring samples enter differentiation.
Its event times are resolution- and tolerance-dependent, not exact failure times.

Every accepted row records time, sample number, stage, cumulative CCW turns,
angles, axial position, radius, coordinate Jacobian, azimuthal rate, winding
fraction, axial fraction, radial rate, and RK4 full-step/two-half-step defect.
That defect measures local integration consistency, not an Einstein equation
residual or proof of physical accuracy. The analytic trajectory supplies an
additional independent position-error check for this prescribed ODE only.

## Isolation contract

`ClosedBoundary.rhs` admits the existing baseline RHS and zero external input.
Any nonzero, nonfinite or malformed external-input vector is recorded and raises
an exception before the attempted step can be accepted. The same restriction
is kept through all five stages. No path correction, projection, clipping, or
forced stage advance is used to conceal a violation.

For this ODE, radius is a fixed material coordinate and its exact derivative is
zero, so relative radial velocity through the moving tube wall is analytically
zero. The timeline also checks sampled radial rate and tube membership. It does
not certify arbitrary between-sample paths or exclude a physical external field.
There are no implemented pressure, electromagnetic, gravitational or other PDE
characteristic boundary channels here. Rejecting a software input is not a
demonstration of physical shielding. The prescribed transition schedule remains
a control assumption, even though no additional forcing is admitted.

## Failure localization

The index stops recording accepted rows at the first detected violation and
retains its time, row and diagnostic label. It detects circulation stall,
clockwise reversal, radial flux, leaving the tube, winding/axial rate mismatch,
nonfinite states, skipped stages and incomplete passes. A toroidal axial pause
is intentionally allowed and is not classified as a circulation stall.

Diagnostic labels identify symptoms. They do not by themselves establish the
physical cause. Controlled reruns of a real field solver must compare equation
residuals, incoming characteristic amplitudes, moving-wall flux, conserved
quantities, boundary residuals and timestep/grid refinements. Different physical
influences cannot be universally excluded by setting radial velocity to zero.
Any physically coupled DEE implementation must specify its governing equations
and admissible boundary data before the no-incoming-influence condition can be
made well posed. Do not overwrite Run 43 acceptance criteria with these checks.

## Actual verification

The completed run passed 11 of 11 instrumentation checks. It generated 1,121
timeline rows. Observed stage entries were at nondimensional times 2.05,
3.9625, 10.05 and 11.9625. Differences from the prescribed joins at 2, 4, 10 and
12 reflect rate thresholds and numerical differentiation.

Checks cover the complete five-stage baseline, positive CCW rate, integration
defect, independent trajectory reference, legitimate axial pause, rejection of
external input during the helix, and known injections of circulation stall,
reversal, radial flux, NaN state and winding/axial mismatch. These are software
verification cases, not measurements of DEE behavior or isolated physical causes.
Numeric evidence and source hash are in `timeline_results.json`.

Run with Python 3.10+, NumPy and Matplotlib:

```bash
python3 geometric_timeline.py
```

This writes `geometric_timeline.csv` and `timeline_results.json` beside the script.
The retained `test_transition.py` contains the original mathematical benchmark;
running it independently writes its own original diagnostic files.

## GitHub automation

The dedicated `geometric-transport-benchmark.yml` workflow runs on changes to
this directory or its workflow, on relevant pull requests, and daily at 03:43
UTC (8:43 p.m. Pacific daylight time; 7:43 p.m. Pacific standard time). It also
supports a manual run. Python 3.12 and the pinned NumPy/Matplotlib versions are
used. The two scripts run independently; a failed command cannot be hidden by
the log-writing pipe. Provenance records the commit, run ID, source hashes and
environment. Diagnostics are uploaded even when tests fail and retained for
14 days. GitHub scheduling can be delayed and runs depend on runner availability,
Actions permissions and account quotas.

This workflow executes only this benchmark. It does not write repository
contents, modify Run 43, launch its workflows, change physical sources or accept
these checks as physical validation. The repository's existing Run 43 automations
retain their own behavior. Source changes are required for substantive new
experiments; a scheduled rerun of unchanged deterministic code is a regression
check, not additional scientific progress.
