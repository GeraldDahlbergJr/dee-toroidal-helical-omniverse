# Reproduce the streamed qualification

Parent evidence checkpoint: `5a18cb1c29e7511ca2a984fa4aafccf04757715b`.
The exact manufactured equations, fixture, mesh, and criteria are specified in
[the pinned protocol](../run43_bspline_qualification/PROTOCOL.md).
The new driver adds the final coupled algebraic defect and the spline
positivity lower bound to the recorded diagnostics. It also stops if either
final coupled relative defect exceeds 1e-10. The linear solve threshold,
boundary and residual-order criteria are inherited.

Use Python 3.12 with NumPy 2.3.5, SciPy 1.17.0, and SymPy 1.14.0.
From the repository root, check the sealed evidence:

```sh
(cd research/run43_bspline_completion && sha256sum -c SHA256SUMS)
```

Run in a disposable directory because the driver writes result files:

```sh
root=$(pwd)
trial=$(mktemp -d)
cp research/run43_bspline_completion/*.py "$trial/"
cd "$trial"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python equivalence.py "$root/research/run43_bspline_qualification/experiment.py" > equivalence.json
python profile_memory.py > profile_memory.jsonl
timeout --signal=TERM --kill-after=5s 1800s python -u run.py > execution.log 2>&1
python compare_results.py "$root/research/run43_bspline_qualification/results.json" results.json > results_equivalence.json
python conditioning.py coefficients_8.npz coefficients_16.npz coefficients_32.npz > conditioning.json
```

The driver and preflights enforce one logical CPU and a 4 GiB address-space
limit before numerical imports. The one full run has a 30-minute wall-clock
limit. The microprofile performs forcing and single operator applications;
it performs no PDE solve. The historical-output comparison needs both
completed small grids. Conditioning can also be run on any subset of the
available archives. Exact replay timing and RSS depend on the environment.

No original-source experiment is part of these commands. In the driver,
`base` retains the original method design checkpoint as historical metadata.
The immediate parent of the streamed change is the evidence checkpoint above.
