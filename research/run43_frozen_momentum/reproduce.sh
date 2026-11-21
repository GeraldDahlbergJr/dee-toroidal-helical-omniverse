#!/usr/bin/env bash
# Run from the repository root. The work directory may contain complete,
# previously generated hybrid and 4N reference fields. diagnose.py verifies them.
set -euo pipefail
work=${1:?provide a work directory}
out=${2:?provide a new output directory}
export OPENBLAS_NUM_THREADS=1
mkdir -p "$work"
if [[ ! -d "$work/hybrid" ]]; then
    python research/run43_hybrid/physical.py "$work/hybrid"
fi
status=0
python research/audit_run43_arrays.py "$work/hybrid" || status=$?
[[ $status -eq 1 ]]
for n in 8 16 32 64; do
    if [[ ! -f "$work/levels/a$((4*n))-r$((8*n))/reduced.npz" ]]; then
        python research/run43_asymptotic/solve_levels.py "$work/levels" --factors 4 --grids "$n"
    fi
done
python research/run43_frozen_momentum/diagnose.py "$work/hybrid" "$work/levels" "$out"
cmp "$out/frozen-acceptance.json" research/run43_asymptotic/results/hybrid-audit.json
