#!/usr/bin/env bash
# Run from the repository root with the versions pinned in REPORT.md.
set -euo pipefail
work=$(realpath -m "${1:?provide a new scratch directory}")
study=research/run43_frozen_boundary
prior=research/run43_frozen_momentum
export OPENBLAS_NUM_THREADS=1
mkdir "$work"
mkdir "$work/coupled"
cp research/run43_coupled_recovery/results/fields/{fields_*,result_*} "$work/coupled/"
status=0
python research/audit_run43_arrays.py "$work/coupled" || status=$?
[[ $status -eq 1 ]]
cmp "$work/coupled/independent-audit.json" research/run43_coupled_recovery/results/fields/independent-audit.json
"$prior/reproduce.sh" "$work/prior" "$work/replay"
diff -r "$prior/results" "$work/replay"
python "$study/investigate.py" "$work" "$work/results"
