#!/usr/bin/env bash
# Fresh names preserve prior evidence, including scientific failures.
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier
tag=${1:?Supply a unique fresh experiment tag}
logs="research/results/$tag-logs"
if [[ -e "$logs" || -e "research/results/$tag-fresh" || -e "research/results/$tag-endpoint" ]]; then
  echo 'Refusing to overwrite existing research evidence' >&2
  exit 1
fi
mkdir -p "$logs"
test "$(git hash-object reproduce/diagnose_toroidal_boundary_convergence.py)" = 5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22 || exit 1
python -m unittest discover -s tests -v > "$logs/regression.log" 2>&1 || exit 1
python research/compare_run43_operators.py "$logs/operator-comparison.json" > "$logs/operator-comparison.log" 2>&1 || exit 1
for variant in fresh endpoint; do
  script=research/run43_consistent_momentum.py
  if [[ "$variant" == endpoint ]]; then script=research/run43_endpoint_continuation.py; fi
  DEE_EXPERIMENT="$tag-$variant" python -u "$script" > "$logs/$variant-solve.log" 2>&1
  status=$?
  echo "$variant solve exit status: $status" >> "$logs/exit-status.log"
  if [[ "$status" == 0 ]]; then
    python research/audit_run43_arrays.py "research/results/$tag-$variant" > "$logs/$variant-frozen.log" 2>&1
    echo "$variant frozen audit exit status: $?" >> "$logs/exit-status.log"
    python research/validate_original_physics.py "research/results/$tag-$variant" > "$logs/$variant-physics.log" 2>&1
    echo "$variant supplementary audit exit status: $?" >> "$logs/exit-status.log"
  fi
done
# A completed script does not imply scientific acceptance. Inspect the twelve
# unrounded orders, source-projection checks, and supplementary solution orders.
