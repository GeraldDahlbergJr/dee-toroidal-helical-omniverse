#!/usr/bin/env bash
# Run from the archived repository root; all source files remain unchanged.
set -uo pipefail
tag=${1:?Supply a fresh output prefix}
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier
mkdir -p "live-evidence/$tag"
for variant in baseline fourth local; do
  case "$variant" in
    baseline) script=research/run43_consistent_momentum.py ;;
    fourth) script=research/run43_fourth_vector.py ;;
    local) script=research/run43_local_boundary.py ;;
  esac
  output="research/results/$tag-$variant"
  if [[ -e "$output" ]]; then
    echo "Refusing to overwrite $output" >&2
    exit 1
  fi
  DEE_EXPERIMENT="$tag-$variant" python -u "$script" > "live-evidence/$tag/$variant-solve.log" 2>&1
  status=$?
  echo "$status" > "live-evidence/$tag/$variant-solve.exit"
  if [[ "$status" == 0 ]]; then
    python research/audit_run43_arrays.py "$output" > "live-evidence/$tag/$variant-audit.log" 2>&1
    echo "$?" > "live-evidence/$tag/$variant-audit.exit"
    python research/validate_original_physics.py "$output" > "live-evidence/$tag/$variant-physics.log" 2>&1
    echo "$?" > "live-evidence/$tag/$variant-physics.exit"
  fi
done
# A zero shell exit is execution completion, not scientific acceptance.
