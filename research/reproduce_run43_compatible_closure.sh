#!/usr/bin/env bash
# A fresh tag is mandatory: never overwrite historical evidence.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${DEE_REPLAY_TAG:?Set a unique DEE_REPLAY_TAG}"
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier
evidence="research/results/${DEE_REPLAY_TAG}-evidence"
mkdir "$evidence"
python -m unittest discover -s tests -v > "$evidence/tests.log" 2>&1
for variant in metric4 coordinate second filter; do
  case "$variant" in
    metric4) script=run43_compatible_closure_metric4.py ;;
    coordinate) script=run43_compatible_closure.py ;;
    second) script=run43_second_hessian.py ;;
    filter) script=run43_radial_filter.py ;;
  esac
  export DEE_EXPERIMENT="${DEE_REPLAY_TAG}-${variant}"
  out="research/results/${DEE_EXPERIMENT}"
  python -u "research/$script" > "$evidence/${variant}-solve.log" 2>&1
  # A scientifically failing gate must not prevent preservation or later audits.
  set +e
  python research/audit_run43_arrays.py "$out" > "$evidence/${variant}-audit.log" 2>&1
  audit_status=$?
  set -e
  test -f "$out/independent-audit.json"
  python research/validate_original_physics.py "$out" > "$evidence/${variant}-physics.log" 2>&1
  set +e
  python research/verify_run43_closure_fields.py "$out" > "$evidence/${variant}-independent-implementation.log" 2>&1
  implementation_status=$?
  set -e
  test -f "$out/independent-implementation-audit.json"
  python research/check_run43_closure_endpoint.py "$out" > "$evidence/${variant}-endpoint.log" 2>&1
  printf '%s frozen=%s implementation=%s\n' "$variant" "$audit_status" "$implementation_status" >> "$evidence/exit-statuses.txt"
done
find "$evidence" "research/results/${DEE_REPLAY_TAG}-metric4" "research/results/${DEE_REPLAY_TAG}-coordinate" "research/results/${DEE_REPLAY_TAG}-second" "research/results/${DEE_REPLAY_TAG}-filter" -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > "$evidence/SHA256SUMS"
