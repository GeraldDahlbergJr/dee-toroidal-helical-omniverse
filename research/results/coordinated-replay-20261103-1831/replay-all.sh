#!/usr/bin/env bash
# Replay genuine original-source solves into new paths; never overwrite evidence.
set -euo pipefail
run43_evidence=$(cd "$(dirname "$0")" && pwd)
cd "$run43_evidence/../../.."
run43_tag=${1:?Supply a unique alphanumeric output tag}
[[ "$run43_tag" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || exit 2
run43_out="research/results/$run43_tag-coordination"
mkdir "$run43_out"
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier
sha256sum -c "$run43_evidence/SHA256SUMS" > "$run43_out/preserved-evidence-check.log"
bash "$run43_evidence/solver/replay.sh" "$run43_tag"
for run43_variant in baseline fourth local; do
  python "$run43_evidence/validation/verify_fresh.py" "$PWD" \
    "research/results/$run43_tag-$run43_variant" \
    "$run43_out/$run43_variant-independent-validation.json" \
    > "$run43_out/$run43_variant-independent-validation.log" 2>&1
done
cp -a "$run43_evidence/operator" "$run43_out/operator"
cp -a "$run43_evidence/boundary" "$run43_out/boundary"
python "$run43_out/operator/verify_operator_independent.py" > "$run43_out/operator/replayed-verification.log" 2>&1
python "$run43_out/boundary/verify_boundary_algebra.py" > "$run43_out/boundary/replayed-algebra.log" 2>&1
python research/certify_run43_first_row.py "$run43_out/coefficient-certificate" > "$run43_out/coefficient-certificate.log" 2>&1
python research/check_run43_physical_row_identity.py "$run43_out/sparse-physical-identity" \
  "research/results/$run43_tag-fourth" > "$run43_out/sparse-physical-identity.log" 2>&1
python "$run43_out/operator/verify_fresh_physical_identity.py" "research/results/$run43_tag-fourth" \
  "$run43_out/array-physical-identity" > "$run43_out/array-physical-identity.log" 2>&1
python -m unittest discover -s tests -v > "$run43_out/regression.log" 2>&1
python - "$run43_out" <<'PY'
import json, pathlib, sys
out=pathlib.Path(sys.argv[1])
reports={name:json.loads((out/f'{name}-independent-validation.json').read_text()) for name in ('baseline','fourth','local')}
summary={name:{'passed':d['passed'],'total':12,'status':d['status']} for name,d in reports.items()}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
# Scientific failure remains a nonzero exit, after retaining all diagnostics.
sys.exit(0 if all(d['status']=='PASS' for d in reports.values()) else 1)
PY
