"""Record provenance and surface benchmark failures without changing results."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import platform
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
results = {}
failed = False
for name in ('results.json', 'timeline_results.json'):
    path = HERE/name
    try:
        data = json.loads(path.read_text())
        checks = data['checks']
        if not checks or any(c.get('passed') is not True for c in checks):
            failed = True
        results[name] = {'passed': sum(c.get('passed') is True for c in checks),
                         'total': len(checks),
                         'failures': [c['name'] for c in checks if c.get('passed') is not True]}
    except (OSError, ValueError, KeyError, TypeError) as exc:
        failed = True
        results[name] = {'error': str(exc)}
sources = [*HERE.glob('*.py'), HERE/'requirements.txt',
           ROOT/'.github/workflows/geometric-transport-benchmark.yml']
report = {
    'scope': 'Prescribed kinematic benchmark and timeline instrumentation; not DEE physical validation',
    'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'python': sys.version, 'platform': platform.platform(),
    'github': {k: os.environ.get(k) for k in
               ('GITHUB_REPOSITORY', 'GITHUB_SHA', 'GITHUB_REF', 'GITHUB_RUN_ID',
                'GITHUB_RUN_ATTEMPT', 'GITHUB_EVENT_NAME')},
    'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in sources if p.exists()},
    'results': results,
    'benchmark_checks_passed': not failed,
}
(HERE/'ci_provenance.json').write_text(json.dumps(report, indent=2)+'\n')
summary = os.environ.get('GITHUB_STEP_SUMMARY')
if summary:
    with open(summary, 'a') as f:
        f.write('## Geometric transport benchmark\n\n')
        f.write('Prescribed geometry and instrumentation only; no DEE physical validation.\n\n')
        for name, result in results.items():
            f.write(f'- {name}: {json.dumps(result)}\n')
        f.write('\nFull diagnostics, source hashes and environment are in the run artifact.\n')
print(json.dumps(report, indent=2))
raise SystemExit(int(failed))
