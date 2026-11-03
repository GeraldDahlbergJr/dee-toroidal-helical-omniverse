"""Record all actual solve results and verify original files were preserved."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
evidence = root / 'live-evidence'
original = json.loads((evidence / 'preexisting-sha256.json').read_text())
preservation = dict(
    checked=len(original),
    missing=[p for p in original if not (root / p).is_file()],
    changed=[p for p, digest in original.items() if (root / p).is_file()
             and hashlib.sha256((root / p).read_bytes()).hexdigest() != digest],
)
assert not preservation['missing'] and not preservation['changed'], preservation
audit = (root / 'reproduce/diagnose_toroidal_boundary_convergence.py').read_bytes()
preservation['frozen_audit_git_blob'] = hashlib.sha1(
    b'blob ' + str(len(audit)).encode() + b'\0' + audit).hexdigest()
assert preservation['frozen_audit_git_blob'] == '5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22'
(evidence / 'preservation-final.json').write_text(json.dumps(preservation, indent=2) + '\n')

summary = {}
new_files = []
for variant in ['baseline', 'fourth', 'local']:
    directory = root / 'research/results' / ('live-20261103-' + variant)
    report = json.loads((directory / 'report.json').read_text())
    assert [r['resolution'][1] for r in report['rows']] == [8, 16, 32, 64]
    assert len(report['checks']) == 12
    summary[variant] = dict(
        status=report['status'], passed=report['passed'], checks=report['checks'],
        solve_exit=int((evidence / f'{variant}-solve.exit').read_text()),
        grids=[dict(n=r['resolution'][1], seconds=r['seconds'],
                    nonlinear_iterations=r['iterations'],
                    maximum_true_vector_residual=max(v['true_residual'] for v in r['gmres']),
                    gmres_infos=[v['info'] for v in r['gmres']],
                    fields_sha256=r['sha256']) for r in report['rows']],
        decision='preserved reference' if variant == 'baseline' else 'REJECTED; not promoted',
    )
    new_files.extend(p for p in directory.rglob('*') if p.is_file())
(evidence / 'solve-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
new_files.extend(p for p in evidence.rglob('*') if p.is_file() and p.name != 'SHA256SUMS.json')
checksums = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted(new_files)}
(evidence / 'SHA256SUMS.json').write_text(json.dumps(checksums, indent=2) + '\n')
print(json.dumps({k: {'status': v['status'], 'passed': v['passed']} for k, v in summary.items()}))
