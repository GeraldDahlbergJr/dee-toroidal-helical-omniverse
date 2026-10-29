"""Audit preserved Run43 fields; a failed boundary gate exits unsuccessfully.

No solver is executed and no input evidence is overwritten. Passing this
boundary gate alone does not establish source equivalence or physical validity.
"""
import argparse
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import zipfile

import numpy as np
from diagnose_toroidal_boundary_convergence import Audit

REGIONS = ('inner_boundary_band', 'outer_boundary_band')
METRICS = ('independent_flux_H', 'independent_tensor_M_order2')
SUPPLEMENTARY = 'independent_tensor_M_order4'
LEVELS = (8, 16, 32, 64)
SOURCE_RUN_ID = 38463009735
SOURCE_ARTIFACT_ID = 11847550871
ARCHIVE_SHA256 = '65606d49c0f9467d40e5ad8ea4856ae572adad24f986a71abd2d7b543101bf02'


def boundary_orders(rows):
    """Apply the original RMS-residual reduction gate to all twelve pairs.

    These are ratios of independently evaluated constraint residual norms,
    not Richardson orders from differences of three solution fields.
    """
    if [r['resolution'] for r in rows] != [[n+1, n, n] for n in LEVELS]:
        raise ValueError('All four comparable resolutions 9/8 through 65/64 are required')
    orders, gates = {}, {}
    for region in REGIONS:
        orders[region] = {}
        for metric in METRICS + (SUPPLEMENTARY,):
            values = [r['independent_audit'][region][metric]['volume_weighted_rms']
                      for r in rows]
            if any(not math.isfinite(v) or v <= 0 for v in values):
                raise ValueError(f'Undefined residual reduction order: {region}/{metric}')
            orders[region][metric] = [math.log2(a/b) for a, b in zip(values, values[1:])]
        gates[region] = all(math.isfinite(p) and p >= 1.5
                            for metric in METRICS for p in orders[region][metric])
    return orders, gates


def audit_field(raw, shape, expected_sha256):
    digest = hashlib.sha256(raw).hexdigest()
    if digest != expected_sha256:
        raise ValueError(f'Field checksum mismatch for {shape}')
    with np.load(io.BytesIO(raw), allow_pickle=False) as fields:
        psi, w = fields['psi'], fields['W_cartesian']
        if psi.shape != shape or w.shape != (3,) + shape:
            raise ValueError(f'Unexpected field shape for {shape}')
        if not (np.isfinite(psi).all() and np.isfinite(w).all() and np.min(psi) > 0):
            raise ValueError(f'Invalid physical fields for {shape}')
        return {'resolution': list(shape), 'independent_audit': Audit(shape).summarize(psi, w),
                'field_sha256': digest}


def audit_preserved_fields(artifact_path, coarse_path):
    archive_sha256 = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
    if archive_sha256 != ARCHIVE_SHA256:
        raise ValueError('Archive checksum does not match the preserved Run43 source artifact')
    coarse = json.loads(coarse_path.read_text())
    if [r['resolution'] for r in coarse['rows']] != [[n+1, n, n] for n in LEVELS[:-1]]:
        raise ValueError('The retained coarse report must identify exactly 8, 16, 32')
    rows = []
    for n in LEVELS[:-1]:
        name = f'toroidal_direct_hessian_fields_{n+1}_{n}.npz'
        rows.append(audit_field((coarse_path.parent/name).read_bytes(), (n+1, n, n),
                                coarse['raw_field_sha256'][name]))
    with zipfile.ZipFile(artifact_path) as archive:
        progress = json.loads(archive.read('toroidal_direct_hessian_candidate_progress.json'))
        if progress['status'] != 'COMPLETE' or progress['completed_resolutions'] != [[65, 64, 64]]:
            raise ValueError('The preserved fine-grid stage is not complete')
        name = 'toroidal_direct_hessian_fields_65_64.npz'
        rows.append(audit_field(archive.read(name), (65, 64, 64),
                                progress['raw_field_sha256'][name]))
    orders, gates = boundary_orders(rows)
    passed = sum(p >= 1.5 for region in REGIONS for metric in METRICS
                 for p in orders[region][metric])
    return {
        'status': 'BOUNDARY_GATE_PASS_VALIDATION_PENDING' if all(gates.values()) else 'BOUNDARY_GATE_NOT_MET',
        'source_run_id': SOURCE_RUN_ID,
        'source_artifact_id': SOURCE_ARTIFACT_ID,
        'source_progress_status': progress['status'],
        'artifact_sha256': archive_sha256,
        'coarse_report_sha256': hashlib.sha256(coarse_path.read_bytes()).hexdigest(),
        'audit_sha256': hashlib.sha256(Path(__file__).with_name('diagnose_toroidal_boundary_convergence.py').read_bytes()).hexdigest(),
        'candidate_source_sha256': {'coarse': coarse['candidate_sha256'], 'fine': progress['candidate_sha256']},
        'source_equivalence': 'REQUIRES_REVIEW',
        'criterion': 'Independent H and second-order M residual RMS reduction order >=1.5 on every refinement pair in both boundary bands',
        'checks_passed': passed, 'checks_total': 12,
        'rows': rows, 'orders': orders, 'boundary_gate': gates,
        'limitations': [
            'All four resolutions are re-audited from checksum-verified saved fields; no solver was run.',
            'Recorded coarse and fine candidate source hashes require equivalence review before joint validation.',
            'Residual RMS reduction is distinct from solution-difference Richardson convergence.',
            'Passing these tests is necessary but not sufficient for independent scientific validation.',
            'Endpoint compatibility and full continuum verification remain separate.',
        ],
    }


def main(argv=None):
    root = Path(__file__).parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact-zip', type=Path, help='Use an existing archive without downloading')
    parser.add_argument('--coarse-report', type=Path, default=root/'toroidal_direct_hessian_candidate.json')
    parser.add_argument('--output', type=Path, default=root/'run43_boundary_validation_20261021.json')
    args = parser.parse_args(argv)
    artifact = int(os.environ.get('RUN43_ARTIFACT_ID', str(SOURCE_ARTIFACT_ID)))
    artifact_path = args.artifact_zip
    if artifact_path is None:
        artifact_path = root/'run43-final-stage.zip'
        with artifact_path.open('wb') as stream:
            subprocess.run(['gh', 'api', f"repos/{os.environ['GITHUB_REPOSITORY']}/actions/artifacts/{artifact}/zip"],
                           stdout=stream, check=True)
    protected = [artifact_path, args.coarse_report]
    protected += [args.coarse_report.parent/f'toroidal_direct_hessian_fields_{n+1}_{n}.npz'
                  for n in LEVELS[:-1]]
    if args.output.resolve() in {p.resolve() for p in protected}:
        parser.error('--output must not overwrite preserved evidence')
    out = audit_preserved_fields(artifact_path, args.coarse_report)
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: out[k] for k in ('status', 'checks_passed', 'checks_total', 'orders', 'boundary_gate')}, indent=2), flush=True)
    return 0 if all(out['boundary_gate'].values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
