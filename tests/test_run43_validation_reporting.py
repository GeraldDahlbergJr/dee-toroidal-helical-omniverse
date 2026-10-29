"""Reporting tests only; fixture values do not represent physical Run43 results."""
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'reproduce'))
import validate_run43_final_boundaries as validation
import solve_toroidal_direct_hessian_candidate as candidate


def reporting_rows():
    return [{'resolution': [n+1, n, n], 'independent_audit': {
        region: {metric: {'volume_weighted_rms': 1/n**2}
                 for metric in validation.METRICS + (validation.SUPPLEMENTARY,)}
        for region in validation.REGIONS}} for n in validation.LEVELS]


class BoundaryReporting(unittest.TestCase):
    def test_all_twelve_pairs_and_supplementary_scope(self):
        rows = reporting_rows()
        for row in rows:
            for region in validation.REGIONS:
                row['independent_audit'][region][validation.SUPPLEMENTARY]['volume_weighted_rms'] = 1.
        orders, gates = validation.boundary_orders(rows)
        self.assertTrue(all(gates.values()))
        self.assertEqual(sum(len(orders[r][m]) for r in validation.REGIONS
                             for m in validation.METRICS), 12)

    def test_coarsest_inner_failure_cannot_be_dropped(self):
        rows = reporting_rows()
        rows[0]['independent_audit'][validation.REGIONS[0]][validation.METRICS[0]]['volume_weighted_rms'] = 1/16**2
        _, gates = validation.boundary_orders(rows)
        self.assertFalse(gates['inner_boundary_band'])
        self.assertTrue(gates['outer_boundary_band'])

    def test_missing_duplicate_or_reordered_grid_is_rejected(self):
        rows = reporting_rows()
        for invalid in (rows[:-1], rows[1:]+rows[:1], rows[:3]+rows[2:3]):
            with self.subTest(resolutions=[r['resolution'] for r in invalid]):
                with self.assertRaises(ValueError):
                    validation.boundary_orders(invalid)

    def test_nonfinite_or_nonpositive_norm_cannot_pass(self):
        for value in (float('nan'), float('inf'), 0., -1.):
            rows = reporting_rows()
            rows[0]['independent_audit'][validation.REGIONS[0]][validation.METRICS[0]]['volume_weighted_rms'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                validation.boundary_orders(rows)

    def test_checksum_mismatch_rejected_before_loading_fields(self):
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            validation.audit_field(b'corrupt evidence', (9, 8, 8), '0'*64)

    def test_another_archive_cannot_inherit_the_original_run_identity(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'other.zip'
            path.write_bytes(b'not the original artifact')
            with self.assertRaisesRegex(ValueError, 'Archive checksum'):
                validation.audit_preserved_fields(path, Path(temp)/'missing.json')

    def test_wrong_shape_and_nonfinite_fields_are_rejected(self):
        for psi in (np.ones((2, 2, 2)), np.full((9, 8, 8), np.nan)):
            buffer = io.BytesIO()
            np.savez_compressed(buffer, psi=psi, W_cartesian=np.zeros((3,)+psi.shape))
            raw = buffer.getvalue()
            with self.assertRaises(ValueError):
                validation.audit_field(raw, (9, 8, 8), hashlib.sha256(raw).hexdigest())

    def test_failed_gate_writes_report_and_returns_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'report.json'
            result = {'status': 'BOUNDARY_GATE_NOT_MET', 'checks_passed': 8, 'checks_total': 12,
                      'orders': {}, 'boundary_gate': {'inner_boundary_band': False,
                                                     'outer_boundary_band': True}}
            with patch.object(validation, 'audit_preserved_fields', return_value=result), patch('sys.stdout', new=io.StringIO()):
                code = validation.main(['--artifact-zip', str(Path(temp)/'stage.zip'), '--output', str(path)])
            self.assertEqual(code, 1)
            self.assertEqual(json.loads(path.read_text())['checks_passed'], 8)

    def test_output_cannot_overwrite_input(self):
        with tempfile.TemporaryDirectory() as temp, patch('sys.stderr', new=io.StringIO()):
            path = str(Path(temp)/'stage.zip')
            with self.assertRaises(SystemExit) as caught:
                validation.main(['--artifact-zip', path, '--output', path])
            self.assertEqual(caught.exception.code, 2)


class CandidateJSONReporting(unittest.TestCase):
    def test_single_and_multiple_resolution_outputs_are_strict_json(self):
        # Mock solves exercise only output serialization in isolated temporary files.
        # The fixture arrays and norms are never used for physical validation.
        class GridFixture:
            def __init__(self, nr, nt, nph):
                self.shape = (nr, nt, nph)
                self.size = nr*nt*nph

            def solve_restartable(self, *args):
                return np.ones(self.size), np.zeros((3, self.size)), 1, True

            def diagnostics(self, *args):
                return {'expanded_momentum_solver_max_abs': 0.}

        class AuditFixture:
            def __init__(self, shape):
                self.n = shape[1]

            def summarize(self, *args):
                return {region: {metric: {'volume_weighted_rms': 1/self.n**2}
                                 for metric in validation.METRICS + (validation.SUPPLEMENTARY,)}
                        for region in validation.REGIONS + ('bulk',)}

        for levels in ('8', '8,16'):
            with self.subTest(levels=levels), tempfile.TemporaryDirectory() as temp:
                source = Path(temp)/'candidate.py'
                source.write_text('serialization fixture\n')
                with patch.object(candidate, '__file__', str(source)), patch.object(candidate, 'DirectHessianGrid', GridFixture), patch.object(candidate, 'Audit', AuditFixture), patch.dict('os.environ', DEE_CANDIDATE_LEVELS=levels), patch('sys.stdout', new=io.StringIO()):
                    candidate.main()
                for name in ('toroidal_direct_hessian_candidate.json', 'toroidal_direct_hessian_candidate_progress.json'):
                    raw = (Path(temp)/name).read_bytes()
                    self.assertTrue(raw.endswith(b'\n'))
                    self.assertFalse(raw.endswith(b'\\n'))
                    json.loads(raw)


if __name__ == '__main__':
    unittest.main()
