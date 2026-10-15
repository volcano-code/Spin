"""Detect lossy request inputs, not geometry or solver acceptance."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest


class RationalRequestViewTests(unittest.TestCase):
    def checker(self):
        import rational_request_view
        return rational_request_view

    def fixture(self, root, vertices):
        path = root/'region.json'
        path.write_text(json.dumps(dict(schema='authoritative-rational-input-boundary-v1',
            vertices=vertices, faces=[[0, 0, 0]])))
        return path, hashlib.sha256(path.read_bytes()).hexdigest()

    def test_thirds_block_authority_reuse(self):
        with tempfile.TemporaryDirectory() as directory:
            path, sha = self.fixture(Path(directory), [['1/3', '1/2', '-2/3']])
            report = self.checker().check_region(path, sha, chunk_size=7)
            self.assertEqual(report['non_binary64_coordinates'], 2)
            self.assertEqual(report['status'], 'NEEDS_BOUNDED_DERIVED_REQUEST_CERTIFICATION')
            self.assertFalse(report['solver_admission_granted'])
            self.assertEqual(report['first_lossy_coordinate']['vertex'], 0)

    def test_lossless_coordinates_do_not_grant_geometry_admission(self):
        with tempfile.TemporaryDirectory() as directory:
            path, sha = self.fixture(Path(directory), [['1/2', '-7/8', '0/1']])
            report = self.checker().check_region(path, sha, chunk_size=3)
            self.assertEqual(report['non_binary64_coordinates'], 0)
            self.assertEqual(report['status'], 'LOSSLESS_BINARY64_COORDINATE_VIEW_ONLY')
            self.assertFalse(report['solver_admission_granted'])
            self.assertFalse(report['geometry_admission_granted'])

    def test_stale_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path, sha = self.fixture(Path(directory), [['0/1']*3])
            path.write_text(path.read_text()+' ')
            with self.assertRaisesRegex(ValueError, 'SOURCE_IDENTITY'):
                self.checker().check_region(path, sha)

    def test_real_boundary_metadata_is_accepted_without_claiming_obj_validity(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = self.fixture(Path(directory), [['0/1']*3])
            record = json.loads(path.read_text())
            record['obj_role'] = 'FLOATING_POINT_DERIVATIVE_NOT_AUTHORITY'
            path.write_text(json.dumps(record))
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            report = self.checker().check_region(path, sha, chunk_size=4)
            self.assertEqual(report['vertex_count'], 1)
            self.assertFalse(report['geometry_admission_granted'])

    def test_overflow_is_loss_not_success(self):
        with tempfile.TemporaryDirectory() as directory:
            path, sha = self.fixture(Path(directory), [[str(10**400)+'/1', '0/1', '0/1']])
            report = self.checker().check_region(path, sha, chunk_size=11)
            self.assertEqual(report['non_binary64_coordinates'], 1)
            self.assertEqual(report['unrepresentable_coordinates'], 1)

    def test_noncanonical_or_malformed_coordinate_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path, sha = self.fixture(Path(directory), [['2/4', '0/1', '0/1']])
            with self.assertRaisesRegex(ValueError, 'CANONICAL'):
                self.checker().check_region(path, sha)

    def test_consumer_gate_raises_on_lossy_admission(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path, sha = self.fixture(root, [['1/3', '0/1', '0/1']])
            for name in ('outer', 'fixed', 'design'):
                (root/(name+'.rational.json')).write_bytes(path.read_bytes())
            (root/'admission.json').write_text(json.dumps(dict(
                problem_digest='fixture', status='INPUT_ADMITTED_FEASIBILITY_UNKNOWN',
                files={name+'.rational.json':sha for name in ('outer', 'fixed', 'design')})))
            with self.assertRaisesRegex(ValueError, 'NEEDS_BOUNDED_DERIVED_REQUEST_CERTIFICATION'):
                self.checker().require_lossless_admission_view(root)

    def test_truncated_document_cannot_count_as_complete_view(self):
        with tempfile.TemporaryDirectory() as directory:
            path, _ = self.fixture(Path(directory), [['0/1']*3])
            path.write_text(path.read_text()[:-1])
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            with self.assertRaises(ValueError):
                self.checker().check_region(path, sha, chunk_size=5)
            path, sha = self.fixture(Path(directory), [['0/1']*2])
            with self.assertRaisesRegex(ValueError, 'VERTEX'):
                self.checker().check_region(path, sha)


if __name__ == '__main__': unittest.main()
