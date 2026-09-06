"""Adversarial local bundle reads must recheck actual exact-native geometry."""
import hashlib
import importlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from test_exact_material_region import box
import test_sealed_spin_artifact as fixtures


class SpinBundleTests(unittest.TestCase):
    problem, claim, seal = fixtures.SpinSealTests.problem, fixtures.SpinSealTests.claim, fixtures.SpinSealTests.seal

    def setUp(self):
        self.a = importlib.import_module('make_it_stand.verification.spin_artifact')
        self.b = importlib.import_module('make_it_stand.verification.spin_bundle')
        self.assertTrue(hasattr(self.b, 'save_spin_artifact'), 'strict bundle writer is missing')
        self.assertTrue(hasattr(self.b, 'load_spin_artifact'), 'strict bundle reload is missing')
        self.n = importlib.import_module('_make_it_stand')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        o = box(); o[0][0, 0] = -0.
        self.p = self.problem(o)
        self.s = self.seal(self.p, o, o, None, o, None)

    def save(self, name='bundle'):
        path = self.root / name
        self.b.save_spin_artifact(self.s, path)
        return path

    def load(self, path, digest=None):
        return self.b.load_spin_artifact(path, expected_problem_digest=self.p.input_digest if digest is None else digest)

    def manifest(self, path):
        return json.loads((path/'manifest.json').read_text())

    def write_manifest(self, path, value):
        # Tampering is test-only; production never overwrites artifact files.
        (path/'manifest.json').write_text(json.dumps(value))

    def test_two_new_directories_reload_same_digest_and_original_bits(self):
        one, two = self.save('one'), self.save('two')
        a, b = self.load(one), self.load(two)
        self.assertEqual(a.scientific_digest, b.scientific_digest)
        self.assertEqual(a.evidence_digest, b.evidence_digest)
        self.assertEqual(a.outer, self.p.original)
        self.assertTrue(np.signbit(a.outer.arrays()[0][0, 0]))
        for path in one.iterdir():
            self.assertEqual(path.read_bytes(), (two/path.name).read_bytes())

    def test_expected_problem_binding_is_required(self):
        path = self.save()
        with self.assertRaises(TypeError): self.b.load_spin_artifact(path)
        with self.assertRaises(ValueError): self.load(path, '0'*64)
        changed = self.problem(wall=.125)
        with self.assertRaises(ValueError): self.load(path, changed.input_digest)

    def test_valid_and_partial_directories_are_never_overwritten(self):
        path = self.save(); before = {p.name: p.read_bytes() for p in path.iterdir()}
        with self.assertRaises(FileExistsError): self.b.save_spin_artifact(self.s, path)
        self.assertEqual(before, {p.name: p.read_bytes() for p in path.iterdir()})
        partial = self.root/'partial'; partial.mkdir(); (partial/'evidence').write_bytes(b'prior partial')
        with self.assertRaises(FileExistsError): self.b.save_spin_artifact(self.s, partial)
        self.assertEqual((partial/'evidence').read_bytes(), b'prior partial')

    def test_duplicate_nonfinite_overflow_and_unknown_manifest_fields(self):
        path = self.save(); original = (path/'manifest.json').read_text()
        cases = ['{"schema":"a","schema":"b"}', '{"x":NaN}', '{"x":Infinity}', '{"x":1e999}',
                 original[:-1]+',"passed":true}', original[:-1]+',"manualpassed":true}']
        for payload in cases:
            (path/'manifest.json').write_text(payload)
            with self.subTest(payload=payload[:35]), self.assertRaises(ValueError): self.load(path)

    def test_unknown_problem_fields_and_duplicate_keys_rejected_even_with_new_hash(self):
        path = self.save(); m = self.manifest(path)
        for payload in ['{"x":1,"x":2}', json.dumps({'passed': True})]:
            data = payload.encode(); (path/'problem.json').write_bytes(data)
            m['problem']['sha256'] = hashlib.sha256(data).hexdigest(); m['problem']['size'] = len(data)
            self.write_manifest(path, m)
            with self.assertRaises(ValueError): self.load(path)

    def test_missing_payload_and_symlink_rejected(self):
        path = self.save(); payload = next(path.glob('*.npy')); backup = payload.read_bytes()
        payload.unlink()
        with self.assertRaises((ValueError, FileNotFoundError)): self.load(path)
        external = self.root/'external.npy'; external.write_bytes(backup); payload.symlink_to(external)
        with self.assertRaises((ValueError, OSError)): self.load(path)
        link = self.root/'link'; link.symlink_to(path, target_is_directory=True)
        with self.assertRaises((ValueError, OSError)): self.load(link)

    def test_payload_traversal_and_unknown_fields_rejected(self):
        path = self.save(); m = self.manifest(path)
        for change in ['../escape.npy', '/tmp/escape.npy']:
            altered = json.loads(json.dumps(m)); altered['final']['outer']['vertices']['name'] = change
            self.write_manifest(path, altered)
            with self.assertRaises(ValueError): self.load(path)
        m['final']['outer']['vertices']['passed'] = True; self.write_manifest(path, m)
        with self.assertRaises(ValueError): self.load(path)

    def test_modified_mass_and_pasted_status_digest_do_not_bypass_recomputation(self):
        path = self.save(); m = self.manifest(path)
        m['facts']['physical']['mass'] = 123.
        m['scientific_digest'] = self.b.digest(m['facts'])
        self.write_manifest(path, m)
        with self.assertRaises(ValueError): self.load(path)

    def test_modified_geometry_with_updated_file_hash_still_fails_partition(self):
        path = self.save(); m = self.manifest(path)
        record = m['final']['material']['vertices']; payload = path/record['name']
        a = np.load(payload, allow_pickle=False); a[:, 0] *= .5
        stream = io.BytesIO(); np.lib.format.write_array(stream, a, version=(1, 0), allow_pickle=False)
        data = stream.getvalue(); payload.write_bytes(data)
        record['size'] = len(data); record['sha256'] = hashlib.sha256(data).hexdigest()
        self.write_manifest(path, m)
        with self.assertRaisesRegex(ValueError, 'partition'): self.load(path)

    def test_stale_native_and_checker_observations_rejected(self):
        path = self.save(); original = self.manifest(path)
        for name in ['native', 'mass_moments', 'spin_artifact']:
            m = json.loads(json.dumps(original)); m['evidence']['files'][name]['sha256'] = '0'*64
            m['evidence_digest'] = self.b.digest({'scientific_digest': m['scientific_digest'], 'observation': m['evidence']})
            self.write_manifest(path, m)
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'observation'): self.load(path)

    def test_target_mode_implementation_is_observed(self):
        self.assertIn('spin_modes', self.s.evidence['files'], 'target-mode validator source is unobserved')

    def test_failed_write_preserves_real_partial_files_without_valid_manifest(self):
        path = self.root/'failed'
        # Only a filesystem failure is injected; actual geometry/native checks run.
        with patch.object(self.b.os, 'fsync', side_effect=OSError('injected sync failure')):
            with self.assertRaises(OSError): self.b.save_spin_artifact(self.s, path)
        self.assertTrue(list(path.glob('*.npy')))
        self.assertFalse((path/'manifest.json').exists())
        with self.assertRaises(FileExistsError): self.b.save_spin_artifact(self.s, path)
        with self.assertRaises(FileNotFoundError): self.load(path)

    def test_npy_dtype_shape_order_object_and_trailing_payload_rejected(self):
        path = self.save(); original = self.manifest(path)
        record = original['final']['material']['vertices']
        payload = path/record['name']; original_array = np.load(payload, allow_pickle=False)
        for array, trailing in [(original_array.astype(np.float32), b''),
                                (np.asfortranarray(original_array), b''),
                                (original_array[:, :2], b''),
                                (original_array.astype(object), b''), (original_array, b'extra')]:
            m = json.loads(json.dumps(original)); rec = m['final']['material']['vertices']
            stream = io.BytesIO(); np.lib.format.write_array(stream, array, version=(1,0), allow_pickle=True)
            data = stream.getvalue()+trailing; payload.write_bytes(data)
            rec['size'] = len(data); rec['sha256'] = hashlib.sha256(data).hexdigest()
            self.write_manifest(path, m)
            with self.subTest(dtype=array.dtype), self.assertRaises(ValueError): self.load(path)

    def test_unsupported_profile_and_extra_files_rejected(self):
        path = self.save()
        with self.assertRaises(ValueError): self.b.load_spin_artifact(path, expected_problem_digest=self.p.input_digest,
                                                                    validation_scope='HYDROSTATIC')
        (path/'passed.json').write_text('{"passed":true}')
        with self.assertRaises(ValueError): self.load(path)


if __name__ == '__main__': unittest.main()
