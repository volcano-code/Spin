"""Adversarial local bundle reads must recheck actual exact-native geometry."""
import hashlib
import importlib
import io
from fractions import Fraction
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

    def test_explicit_reaudit_rechecks_old_payload_without_upgrading_old_evidence(self):
        reaudit = getattr(self.b, 'reaudit_spin_artifact', None)
        self.assertTrue(callable(reaudit), 'saved-material reaudit entry is missing')
        path = self.save()
        manifest = self.manifest(path)
        # Simulate an old observed native hash; this is not build attestation.
        manifest['evidence']['files']['native']['sha256'] = '0'*64
        manifest['evidence_digest'] = self.b.digest({
            'scientific_digest':manifest['scientific_digest'], 'observation':manifest['evidence']})
        self.write_manifest(path, manifest)
        before = {p.name:p.read_bytes() for p in path.iterdir()}
        with self.assertRaisesRegex(ValueError, 'stale implementation'):
            self.load(path)
        fresh, receipt = reaudit(path, expected_problem_digest=self.p.input_digest,
                                expected_scientific_digest=self.s.scientific_digest)
        self.assertEqual(fresh.scientific_digest, self.s.scientific_digest)
        self.assertNotEqual(fresh.evidence_digest, manifest['evidence_digest'])
        self.assertEqual(receipt['source_evidence_digest'], manifest['evidence_digest'])
        self.assertEqual(receipt['current_evidence_digest'], fresh.evidence_digest)
        self.assertEqual(receipt['solver_freshness'], 'NOT_REASSESSED')
        self.assertEqual(receipt['status'], 'MATERIAL_REAUDIT_IDENTICAL_FACTS')
        self.assertEqual(before, {p.name:p.read_bytes() for p in path.iterdir()})
        with self.assertRaisesRegex(ValueError, 'stale implementation'):
            self.load(path)
        with self.assertRaisesRegex(ValueError, 'scientific.*binding'):
            reaudit(path, expected_problem_digest=self.p.input_digest,
                    expected_scientific_digest='0'*64)
        manifest['facts']['physical']['mass'] += 1.
        self.write_manifest(path, manifest)
        with self.assertRaisesRegex(ValueError, 'scientific facts mismatch'):
            reaudit(path, expected_problem_digest=self.p.input_digest,
                    expected_scientific_digest=self.s.scientific_digest)

    def test_reaudit_still_rejects_payload_and_historical_digest_tampering(self):
        reaudit = getattr(self.b, 'reaudit_spin_artifact', None)
        self.assertTrue(callable(reaudit), 'saved-material reaudit entry is missing')
        path = self.save(); m = self.manifest(path)
        m['evidence_digest'] = '0'*64; self.write_manifest(path, m)
        with self.assertRaisesRegex(ValueError, 'source evidence digest'):
            reaudit(path, expected_problem_digest=self.p.input_digest,
                    expected_scientific_digest=self.s.scientific_digest)
        m['evidence_digest'] = self.s.evidence_digest; self.write_manifest(path, m)
        p = next(path.glob('*.npy')); p.write_bytes(p.read_bytes()+b'bad')
        with self.assertRaises(ValueError):
            reaudit(path, expected_problem_digest=self.p.input_digest,
                    expected_scientific_digest=self.s.scientific_digest)

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


    def test_rational_material_is_authoritative_through_seal_and_bundle_reload(self):
        from make_it_stand.verification.mass_moments import integrate_exact_rational_raw_moments
        from make_it_stand.verification.spin_problem import MeshSnapshot

        # A cube whose side length is exactly 1/3 cannot be represented by
        # binary64 vertices. The rational mesh, not its display/export copy,
        # must remain the geometry used by the sealed artifact.
        unit = box()
        exact_vertices = [[('1/3' if coordinate else '0/1') for coordinate in vertex]
                          for vertex in unit[0].astype(int)]
        exact = MeshSnapshot.from_rationals(exact_vertices, unit[1])
        self.assertNotEqual(MeshSnapshot(exact.arrays()), exact)
        self.assertEqual(exact.rational_vertices, tuple(tuple(row) for row in exact_vertices))
        exact_moments = integrate_exact_rational_raw_moments(exact.rational_vertices,
                                                              exact.arrays()[1].tolist())
        self.assertEqual(exact_moments[0], Fraction(1, 27))
        self.assertEqual(exact_moments[1], Fraction(1, 162))

        p = self.problem(outer=exact, design=[('d', exact, 2.)], wall=.25)
        normalized_claim = {'frame': 'INPUT_NORMALIZED_V1', 'order': [
            'M', 'LX', 'LY', 'LZ', 'QXX', 'QYY', 'QZZ', 'QYZ', 'QXZ', 'QXY'],
            'moments': [1., 0., 0., 0., 1/12, 1/12, 1/12, 0., 0., 0.]}
        sealed = self.a.SealedSpinArtifact.create(p, outer=exact,
            design_results=[('d', exact, None)], final_material=exact,
            newly_carved=None, solver_claim=normalized_claim)
        region = sealed.facts['regions'][0]
        self.assertEqual(region['physical']['exact_rational_moments'][0], '2/27')

        path = self.root / 'exact-bundle'
        self.b.save_spin_artifact(sealed, path)
        reloaded = self.b.load_spin_artifact(path, expected_problem_digest=p.input_digest)
        self.assertEqual(reloaded.outer.rational_vertices, exact.rational_vertices)
        self.assertEqual(reloaded.scientific_digest, sealed.scientific_digest)

        tampered = self.root / 'tampered-rational-bundle'
        self.b.save_spin_artifact(sealed, tampered)
        manifest = self.manifest(tampered)
        record = manifest['final']['outer']['rational_vertices']
        payload = tampered / record['name']
        rows = json.loads(payload.read_text())
        rows[0][0] = '1/9'
        data = self.b.canonical(rows).encode('ascii')
        payload.write_bytes(data)
        record['size'] = len(data)
        record['sha256'] = hashlib.sha256(data).hexdigest()
        self.write_manifest(tampered, manifest)
        with self.assertRaisesRegex(ValueError, 'derived binary64'):
            self.b.load_spin_artifact(tampered, expected_problem_digest=p.input_digest)


    def test_rational_nonempty_carve_seals_with_exact_wall_and_bundle_reload(self):
        from make_it_stand.verification.mass_moments import integrate_exact_rational_raw_moments
        from make_it_stand.verification.spin_problem import MeshSnapshot

        _, faces = box()
        outer_rows = [[('1/1' if coordinate else '0/1') for coordinate in vertex]
                      for vertex in box()[0].astype(int)]
        cavity_rows = [[('2/3' if coordinate else '1/3') for coordinate in vertex]
                       for vertex in box()[0].astype(int)]
        outer = MeshSnapshot.from_rationals(outer_rows, faces)
        cavity = MeshSnapshot.from_rationals(cavity_rows, faces)
        retained_rows = outer_rows + cavity_rows
        retained_faces = np.vstack([faces, faces[:, ::-1] + len(outer_rows)])
        retained = MeshSnapshot.from_rationals(retained_rows, retained_faces)
        problem = self.problem(outer=outer, design=[('d', outer, 1.)], wall=.3)
        frame = problem.normalization
        origin = [Fraction(float(value)) for value in frame['origin']]
        length = Fraction(float(frame['L']))
        normalized_rows = [[str((Fraction(value) - origin[axis]) / length)
                            for axis, value in enumerate(row)] for row in retained_rows]
        moments = integrate_exact_rational_raw_moments(normalized_rows, retained_faces.tolist())
        moments = [value / Fraction(float(frame['R'])) for value in moments]
        claim = {'order': ['M', 'LX', 'LY', 'LZ', 'QXX', 'QYY', 'QZZ', 'QYZ', 'QXZ', 'QXY'],
                 'frame': 'INPUT_NORMALIZED_V1', 'moments': [float(value) for value in moments]}
        sealed = self.a.SealedSpinArtifact.create(problem, outer=outer,
            design_results=[('d', retained, cavity)], final_material=retained,
            newly_carved=cavity, solver_claim=claim)
        self.assertEqual(sealed.facts['clearance']['squared_distance_exact'], '1/9')
        self.assertTrue(sealed.facts['clearance']['at_least_minimum'])
        self.assertEqual(sealed.facts['newly_carved_components'], 1)
        path = self.root / 'rational-carve-bundle'
        self.b.save_spin_artifact(sealed, path)
        exact_constructor = self.n.make_exact_material_region_from_rationals
        with patch.object(self.n, 'make_exact_material_region_from_rationals',
                          wraps=exact_constructor) as exact_constructor_spy:
            loaded = self.b.load_spin_artifact(path, expected_problem_digest=problem.input_digest)
        # The bundle stores this exact cavity both as the top-level newly-carved
        # set and as the design-result carved set. Both serialized payloads must
        # still be read and checked, but the expensive exact solid is identical.
        self.assertEqual(exact_constructor_spy.call_count, 3)
        self.assertIs(loaded.newly_carved.exact_region(),
                      loaded.design_results[0][2].exact_region())
        self.assertEqual(loaded.outer.rational_vertices, outer.rational_vertices)
        self.assertEqual(loaded.newly_carved.rational_vertices, cavity.rational_vertices)
        self.assertEqual(loaded.scientific_digest, sealed.scientific_digest)

    def test_legacy_v1_bundle_without_rational_sidecars_remains_readable(self):
        path = self.save('legacy-v1')
        manifest = self.manifest(path)
        manifest['schema'] = 'sealed-spin-bundle-v1'

        def strip_v2_mesh_fields(value):
            if isinstance(value, dict):
                value.pop('rational_vertices', None)
                for child in value.values():
                    strip_v2_mesh_fields(child)
            elif isinstance(value, list):
                for child in value:
                    strip_v2_mesh_fields(child)

        # problem.json is immutable; v1 mesh descriptors live in that payload,
        # so rebuild its canonical bytes and update only its integrity record.
        problem_path = path / 'problem.json'
        problem = json.loads(problem_path.read_text())
        strip_v2_mesh_fields(problem)
        problem_data = self.b.canonical(problem).encode('ascii')
        problem_path.write_bytes(problem_data)
        manifest['problem']['size'] = len(problem_data)
        manifest['problem']['sha256'] = hashlib.sha256(problem_data).hexdigest()
        strip_v2_mesh_fields(manifest['final'])
        self.write_manifest(path, manifest)
        loaded = self.load(path)
        self.assertEqual(loaded.scientific_digest, self.s.scientific_digest)


if __name__ == '__main__': unittest.main()
