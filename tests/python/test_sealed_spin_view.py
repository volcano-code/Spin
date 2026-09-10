"""Scene admission uses real independent integration and exact geometry."""
import importlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from test_exact_material_region import box, combine, l_prism
from test_sealed_spin_artifact import SpinSealTests, ORDER


class SceneFixture(unittest.TestCase):
    problem, claim, seal = SpinSealTests.problem, SpinSealTests.claim, SpinSealTests.seal

    def setUp(self):
        self.a = importlib.import_module('make_it_stand.verification.spin_artifact')
        self.n = importlib.import_module('_make_it_stand')
        self.b = importlib.import_module('make_it_stand.verification.spin_bundle')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def fixture(self, empty=False):
        outer = box(); outer[0][0, 0] = -0.
        design = box((.125,)*3, (.875,)*3)
        fixed = combine(outer, (design[0], design[1][:, ::-1]))
        cavity = None if empty else box((.25,)*3, (.75,)*3)
        retained = design if empty else combine(design, (cavity[0], cavity[1][:, ::-1]))
        material = outer if empty else combine(outer, (cavity[0], cavity[1][:, ::-1]))
        p = self.problem(outer, fixed=[('shell', fixed, 2.)], design=[('d', design, 3.)])
        seal = self.seal(p, outer, retained, cavity, material, cavity,
                         self.claim(p, [(fixed, 2.), (retained, 3.)]))
        bundle = self.root/'bundle'
        self.b.save_spin_artifact(seal, bundle)
        return bundle, p, seal


class SealedSpinViewTests(SceneFixture):
    def load(self, path, digest):
        name = 'make_it_stand.verification.spin_view'
        self.assertIsNotNone(importlib.util.find_spec(name), 'sealed scene admission is missing')
        return importlib.import_module(name).load_spin_scene(path, expected_problem_digest=digest)

    def test_actual_layers_preserve_signed_boundary_bytes_and_densities(self):
        path, p, seal = self.fixture()
        scene = self.load(path, p.input_digest)
        meshes = {m['key']: m for m in scene['meshes']}
        expected = {'original': p.original, 'outer': seal.outer,
                    'fixed:shell': p.fixed_regions[0].mesh, 'design:d': p.design_regions[0].mesh,
                    'retained:d': seal.design_results[0][1], 'material': seal.final_material,
                    'cavity': seal.newly_carved}
        self.assertEqual(set(meshes), set(expected))
        for key, mesh in expected.items():
            v, f = mesh.arrays()
            self.assertEqual(meshes[key]['vertices'].tobytes(), v.tobytes())
            self.assertEqual(meshes[key]['faces'].tobytes(), f.tobytes())
            self.assertEqual(meshes[key]['identity'], mesh.identity)
        self.assertEqual(meshes['fixed:shell']['density'], 2.)
        self.assertEqual(meshes['retained:d']['density'], 3.)
        self.assertEqual(meshes['fixed:shell']['faces'].shape, (24, 3))
        self.assertEqual(meshes['outer']['faces'].shape, (12, 3))
        self.assertAlmostEqual(scene['physical']['mass'], 131/64)
        self.assertEqual(scene['artifact_digest'], seal.scientific_digest)
        self.assertEqual(scene['geometry_status'], 'MESH_MATERIAL_VALIDATED')
        meshes['material']['vertices'][:] = 17
        scene['physical']['center_of_mass'][0] = 99
        again = self.load(path, p.input_digest)
        np.testing.assert_allclose(again['physical']['center_of_mass'], [.5]*3)
        self.assertEqual(again['meshes'][1]['vertices'].tobytes(), seal.outer.arrays()[0].tobytes())

    def test_empty_cavity_has_no_dummy_mesh(self):
        path, p, _ = self.fixture(empty=True)
        cavity = next(m for m in self.load(path, p.input_digest)['meshes'] if m['key'] == 'cavity')
        self.assertTrue(cavity['identity']['explicit_empty'])
        self.assertEqual(cavity['vertices'].shape, (0, 3))
        self.assertEqual(cavity['faces'].shape, (0, 3))

    def test_repeated_eigenvalues_suppress_unique_axes(self):
        path, p, _ = self.fixture()
        axes = self.load(path, p.input_digest)['principal_axes']
        self.assertEqual(axes['unique_axes'], [])
        self.assertEqual(axes['multiplicities'], [3])

    def test_full_nondiagonal_world_inertia_defines_axes(self):
        o = l_prism(); p = self.problem(o)
        s = self.seal(p, o, o, None, o, None)
        path = self.root/'l'; self.b.save_spin_artifact(s, path)
        scene = self.load(path, p.input_digest)
        np.testing.assert_allclose(scene['physical']['central_inertia'],
            [[121/30, 9/5, 0], [9/5, 121/30, 0], [0, 0, 217/30]], atol=1e-10)
        axes = scene['principal_axes']
        np.testing.assert_allclose(axes['eigenvalues'], [67/30, 175/30, 217/30], atol=1e-10)
        self.assertEqual(axes['multiplicities'], [1, 1, 1])
        np.testing.assert_allclose(np.abs(axes['unique_axes'][0]['direction']),
                                   [2**-.5, 2**-.5, 0], atol=1e-10)

    def test_two_dimensional_eigenspace_keeps_only_unique_z(self):
        o = box(hi=(1., 1., 2.)); p = self.problem(o)
        s = self.seal(p, o, o, None, o, None)
        path = self.root/'tall'; self.b.save_spin_artifact(s, path)
        axes = self.load(path, p.input_digest)['principal_axes']
        self.assertEqual(axes['multiplicities'], [1, 2])
        self.assertEqual(len(axes['unique_axes']), 1)
        np.testing.assert_allclose(np.abs(axes['unique_axes'][0]['direction']), [0, 0, 1])

    def test_small_model_scale_preserves_distinct_nondiagonal_axes(self):
        v, f = l_prism(); o = (v*.01, f); p = self.problem(o, wall=.0025)
        # Uniform scaling cancels from input-normalized moments. The generic L
        # fixture claim helper is hardcoded to the unit L and cannot scale it.
        claim = {'frame': 'INPUT_NORMALIZED_V1', 'order': ORDER,
                 'moments': [5/27, -2/81, -2/81, 0, 53/2916, 53/2916, 5/2916, 0, 0, -1/243]}
        s = self.seal(p, o, o, None, o, None, claim)
        path = self.root/'small'; self.b.save_spin_artifact(s, path)
        axes = self.load(path, p.input_digest)['principal_axes']
        self.assertEqual(axes['multiplicities'], [1, 1, 1])
        np.testing.assert_allclose(axes['eigenvalues'], np.array([67/30, 175/30, 217/30])*1e-10,
                                   atol=1e-20, rtol=0)

    def test_no_support_square_or_solver_success_from_geometry(self):
        path, p, _ = self.fixture()
        scene = self.load(path, p.input_digest)
        self.assertEqual(scene['target']['support_points'], [])
        self.assertEqual(scene['support_display'], 'TARGET_POINTS_ONLY_NO_CERTIFIED_HULL')
        self.assertEqual(scene['disclosures']['project_status'], 'PARTIAL_BLOCKED')
        self.assertEqual(scene['disclosures']['solver_target_feasibility'], 'UNKNOWN')
        self.assertEqual(scene['disclosures']['numerical_solver_trace_convergence'],
                         'NOT_AUTHENTICATED_BY_GEOMETRY_BUNDLE')

    def test_actual_target_point_and_nonunit_axis_are_preserved(self):
        o = box()
        target = {'mode': 'STATIC_BALANCE', 'axis': [0., 0., 2.], 'origin': [.2, .3, 0.],
                  'support_points': [[.2, .3, 0.]]}
        p = self.a.SpinArtifactProblem(original=o, fixed_regions=[], design_regions=[('d', o, 1.)],
            length_unit='m', mass_unit='kg', minimum_outer_wall=.25, target=target)
        s = self.seal(p, o, o, None, o, None)
        path = self.root/'support'; self.b.save_spin_artifact(s, path)
        scene = self.load(path, p.input_digest)
        self.assertEqual(scene['target'], target)
        self.assertEqual(scene['support_display'], 'TARGET_POINTS_ONLY_NO_CERTIFIED_HULL')

    def test_wrong_binding_mutation_staleness_and_fake_passed_rejected(self):
        path, p, _ = self.fixture()
        with self.assertRaises(ValueError): self.load(path, '0'*64)
        mod = importlib.import_module('make_it_stand.verification.spin_view')
        with self.assertRaises(TypeError): mod.load_spin_scene(path)
        original = (path/'manifest.json').read_text()
        for change in ['passed', 'stale', 'scope']:
            manifest = json.loads(original)
            if change == 'passed': manifest['geometry_sealed'] = True
            elif change == 'stale': manifest['evidence']['files']['native']['sha256'] = '0'*64
            else: manifest['validation_scope'] = 'HYDROSTATIC'
            (path/'manifest.json').write_text(json.dumps(manifest))
            with self.subTest(change=change), self.assertRaises(ValueError): self.load(path, p.input_digest)
        (path/'manifest.json').write_text(original)
        payload = next(path.glob('*.npy')); payload.write_bytes(payload.read_bytes()+b'altered')
        with self.assertRaises(ValueError): self.load(path, p.input_digest)


class NativeSealedSpinViewTests(SceneFixture):
    def run_viewer(self, *args, env=None):
        viewer = os.environ['MAKE_IT_STAND_SEALED_VIEWER']
        self.assertTrue(Path(viewer).is_file(), 'native sealed viewer entry is unavailable')
        env = dict(os.environ if env is None else env)
        env.pop('LD_PRELOAD', None)  # The embedded executable must start independently.
        env['DISPLAY'] = ''
        return subprocess.run([viewer, *args], capture_output=True, text=True, timeout=120, env=env)

    def test_native_same_admission_conversion_and_failures_without_display(self):
        path, p, seal = self.fixture()
        before = {f.name: f.read_bytes() for f in path.iterdir()}
        env = dict(os.environ, DISPLAY='', PYTHONPATH='/deliberately/unrelated')
        result = self.run_viewer('--check-sealed-spin-bundle', str(path), p.input_digest, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(result.stdout.startswith('{'), 'native check printed non-summary output: ' + result.stdout[:100])
        summary = json.loads(result.stdout)
        self.assertEqual(summary['artifact_digest'], seal.scientific_digest)
        self.assertEqual(summary['geometry_status'], 'MESH_MATERIAL_VALIDATED')
        self.assertEqual(summary['layers']['fixed:shell']['faces'], 24)
        self.assertEqual(summary['layers']['outer']['faces'], 12)
        self.assertEqual(summary['layers']['cavity']['faces'], 12)
        self.assertEqual(summary['layers']['material']['faces'], 24)
        self.assertEqual(summary['layers']['retained:d']['faces'], 24)
        self.assertEqual(summary['layers']['design:d']['faces'], 12)
        self.assertEqual(summary['layers']['original']['vertices_sha256'], p.original.identity['vertices']['sha256'])
        self.assertEqual(summary['disclosures']['solver_target_feasibility'], 'UNKNOWN')
        self.assertEqual(Path(summary['loaded_paths']['native']).parent,
                         Path(os.environ['MAKE_IT_STAND_SEALED_VIEWER']).parent/'python')
        self.assertEqual({f.name: f.read_bytes() for f in path.iterdir()}, before)
        for args in [('--check-sealed-spin-bundle', str(path), '0'*64),
                     ('--sealed-spin-bundle', str(path)), ('--check-sealed-spin-bundle',),
                     ('--sealed-spin-bundle', str(path), p.input_digest, '--frames', '-1'),
                     ('--sealed-spin-bundle', str(path), p.input_digest, '--frames', '1x'),
                     ('--check-sealed-spin-bundle', str(path), p.input_digest, '--passed')]:
            with self.subTest(args=args):
                bad = self.run_viewer(*args, env=env)
                self.assertNotEqual(bad.returncode, 0)
                self.assertIn('sealed-spin', bad.stderr)
        manifest = json.loads((path/'manifest.json').read_text()); manifest['passed'] = True
        (path/'manifest.json').write_text(json.dumps(manifest))
        bad = self.run_viewer('--check-sealed-spin-bundle', str(path), p.input_digest, env=env)
        self.assertNotEqual(bad.returncode, 0)

    def test_native_unique_axis_and_actual_support_point_conversion(self):
        o = box(hi=(1., 1., 2.))
        target = {'mode': 'STATIC_BALANCE', 'axis': [0., 0., 2.], 'origin': [.2, .3, 0.],
                  'support_points': [[.2, .3, 0.]]}
        p = self.a.SpinArtifactProblem(original=o, fixed_regions=[], design_regions=[('d', o, 1.)],
            length_unit='m', mass_unit='kg', minimum_outer_wall=.25, target=target)
        s = self.seal(p, o, o, None, o, None)
        path = self.root/'tall'; self.b.save_spin_artifact(s, path)
        result = self.run_viewer('--check-sealed-spin-bundle', str(path), p.input_digest)
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads(result.stdout)
        self.assertEqual(summary['target'], target)
        self.assertEqual(summary['principal_axes']['multiplicities'], [1, 2])
        self.assertEqual(len(summary['principal_axes']['unique_axes']), 1)

    def test_native_empty_cavity_and_stale_or_mutated_payload(self):
        path, p, _ = self.fixture(empty=True)
        result = self.run_viewer('--check-sealed-spin-bundle', str(path), p.input_digest)
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads(result.stdout)
        self.assertEqual(summary['layers']['cavity']['faces'], 0)
        self.assertEqual(summary['layers']['cavity']['vertices'], 0)
        original = (path/'manifest.json').read_text()
        manifest = json.loads(original)
        manifest['evidence']['files']['native']['sha256'] = '0'*64
        (path/'manifest.json').write_text(json.dumps(manifest))
        bad = self.run_viewer('--check-sealed-spin-bundle', str(path), p.input_digest)
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn('stale implementation', bad.stderr)
        (path/'manifest.json').write_text(original)
        payload = next(path.glob('*.npy')); payload.write_bytes(payload.read_bytes()+b'changed')
        bad = self.run_viewer('--check-sealed-spin-bundle', str(path), p.input_digest)
        self.assertNotEqual(bad.returncode, 0)


if __name__ == '__main__': unittest.main()
