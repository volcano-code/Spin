"""Real sealed analytic meshes must not bypass incomplete material admission."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import unittest

from test_sealed_spin_view import SceneFixture

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    'g7_dynamics_subject', ROOT/'experiments/teapot_spin2014/reconstructed_g7_dynamics.py')
SUBJECT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUBJECT)


class DynamicsAdmissionTests(SceneFixture):
    def load(self, path, problem, seal):
        self.assertTrue(hasattr(SUBJECT, 'load_qualified_dynamics_bundle'),
                        'physics entry must independently enforce the complete material gate')
        return SUBJECT.load_qualified_dynamics_bundle(path,
            expected_problem_digest=problem.input_digest,
            expected_scientific_digest=seal.scientific_digest)

    def test_real_hollow_seal_unknown_feature_thickness_cannot_simulate(self):
        path, problem, seal = self.fixture()
        self.assertEqual(seal.facts['geometry_status'], 'MESH_MATERIAL_VALIDATED')
        self.assertEqual(seal.facts['wall_thickness_audit']['overall'], 'UNKNOWN')
        with self.assertRaisesRegex(ValueError, 'COMPLETE_MATERIAL_GEOMETRY_REQUIRED'):
            self.load(path, problem, seal)

    def test_empty_carving_does_not_exempt_original_same_component_features(self):
        path, problem, seal = self.fixture(empty=True)
        with self.assertRaisesRegex(ValueError, 'COMPLETE_MATERIAL_GEOMETRY_REQUIRED'):
            self.load(path, problem, seal)

    def test_hand_promoted_wall_facts_and_rehashed_json_are_rejected(self):
        from make_it_stand.verification.spin_common import digest
        path, problem, seal = self.fixture()
        file = path/'manifest.json'
        value = json.loads(file.read_text())
        forged = copy.deepcopy(value)
        forged['facts']['wall_thickness_audit']['overall'] = 'VERIFIED'
        forged['facts']['wall_thickness_audit']['same_component_feature'] = 'VERIFIED'
        forged['facts']['same_component_feature_thickness'] = 'VERIFIED'
        forged['facts']['manufacturing_admission'] = 'VERIFIED'
        forged['scientific_digest'] = digest(forged['facts'])
        forged['evidence_digest'] = digest({'scientific_digest': forged['scientific_digest'],
                                          'observation': forged['evidence']})
        file.write_text(json.dumps(forged))
        self.assertTrue(hasattr(SUBJECT, 'load_qualified_dynamics_bundle'))
        with self.assertRaisesRegex(ValueError, 'recomputed scientific facts mismatch'):
            SUBJECT.load_qualified_dynamics_bundle(path,
                expected_problem_digest=problem.input_digest,
                expected_scientific_digest=forged['scientific_digest'])

    def test_clean_worker_environment_loads_exact_sealed_native_extension(self):
        path, _, seal = self.fixture(empty=True)
        self.assertTrue(hasattr(SUBJECT, 'sealed_dynamics_runtime'),
                        'isolated worker must receive its source-bound native dependency')
        runtime = SUBJECT.sealed_dynamics_runtime(path, native_extension=self.n.__file__)
        command = [*runtime['command_prefix'], sys.executable, '-B', '-c',
                   'import _make_it_stand; print(_make_it_stand.__file__)']
        result = subprocess.run(command, capture_output=True, text=True, timeout=30,
            env={'PATH': '/usr/bin:/bin', **runtime['environment']})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(Path(result.stdout.strip()).resolve(),
                         Path(self.n.__file__).resolve())
        self.assertIn('sealed_native_extension', runtime['identity_files'])
        self.assertIn('native_cxx_runtime', runtime['identity_files'])

    def test_runtime_rejects_wrong_native_digest_before_import(self):
        path, _, _ = self.fixture(empty=True)
        manifest = path/'manifest.json'
        record = json.loads(manifest.read_text())
        record['evidence']['files']['native']['sha256'] = '0'*64
        manifest.write_text(json.dumps(record))
        self.assertTrue(hasattr(SUBJECT, 'sealed_dynamics_runtime'))
        with self.assertRaisesRegex(ValueError, 'SEALED_NATIVE_RUNTIME_MISMATCH'):
            SUBJECT.sealed_dynamics_runtime(path, native_extension=self.n.__file__)


if __name__ == '__main__':
    unittest.main()
