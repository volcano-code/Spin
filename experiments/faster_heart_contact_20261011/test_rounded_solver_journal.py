"""Small actual CGAL mapping tests, not rounded Heart optimization results."""
import importlib.util
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import trimesh
from make_it_stand.spin import SpinItFasterRequest, SpinSolverConfig, SpinEvaluatorConfig


class DurableMappingTests(unittest.TestCase):
    def implementation(self):
        path = Path(__file__).with_name('rounded_solver_journal.py')
        self.assertTrue(path.exists(), 'durable actual full-rank adapter is missing')
        spec = importlib.util.spec_from_file_location('rounded_solver_journal', path)
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        return module

    def request(self):
        outer = trimesh.creation.box(extents=[4, 4, 4])
        domain = trimesh.creation.box(extents=[2, 2, 2])
        return SpinItFasterRequest(outer_vertices=outer.vertices.tolist(), outer_faces=outer.faces.tolist(),
            design_vertices=domain.vertices.tolist(), design_faces=domain.faces.tolist(),
            mode='SPINNING_TOP', center_of_mass_weight=1.,
            spinning_top_objective='INERTIA_AND_LOWERING',
            settings=SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1',
                maximum_restarts=1, maximum_function_evaluations_per_restart=1,
                capture_evaluation_payload=True))

    def settings(self):
        return SpinEvaluatorConfig(icosphere_subdivisions=0,
            derivative_method='DISCRETE_FULL_RANK_QUADRIC')

    def test_actual_evaluation_is_saved_once_with_matching_jacobian(self):
        module = self.implementation()
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)/'mapping'
            with module.DurableFullRankMap(self.request(), self.settings(), folder) as mapping:
                response = mapping(np.array([-.25, 0, 0, 0, 1, 1.2, 1.4, 0, 0, 0]))
            saved = json.loads((folder/'evaluation-000001.json').read_text())
            self.assertEqual(saved['status'], 'RETURNED_ACTUAL_MAPPING')
            self.assertEqual(saved['evaluation']['total_raw_moments'], list(response['raw_moments']))
            np.testing.assert_array_equal(saved['evaluation']['coefficient_jacobian'], response['moment_jacobian'])
            self.assertEqual(len(list(folder.glob('evaluation-*.json'))), 1)
            self.assertGreater(len(saved['evaluation']['material_faces']), 0)
            self.assertEqual(saved['evaluation']['geometry_hash'], response['geometry_hash'])

    def test_native_failure_retains_started_prefix_and_exception(self):
        module = self.implementation()
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)/'mapping'
            with module.DurableFullRankMap(self.request(), self.settings(), folder) as mapping:
                with self.assertRaises((ValueError, RuntimeError)):
                    mapping(np.zeros(10))
            rows = [json.loads(x) for x in (folder/'mapping-events.jsonl').read_text().splitlines()]
            self.assertEqual(rows[0]['status'], 'STARTED')
            self.assertEqual(rows[-1]['status'], 'FAILED')
            self.assertTrue(rows[-1]['exception_type'])
            self.assertEqual(rows[-1]['coefficients'], [0.]*10)

    def test_config_binds_the_actually_loaded_native_module(self):
        module = self.implementation()
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)/'mapping'
            with module.DurableFullRankMap(self.request(), self.settings(), folder):
                pass
            config = json.loads((folder/'mapping-config.json').read_text())
            self.assertIn('native_module', config, 'loaded binary binding is missing')
            actual = config['native_module']
            self.assertEqual(actual['sha256'], hashlib.sha256(Path(actual['path']).read_bytes()).hexdigest())

    def test_actual_newton_budget_events_are_durable_and_not_sealed(self):
        module = self.implementation()
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)/'run'
            result = module.run_durable_fullrank(self.request(), self.settings(), folder)
            receipt = json.loads((folder/'numerical-result.json').read_text())
            rows = [json.loads(x) for x in (folder/'newton-events.jsonl').read_text().splitlines()]
            self.assertEqual(rows, module.json_value(result['events']))
            self.assertEqual(receipt['result'], module.json_value(result))
            self.assertEqual(receipt['material_acceptance'], 'NOT_CHECKED')
            self.assertFalse((folder/'sealed.json').exists())
            self.assertLessEqual(result['function_evaluations'], 1)
            self.assertEqual(receipt['initialization'], 'REQUEST_SEED_ONLY_NO_EXPLICIT_QUADRICS')
            config = json.loads((folder/'mapping-config.json').read_text())
            self.assertEqual(config['request']['spinning_top_objective'], 'INERTIA_AND_LOWERING')
            self.assertEqual(config['request']['center_of_mass_weight'], 1.)


if __name__ == '__main__':
    unittest.main()
