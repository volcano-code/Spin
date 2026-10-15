import hashlib
import json
from pathlib import Path
import unittest
import numpy as np
import rounded_solver_journal as journal


class BackendIdentityTests(unittest.TestCase):
    def test_native_provenance_is_bound_without_newton_delimiters(self):
        encode = getattr(journal, 'newton_backend_identity', None)
        self.assertTrue(callable(encode), 'adapter must encode native provenance for Newton')
        provenance = 'libigl::copyleft::cgal::mesh_boolean;revision=f06fb130'
        identity = encode(provenance)
        self.assertNotIn(';', identity)
        self.assertEqual(identity, 'backend-sha256-' + hashlib.sha256(provenance.encode()).hexdigest())
        self.assertNotEqual(identity, encode(provenance + ',tree=changed'))
        with self.assertRaises(ValueError):
            encode('')

    def test_saved_actual_mapping_reaches_native_newton_without_boolean_recompute(self):
        from make_it_stand.spin import SpinItFasterRequest, SpinSolverConfig, solve_spin_it_faster_with_external_moment_evaluator
        here = Path(__file__).resolve().parent
        frozen = json.loads((here/'rounded-clearance-v3-faster-frozen-v2/request.json').read_text())
        saved = json.loads((here/'rounded-clearance-v3-core-faster-v2/seed-0-job-01/output/run/evaluation-000001.json').read_text())
        selected = dict(frozen['requests'][0])
        settings = dict(selected['settings'])
        settings.update(maximum_restarts=1, maximum_function_evaluations_per_restart=1)
        selected['settings'] = SpinSolverConfig(**settings)
        evaluation = saved['evaluation']
        calls = []
        def replay(q):
            np.testing.assert_array_equal(q, saved['coefficients'])
            calls.append(np.asarray(q).copy())
            return dict(raw_moments=evaluation['total_raw_moments'],
                variable_raw_moments=evaluation['variable_raw_moments'],
                derivative_defined=evaluation['derivative_defined'],
                derivative_diagnostic=evaluation['derivative_diagnostic'],
                moment_jacobian=np.asarray(evaluation['coefficient_jacobian']),
                geometry_hash=evaluation['geometry_hash'], topology_hash='saved-actual-topology',
                evaluation_map_identity=saved['mapping_identity'],
                boolean_backend_identity=journal.newton_backend_identity(evaluation['boolean_backend_identity']),
                derivative_identity=evaluation['scientific_hash'])
        result = solve_spin_it_faster_with_external_moment_evaluator(SpinItFasterRequest(**selected), replay)
        self.assertEqual(len(calls), 1, result['events'])
        self.assertFalse(any(row.get('reason') == 'evaluation_exception' for row in result['events']))


if __name__ == '__main__':
    unittest.main()
