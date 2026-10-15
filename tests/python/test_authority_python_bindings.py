import importlib
import os
import json
import copy
import tempfile
import sys
from pathlib import Path
import unittest
from fractions import Fraction
import numpy as np
import _make_it_stand as native
from test_spin2014_authority_numerical_checker import checker

api = importlib.import_module(os.environ.get('AUTHORITY_BINDING_TEST_MODULE', '_make_it_stand'))

class AuthorityBindingTests(unittest.TestCase):
    def test_saved_authority_geodesic_replay_keeps_mpfr_state(self):
        # Break caught: reconstructing an archived trial through binary64 or
        # silently renormalizing an inconsistent saved geometry/dual pair.
        import gmpy2
        initial = ['1/1']+['0/1']*9
        direction = ['0/1','1/1']+['0/1']*8
        result = api.update_spin2014_authority_coefficients(initial,initial,direction,'1/4')
        self.assertEqual(result['scope'],'COEFFICIENT_UPDATE_NOT_MATERIAL_OR_CONVERGENCE')
        self.assertEqual(result['precision_bits'],1024)
        with gmpy2.context(precision=2048):
            first = gmpy2.mpfr(gmpy2.mpq(result['dual_authority'][0]))
            second = gmpy2.mpfr(gmpy2.mpq(result['dual_authority'][1]))
            arc = gmpy2.mpfr(gmpy2.mpq(1,4))
            self.assertLess(abs(first-gmpy2.cos(arc)),gmpy2.mpfr(2)**-1000)
            self.assertLess(abs(second-gmpy2.sin(arc)),gmpy2.mpfr(2)**-1000)
        for value in result['geometry_authority']+result['dual_authority']:
            self.assertEqual(str(Fraction(value).numerator)+'/'+str(Fraction(value).denominator),value)
        self.assertEqual(api.update_spin2014_authority_coefficients(initial,initial,direction,'0/1')['geometry_authority'],initial)
        broken = initial.copy();broken[0]='2/1'
        with self.assertRaises(ValueError):
            api.update_spin2014_authority_coefficients(broken,initial,direction,'1/4')

    def test_restart_binding_keeps_both_failed_attempts(self):
        vertices = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], dtype=float)
        faces = np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[3,7,6],[3,6,2],[0,4,7],[0,7,3],[1,2,6],[1,6,5]], dtype=np.int32)
        outer, design = native.SpinTriangleMesh(), native.SpinTriangleMesh()
        outer.vertices, outer.faces = vertices * 3, faces
        design.vertices, design.faces = vertices * 2, faces
        request = native.SpinDesignRequest()
        request.outer_surface, request.design_domain = outer, design
        request.material_density = 3
        request.mode = native.SpinMode.SPINNING_TOP
        solver = request.settings
        solver.maximum_restarts = 2
        request.settings = solver
        evaluation = native.SpinEvaluatorSettings()
        tri = evaluation.triangulation
        tri.icosphere_subdivisions = 0
        tri.rational_source_coordinates = True
        tri.compact_rational_source_coordinates = True
        tri.truncation_radius = 8
        tri.maximum_radius_doublings = 0
        evaluation.triangulation = tri
        evaluation.derivative_method = native.SpinDerivativeMethod.DISCRETE_FULL_RANK_QUADRIC
        evaluation.exact_binary64_source_boundary = True
        evaluation.exact_binary64_material_boundary = True
        result = api.solve_spin2014_authority_quadric_hollowing(request,
            native.Spin2014QuadricObjectiveSettings(), 2, evaluation,
            evaluation_budget_per_attempt=1)
        self.assertEqual(result['function_evaluations'], 2)
        self.assertEqual(len(result['attempts']), 2)
        self.assertIsNone(result['stationary_attempt_index'])
        self.assertIsNone(result['stationary_point'])
        self.assertEqual(result['protocol_scope'], 'DIAGNOSTIC_REDUCED_ATTEMPT_BUDGET')
        self.assertNotEqual(result['attempts'][0]['initial_geometry_authority'], result['attempts'][1]['initial_geometry_authority'])
        self.assertEqual(json.loads(json.dumps(result, allow_nan=False))['attempts'], result['attempts'])
        expected = api.spin2014_authority_initial_sequence(request.mode,request.settings)
        checked = checker().check_authority_restart_ledger(result,
            objective_settings={'objective':'TOP_MASS_AND_LOWERING'},
            expected_initial_sequence=expected,evaluation_budget_per_attempt=1,
            expected_request_identity=result['request_identity'],
            expected_context_identity=result['context_identity'])
        self.assertEqual(checked['function_evaluations'],2)
        self.assertFalse(checked['full_run_candidate'])
        for change in ('count','missing','start'):
            changed = copy.deepcopy(result)
            if change == 'count':
                changed['function_evaluations'] += 1
            elif change == 'missing':
                changed['attempts'].pop()
                changed['function_evaluations'] -= 1
            else:
                changed['attempts'][0]['initial_geometry_authority'][0] = '0/1'
            with self.assertRaises(ValueError):
                checker().check_authority_restart_ledger(changed,
                    objective_settings={'objective':'TOP_MASS_AND_LOWERING'},
                    expected_initial_sequence=expected,evaluation_budget_per_attempt=1,
                    expected_request_identity=result['request_identity'],
                    expected_context_identity=result['context_identity'])

    def test_actual_material_point_preserves_authority_and_density(self):
        def cube(radius):
            mesh = native.SpinTriangleMesh()
            mesh.vertices = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], dtype=float) * radius
            mesh.faces = np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[3,7,6],[3,6,2],[0,4,7],[0,7,3],[1,2,6],[1,6,5]], dtype=np.int32)
            return mesh
        request = native.SpinDesignRequest()
        request.outer_surface = cube(3)
        request.design_domain = cube(2)
        request.material_density = 3
        request.mode = native.SpinMode.SPINNING_TOP
        settings = native.Spin2014QuadricObjectiveSettings()
        evaluation = native.SpinEvaluatorSettings()
        tri = evaluation.triangulation
        tri.icosphere_subdivisions = 0
        tri.rational_source_coordinates = True
        tri.compact_rational_source_coordinates = True
        tri.truncation_radius = 8
        tri.maximum_radius_doublings = 0
        evaluation.triangulation = tri
        evaluation.derivative_method = native.SpinDerivativeMethod.DISCRETE_FULL_RANK_QUADRIC
        evaluation.exact_binary64_source_boundary = True
        evaluation.exact_binary64_material_boundary = True
        authority = api.normalize_spin2014_authority_coefficients(['-1/2','1/10','0/1','0/1','1/2','13/20','2/5','0/1','0/1','0/1'])
        point = api.evaluate_spin2014_authority_quadric_point(request, settings, 2,
            authority['geometry_authority'], authority['dual_authority'], evaluation)
        mapping = api.evaluate_spin2014_authority_material_mapping(request,2,
            authority['geometry_authority'],authority['dual_authority'],evaluation)
        self.assertEqual(mapping['schema'],'spin2014-authority-material-map-v1')
        self.assertEqual(mapping['scope'],'MOMENT_MAPPING_NOT_SOLVER_OR_FINAL_SEAL')
        self.assertNotIn('raw_equation_residual',mapping)
        self.assertEqual(len(mapping['provenance_topology_hash']),64)
        independently=checker().check_independent_authority_material_mapping(mapping,
            expected_request_identity=mapping['request_identity'],
            expected_context_identity=mapping['context_identity'])
        self.assertEqual(independently['independent_total_moments_target'],mapping['total_moments_target'])
        self.assertNotIn('numerical_stationarity_qualified',independently)
        # Catch the protocol adapter accidentally selecting variable-only
        # moments or the binary64 coefficient view instead of actual authority.
        import importlib.util
        runner_path=Path(__file__).resolve().parents[2]/'.superpowers/sdd/2026-10-07-spin2014-objectives-quadric/validate_production_five_point_map.py'
        spec=importlib.util.spec_from_file_location('production_five_point_runner',runner_path)
        runner=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        observed=runner.independent_authority({'evaluation':mapping},
            expected_request_identity=mapping['request_identity'],
            expected_context_identity=mapping['context_identity'])
        self.assertEqual(observed,tuple(Fraction(x) for x in mapping['total_moments_target']))
        self.assertEqual(observed[0]-Fraction(mapping['variable_moments_target'][0]),304)
        for change in ('density','moment','affine'):
            damaged=copy.deepcopy(mapping)
            if change=='density':damaged['fixed_material_density']=4
            elif change=='moment':damaged['variable_moments_target'][0]='0/1'
            else:damaged['model_to_target_authority'][0][3]='1/1'
            with self.assertRaises(ValueError):
                checker().check_independent_authority_material_mapping(damaged,
                    expected_request_identity=mapping['request_identity'],
                    expected_context_identity=mapping['context_identity'])
            with self.assertRaises(ValueError):
                runner.independent_authority({'evaluation':damaged},
                    expected_request_identity=mapping['request_identity'],
                    expected_context_identity=mapping['context_identity'])
        for key in ('geometry_authority','dual_authority','context_identity','request_identity',
                    'material_vertices_target_authority','material_faces','material_birth_faces',
                    'source_vertices_target_authority','source_faces','source_construction_identity',
                    'fixed_moments_target','variable_moments_target','total_moments_target',
                    'native_moment_jacobian','material_scientific_hash',
                    'provenance_topology_hash','material_topology_hash',
                    'material_vertex_provenance','surface_type','effective_normal_subdivisions',
                    'source_truncation_radius','source_interface_face_count','source_cap_face_count'):
            self.assertEqual(mapping[key],point[key],key)
        self.assertEqual(point['geometry_authority'], authority['geometry_authority'])
        self.assertEqual(Fraction(point['fixed_moments_target'][0]), 304)
        self.assertEqual(point['evaluation_map_identity'], point['material_scientific_hash'])
        reloaded = json.loads(json.dumps(point, allow_nan=False))
        self.assertEqual(len(point['request_identity']), 64)
        self.assertEqual(point['outer_vertices_model_authority'][0], ['-3/1','-3/1','-3/1'])
        self.assertEqual(point['design_vertices_model_authority'][0], ['-2/1','-2/1','-2/1'])
        self.assertEqual(point['outer_faces'], request.outer_surface.faces.tolist())
        self.assertEqual(point['design_faces'], request.design_domain.faces.tolist())
        self.assertEqual(len(point['objective_gradient']), 10)
        self.assertEqual(len(point['objective_hessian']), 10)
        self.assertEqual(len(point['constraint_jacobian']), len(point['constraints']))
        self.assertEqual(reloaded['objective_gradient'], point['objective_gradient'])
        self.assertEqual(reloaded['geometry_authority'], authority['geometry_authority'])
        self.assertEqual(reloaded['native_moment_jacobian'], point['native_moment_jacobian'])
        vertices = [[Fraction(c) for c in v] for v in point['material_vertices_target_authority']]
        volume = Fraction(0)
        for ia, ib, ic in point['material_faces']:
            a, b, c = vertices[ia], vertices[ib], vertices[ic]
            volume += (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0])) / 6
        self.assertEqual(Fraction(point['variable_moments_target'][0]), volume * 3)
        self.assertEqual(Fraction(point['total_moments_target'][0]), 304 + volume * 3)
        self.assertGreater(len(point['source_vertices_target_authority']), 0)
        self.assertEqual(len(point['source_coordinate_lower']), len(point['source_vertices_target_authority']))
        self.assertEqual(len(point['source_coordinate_upper']), len(point['source_vertices_target_authority']))
        self.assertEqual(len(point['native_moment_jacobian']), 10)
        wrong_dual = list(authority['dual_authority'])
        wrong_dual[0] = '0/1'
        with self.assertRaises(ValueError):
            api.evaluate_spin2014_authority_quadric_point(request, settings, 2,
                authority['geometry_authority'], wrong_dual, evaluation)
        for lower, upper in zip(point['source_coordinate_lower'], point['source_coordinate_upper']):
            for lo, hi in zip(lower, upper):
                self.assertLessEqual(Fraction(hi)-Fraction(lo), Fraction(1, 2**384))
        checked = checker().check_independent_authority_point(reloaded,
            {'objective':'TOP_MASS_AND_LOWERING'},
            expected_request_identity=point['request_identity'],
            expected_context_identity=point['context_identity'])
        self.assertEqual(checked['independent_total_moments_target'], point['total_moments_target'])
        self.assertFalse(checked['numerical_stationarity_qualified'])
        regions, audit = checker().prepare_authority_material_boundaries(reloaded,
            {'objective':'TOP_MASS_AND_LOWERING'},
            expected_request_identity=point['request_identity'],
            expected_context_identity=point['context_identity'])
        self.assertEqual(set(regions), {'outer','fixed','design','retained','cavity','material'})
        self.assertEqual(regions['outer'].arrays()[0].tolist(), request.outer_surface.vertices.tolist())
        self.assertEqual(regions['outer'].arrays()[1].tolist(), request.outer_surface.faces.tolist())
        for field in ('design_partition','whole_partition','material_partition'):
            self.assertTrue(audit[field]['complete_partition'])
            self.assertTrue(audit[field]['no_positive_volume_overlap'])
        self.assertFalse(audit['manufacturing_approved'])
        self.assertFalse(audit['numerical_stationarity_qualified'])
        self.assertEqual(sorted(regions['fixed'].exact_region().facts['original_signs']), [-1,1])
        self.assertEqual(regions['fixed'].exact_region().facts['bounded_void_components'], 1)
        bundle_path = Path(__file__).parents[2]/'python/make_it_stand/verification/spin_bundle.py'
        spec = importlib.util.spec_from_file_location('make_it_stand.verification._candidate_bundle',bundle_path)
        bundle = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(bundle)
        sys.modules['make_it_stand.verification._spin2014_authority'] = checker()
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)/'authority'
            saved_regions,saved_audit = bundle._save_spin2014_authority_point_attachment(
                reloaded,{'objective':'TOP_MASS_AND_LOWERING'},destination,
                expected_request_identity=point['request_identity'],
                expected_context_identity=point['context_identity'])
            self.assertEqual(saved_audit,audit)
            for name in regions:
                self.assertEqual(saved_regions[name].identity,regions[name].identity)
            self.assertFalse(saved_audit['sealed_artifact'])
            with self.assertRaises(FileExistsError):
                bundle._save_spin2014_authority_point_attachment(
                    reloaded,{'objective':'TOP_MASS_AND_LOWERING'},destination,
                    expected_request_identity=point['request_identity'],
                    expected_context_identity=point['context_identity'])
            problem = destination/'problem.json'
            problem.write_bytes(problem.read_bytes()+b' ')
            with self.assertRaises(ValueError):
                bundle._load_spin2014_authority_point_attachment(destination,
                    expected_point_digest=audit['source_point_digest'],
                    expected_request_identity=point['request_identity'],
                    expected_context_identity=point['context_identity'])
        changed = copy.deepcopy(reloaded)
        changed['source_coordinate_lower'][0][0] = '100/1'
        with self.assertRaises(ValueError):
            checker().check_independent_authority_point(changed, {'objective':'TOP_MASS_AND_LOWERING'},
                expected_request_identity=point['request_identity'],
                expected_context_identity=point['context_identity'])
        changed = copy.deepcopy(reloaded)
        changed['source_coordinate_upper'][0][0] = str(
            Fraction(changed['source_coordinate_lower'][0][0]) + Fraction(1,2**383))
        with self.assertRaises(ValueError):
            checker().check_independent_authority_point(changed, {'objective':'TOP_MASS_AND_LOWERING'},
                expected_request_identity=point['request_identity'],
                expected_context_identity=point['context_identity'])
        for field in ('fixed_moments_target','total_moments_target','objective_gradient','geometry_authority'):
            changed = copy.deepcopy(reloaded)
            changed[field][0] = '0/1'
            with self.assertRaises(ValueError):
                checker().check_independent_authority_point(changed, {'objective':'TOP_MASS_AND_LOWERING'},
                    expected_request_identity=point['request_identity'],
                    expected_context_identity=point['context_identity'])
        changed = copy.deepcopy(reloaded)
        changed['fixed_material_density'] = 4
        with self.assertRaises(ValueError):
            checker().check_independent_authority_point(changed, {'objective':'TOP_MASS_AND_LOWERING'},
                expected_request_identity=point['request_identity'],
                expected_context_identity=point['context_identity'])

    def test_five_objectives_use_authority_moments(self):
        moments = ['2/1', '0/1', '0/1', '1/1', '3/1', '4/1', '5/1', '0/1', '0/1', '0/1']
        names = ['YOYO', 'TOP_INERTIA_ONLY', 'TOP_LOWERING_ONLY', 'TOP_MASS_AND_LOWERING', 'STATIC_LOWERING']
        expected = [Fraction(257, 98), Fraction(257, 98), Fraction(563, 196), Fraction(355, 98), Fraction(1, 4)]
        for name, want in zip(names, expected):
            settings = native.Spin2014QuadricObjectiveSettings()
            settings.objective = getattr(native.Spin2014Objective, name)
            result = api.evaluate_spin2014_authority_moment_problem(settings, moments)
            self.assertLess(abs(Fraction(result['objective']) - want), Fraction(1, 10**100))
            self.assertEqual(result['raw_moments_authority'], moments)
            self.assertEqual(len(result['objective_gradient']), 10)
            non_diagonal = ['2/1','1/4','-1/3','2/5','3/1','4/1','5/1','1/7','-1/9','2/11']
            settings.inertia_weight = 0.75
            settings.center_weight = 1.25
            settings.solver_mass_unit_in_declared_units = 3.125
            settings.solver_length_unit_in_declared_units = 0.75
            independent = checker().independent_authority_physics(non_diagonal,
                dict(objective=name,inertia_weight=0.75,center_weight=1.25,
                     solver_mass_unit_in_declared_units=3.125,solver_length_unit_in_declared_units=0.75))
            observed = api.evaluate_spin2014_authority_moment_problem(settings,non_diagonal)
            for field in ('objective_gradient','constraints'):
                for actual,wanted in zip(observed[field],independent[field]):
                    self.assertLess(abs(Fraction(actual)-Fraction(wanted)),Fraction(1,10**100))
            for actual_row,wanted_row in zip(observed['constraint_jacobian'],independent['constraint_jacobian']):
                for actual,wanted in zip(actual_row,wanted_row):
                    self.assertLess(abs(Fraction(actual)-Fraction(wanted)),Fraction(1,10**100))

    def test_noncanonical_moments_are_rejected(self):
        settings = native.Spin2014QuadricObjectiveSettings()
        moments = ['2/2'] + ['0/1'] * 9
        with self.assertRaises(ValueError):
            api.evaluate_spin2014_authority_moment_problem(settings, moments)

    def test_seed_sequence_preserves_rational_states(self):
        settings = native.SpinSolverSettings()
        first = api.spin2014_authority_initial_sequence(native.SpinMode.SPINNING_TOP, settings)
        again = api.spin2014_authority_initial_sequence(native.SpinMode.SPINNING_TOP, settings)
        self.assertEqual(first, again)
        self.assertEqual(len(first), 8)
        settings.deterministic_seed = 1
        other = api.spin2014_authority_initial_sequence(native.SpinMode.SPINNING_TOP, settings)
        self.assertNotEqual(first[0]['geometry_authority'], other[0]['geometry_authority'])
        for state in first:
            self.assertEqual(len(state['geometry_authority']), 10)
            self.assertTrue(all(isinstance(value, str) for value in state['geometry_authority']))
