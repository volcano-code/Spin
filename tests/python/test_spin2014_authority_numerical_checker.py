import importlib.util
import subprocess
import sys
import unittest
from unittest.mock import patch
from fractions import Fraction
from pathlib import Path

import make_it_stand.verification


def checker():
    source = Path(__file__).resolve().parents[2] / 'python/make_it_stand/verification/_spin2014_authority.py'
    spec = importlib.util.spec_from_file_location('make_it_stand.verification._spin2014_authority', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def text(value):
    value = Fraction(value)
    return str(value.numerator) + '/' + str(value.denominator)


class AuthorityNumericalCheckerTests(unittest.TestCase):
    moments = ['2/1', '0/1', '0/1', '1/1', '3/1', '4/1', '5/1', '0/1', '0/1', '0/1']

    def test_actual_material_mapping_accepts_explicit_exact_gmp_backend(self):
        # A literal signed-tetrahedron mapping fixture, not a quadric or seal.
        api=checker()
        faces=[[0,2,1],[0,1,3],[0,3,2],[1,2,3]]
        outer=[['0/1','0/1','0/1'],['1/1','0/1','0/1'],
               ['0/1','1/1','0/1'],['0/1','0/1','1/1']]
        material=[list(row) for row in outer];material[1][0]='1/3'
        expected=['1/6','1/72','1/24','1/24','1/540','1/60','1/60','1/120','1/360','1/360']
        affine=[['1/1','0/1','0/1','0/1'],['0/1','1/1','0/1','0/1'],['0/1','0/1','1/1','0/1']]
        point=dict(schema='spin2014-authority-material-map-v1',
            scope='MOMENT_MAPPING_NOT_SOLVER_OR_FINAL_SEAL',request_identity='a'*64,context_identity='b'*64,
            coefficient_policy='COEFFICIENT_MPFR1024_V1',precision_bits=1024,
            geometry_authority=['1/1']+['0/1']*9,dual_authority=['1/1']+['0/1']*9,
            model_to_target_authority=affine,target_to_model_authority=affine,
            source_vertices_target_authority=material,source_coordinate_lower=material,source_coordinate_upper=material,
            outer_vertices_model_authority=outer,outer_faces=faces,
            design_vertices_model_authority=outer,design_vertices_target_authority=outer,design_faces=faces,
            material_vertices_target_authority=material,material_faces=faces,
            fixed_material_density=2,design_material_density=3,fixed_moments_target=['0/1']*10,
            variable_moments_target=expected,total_moments_target=expected)
        options=dict(expected_request_identity='a'*64,expected_context_identity='b'*64)
        reference=api.check_independent_authority_material_mapping(point,**options)
        self.assertNotIn('arithmetic_identity',reference,
            'legacy reports must not gain metadata that invalidates saved audits')
        actual=api.check_independent_authority_material_mapping(point,**options,
            arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
        self.assertEqual(actual['independent_total_moments_target'],expected)
        for key in ('independent_fixed_moments_target','independent_variable_moments_target','independent_total_moments_target'):
            self.assertEqual(actual[key],reference[key])
        self.assertEqual(actual['arithmetic_identity']['gmpy2'],'2.2.1')
        bad=dict(point,variable_moments_target=['1/1']+expected[1:])
        with self.assertRaisesRegex(ValueError,'moments mismatch'):
            api.check_independent_authority_material_mapping(bad,**options,arithmetic_backend='GMPY2_2_2_1_EXACT_V1')

    def test_independent_initial_screen_rejects_nonpositive_mu_and_zero_sensitivity(self):
        api=checker()
        point=dict(multiplier_mu='1/1',native_moment_jacobian=[['0/1']*10 for _ in range(10)])
        with self.assertRaisesRegex(ValueError,'zero initial moment sensitivity'):
            api.require_authority_initial_screen(point)
        point['native_moment_jacobian'][0][0]='1/1'
        api.require_authority_initial_screen(point)
        point['native_moment_jacobian'][9][9]='NaN'
        with self.assertRaises(ValueError):
            api.require_authority_initial_screen(point)
        point['native_moment_jacobian'][9][9]='0/1'
        for mu in ('0/1','-1/1'):
            point['multiplier_mu']=mu
            with self.assertRaisesRegex(ValueError,'nonpositive initial multiplier'):
                api.require_authority_initial_screen(point)

    def test_full_run_admission_rejects_diagnostic_and_missing_initial_material(self):
        api=checker()
        args=dict(settings={},expected_initial_sequence=[],journal_path='missing.jsonl',
            expected_request_identity='a'*64,expected_context_identity='b'*64,
            expected_objective_identity='c'*64,initial_point=None)
        with self.assertRaisesRegex(ValueError,'complete full-budget run required'):
            api.check_authority_complete_run_evidence(dict(schema='spin2014-authority-restarts-v1',
                protocol_scope='DIAGNOSTIC_REDUCED_ATTEMPT_BUDGET'),**args)
        with self.assertRaisesRegex(ValueError,'initial actual material required'):
            api.check_authority_complete_run_evidence(dict(schema='spin2014-authority-restarts-v1',
                protocol_scope='CONFIGURED_FULL_ATTEMPT_BUDGET'),**args)

    def test_material_preparation_rejects_unsupported_arithmetic_before_geometry(self):
        with self.assertRaisesRegex(ValueError,'unsupported authority arithmetic backend'):
            checker().prepare_authority_material_boundaries({}, {},
                expected_request_identity='a'*64,expected_context_identity='b'*64,
                arithmetic_backend='FLOAT_APPROXIMATION')

    def test_attachment_rejects_unsupported_arithmetic_without_creating_output(self):
        import tempfile
        from make_it_stand.verification.spin_bundle import _save_spin2014_authority_point_attachment
        with tempfile.TemporaryDirectory() as directory:
            destination=Path(directory)/'not-created'
            with self.assertRaisesRegex(ValueError,'unsupported authority arithmetic backend'):
                _save_spin2014_authority_point_attachment({}, {},destination,
                    expected_request_identity='a'*64,expected_context_identity='b'*64,
                    arithmetic_backend='FLOAT_APPROXIMATION')
            self.assertFalse(destination.exists())

    def test_material_change_requires_moments_not_only_changed_hash(self):
        api=checker()
        attempt=dict(accepted_updates=1,initial_total_moments=self.moments,
                     initial_material_identity='a'*64)
        point=dict(total_moments_target=list(self.moments),material_scientific_hash='b'*64)
        report=api.check_authority_material_change(attempt,point)
        self.assertFalse(report['material_change_witness'])
        self.assertFalse(report['sealed_artifact'])
        point['total_moments_target'][0]='3/1'
        self.assertTrue(api.check_authority_material_change(attempt,point)['material_change_witness'])
        attempt['accepted_updates']=0
        self.assertFalse(api.check_authority_material_change(attempt,point)['material_change_witness'])

    def test_structural_ledger_accepts_native_relinearized_topology_reason(self):
        # Break caught: rejecting the native accepted transition label although
        # the current geometry was rebuilt and all accounting is consistent.
        # Literal metadata tests the structural parser, NOT an optimization pass.
        q = ['1/1']+['0/1']*9
        attempt = dict(initial_geometry_authority=q,initial_dual_authority=q,
            function_evaluations=2,initial_state_identity='a'*64,initial_material_identity='f'*64,
            initial_total_moments=['1/1','0/1','0/1','0/1','1/1','1/1','1/1','0/1','0/1','0/1'],
            initial_objective='2/1',trials=[dict(accepted=True,arc_length='1/4',
                reason='topology_transition_relinearized',raw_equation_residual='1/1000000')],
            accepted_updates=1,accepted_state_identities=['b'*64],step_identities=['c'*64],
            final_state_identity='b'*64,status='evaluation_budget_exhausted')
        run = dict(schema='spin2014-authority-restarts-v1',request_identity='d'*64,
            context_identity='e'*64,protocol_scope='DIAGNOSTIC_REDUCED_ATTEMPT_BUDGET',
            attempts=[attempt],stationary_attempt_index=None,stationary_point=None,
            status='restart_budget_exhausted',function_evaluations=2)
        report = checker().check_authority_restart_ledger(run,
            objective_settings=dict(objective='TOP_INERTIA_ONLY'),
            expected_initial_sequence=[dict(geometry_authority=q,dual_authority=q)],
            evaluation_budget_per_attempt=2,expected_request_identity='d'*64,
            expected_context_identity='e'*64)
        self.assertFalse(report['full_run_candidate'])
        self.assertFalse(report['sealed_artifact'])

    def test_gmp_exact_backend_matches_fraction_on_all_objectives_and_cross_moments(self):
        api = checker()
        moments = ['2/1','1/5','-1/7','1/1','3/1','4/1','5/1','-1/9','1/11','1/13']
        for name in ('YOYO','TOP_INERTIA_ONLY','TOP_LOWERING_ONLY',
                     'TOP_MASS_AND_LOWERING','STATIC_LOWERING'):
            settings = dict(objective=name,inertia_weight=2,center_weight=3,
                            solver_mass_unit_in_declared_units=3,
                            solver_length_unit_in_declared_units=2)
            reference = api.independent_authority_physics(moments,settings)
            actual = api.independent_authority_physics(moments,settings,
                arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
            self.assertEqual(actual,reference)

    def test_gmp_stationarity_keeps_exact_threshold_and_canonical_rejection(self):
        api = checker()
        moments = ['1/1','0/1','0/1','0/1','1/1','1/1','1/1','0/1','0/1','0/1']
        dual = ['1/1']+['0/1']*9
        settings = dict(objective='STATIC_LOWERING',center_weight=0)
        for mu in (Fraction(1,10**12),Fraction(1,10**12)+Fraction(1,10**50),Fraction(0)):
            reference = api.check_independent_authority_stationarity(
                moments,settings,dual,['0/1','0/1'],text(mu))
            actual = api.check_independent_authority_stationarity(
                moments,settings,dual,['0/1','0/1'],text(mu),
                arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
            self.assertEqual(actual,reference)
        for bad in ('-0/1','2/2','0/2','1/-2','01/1','1/01'):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                api.independent_authority_physics([bad]+moments[1:],settings,
                    arithmetic_backend='GMPY2_2_2_1_EXACT_V1')

    def test_large_canonical_text_conversion_has_bounded_cost_and_exact_value(self):
        # Resource regression for the measured quadratic decimal conversion.
        # Expected digits and integer are independently specified, with zero
        # runs on both sides of the split and a non-power-of-ten remainder.
        code = '''
import sys
from fractions import Fraction
from make_it_stand.verification.rational_text import fraction_from_rational_text, fraction_to_rational_text
guard = sys.get_int_max_str_digits()
value = Fraction(-(10**800000 + 10**400000 + 5), 3)
expected = '-1' + '0'*399999 + '1' + '0'*399999 + '5/3'
assert fraction_to_rational_text(value) == expected
assert fraction_from_rational_text(expected, require_canonical=True) == value
assert sys.get_int_max_str_digits() == guard
'''
        subprocess.run([sys.executable, '-B', '-c', code], timeout=10,
                       check=True, capture_output=True, text=True)

    def test_canonical_authority_read_does_not_render_large_integers(self):
        # Rendering is the measured slow pure operation, not the thing under
        # test: canonicality must be decided from syntax and reduced integers.
        api = checker()
        with patch.object(api, 'fraction_to_rational_text',
                          side_effect=AssertionError('unexpected decimal round trip')):
            self.assertEqual(api._read('1/' + '1' + '0'*5000), Fraction(1, 10**5000))
            self.assertEqual(api._read('-7/3'), Fraction(-7, 3))
            self.assertEqual(api._read('0/1'), Fraction(0))
            for token in ('-0/1', '0/2', '2/2', '01/1', '1/01', '+1/1',
                          '1/-2', '1/0', ' 1/1', '1/1\n'):
                with self.subTest(token=token), self.assertRaises(ValueError):
                    api._read(token)

    def test_five_objectives_and_mass_weighted_gradient(self):
        api = checker()
        names = ['YOYO', 'TOP_INERTIA_ONLY', 'TOP_LOWERING_ONLY', 'TOP_MASS_AND_LOWERING', 'STATIC_LOWERING']
        expected = [Fraction(257,98), Fraction(257,98), Fraction(563,196), Fraction(355,98), Fraction(1,4)]
        results = [api.independent_authority_physics(self.moments, {'objective': name}) for name in names]
        for row, wanted in zip(results, expected):
            self.assertEqual(Fraction(row['objective']), wanted)
        self.assertEqual([Fraction(x) for x in results[4]['objective_gradient']],
                         [Fraction(-1,4), 0, 0, Fraction(1,2), 0, 0, 0, 0, 0, 0])
        self.assertEqual([Fraction(a)-Fraction(b) for a,b in zip(results[3]['objective_gradient'],results[1]['objective_gradient'])],
                         [0, 0, 0, 2, 0, 0, 0, 0, 0, 0])
        self.assertEqual(len(results[0]['constraints']), 5)
        self.assertEqual(len(results[4]['constraints']), 2)

    def test_units_change_height_and_mass_weighted_height_differently(self):
        api = checker()
        expected = {'TOP_INERTIA_ONLY': Fraction(257,98), 'TOP_LOWERING_ONLY': Fraction(355,98),
                    'TOP_MASS_AND_LOWERING': Fraction(3785,98), 'STATIC_LOWERING': Fraction(1)}
        for name,wanted in expected.items():
            settings = dict(objective=name,solver_mass_unit_in_declared_units=3,
                            solver_length_unit_in_declared_units=2)
            self.assertEqual(Fraction(api.independent_authority_physics(self.moments,settings)['objective']), wanted)

    def test_raw_threshold_does_not_round_to_binary64(self):
        api = checker()
        moments = ['1/1','0/1','0/1','0/1','1/1','1/1','1/1','0/1','0/1','0/1']
        dual = ['1/1'] + ['0/1']*9
        settings = dict(objective='STATIC_LOWERING',center_weight=0)
        limit = Fraction(1,10**12)
        edge = api.check_independent_authority_stationarity(moments,settings,dual,['0/1','0/1'],text(limit))
        over = api.check_independent_authority_stationarity(moments,settings,dual,['0/1','0/1'],text(limit+Fraction(1,10**50)))
        self.assertTrue(edge['raw_equation_within_tolerance'])
        self.assertFalse(over['raw_equation_within_tolerance'])
        self.assertFalse(over['numerical_stationarity_qualified'])
        self.assertEqual(Fraction(over['raw_equation_residual_squared']), (limit+Fraction(1,10**50))**2)
        nonpositive = api.check_independent_authority_stationarity(moments,settings,dual,['0/1','0/1'],'0/1')
        self.assertFalse(nonpositive['numerical_stationarity_qualified'])

    def test_invalid_authority_and_mass_are_rejected(self):
        api = checker()
        with self.assertRaises(ValueError):
            api.independent_authority_physics(['2/2'] + self.moments[1:], {'objective':'YOYO'})
        with self.assertRaises(ValueError):
            api.independent_authority_physics(['-1/1'] + self.moments[1:], {'objective':'YOYO'})
        with self.assertRaises(ValueError):
            api.independent_authority_physics(self.moments, {'objective':'TOP_MASS_AND_LOWERING','center_weight':float('nan')})

    def test_initial_objective_is_recomputed_from_actual_recorded_moments(self):
        api = checker()
        initial = dict(initial_total_moments=self.moments,
                       initial_objective='355/98',initial_material_identity='a'*64)
        checked = api.check_authority_initial_objective(initial,{'objective':'TOP_MASS_AND_LOWERING'})
        self.assertEqual(checked['independent_initial_objective'],'355/98')
        changed = dict(initial,initial_objective='257/98')
        with self.assertRaises(ValueError):
            api.check_authority_initial_objective(changed,{'objective':'TOP_MASS_AND_LOWERING'})
        with self.assertRaises(ValueError):
            api.check_authority_initial_objective(dict(initial,initial_total_moments=None),
                                                  {'objective':'TOP_MASS_AND_LOWERING'})


if __name__ == '__main__':
    unittest.main()
