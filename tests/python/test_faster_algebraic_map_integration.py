"""Explicit finite-map selection must reach geometry, derivative and evidence."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from make_it_stand import spin
from test_faster_cherchi_newton_integration_runner import RUNNER, _Extension
from test_faster_cherchi_table_dispatch import locked_case
from test_spin_it_faster_table1_runner import load_runner
import test_faster_cherchi_basis as basis_tests

ROOT = Path(__file__).resolve().parents[2]


class AlgebraicMapIntegrationTests(unittest.TestCase):
    def test_explicit_selection_reaches_actual_source_and_branch_identity(self):
        self.assertIn('algebraic_finite_map', spin.SpinEvaluatorConfig.__dataclass_fields__,
                      'Python evaluator cannot select the validated finite map')
        native = spin._native_module()
        self.assertTrue(hasattr(native.SpinTriangulationSettings(), 'algebraic_finite_map'),
                        'loaded binary does not expose the finite-map policy')
        parser = RUNNER._argument_parser()
        self.assertFalse(parser.parse_args(['--output', 'unused']).algebraic_finite_map)
        self.assertTrue(parser.parse_args(['--output', 'unused', '--algebraic-finite-map']).algebraic_finite_map)
        settings, config = RUNNER._make_triangulation_configuration(native,
            icosphere_subdivisions=1, algebraic_finite_map=True)
        q = (.36,.04,-.03,.02,-1,-1.4,-1.8,.01,-.02,.03)
        domain = ((-1.,-1.,-1.), (1.,1.,1.))
        source = native.evaluate_spin_it_faster_source_motion(q,domain,settings)
        self.assertTrue(source['defined'], source['diagnostic'])
        self.assertIs(source['algebraic_finite_map'], True)
        via_config = native.evaluate_spin_it_faster_source_motion(q,domain,
            spin._native_evaluator_settings(native,config).triangulation)
        for key in ('vertices','faces','vertex_velocity'):
            np.testing.assert_array_equal(source[key],via_config[key])
        branch,_ = RUNNER._source_derivative_branch(source,algebraic_finite_map=True)
        old_settings,_ = RUNNER._make_triangulation_configuration(native,icosphere_subdivisions=1)
        old = native.evaluate_spin_it_faster_source_motion(q,domain,old_settings)
        self.assertNotEqual(source['triangulation_hash'], old['triangulation_hash'])
        old_branch,_ = RUNNER._source_derivative_branch(old)
        self.assertNotEqual(RUNNER._derivative_stratum_identity('intersection','b'*64,branch),
                            RUNNER._derivative_stratum_identity('intersection','b'*64,old_branch))
        for bad in (dict(source,algebraic_finite_map=False),
                    {k:v for k,v in source.items() if k!='algebraic_finite_map'}):
            with self.assertRaisesRegex(ValueError,'algebraic'):
                RUNNER._source_derivative_branch(bad,algebraic_finite_map=True)

    def test_bad_flags_and_old_native_are_fail_closed_but_legacy_remains_usable(self):
        self.assertIn('algebraic_finite_map', spin.SpinEvaluatorConfig.__dataclass_fields__)
        for bad in (1,0,None,'true'):
            with self.subTest(bad=bad),self.assertRaises(ValueError):
                spin.SpinEvaluatorConfig(algebraic_finite_map=bad)
            with self.subTest(runner=bad),self.assertRaises(ValueError):
                RUNNER._make_triangulation_configuration(_Extension,
                    icosphere_subdivisions=1,algebraic_finite_map=bad)
        with self.assertRaisesRegex(RuntimeError,'unavailable'):
            RUNNER._make_triangulation_configuration(_Extension,
                icosphere_subdivisions=1,algebraic_finite_map=True)
        settings,config=RUNNER._make_triangulation_configuration(_Extension,icosphere_subdivisions=1)
        self.assertFalse(config.algebraic_finite_map)
        self.assertFalse(hasattr(settings,'algebraic_finite_map'))

    def test_basis_adapter_binds_policy_and_rejects_unsupported_paraboloid(self):
        native=spin._native_module();runner=basis_tests.CherchiBasisTests().module()
        fixture=ROOT/'experiments/spin_it_reproduction/faster/fixtures/cherchi-full-basis-method-v1/problem.json'
        p=json.loads(fixture.read_text());p['algebraic_finite_map']=True
        settings=runner._settings(native,p)
        self.assertTrue(settings.algebraic_finite_map)
        p['schema']='cherchi-paraboloid-basis-method-v1'
        with self.assertRaisesRegex(ValueError,'algebraic|rank'):
            runner._settings(native,p)

    def test_table_plan_preserves_opt_in_and_rejects_non_boolean(self):
        runner=load_runner()
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);case=locked_case(root)
            old=runner.cherchi_case_plan(case,1,1,asset_root=root)
            self.assertNotIn('algebraic_finite_map',old['settings'])
            for value in (True,False):
                c=copy.deepcopy(case);c['cherchi_settings']['algebraic_finite_map']=value
                self.assertIs(runner.cherchi_case_plan(c,1,1,asset_root=root)['settings']['algebraic_finite_map'],value)
            for value in (1,None,'true'):
                c=copy.deepcopy(case);c['cherchi_settings']['algebraic_finite_map']=value
                with self.assertRaises(ValueError):runner.cherchi_case_plan(c,1,1,asset_root=root)

    def test_actual_budgeted_table_call_binds_policy_without_claiming_success(self):
        from make_it_stand.verification.cherchi_application import audit_saved_cherchi_case
        native=spin._native_module();runner=load_runner()
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);case=locked_case(root/'inputs')
            case['cherchi_settings']['algebraic_finite_map']=True
            result=runner.run_case(case,1,1,bundle_root=root/'runs',asset_root=root/'inputs',
                cherchi_build=Path(native.__file__).resolve().parent.parent,
                cherchi_probe=ROOT/'build/cherchi-reference-singlethread-20260917/cherchi_probe')
            inputs=json.loads((root/'runs/locked-method/run-000000/inputs.json').read_text())
            self.assertIs(inputs['algebraic_finite_map'],True)
            self.assertFalse(result['passed'])
            self.assertEqual(result['summary']['total_function_evaluations_all_attempts'],1)
            record=json.loads((root/'runs/locked-method/run-000000/result.json').read_text())
            successful=[e for e in record['evaluations'] if 'moment_map' in e]
            self.assertGreater(len(successful),0,'budget consumption hid a geometry integration failure')
            self.assertTrue(all(e['derivative_branch']['finite_point_map']=='ALGEBRAIC_V1' for e in successful))
            audit_saved_cherchi_case(result,bundle_root=root/'runs')
            for remove in (False,True):
                bad=copy.deepcopy(result);s=bad['initialization_plan']['settings']
                if remove:s.pop('algebraic_finite_map')
                else:s['algebraic_finite_map']=False
                bad['scientific_hash']=runner.scientific_hash({k:v for k,v in bad.items() if k!='scientific_hash'})
                with self.assertRaisesRegex(ValueError,'algebraic'):
                    audit_saved_cherchi_case(bad,bundle_root=root/'runs')


if __name__=='__main__': unittest.main()
