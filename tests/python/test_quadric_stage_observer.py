"""Diagnostic callback must observe real production stages without changing science."""
import unittest
import numpy as np

from make_it_stand.spin import (SpinEvaluatorConfig, SpinItFasterRequest,
    _native_request, _native_evaluator_settings, evaluate_spin_it_faster_fullrank_quadric)
from make_it_stand.verification.faster_material import _json_value
from test_spin2014_quadric_interface import box


class QuadricStageObserver(unittest.TestCase):
    def test_public_observer_is_not_part_of_scientific_configuration(self):
        outer, faces=box(-2.,2.)
        design, df=box(-1.,1.)
        request=SpinItFasterRequest(outer,faces,design,df,material_density=2.7)
        config=SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
            icosphere_subdivisions=1,rational_source_coordinates=True)
        q=(-.25,0.,0.,0.,1.,1.,1.,0.,0.,0.)
        identity=config.identity_dict()
        stages=[]
        baseline=evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=config)
        measured=evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=config,
            diagnostic_stage_observer=stages.append)
        self.assertEqual(baseline,measured)
        self.assertEqual(identity,config.identity_dict())
        self.assertIn('boolean-end',stages)
        with self.assertRaisesRegex(ValueError,'callable'):
            evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=config,
                diagnostic_stage_observer=7)

    def test_real_callback_preserves_complete_production_packet(self):
        outer, faces=box(-2.,2.)
        design, df=box(-1.,1.)
        request=SpinItFasterRequest(outer,faces,design,df,material_density=2.7)
        native,nr=_native_request(request)
        config=SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
            icosphere_subdivisions=1,rational_source_coordinates=True)
        settings=_native_evaluator_settings(native,config)
        # Removing the callback binding or skipping a real stage must fail;
        # equality covers actual material, moments, J and scientific identity.
        self.assertIsNone(settings.diagnostic_stage_observer)
        q=np.array([-.25,0.,0.,0.,1.,1.,1.,0.,0.,0.])
        baseline=native.evaluate_spin_it_faster_fullrank_quadric(nr,q,settings)
        stages=[]
        settings.diagnostic_stage_observer=stages.append
        measured=native.evaluate_spin_it_faster_fullrank_quadric(nr,q,settings)
        self.assertEqual(_json_value(baseline),_json_value(measured))
        self.assertTrue(measured['derivative_defined'])
        for name in ('source-construction','boolean','material-authority-audit',
                     'moments','source-motion','material-transport'):
            self.assertIn(name+'-begin',stages)
            self.assertIn(name+'-end',stages)
            self.assertLess(stages.index(name+'-begin'),stages.index(name+'-end'))
        settings.diagnostic_stage_observer=None
        self.assertIsNone(settings.diagnostic_stage_observer)

    def test_observer_error_is_not_swallowed_as_a_scientific_result(self):
        outer, faces=box(-2.,2.)
        design, df=box(-1.,1.)
        native,nr=_native_request(SpinItFasterRequest(outer,faces,design,df))
        settings=_native_evaluator_settings(native,SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC',
            icosphere_subdivisions=1,rational_source_coordinates=True))
        def fail(stage):
            raise LookupError('diagnostic sink failed at '+stage)
        settings.diagnostic_stage_observer=fail
        with self.assertRaisesRegex(LookupError,'diagnostic sink failed'):
            native.evaluate_spin_it_faster_fullrank_quadric(nr,
                np.array([-.25,0.,0.,0.,1.,1.,1.,0.,0.,0.]),settings)


if __name__=='__main__':
    unittest.main()
