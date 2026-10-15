"""Catch successful Newton solves silently omitting the requested payload."""
import json
from pathlib import Path
import unittest
import numpy as np
from make_it_stand import spin
from cherchi_core_map import CherchiCoreMap


class SuccessfulNewtonCaptureTests(unittest.TestCase):
    def test_regular_linear_solve_emits_matrix_without_extra_map_evaluation(self):
        here=Path(__file__).parent
        frozen=json.loads((here/'same-input-job-03/output/frozen-inputs.json').read_text())
        payload=dict(frozen['request'])
        payload.update(mode='SPINNING_TOP',target_origin=(0.,0.,-1.1),center_of_mass_weight=1.)
        config=dict(payload['settings']);config.update(equation_profile='DUAL_MONOMIAL_RAW_V1',
            deterministic_seed=1,maximum_iterations=3,maximum_restarts=1,
            maximum_function_evaluations_per_restart=3,capture_evaluation_payload=True)
        payload['settings']=spin.SpinSolverConfig(**config)
        request=spin.SpinItFasterRequest(**payload)
        mapping=CherchiCoreMap(request,spin.SpinEvaluatorConfig(**frozen['settings']))
        events=[]
        result=spin.solve_spin_it_faster_with_external_moment_evaluator(request,mapping,event_sink=events.append)
        captured=[r for r in events if r.get('newton_system_payload')]
        self.assertGreaterEqual(len(captured),1,'Successful regular Newton matrix was not captured')
        self.assertEqual(result['function_evaluations'],3)
        self.assertEqual(len([r for r in events if r['type']=='evaluation']),3)
        first=captured[0]['newton_system_payload']
        self.assertEqual(np.asarray(first['matrix']).shape,(15,15))
        self.assertEqual(np.asarray(first['right_hand_side']).shape,(15,))
        self.assertEqual(first['computed_rank'],15)
        self.assertGreater(first['multiplier_mu'],0)


if __name__=='__main__':unittest.main()
