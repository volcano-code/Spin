import json
from pathlib import Path
import unittest
import numpy as np
from make_it_stand import spin
from cherchi_core_map import CherchiCoreMap


class ActualCoreMapTests(unittest.TestCase):
    def test_same_call_retains_original_construction_evidence(self):
        folder=Path(__file__).parent/'same-input-job-03/output'
        packet=json.loads((folder/'frozen-inputs.json').read_text())
        payload=dict(packet['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
        mapping=CherchiCoreMap(spin.SpinItFasterRequest(**payload),spin.SpinEvaluatorConfig(**packet['settings']))
        q=next(iter(packet['states'].values()));mapping(q)
        records=mapping.last.get('original_construction_records',[])
        self.assertEqual(len(records),len(mapping.last['vertices']))
        self.assertEqual(sorted(r['vertex_index'] for r in records),list(range(len(records))))
        self.assertTrue(all('operand_coordinates_exact' in r and 'point_type' in r for r in records))

    def test_feature_branch_profile_does_not_change_actual_moments_or_derivatives(self):
        here=Path(__file__).parent;folder=here/'same-input-job-03/output'
        packet=json.loads((folder/'frozen-inputs.json').read_text())
        payload=dict(packet['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
        request=spin.SpinItFasterRequest(**payload);settings=spin.SpinEvaluatorConfig(**packet['settings'])
        legacy=CherchiCoreMap(request,settings)
        featured=CherchiCoreMap(request,settings,branch_identity_policy='SOURCE_FEATURE_BRANCH_V2')
        for name,q in packet['states'].items():
            with self.subTest(state=name):
                a=legacy(q);b=featured(q)
                np.testing.assert_array_equal(a['raw_moments'],b['raw_moments'])
                np.testing.assert_array_equal(a['moment_jacobian'],b['moment_jacobian'])
                self.assertEqual(a['geometry_hash'],b['geometry_hash'])
                self.assertNotEqual(a['evaluation_map_identity'],b['evaluation_map_identity'])

    def test_three_same_operand_maps_match_saved_production_moments_and_jacobians(self):
        here=Path(__file__).parent;folder=here/'same-input-job-03/output'
        packet=json.loads((folder/'frozen-inputs.json').read_text())
        payload=dict(packet['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
        mapping=CherchiCoreMap(spin.SpinItFasterRequest(**payload),spin.SpinEvaluatorConfig(**packet['settings']))
        for name,q in packet['states'].items():
            with self.subTest(state=name):
                result=mapping(q)
                target=json.loads((folder/(name+'-fullrank-evaluation.json')).read_text())
                np.testing.assert_allclose(result['raw_moments'],target['total_raw_moments'],rtol=0,atol=1e-10)
                np.testing.assert_allclose(result['moment_jacobian'],target['coefficient_jacobian'],rtol=0,atol=1e-6)


if __name__=='__main__':unittest.main()
