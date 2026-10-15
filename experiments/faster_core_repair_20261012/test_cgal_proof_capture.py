"""Catch loss of same-call exact coordinates/birth facets at the Python boundary."""
import json,unittest
from pathlib import Path
from fractions import Fraction
import numpy as np
from make_it_stand import spin
class CgalProofCapture(unittest.TestCase):
    def test_opt_in_attachment_preserves_actual_default_result(self):
        frozen=json.loads((Path(__file__).parent/'same-input-job-03/output/frozen-inputs.json').read_text())
        payload=dict(frozen['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
        native,request=spin._native_request(spin.SpinItFasterRequest(**payload))
        config=dict(frozen['settings']);config['exact_binary64_material_boundary']=True
        settings=spin._native_evaluator_settings(native,spin.SpinEvaluatorConfig(**config))
        q=np.asarray(frozen['states']['ellipsoid'])
        default=native.evaluate_spin_it_faster_fullrank_quadric(request,q,settings)
        try:captured=native.evaluate_spin_it_faster_fullrank_quadric(request,q,settings,capture_material_proof=True)
        except TypeError as exc:self.fail('same-call CGAL proof capture unavailable: '+str(exc).splitlines()[0])
        proof=captured.pop('material_proof_attachment')
        self.assertNotIn('material_proof_attachment',default)
        for key in default:
            if isinstance(default[key],np.ndarray):np.testing.assert_array_equal(default[key],captured[key])
            else:self.assertEqual(default[key],captured[key])
        self.assertEqual(len(proof['exact_material_vertices']),len(default['material_vertices']))
        self.assertEqual(len(proof['material_birth_faces']),len(default['material_faces']))
        self.assertGreater(len(proof['source_faces']),0)
        for i,row in enumerate(proof['exact_material_vertices']):
            for k,x in enumerate(row):self.assertAlmostEqual(float(Fraction(x)),default['material_vertices'][i,k],places=14)
        self.assertTrue(all(0<=f<len(payload['design_faces'])+len(proof['source_faces']) for f in proof['material_birth_faces']))
if __name__=='__main__':unittest.main()
