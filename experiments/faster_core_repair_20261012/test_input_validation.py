from dataclasses import replace
import json
from pathlib import Path
import unittest
from make_it_stand import spin


class WarmInputValidationTests(unittest.TestCase):
    def test_warm_geometry_does_not_admit_open_boundary_or_invalid_request(self):
        frozen=json.loads((Path(__file__).parent/'same-input-job-03/output/frozen-inputs.json').read_text())
        payload=dict(frozen['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
        request=spin.SpinItFasterRequest(**payload);settings=spin.SpinEvaluatorConfig(**frozen['settings'])
        q=frozen['states']['ellipsoid']
        spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=settings)
        with self.assertRaises((ValueError,RuntimeError)):
            spin.evaluate_spin_it_faster_fullrank_quadric(replace(request,design_faces=request.design_faces[:-1]),q,
                evaluator_settings=settings)
        with self.assertRaises((ValueError,RuntimeError)):
            spin.evaluate_spin_it_faster_fullrank_quadric(replace(request,material_density=-1.),q,
                evaluator_settings=settings)
        with self.assertRaises((ValueError,RuntimeError)):
            spin.evaluate_spin_it_faster_fullrank_quadric(replace(request,target_axes=((1.,0.,0.),(0.,1.,0.),(0.,0.,2.))),q,
                evaluator_settings=settings)


if __name__=='__main__':unittest.main()
