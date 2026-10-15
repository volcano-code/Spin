"""Real native regressions: reuse only unchanged fixed design geometry."""
from dataclasses import replace
import json
from pathlib import Path
import unittest
from make_it_stand import spin


class FixedDesignCacheTests(unittest.TestCase):
    def test_fixed_proof_reused_and_changed_vertex_invalidates_it(self):
        frozen=json.loads((Path(__file__).parent/'same-input-job-03/output/frozen-inputs.json').read_text())
        payload=dict(frozen['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
        request=spin.SpinItFasterRequest(**payload)
        settings=spin.SpinEvaluatorConfig(**frozen['settings'])
        q=frozen['states']['ellipsoid']
        first=[];second=[];changed=[]
        a=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=settings,diagnostic_stage_observer=first.append)
        b=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=settings,diagnostic_stage_observer=second.append)
        self.assertEqual(a,b)
        self.assertIn('transport-design-proof-cache-hit',second)
        vertices=[list(row) for row in request.design_vertices];vertices[0][0]-=.01
        spin.evaluate_spin_it_faster_fullrank_quadric(replace(request,design_vertices=vertices),q,
            evaluator_settings=settings,diagnostic_stage_observer=changed.append)
        self.assertNotIn('transport-design-proof-cache-hit',changed)
        policy_changed=[]
        spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=settings)
        spin.evaluate_spin_it_faster_fullrank_quadric(request,q,
            evaluator_settings=replace(settings,exact_binary64_source_boundary=False),
            diagnostic_stage_observer=policy_changed.append)
        self.assertNotIn('transport-design-proof-cache-hit',policy_changed)


if __name__=='__main__':unittest.main()
