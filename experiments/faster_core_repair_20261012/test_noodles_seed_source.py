"""Actual first seeded state: no author answer or successful terminal reused."""
import json
from pathlib import Path
import unittest
import numpy as np
import _make_it_stand as native

class NoodlesSeedSource(unittest.TestCase):
    def test_first_frozen_seed_has_supported_projective_source_and_motion(self):
        here=Path(__file__).resolve().parent
        run=here/'noodles-bound-seed0-production-probe-job-01/output/run'
        config=json.loads((run/'mapping-config.json').read_text())
        rows=[json.loads(x) for x in (run/'mapping-events.jsonl').read_text().splitlines()]
        first=next(row for row in rows if row['status']=='STARTED')
        request=config['request'];evaluator=config['evaluator']
        design=(np.asarray(request['design_vertices'])-np.asarray(request['target_origin']))@np.asarray(request['target_axes'])
        settings=native.SpinTriangulationSettings()
        for name in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
            'retain_fixed_normal_provenance','truncation_radius','paper_projective_hyperboloid','algebraic_finite_map'):
            setattr(settings,name,evaluator[name])
        try:
            source=native.evaluate_spin_it_faster_source_motion(np.asarray(first['coefficients']),design,settings)
        except RuntimeError as exc:
            self.fail('frozen first seed projective source unsupported: '+str(exc))
        self.assertTrue(source['defined'],source['diagnostic'])
        self.assertTrue(source['side_selection_defined'])
        self.assertEqual(np.asarray(source['vertex_velocity']).size,30*len(source['vertices']))

if __name__=='__main__':unittest.main()
