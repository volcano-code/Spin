"""Catch omission of connected finite-normal patches in actual Heart states."""
import json
from pathlib import Path
import unittest
import numpy as np
from make_it_stand import spin

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


class SavedHeartProjectiveTests(unittest.TestCase):
    def test_saved_failed_states_have_closed_source_and_actual_motion(self):
        packet=json.loads((ROOT/'experiments/faster_heart_contact_20261011'
            /'rounded-clearance-v3-faster-frozen-v2/request.json').read_text())
        request=packet['requests'][0]
        rotation=np.asarray(request['target_axes']).T
        design=(np.asarray(request['design_vertices'])-request['target_origin'])@rotation.T
        states=json.loads((HERE/'frozen-heart-states.json').read_text())['records']
        native=spin._native_module()
        for name in ('heart-first-cap-failure','heart-first-small-face-rejection'):
            with self.subTest(state=name):
                settings=native.SpinTriangulationSettings()
                for attr in ('icosphere_subdivisions','maximum_icosphere_subdivisions',
                    'maximum_radius_doublings','retain_fixed_normal_provenance',
                    'truncation_radius','algebraic_finite_map'):
                    setattr(settings,attr,packet['evaluator_settings'][attr])
                settings.paper_projective_hyperboloid=True
                source=native.evaluate_spin_it_faster_source_motion(
                    np.asarray(states[name]['coefficients']),design,settings)
                self.assertTrue(source['defined'],source['diagnostic'])
                self.assertTrue(source['side_selection_defined'])
                self.assertEqual(source['effective_icosphere_subdivisions'],2)
                self.assertTrue(np.isfinite(source['vertex_velocity']).all())
                audit=native.audit_exact_triangle_boundary(source['vertices'],source['faces'],'OutwardPositive')
                self.assertTrue(audit['boundary_valid'],str(audit))


if __name__=='__main__':unittest.main()
