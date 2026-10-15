"""The diagnostic observer must not change evaluation or leak between calls."""
import json
import unittest
from pathlib import Path
import numpy as np
from make_it_stand import spin


class NativeFailureTape(unittest.TestCase):
    def test_same_call_tape_preserves_material_and_does_not_leak(self):
        frozen = json.loads((Path(__file__).parent / 'same-input-job-03/output/frozen-inputs.json').read_text())
        payload = dict(frozen['request'])
        payload['settings'] = spin.SpinSolverConfig(**payload['settings'])
        native, request = spin._native_request(spin.SpinItFasterRequest(**payload))
        config = dict(frozen['settings'])
        config['exact_binary64_material_boundary'] = True
        settings = spin._native_evaluator_settings(native, spin.SpinEvaluatorConfig(**config))
        q = np.asarray(frozen['states']['ellipsoid'])
        normal = native.evaluate_spin_it_faster_fullrank_quadric(request, q, settings, True)
        diagnostic = native.capture_spin_it_faster_fullrank_geometry(request, q, settings)
        self.assertTrue(diagnostic['evaluation_completed'])
        self.assertFalse(diagnostic['production_admission'])
        self.assertEqual(diagnostic['error'], '')
        self.assertEqual([s['stage'] for s in diagnostic['snapshots']],
                         ['source-before-side', 'boolean-before-normalization', 'boolean-after-normalization'])
        final = diagnostic['snapshots'][-1]
        np.testing.assert_array_equal(final['material_vertices'], normal['material_vertices'])
        np.testing.assert_array_equal(final['material_faces'], normal['material_faces'])
        self.assertEqual(final['exact_material_vertices'],
                         normal['material_proof_attachment']['exact_material_vertices'])
        np.testing.assert_array_equal(final['material_birth_faces'],
                                      normal['material_proof_attachment']['material_birth_faces'])
        again = native.capture_spin_it_faster_fullrank_geometry(request, q, settings)
        self.assertEqual(len(again['snapshots']), 3)
        np.testing.assert_array_equal(again['snapshots'][-1]['material_faces'], final['material_faces'])


if __name__ == '__main__':
    unittest.main()
