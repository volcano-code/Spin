"""Public consumers must retain same-call proof without changing legacy output."""
from dataclasses import asdict
import json
from pathlib import Path
import unittest
from make_it_stand import spin


class PublicMaterialProof(unittest.TestCase):
    def test_public_attachment_survives_serialization_without_changing_science(self):
        frozen = json.loads((Path(__file__).parent / 'same-input-job-03/output/frozen-inputs.json').read_text())
        payload = dict(frozen['request'])
        payload['settings'] = spin.SpinSolverConfig(**payload['settings'])
        request = spin.SpinItFasterRequest(**payload)
        config = dict(frozen['settings'], exact_binary64_material_boundary=True)
        settings = spin.SpinEvaluatorConfig(**config)
        old = spin.evaluate_spin_it_faster_fullrank_quadric(request, frozen['states']['ellipsoid'], evaluator_settings=settings)
        try:
            captured = spin.evaluate_spin_it_faster_fullrank_quadric(request, frozen['states']['ellipsoid'], evaluator_settings=settings, capture_material_proof=True)
        except TypeError as exc:
            self.fail('public same-call material proof unavailable: ' + str(exc))
        packet = asdict(captured)
        proof = json.loads(packet.pop('material_proof_json'))
        self.assertEqual(packet, asdict(old))
        self.assertEqual(proof['geometry_hash'], old.geometry_hash)
        self.assertEqual(len(proof['exact_material_vertices']), len(old.material_vertices))
        self.assertEqual(len(proof['material_birth_faces']), len(old.material_faces))
        captured.material_proof['geometry_hash'] = 'caller mutation'
        self.assertEqual(captured.material_proof['geometry_hash'], old.geometry_hash)
        self.assertEqual(json.loads(json.dumps(asdict(captured)))['material_proof_json'], captured.material_proof_json)

    def test_capture_rejects_policy_without_exact_material_tape(self):
        frozen = json.loads((Path(__file__).parent / 'same-input-job-03/output/frozen-inputs.json').read_text())
        payload = dict(frozen['request'])
        payload['settings'] = spin.SpinSolverConfig(**payload['settings'])
        settings = spin.SpinEvaluatorConfig(**dict(frozen['settings'], exact_binary64_material_boundary=False))
        with self.assertRaisesRegex(ValueError, 'exact.*material'):
            spin.evaluate_spin_it_faster_fullrank_quadric(spin.SpinItFasterRequest(**payload), frozen['states']['ellipsoid'], evaluator_settings=settings, capture_material_proof=True)


if __name__ == '__main__':
    unittest.main()
