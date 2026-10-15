import unittest
from pathlib import Path
import numpy as np
import pybullet as p
import flower_table_scene as scene


class TableSceneTests(unittest.TestCase):
    def test_material_inertia_survives_pybullet_joint_configuration(self):
        bundle = Path(__file__).parent/'paper-flower-method-v2/table-scene-v1'
        packet = scene.load_bundle(bundle)
        client, bodies, diagnostics = scene.build_world(bundle, packet, .0005, visuals=False)
        try:
            for name, body in bodies.items():
                inertia = np.array(packet['properties'][name]['inertia_kg_m2'])
                self.assertGreater(abs(inertia[0, 1]), 5e-9)
                np.testing.assert_allclose(diagnostics[name]['loaded_inertia_kg_m2'], inertia, rtol=1e-12, atol=1e-18)
            for _ in range(200):
                p.stepSimulation(physicsClientId=client)
            for name, body in bodies.items():
                angle, omega, reaction, torque = p.getJointState(body, 0, physicsClientId=client)
                self.assertEqual(torque, 0)
                self.assertAlmostEqual(omega, 2*np.pi*35, places=8)
                self.assertLess(abs(np.linalg.norm(reaction[:2])-diagnostics[name]['analytic_imbalance_force_N']), 1e-6)
        finally:
            p.disconnect(client)

    def test_complete_assembly_preserves_all_nine_parts(self):
        import flower_table_delivery as delivery
        bundle = Path(__file__).parent/'paper-flower-method-v2/table-scene-v1'
        import json
        checks = json.loads((bundle/'assembly-check.json').read_text())
        check = checks['flower-complete.obj']
        self.assertEqual(len(check['parts']), 9)
        self.assertNotIn('cavity', check['parts'])
        self.assertNotIn('design', check['parts'])
        self.assertEqual(delivery.sha(bundle/'flower-complete.obj'), check['sha256'])


if __name__ == '__main__':
    unittest.main()
