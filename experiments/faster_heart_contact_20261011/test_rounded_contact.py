"""Synthetic backend validation, separate from reconstructed Heart outcomes."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation

from contact_matrix import scenarios
from rounded_contact import simulate_scene
from rounded_input import write_mesh


class ActualContactBackendTests(unittest.TestCase):
    def test_records_real_motion_using_non_diagonal_inertia(self):
        protocol = json.loads((Path(__file__).parent/'rounded-asymmetric-v2/contact-protocol.json').read_text())
        scene = scenarios(protocol)[0]
        dimensions = np.array([.02, .03, .04])
        center = np.array([.001, -.002, .02])
        axes = Rotation.from_euler('xyz', [12, 23, 34], degrees=True).as_matrix()
        mass = .02
        diagonal = mass/12*(np.sum(dimensions**2)-dimensions**2)
        # eigh ordering is not mandatory; the proper orthogonal principal frame
        # and the corresponding diagonal must reconstruct the actual inertia.
        inertia = axes@np.diag(diagonal)@axes.T
        mesh = trimesh.creation.box(extents=dimensions)
        mesh.vertices = mesh.vertices@axes.T+center
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'analytic-box.obj'
            write_mesh(path, mesh)
            props = mass, center, inertia, diagonal, axes
            result = simulate_scene(path, [path], props, mesh.vertices, scene, protocol)
        self.assertEqual(result['actual_steps'], 16000)
        self.assertEqual(result['records'][-1]['time_s'], 8.)
        self.assertEqual(result['world_constraint_count'], 0)
        self.assertFalse(result['persistent_drive_torque'])
        self.assertEqual(result['loaded_mass_kg'], mass)
        np.testing.assert_array_equal(result['loaded_inertia_kg_m2'], inertia)
        self.assertTrue(any(x['contacts'] for x in result['records']))
        self.assertNotEqual(result['records'][0]['model_quaternion'], result['records'][-1]['model_quaternion'])


if __name__ == '__main__':
    unittest.main()
