"""Independent contracts for the reconstructed heart experiment, not paper credit."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import trimesh

PATH = Path(__file__).resolve().parents[2] / 'experiments/reconstructed_heart_top/heart_top.py'
module = None
if PATH.exists():
    spec = importlib.util.spec_from_file_location('heart_top', PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)


class HeartTopTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(module, 'real heart geometry and physical experiment are not implemented')

    def test_closed_heart_has_two_lobes_notch_and_one_bottom_tip(self):
        v, f = module.heart_mesh()
        m = trimesh.Trimesh(v, f, process=False)
        self.assertTrue(m.is_watertight and m.is_winding_consistent)
        self.assertEqual(m.euler_number, 2)
        self.assertGreater(m.volume, .1)
        self.assertLess(abs(v[:, 2].min()), 1e-14)
        self.assertEqual(np.count_nonzero(np.linalg.norm(v, axis=1) < 1e-12), 1)
        front = v[np.abs(v[:, 1]) < 1e-12]
        self.assertGreater(front[front[:, 0] > .25, 2].max(),
                           front[np.abs(front[:, 0]) < 1e-12, 2].max() + .2)
        np.testing.assert_allclose(module.integrate(v, f)[1:3], 0, atol=1e-13)

    def test_oriented_tetra_moments_match_translated_cube_analytic_values(self):
        m = trimesh.creation.box(extents=(2,2,2))
        v = m.vertices + [.4,-.2,1.3]
        s = module.integrate(v, m.faces)
        expected = [8,3.2,-1.6,10.4,8/3+1.28,8/3+.32,8/3+13.52,-2.08,4.16,-.64]
        np.testing.assert_allclose(s, expected, atol=1e-13)
        mass, c, I = module.properties(s)
        self.assertAlmostEqual(mass, 8, delta=1e-13)
        np.testing.assert_allclose(c,[.4,-.2,1.3],atol=1e-14)
        np.testing.assert_allclose(I, np.eye(3)*16/3 + 8*((c@c)*np.eye(3)-np.outer(c,c)),atol=1e-13)

    def test_torque_free_rotation_conserves_energy_and_world_momentum(self):
        r = module.simulate(1., np.zeros(3), np.diag([.2,.3,.4]),
                            [1.1,.7,5.], duration=1., gravity=0., tilt_degrees=2.)
        self.assertLess(r['relative_energy_drift'], 1e-8)
        self.assertLess(r['relative_momentum_drift'], 1e-8)

    def test_real_gravity_destabilizes_slow_top_but_not_fast_symmetric_top(self):
        args = (.1,np.array([0.,0.,.05]),np.diag([.002,.002,.001]))
        slow = module.simulate(*args,[0.,0.,5.],duration=3.)
        fast = module.simulate(*args,[0.,0.,40.],duration=3.)
        self.assertGreater(slow['max_tilt_degrees'], 30.)
        self.assertLess(fast['max_tilt_degrees'], 5.)
        self.assertLess(fast['relative_energy_drift'], 1e-7)

    def test_invalid_inertia_cannot_be_simulated(self):
        with self.assertRaises(ValueError):
            module.simulate(1.,np.zeros(3),np.diag([1.,0.,-1.]),[0.,0.,1.])


if __name__ == '__main__': unittest.main()
