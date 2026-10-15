import importlib.util
from pathlib import Path
import unittest
import numpy as np

class BearingTests(unittest.TestCase):
    def load(self):
        path=Path(__file__).with_name('bearing_physics.py')
        self.assertTrue(path.exists(), 'Missing actual bearing simulation implementation')
        spec=importlib.util.spec_from_file_location('bearing',path)
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        return mod

    def test_eccentric_rotor_reaction_and_unpowered_speed(self):
        # Removing the COM offset, leaving the default motor enabled or loading
        # a default inertia must break this real Bullet dynamics benchmark.
        mod=self.load()
        out=mod.simulate(.02,[.003,0,0],np.diag([1e-6,2e-6,2.5e-6]),dt=.00025,duration=.05,hz=10)
        self.assertAlmostEqual(out['axis_inertia'],2.68e-6,places=14)
        self.assertLess(out['maximum_relative_energy_drift'],1e-8)
        self.assertLess(abs(out['records'][-1]['omega']-20*np.pi),1e-8)
        force=np.array([r['reaction'][:2] for r in out['records'][1:]])
        self.assertAlmostEqual(np.linalg.norm(force,axis=1).mean(),.236870505626,delta=.001)

    def test_off_diagonal_inertia_is_not_discarded(self):
        mod=self.load();I=np.array([[2e-6,.3e-6,0],[.3e-6,1e-6,0],[0,0,2.5e-6]])
        out=mod.simulate(.02,[.002,.001,0],I,dt=.0005,duration=.01,hz=5)
        np.testing.assert_allclose(out['loaded_inertia'],I,rtol=1e-12,atol=1e-18)

    def test_registered_35_hz_is_not_clipped_by_engine_limit(self):
        out=self.load().simulate(.02,[.003,0,0],np.diag([1e-6,2e-6,2.5e-6]),dt=.00025,duration=.01,hz=35)
        self.assertAlmostEqual(out['records'][1]['omega'],219.9114857512855,places=8)

if __name__=='__main__':unittest.main()
