"""Guard source authentication and exact-material-to-physics conversion."""
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path
import numpy as np

class DemoTests(unittest.TestCase):
    def module(self):
        path=Path(__file__).with_name('demo.py')
        self.assertTrue(path.exists(), 'demonstration adapter missing')
        spec=importlib.util.spec_from_file_location('demo',path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        return module

    def test_modified_saved_payload_is_rejected(self):
        mod=self.module()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'payload';p.write_bytes(b'original')
            record={'name':'payload','size':8,'sha256':hashlib.sha256(b'original').hexdigest()}
            self.assertEqual(mod.authenticated_bytes(Path(d),record),b'original')
            p.write_bytes(b'modified')
            with self.assertRaises(ValueError):mod.authenticated_bytes(Path(d),record)

    def test_principal_frame_and_si_conversion_preserve_non_diagonal_inertia(self):
        mod=self.module()
        # Unit mass at COM zero: Q diag 2,3,4 with Qxy=.25.
        props=mod.to_physics([1,0,0,0,2,3,4,0,0,.25])
        self.assertAlmostEqual(props[0],.01875)
        expected=np.array([[7,-.25,0],[-.25,6,0],[0,0,5]])*.01875*.025**2
        np.testing.assert_allclose(props[2],expected,rtol=1e-14,atol=1e-20)
        np.testing.assert_allclose(props[4]@np.diag(props[3])@props[4].T,expected,rtol=1e-14,atol=1e-20)

if __name__=='__main__':unittest.main()
