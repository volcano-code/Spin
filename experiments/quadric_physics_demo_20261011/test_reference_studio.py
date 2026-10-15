"""Guard actual source identities and material roles, not image styling."""
import hashlib
from pathlib import Path
import tempfile
import unittest

class StudioSourceTests(unittest.TestCase):
    def test_modified_mesh_is_rejected_before_render(self):
        from reference_studio import obj_part
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'mesh.obj';p.write_text('v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n')
            digest=hashlib.sha256(p.read_bytes()).hexdigest()
            p.write_text(p.read_text()+'v 9 9 9\n')
            with self.assertRaisesRegex(ValueError,'SOURCE_MESH_IDENTITY_MISMATCH'):
                obj_part(p,digest,'fixed',1)

    def test_fixed_tubes_keep_fixed_role_and_original_geometry(self):
        from reference_studio import obj_part
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'mesh.obj';p.write_text('v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n')
            d=obj_part(p,hashlib.sha256(p.read_bytes()).hexdigest(),'fixed',0.025)
            self.assertEqual(d['role'],'fixed')
            self.assertEqual(d['vertices'],[[0.,0.,0.],[.025,0.,0.],[0.,.025,0.]])
            self.assertEqual(d['faces'],[[0,1,2]])

if __name__=='__main__':unittest.main()
