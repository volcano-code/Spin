"""CAD import must retain hollow shells, islands and disconnected material."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from test_exact_material_region import box, combine
from make_it_stand.verification.mass_moments import TriangleMesh, compute_mass_properties

ROOT=Path(__file__).resolve().parents[2]
PROBE=Path(os.environ.get('MAKE_IT_STAND_OCCT_STEP_PROBE',
    ROOT/'build/occt-buoyant-chain-20261001/occt_step_roundtrip'))


class NestedMaterialCadTests(unittest.TestCase):
    def run_mesh(self,mesh,root):
        v,f=mesh
        path=root/'input.obj'
        path.write_text(''.join('v %.17g %.17g %.17g\n'%tuple(p) for p in v)+
                        ''.join('f %d %d %d\n'%tuple(t+1) for t in f))
        return subprocess.run([str(PROBE),'--mesh',str(path),'1',str(root/'out.brep'),str(root/'out.step')],
                              capture_output=True,text=True,timeout=60)

    def test_nested_cavity_island_and_disconnected_solid_keep_signed_material(self):
        mesh=combine(box((-2,)*3,(2,)*3),box((-1,)*3,(1,)*3,inward=True),
                     box((-.25,)*3,(.25,)*3),box((4,0,0),(5,1,1)))
        props=compute_mass_properties(TriangleMesh(mesh[0].tolist(),mesh[1].tolist()),1.)
        q=props.second_moment
        expected=np.array([props.mass,*props.first_moment,q[0][0],q[1][1],q[2][2],q[1][2],q[0][2],q[0][1]])
        self.assertAlmostEqual(props.mass,57.125,places=12)
        with tempfile.TemporaryDirectory() as d:
            done=self.run_mesh(mesh,Path(d))
            self.assertEqual(done.returncode,0,done.stderr)
            record=json.loads(done.stdout)
            for key in ('original','brep','step'):
                np.testing.assert_allclose(record['raw_moments'][key],expected,rtol=0.,atol=1e-10)

    def test_nested_wrong_orientation_or_intersection_never_becomes_material(self):
        variants=(combine(box((-2,)*3,(2,)*3),box((-1,)*3,(1,)*3)),
                  combine(box((0,0,0),(2,2,2)),box((1,1,1),(3,3,3))),
                  combine(box((0,0,0),(1,1,1)),box((1,1,1),(2,2,2))))
        for mesh in variants:
            with tempfile.TemporaryDirectory() as d:
                root=Path(d);done=self.run_mesh(mesh,root)
                self.assertNotEqual(done.returncode,0)
                self.assertFalse((root/'out.step').exists())

    def test_stray_open_face_cannot_disappear_beside_closed_shell(self):
        mesh=combine(box(),(np.array([[3.,0.,0.],[4.,0.,0.],[3.,1.,0.]]),np.array([[0,1,2]])))
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);done=self.run_mesh(mesh,root)
            self.assertNotEqual(done.returncode,0,'an open component was silently discarded')
            self.assertFalse((root/'out.step').exists())


if __name__=='__main__':unittest.main()
