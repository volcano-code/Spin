"""Portable synthetic contract tests; never an original-paper experiment.

Run directly with Python. Missing pybullet/trimesh/scipy is an error, not a skip.
Set CONSTRUCTOR_MODULE=bullet_contact to exercise the legacy regression witness.
"""
import importlib
import os
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pybullet as p
import trimesh
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/heart_stem_contact'))
create_body=importlib.import_module(os.environ.get('CONSTRUCTOR_MODULE','bullet_contact_common_frame')).create_body

class CommonFrameContract(unittest.TestCase):
    def test_collision_vertices_do_not_depend_on_inertial_frame(self):
        # Changing COM/principal axes must not change the collision outer shape.
        with tempfile.TemporaryDirectory() as directory:
            mesh=trimesh.creation.icosphere(subdivisions=2,radius=.02)
            mesh.vertices[:,0]*=1.3
            path=Path(directory)/'convex.obj';mesh.export(path)
            R=Rotation.from_euler('xyz',[5,7,19],degrees=True).as_matrix()
            origin=np.array([.03,-.02,.06]);worlds=[]
            for angles,c in [([0,0,0],[0,0,0]),([13,29,47],[.001,-.002,.003]),([-27,11,63],[-.002,.001,.004])]:
                U=Rotation.from_euler('xyz',angles,degrees=True).as_matrix();c=np.array(c)
                d=np.array([3e-6,4e-6,5e-6]);I=U@np.diag(d)@U.T;m=.025
                client=p.connect(p.DIRECT)
                try:
                    body=create_body(p,client,path,[path],(m,c,I,d,U),R,origin,20)
                    pos,q=p.getBasePositionAndOrientation(body);Q=np.array(p.getMatrixFromQuaternion(q)).reshape(3,3)
                    points=np.array(p.getMeshData(body,-1,flags=p.MESH_DATA_SIMULATION_MESH)[1])
                    worlds.append(points@Q.T+pos)
                    info=p.getDynamicsInfo(body,-1)
                    self.assertEqual(info[0],m)
                    np.testing.assert_allclose(Q@np.diag(info[2])@Q.T,R@I@R.T,rtol=1e-12,atol=1e-16)
                    np.testing.assert_allclose(pos,origin+R@c,rtol=0,atol=1e-12)
                    linear,angular=p.getBaseVelocity(body);omega=R@np.array([0,0,40*np.pi])
                    np.testing.assert_allclose(angular,omega,rtol=0,atol=1e-12)
                    np.testing.assert_allclose(linear,np.cross(omega,R@c),rtol=0,atol=1e-12)
                    self.assertEqual(p.getNumConstraints(),0)
                finally:p.disconnect(client)
            for points in worlds[1:]:
                error=max(max(cKDTree(worlds[0]).query(points)[0]),max(cKDTree(points).query(worlds[0])[0]))
                self.assertLess(error,1e-9)

    def test_empty_collision_parts_rejected(self):
        with self.assertRaisesRegex(ValueError,'missing convex collision parts'):
            create_body(p,0,'unused',[],(1.,np.zeros(3),np.eye(3),np.ones(3),np.eye(3)),np.eye(3),np.zeros(3),0)

if __name__=='__main__':unittest.main()
