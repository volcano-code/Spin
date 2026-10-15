"""Collision-query surface must match returned hull, including fixed links."""
import importlib,os,sys,tempfile,unittest
from pathlib import Path
import numpy as np
import pybullet as p
import trimesh
from scipy.spatial.transform import Rotation
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/heart_stem_contact'))
create_body=importlib.import_module(os.environ.get('CONSTRUCTOR_MODULE','bullet_contact_inertial_frame')).create_body

class LinkSupport(unittest.TestCase):
    def test_all_link_support_queries_match_engine_hulls(self):
        with tempfile.TemporaryDirectory() as directory:
            paths=[]
            for i in range(3):
                mesh=trimesh.creation.icosphere(subdivisions=1,radius=.008)
                mesh.apply_translation([i*.015,0,i*.005])
                path=Path(directory)/f'part{i}.obj';mesh.export(path);paths.append(path)
            for angles in [[0,0,0],[13,29,47]]:
                U=Rotation.from_euler('xyz',angles,degrees=True).as_matrix()
                R=Rotation.from_euler('xyz',[11,17,23],degrees=True).as_matrix()
                c=np.array([.001,-.002,.003]);d=np.array([3e-6,4e-6,5e-6])
                client=p.connect(p.DIRECT)
                try:
                    body=create_body(p,client,paths[0],paths,(.025,c,U@np.diag(d)@U.T,d,U),R,np.array([0,0,.05]),0)
                    floor=p.createMultiBody(0,p.createCollisionShape(p.GEOM_PLANE))
                    for link in range(-1,p.getNumJoints(body)):
                        with self.subTest(angles=angles,link=link):
                            vertices=np.asarray(p.getMeshData(body,link,flags=p.MESH_DATA_SIMULATION_MESH)[1])
                            pos,q=p.getBasePositionAndOrientation(body) if link==-1 else p.getLinkState(body,link,computeForwardKinematics=1)[:2]
                            Q=Rotation.from_quat(q).as_matrix()
                            support=float(np.min(vertices@Q[2,:]+pos[2]))
                            contacts=p.getClosestPoints(body,floor,.2,linkIndexA=link)
                            self.assertTrue(contacts)
                            self.assertLess(abs(support-min(cp[8] for cp in contacts)),1e-9)
                finally:p.disconnect(client)

if __name__=='__main__':unittest.main()
