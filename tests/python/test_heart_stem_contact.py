import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'experiments/heart_stem_contact'))


class StemTests(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec('stem_geometry'),'missing stem geometry implementation')
        import stem_geometry
        return stem_geometry

    def test_boolean_union_does_not_double_count_overlap(self):
        g=self.module(); a=trimesh.creation.box(extents=[2,2,2]); b=a.copy(); b.apply_translation([1,0,0])
        with tempfile.TemporaryDirectory() as d:
            u=g.boolean(a,b,'union',Path(d))
        self.assertTrue(u.is_watertight and u.is_winding_consistent)
        self.assertAlmostEqual(u.volume,12,places=10)

    def test_stem_is_connected_and_preserved_outside_design(self):
        g=self.module()
        with tempfile.TemporaryDirectory() as d:
            shapes=g.make_input(Path(d))
            outer,design,fixed,rod=[shapes[k] for k in ['outer','design','fixed','rod']]
            for mesh in [outer,design,fixed,rod]:
                self.assertTrue(mesh.is_watertight and mesh.is_winding_consistent)
            self.assertEqual(len(outer.split()),1)
            self.assertLess(outer.bounds[0,2],-.1)
            self.assertGreater(outer.bounds[1,2],2.0)
            self.assertAlmostEqual(outer.volume,design.volume+fixed.volume,places=9)
            intersection=g.boolean(design,rod,'intersect',Path(d))
            # Closed design and fixed shaft share a boundary. Triangle Boolean
            # may retain a zero-volume contact sheet; it must not retain volume.
            self.assertAlmostEqual(intersection.volume,0,delta=1e-12)

    def test_central_inertia_is_translated_and_diagonalized_correctly(self):
        self.assertIsNotNone(importlib.util.find_spec('bullet_contact'),'missing contact implementation')
        from bullet_contact import material_properties
        box=trimesh.creation.box(extents=[.02,.04,.06])
        R=Rotation.from_euler('xyz',[.4,.7,-.2]).as_matrix()
        box.vertices=box.vertices@R.T+[.01,-.02,.04]
        m,c,I,d,U=material_properties(np.asarray(box.vertices),np.asarray(box.faces),density=1200)
        self.assertAlmostEqual(m,.0576,places=12)
        np.testing.assert_allclose(c,[.01,-.02,.04],atol=1e-12)
        expected=R@np.diag([.0576*(.04**2+.06**2)/12,.0576*(.02**2+.06**2)/12,.0576*(.02**2+.04**2)/12])@R.T
        np.testing.assert_allclose(I,expected,rtol=1e-11,atol=1e-14)
        np.testing.assert_allclose(U@np.diag(d)@U.T,expected,rtol=1e-11,atol=1e-14)
        self.assertGreater(np.linalg.det(U),.999999)

    def test_bullet_frame_and_unforced_spin_are_not_fixed_or_clamped(self):
        import bullet_contact as b
        self.assertTrue(hasattr(b,'create_body'),'missing free-contact body setup')
        import pybullet as p
        box=trimesh.creation.box(extents=[.02,.04,.06]); box.apply_translation([.01,-.02,.04])
        R=Rotation.from_euler('xyz',[.2,-.3,.4]).as_matrix(); origin=np.array([.1,.2,.3])
        props=b.material_properties(box.vertices,box.faces,1200)
        client=p.connect(p.DIRECT)
        try:
            p.setGravity(0,0,0,physicsClientId=client); p.setTimeStep(1/4000,physicsClientId=client)
            with tempfile.TemporaryDirectory() as d:
                file=Path(d)/'cube.obj'; box.export(file)
                body=b.create_body(p,client,file,[file,file],props,R,origin,35.)
                for visual in p.getVisualShapeData(body,physicsClientId=client):
                    if visual[1]>=0: self.assertEqual(visual[7][3],0,'proxy must not cover actual outer visual')
                x,q=p.getBasePositionAndOrientation(body,physicsClientId=client)
                self.assertAlmostEqual(p.getDynamicsInfo(body,-1,physicsClientId=client)[0],props[0],places=12)
                np.testing.assert_allclose(x,origin+R@props[1],atol=1e-12)
                Q=Rotation.from_quat(q).as_matrix()
                np.testing.assert_allclose(Q@np.diag(p.getDynamicsInfo(body,-1,physicsClientId=client)[2])@Q.T,R@props[2]@R.T,atol=1e-12)
                self.assertEqual(p.getNumConstraints(physicsClientId=client),0)
                for _ in range(40): p.stepSimulation(physicsClientId=client)
                velocity=p.getBaseVelocity(body,physicsClientId=client)[1]
                self.assertAlmostEqual(np.linalg.norm(velocity),2*np.pi*35,delta=.01)
        finally: p.disconnect(client)

    def test_free_drop_produces_table_contact_without_a_constraint(self):
        import bullet_contact as b
        self.assertTrue(hasattr(b,'simulate_contact'),'missing contact simulation')
        box=trimesh.creation.box(extents=[.02,.04,.06]); box.apply_translation([0,0,.03])
        props=b.material_properties(box.vertices,box.faces,1200)
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'box.obj'; box.export(file)
            r=b.simulate_contact(file,[file],props,np.asarray(box.vertices),hz=0,tilt=0,azimuth=0,
                duration=.25,dt=1/2000,clearance=.005)
            self.assertGreater(r['maximum_contact_count'],0)
            self.assertEqual(r['constraint_count'],0)
            self.assertLess(r['maximum_penetration_m'],.0002)
            self.assertAlmostEqual(r['samples'][-1]['position'][2],.03,delta=.0002)

    def test_lowest_child_margin_cannot_be_reported_as_valid_support(self):
        import bullet_contact as b
        lower=trimesh.creation.box(extents=[.02,.04,.06]); lower.apply_translation([0,0,.03])
        upper=lower.copy(); upper.apply_translation([0,0,.04])
        props=b.material_properties(lower.vertices,lower.faces,1200)
        with tempfile.TemporaryDirectory() as d:
            lo=Path(d)/'lower.obj'; hi=Path(d)/'upper.obj'
            lower.export(lo); upper.export(hi)
            with self.assertRaisesRegex(ValueError,'COLLISION_MODEL_BLOCKER'):
                b.simulate_contact(lo,[hi,lo],props,np.asarray(lower.vertices),
                    hz=0,tilt=0,azimuth=0,duration=.001,clearance=.0002)

    def test_tip_base_compound_has_correct_runtime_table_gap(self):
        import bullet_contact as b
        lower=trimesh.creation.box(extents=[.02,.04,.06]); lower.apply_translation([0,0,.03])
        upper=lower.copy(); upper.apply_translation([0,0,.04])
        props=b.material_properties(lower.vertices,lower.faces,1200)
        with tempfile.TemporaryDirectory() as d:
            lo=Path(d)/'lower.obj'; hi=Path(d)/'upper.obj'
            lower.export(lo); upper.export(hi)
            r=b.simulate_contact(lo,[lo,hi],props,np.asarray(lower.vertices),
                hz=0,tilt=0,azimuth=0,duration=.01,clearance=.0002)
            self.assertTrue(r['collision_initial_support_passed'])
            self.assertAlmostEqual(r['initial_gap_m'],.0002,delta=.0001)


if __name__=='__main__': unittest.main()
