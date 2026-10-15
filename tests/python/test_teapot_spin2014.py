"""Protect planar closure, strategy dispatch, and stale-input admission."""
import sys
import unittest
from pathlib import Path
import numpy as np
import trimesh

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'experiments/teapot_spin2014'))
import teapot
import video


class TeapotTests(unittest.TestCase):
    def test_single_optimized_trajectory_is_not_mislabeled_original(self):
        self.assertTrue(callable(getattr(video,'comparison_scope',None)), 'missing comparison scope')
        self.assertEqual(video.comparison_scope(['original']),'ORIGINAL_BASELINE_ONLY')
        self.assertEqual(video.comparison_scope(['inertia']),'PARTIAL_COMPARISON')
        self.assertEqual(video.comparison_scope(['original','inertia','lowering','mass_lowering']),'FOUR_WAY_RECONSTRUCTED_COMPARISON')

    def test_contact_pose_removes_rotated_com_offset(self):
        self.assertTrue(callable(getattr(video,'trajectory_pose',None)), 'missing trajectory replay')
        q=[0,0,np.sqrt(.5),np.sqrt(.5)]
        record={'samples':[dict(time=t,model_quaternion=q,position=[1,2,3]) for t in [0.,1.]]}
        position,rotation=video.trajectory_pose(record,.5,[.1,.2,.3])
        np.testing.assert_allclose(position,[1.2,1.9,2.7],atol=1e-14)
        np.testing.assert_allclose(rotation.apply([1,0,0]),[0,1,0],atol=1e-14)
        with self.assertRaises(ValueError):video.trajectory_pose(record,2.,[.1,.2,.3])

    def test_replay_rejects_nonfinite_trajectory(self):
        self.assertTrue(callable(getattr(video,'trajectory_pose',None)), 'missing trajectory replay')
        record={'time':[0.,1.],'quaternions':[[0,0,0,1],[0,0,float('nan'),1]]}
        with self.assertRaises(ValueError):video.trajectory_pose(record,.5,[0,0,1])

    def test_exact_index_cleanup_removes_only_collapsed_faces(self):
        self.assertTrue(callable(getattr(teapot, 'canonical_indices', None)), 'missing exact index cleanup')
        cube=trimesh.creation.box(extents=[1e-14,1,1])
        v=np.vstack([cube.vertices,cube.vertices[0]])
        f=np.vstack([cube.faces,[0,8,1]])
        actual,audit=teapot.canonical_indices(trimesh.Trimesh(v,f,process=False))
        self.assertEqual(len(actual.vertices),8)
        self.assertEqual(len(actual.faces),12)
        self.assertEqual(audit['collapsed_faces_removed'],1)
        self.assertTrue(actual.is_watertight)
        self.assertAlmostEqual(actual.volume/1e-14,1.)

    def test_caps_open_box_without_moving_source_vertices(self):
        self.assertTrue(callable(getattr(teapot, 'cap_planar_boundaries', None)), 'missing planar closure')
        box = trimesh.creation.box()
        keep = box.face_normals[:, 2] < .9
        open_box = trimesh.Trimesh(box.vertices.copy(), box.faces[keep], process=False)
        closed, audit = teapot.cap_planar_boundaries(open_box)
        np.testing.assert_array_equal(closed.vertices[:len(box.vertices)], box.vertices)
        self.assertTrue(closed.is_watertight)
        self.assertTrue(closed.is_winding_consistent)
        self.assertAlmostEqual(closed.volume, 1.)
        self.assertEqual(audit['capped_loops'], 1)
        self.assertEqual(audit['added_triangles'], 4)

    def test_rejects_nonplanar_boundary(self):
        self.assertTrue(callable(getattr(teapot, 'cap_planar_boundaries', None)), 'missing planar closure')
        box = trimesh.creation.box()
        keep = box.face_normals[:, 2] < .9
        mesh = trimesh.Trimesh(box.vertices.copy(), box.faces[keep], process=False)
        mesh.vertices[7, 2] += .1
        with self.assertRaisesRegex(ValueError, 'NONPLANAR'):
            teapot.cap_planar_boundaries(mesh)

    def test_three_strategies_dispatch_to_distinct_native_objectives(self):
        self.assertTrue(callable(getattr(teapot, 'strategy_settings', None)), 'missing strategy dispatch')
        expected = {'inertia': 'TOP_INERTIA_ONLY', 'lowering': 'TOP_LOWERING_ONLY',
                    'mass_lowering': 'TOP_MASS_AND_LOWERING'}
        actual = {k: teapot.strategy_settings(k).objective for k in expected}
        self.assertEqual(actual, expected)
        with self.assertRaises(ValueError):
            teapot.strategy_settings('invented')

    def test_rounded_shaft_is_closed_and_tip_is_origin(self):
        self.assertTrue(callable(getattr(teapot, 'shaft_mesh', None)), 'missing shaft')
        shaft = teapot.shaft_mesh()
        self.assertTrue(shaft.is_watertight)
        self.assertTrue(shaft.is_winding_consistent)
        self.assertGreater(shaft.volume, 0.)
        self.assertAlmostEqual(shaft.vertices[:,2].min(), 0.)
        self.assertGreater(shaft.vertices[:,2].max(), 1.41)

    def test_rejects_changed_source_instead_of_accepting_old_manifest(self):
        import tempfile
        self.assertTrue(callable(getattr(teapot, 'check_files', None)), 'missing source gate')
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'source'; p.write_bytes(b'original')
            files = teapot.file_hashes([p]); teapot.check_files(files)
            p.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'STALE_RESULT'):
                teapot.check_files(files)


if __name__ == '__main__': unittest.main()
