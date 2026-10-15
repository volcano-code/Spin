import unittest
import numpy as np
import trimesh
from paper_flower_saved_check import integrate_saved_mesh,verify_native_identity,validate_run_identity


class SavedChecks(unittest.TestCase):
    def test_shared_wall_rejects_twisted_top_and_unpaired_half(self):
        import paper_flower_saved_check as check
        proof={'allowed_directed_wall_edges':{'retained':{((0.,0.),(1.,0.))}}}
        # Old checker examined only the valid bottom edge and skipped the
        # one-bottom/two-top half, accepting a twisted wall.
        mesh=trimesh.Trimesh(vertices=[[0,0,0],[1,0,0],[1,1,1],[0,1,1]],
                              faces=[[0,1,2],[0,2,3]],process=False)
        with self.assertRaisesRegex(ValueError,'SHARED_WALL'):
            check.check_shared_wall_sources(mesh,'retained',proof,[0,1],1.)
        half=trimesh.Trimesh(vertices=[[0,0,0],[1,0,0],[1,0,1]],faces=[[0,1,2]],process=False)
        with self.assertRaisesRegex(ValueError,'SHARED_WALL'):
            check.check_shared_wall_sources(half,'retained',proof,[0,1],1.)

    def test_shared_wall_accepts_both_oriented_vertical_diagonals(self):
        import paper_flower_saved_check as check
        proof={'allowed_directed_wall_edges':{'retained':{((0.,0.),(1.,0.))}}}
        points=[[0,0,0],[1,0,0],[1,0,1],[0,0,1]]
        for faces in [[[0,1,2],[0,2,3]],[[0,1,3],[1,2,3]]]:
            mesh=trimesh.Trimesh(vertices=points,faces=faces,process=False)
            check.check_shared_wall_sources(mesh,'retained',proof,[0,1],1.)

    def test_independent_tetrahedra_non_diagonal_box(self):
        mesh=trimesh.creation.box(extents=[2.,4.,6.]);mesh.apply_translation([1.,2.,3.])
        expected=2*np.array([48.,48.,96.,144.,64.,256.,576.,288.,144.,96.])
        np.testing.assert_allclose(integrate_saved_mesh(mesh,2.),expected,rtol=1e-14)

    def test_wrong_native_module_sha_is_rejected(self):
        from types import SimpleNamespace
        with self.assertRaisesRegex(ValueError,'NATIVE_BINARY_IDENTITY'):
            verify_native_identity(SimpleNamespace(__file__=__file__),'0'*64)

    def test_diagnostic_or_no_update_cannot_be_success(self):
        prep=dict(circle_center=[0.,0.],circle_radius=1.,selected_seed=0)
        case=dict(planar_problem={'x':1},attempts=[dict(seed=0,status='CONVERGED',solver_result=dict(accepted_steps=0,center=[0.,0.],radius=1.))])
        with self.assertRaisesRegex(ValueError,'NO_ACTUAL_ACCEPTED_UPDATE'):
            validate_run_identity(case,{'x':1},prep)

    def test_changed_circle_is_rejected_even_with_hand_passed_flag(self):
        prep=dict(circle_center=[1.,0.],circle_radius=1.,selected_seed=0,moment_gate_pass=True)
        case=dict(planar_problem={'x':1},attempts=[dict(seed=0,status='CONVERGED',solver_result=dict(accepted_steps=1,center=[0.,0.],radius=1.))])
        with self.assertRaisesRegex(ValueError,'SAVED_CIRCLE_IDENTITY'):
            validate_run_identity(case,{'x':1},prep)

    def test_actual_saved_preflight_never_grants_geometry_admission(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path
        from types import SimpleNamespace
        from unittest.mock import patch
        import json
        import paper_flower_saved_check as check
        root=Path(__file__).resolve().parent
        folder=root/'paper-flower-method-v1'
        source=folder/'derived-output-v5'
        application=folder/'solve-job-01/output/application.json'
        self.assertTrue(application.is_file())
        with TemporaryDirectory() as directory:
            derived=Path(directory)/'derived';derived.mkdir()
            prep=json.loads((source/'output-preparation.json').read_text())
            prep['moment_gate_pass']=False  # This caller-owned field is ignored.
            (derived/'output-preparation.json').write_text(json.dumps(prep))
            for info in prep['files'].values():(derived/info['file']).symlink_to(source/info['file'])
            native=SimpleNamespace(__file__=__file__)
            with patch.dict('sys.modules',{'_make_it_stand':native}):
                check.run_check(folder,application,derived,root/'flower-grip-preview-v1/assembly.json',check.sha(__file__),Path(directory)/'report',preflight_only=True)
            result=json.loads((Path(directory)/'report/preflight-check.json').read_text())
            self.assertFalse(result['simulation_admission'])
            self.assertEqual(result['geometry_check'],'NOT_EXECUTED_IN_PREFLIGHT')
            self.assertLess(max(result['moment_errors']['optimized']),1e-6)


if __name__=='__main__':unittest.main()
