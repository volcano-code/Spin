"""Fixed-normal source must survive request, Newton and saved-map admission."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from fractions import Fraction as Q

import numpy as np
from make_it_stand import spin
from test_faster_cherchi_newton_integration_runner import RUNNER, ROOT


class FixedNormalRunnerTests(unittest.TestCase):
    def test_explicit_config_drives_native_source_and_rejects_mixed_profiles(self):
        self.assertIn('fixed_normal_paraboloid', spin.SpinEvaluatorConfig.__dataclass_fields__)
        config=spin.SpinEvaluatorConfig(derivative_method='DISCRETE_REDUCED_STRATUM',
            fixed_normal_paraboloid=True,icosphere_subdivisions=2,truncation_radius=4.)
        native=spin._native_module()
        settings=spin._native_evaluator_settings(native,config).triangulation
        self.assertTrue(settings.fixed_normal_paraboloid)
        q=np.array([.125,.3,.7,-1,0,0,0,.2,.4,0.]);q/=np.linalg.norm(q)
        source=RUNNER._pure_height_source_motion(native,q,np.array([[-.2]*3,[.2]*3]),settings)
        branch,features=RUNNER._source_derivative_branch(source)
        self.assertEqual(branch['source_map'],'BOUNDED_FINITE_NORMAL_PARABOLOID_V1')
        self.assertEqual(branch['initial_icosphere_subdivisions'],2)
        self.assertEqual(np.asarray(source['vertex_velocity']).shape[1],5)
        self.assertEqual(len(RUNNER._derivative_stratum_identity('intersection','a'*64,branch)),64)
        self.assertEqual(len(features['normals']),162)
        self.assertEqual(RUNNER._finite_normal_feature_digest(features),branch['finite_normal_features_sha256'])
        changed=copy.deepcopy(features);changed['normal_vertices'][0]+=1
        with self.assertRaises(ValueError):RUNNER._finite_normal_feature_digest(changed)
        changed=copy.deepcopy(features);changed['normals'][0][0]+=.01
        with self.assertRaises(ValueError):RUNNER._finite_normal_feature_digest(changed)
        # Independent exact face-normal separation for represented slanted
        # closures, including faces whose axis-aligned boxes overlap the domain.
        vertices=np.asarray(source['vertices']);faces=np.asarray(source['faces'])
        corners=np.array([[.2 if mask&(1<<j) else -.2 for j in range(3)] for mask in range(8)])
        missed_by_aabb=0
        for face in faces[np.asarray(source['face_is_cap'],dtype=bool)]:
            tri=vertices[face]
            if np.any(np.min(tri,axis=0)>.2) or np.any(np.max(tri,axis=0)<-.2):continue
            missed_by_aabb+=1
            points=[[Q(float(x)) for x in row] for row in tri]
            u=[points[1][j]-points[0][j] for j in range(3)]
            v=[points[2][j]-points[0][j] for j in range(3)]
            normal=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
            project=lambda row:sum(a*Q(float(b)) if not isinstance(b,Q) else a*b for a,b in zip(normal,row))
            a=list(map(project,points));b=list(map(project,corners))
            self.assertTrue(max(a)<min(b) or max(b)<min(a))
        self.assertGreater(missed_by_aabb,0)
        for changes in ({'paper_projective_hyperboloid':True},
                        {'algebraic_finite_map':True}, {'maximum_icosphere_subdivisions':3},
                        {'derivative_method':'DISCRETE_FULL_RANK_QUADRIC'},
                        {'fixed_normal_paraboloid':1}):
            with self.subTest(changes=changes),self.assertRaises(ValueError):
                spin.SpinEvaluatorConfig(**{**config.__dict__,**changes})

    def test_real_seeded_run_records_valid_map_without_promoting_budget_exit(self):
        self.assertIn('fixed_normal_paraboloid',spin.SpinEvaluatorConfig.__dataclass_fields__)
        native=spin._native_module()
        build=Path(native.__file__).resolve().parent.parent
        fixture=ROOT/'experiments/spin_it_reproduction/faster/fixtures/cherchi-centered-axis-pure-height-method-v1/problem.json'
        probe=ROOT/'build/cherchi-reference-singlethread-20260917/cherchi_probe'
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'run'
            RUNNER.run(output,build=build,probe=probe,problem_manifest=fixture,
                maximum_function_evaluations=1,maximum_restarts=1,maximum_radius_doublings=4,
                icosphere_subdivisions=2,fixed_normal_paraboloid=True,
                spinning_top_objective='LOWER_CENTER_OF_MASS',
                directional_profile='PURE_HEIGHT_ACTIVE_STRATUM_V1',
                initialization_mode='supported_mode_seed',deterministic_seed=0,defer_finalization=True)
            inputs=json.loads((output/'inputs.json').read_text())
            self.assertEqual(inputs['paper_source_status'],'METHOD_BOUNDED_FINITE_NORMAL_NOT_COMPLETE_IDEAL_MAP')
            np.testing.assert_array_equal(inputs['initial_quadrics'],RUNNER._pure_height_initial_quadrics(0,1))
            saved=RUNNER._load_terminal_audit_checkpoint(output/'terminal-audit-checkpoint.json')
            result=saved['run_finalization']['native_result']
            self.assertEqual(result['function_evaluations'],1)
            maps=saved['evaluations']
            self.assertEqual(len(maps),1)
            # Exact local cap certification removes the former false rejection.
            # The same frozen start now yields a real map, NOT a solved object.
            self.assertNotEqual(maps[0].get('outcome'),'failed')
            self.assertEqual(maps[0]['jacobian_shape'],[10,5])
            self.assertLessEqual(maps[0]['independent_moment_audit']['normalized_inf_residual'],1e-6)
            self.assertEqual(maps[0]['solver_evaluation_index'],1)
            self.assertIsNotNone(maps[0]['geometry_hash'])
            self.assertNotEqual(result['status'],'CONVERGED')
            self.assertFalse(saved['final_geometry_verified'])
            self.assertEqual(result['supported_space']['evaluations'][0]['error'],'')

    def test_saved_genuinely_uncovered_source_remains_rejected(self):
        # Seed 1's first actual uncovered trial (zero-based evaluation 18).
        # Independent rational clipping proves nonzero missing area; do not
        # confuse this with the repaired cap or polygon-kernel false negatives.
        p=[.10945843568001104,.4637369557466022,-.11417428584951093,-.4664538241279329,
           0,0,0,.4893755887946585,-.5503303391894924,0]
        domain=np.array([[.25 if i&(1<<j) else -.25 for j in range(3)] for i in range(8)])
        native=spin._native_module();settings=native.SpinTriangulationSettings()
        settings.fixed_normal_paraboloid=True;settings.icosphere_subdivisions=2
        settings.truncation_radius=4.;settings.maximum_radius_doublings=4
        with self.assertRaisesRegex(RuntimeError,'finite patch does not cover design hull'):
            RUNNER._pure_height_source_motion(native,p,domain,settings)

    def test_fixed_normal_mode_cannot_silently_select_full_rank_or_adaptive_normals(self):
        self.assertIn('fixed_normal_paraboloid',spin.SpinEvaluatorConfig.__dataclass_fields__)
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'absent'
            with self.assertRaises(ValueError):
                RUNNER.run(output,build=Path('/absent'),probe=Path('/absent'),fixed_normal_paraboloid=True)
            self.assertFalse(output.exists())

    def test_source_and_saved_branch_reject_incompatible_map_flags(self):
        native=spin._native_module();settings=native.SpinTriangulationSettings()
        settings.fixed_normal_paraboloid=True;settings.icosphere_subdivisions=2
        settings.truncation_radius=16.
        p=json.loads((ROOT/'experiments/spin_it_reproduction/faster/fixtures/cherchi-fixed-normal-paraboloid-method-v1/problem.json').read_text())
        source=dict(native.evaluate_spin_it_faster_stratum_source_motion(p['coefficients'],p['design_vertices'],settings))
        source.update(directional_profile='PURE_HEIGHT_ACTIVE_STRATUM_V1',paraboloid_grid_resolution=8)
        branch,_=RUNNER._source_derivative_branch(source)
        for name in ('paper_projective_hyperboloid','algebraic_finite_map'):
            with self.subTest(name=name),self.assertRaises(ValueError):
                RUNNER._source_derivative_branch({**source,name:True})
        for mix in ({'finite_point_map':'ALGEBRAIC_V1'},{'paper_source_features_sha256':'a'*64}):
            with self.subTest(mix=mix),self.assertRaises(ValueError):
                RUNNER._derivative_stratum_identity('intersection','b'*64,{**branch,**mix})
