"""Actual prescribed-waterline dispatch and sealed A-to-B method regression."""
import copy
import json
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest
import numpy as np
from test_exact_material_region import box,combine
from test_spin_it_faster_table1_runner import load_runner
from make_it_stand import spin

ROOT=Path(__file__).resolve().parents[2]
def mesh(value):return {"vertices":value[0].tolist(),"faces":value[1].tolist()}
def buoyant_case():
    outer=box((-2.,)*3,(2.,)*3)
    design=box((-1.9,)*3,(1.9,)*3)
    shell=combine(outer,(design[0],design[1][:,::-1]))
    # Input-first method fixture: same fluid, A draft0 and B draft1.2.
    # Neither liquid density nor initial planes are inferred from a solution.
    return {"id":"buoyant-method","mode":"BUOYANT_STABILITY","asset_class":"RECONSTRUCTED",
        "buoyant_problem":{"outer":mesh(outer),"design":{**mesh(design),"density":1.},
            "fixed_regions":[{"id":"shell",**mesh(shell),"density":.1}],
            "displacement":{"semantics":"SEALED_OUTER_HULL",**mesh(outer)},
            "liquid_density":.6,"waterline":{"origin":[0.,0.,0.],"normal":[0.,0.,1.]},
            "target_origin":[0.,0.,0.],"target_axes":np.eye(3).tolist(),
            "length_unit":"m","mass_unit":"kg","minimum_outer_wall":.05,
            "controls":{"maximum_iterations":60,"maximum_line_search_steps":24,"constraint_tolerance":1e-10},
            "initial_parameters":[0.,0.,0.],
            "stage_b":{"design_from":"A_CAVITY","waterline":{"origin":[0.,0.,1.2],"normal":[0.,0.,1.]},
                       "initial_parameters":[.8,0.,0.]}}}

class BuoyantTableDispatchTests(unittest.TestCase):
    def test_missing_actual_displacement_blocks_legacy_dispatch(self):
        runner=load_runner();case={"id":"missing","mode":"BUOYANT_STABILITY"}
        planned=runner.plan_payload([case],1)["cases"][0]
        self.assertEqual(planned["status"],"BLOCKED_MISSING_BUOYANT_GEOMETRY")
        with self.assertRaisesRegex(ValueError,"buoyant_problem"):runner.run_case(case,1,1)

    def test_initial_failed_native_geometry_is_present_in_evaluation_trace(self):
        case=buoyant_case()["buoyant_problem"]
        request=spin.SpinItFasterRequest(case["outer"]["vertices"],case["outer"]["faces"],
            case["design"]["vertices"],case["design"]["faces"])
        fixed=case["fixed_regions"][0]
        request=spin.with_measured_buoyancy_geometry(request,case["outer"]["vertices"],case["outer"]["faces"],
            fixed["vertices"],fixed["faces"],waterline_origin=(0.,0.,0.),waterline_normal=(0.,0.,1.),
            liquid_density=.6,fixed_shell_density=.1)
        result=spin.solve_spin_it_faster_buoyant_plane(request,initial_parameters=(1e308,1e308,1e308))
        self.assertEqual(result["status"],"GEOMETRY_EVALUATION_FAILED")
        self.assertEqual(len(result["trace"]),1)
        self.assertFalse(result["trace"][0]["evaluation_succeeded"])
        self.assertTrue(result["trace"][0]["diagnostic"])

    def test_real_a_then_b_saved_material_and_same_liquid(self):
        runner=load_runner();case=buoyant_case()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            result=runner.run_case(case,1,1,bundle_root=root)
            self.assertEqual(result["solver_family"],"BUOYANT_PLANE_GEOMETRY_V1")
            self.assertTrue(result["passed"],result)
            stages=result["runs"][0]["stages"]
            self.assertEqual([s["stage"] for s in stages],["A","B"])
            self.assertGreater(stages[0]["accepted_steps"],0)
            self.assertGreater(stages[1]["accepted_steps"],0)
            self.assertEqual(stages[1]["source_a_scientific_digest"],stages[0]["artifact_scientific_digest"])
            self.assertAlmostEqual(stages[0]["total_raw_moments"][0],19.2,places=8)
            self.assertAlmostEqual(stages[1]["total_raw_moments"][0],30.72,places=8)
            for s in stages:
                self.assertEqual(s["liquid_density"],.6)
                self.assertEqual(s["function_evaluations"],len(s["trace"]))
                self.assertTrue((root/s["relative_bundle"]/"manifest.json").is_file())
            from make_it_stand.verification.buoyant_material import audit_saved_buoyant_case
            audit=audit_saved_buoyant_case(result,bundle_root=root)
            self.assertEqual(audit["qualified_runs"],1)
            self.assertLessEqual(audit["maximum_normalized_moment_error"],1e-6)
            application=root/'application.json'
            application.write_text(json.dumps({'cases':[result]}))
            checked=subprocess.run([sys.executable,str(ROOT/'experiments/spin_it_reproduction/faster/check_saved_buoyant_application.py'),
                '--application',str(application),'--bundle-root',str(root),'--output',str(root/'checked')],
                capture_output=True,text=True,timeout=60)
            self.assertEqual(checked.returncode,0,checked.stderr)
            self.assertEqual(json.loads((root/'checked/scientific.json').read_text())['cases'][0]['qualified_runs'],1)
            for bad in ("density","source","geometry","trace","plane","eq7"):
                forged=copy.deepcopy(result)
                if bad=="density":forged["buoyant_problem"]["liquid_density"]=.7
                elif bad=="source":forged["runs"][0]["stages"][1]["source_a_scientific_digest"]="0"*64
                elif bad=="geometry":forged["runs"][0]["stages"][1]["total_raw_moments"][0]+=1.
                else:
                    a=forged['runs'][0]['stages'][0]
                    if bad=='trace':
                        a['trace']=a['trace'][:1];a['solver_result']['trace']=a['trace']
                        a['accepted_steps']=a['solver_result']['accepted_steps']=0
                        a['function_evaluations']=1
                    elif bad=='plane':a['solver_result']['plane_coefficients']=[0.,0.,0.,-1.]
                    else:a['solver_result']['eq7_numerical_conditions_verified']=False
                from make_it_stand.verification.spin_common import digest
                forged['buoyant_problem_sha256']=digest(forged['buoyant_problem'])
                forged['scientific_hash']=digest({k:v for k,v in forged.items() if k!='scientific_hash'})
                with self.subTest(bad=bad),self.assertRaises(ValueError):audit_saved_buoyant_case(forged,bundle_root=root)

    def test_no_fake_random_restarts_or_flooding(self):
        from make_it_stand.verification.buoyant_application import buoyant_case_plan
        with self.assertRaisesRegex(ValueError,'one attempt'):
            buoyant_case_plan(buoyant_case(),1,8,[0])
        case=buoyant_case();case['buoyant_problem']['displacement']['semantics']='FLOODABLE'
        with self.assertRaisesRegex(ValueError,'semantics'):
            buoyant_case_plan(case,1,1,[0])

    def test_failed_initial_retained_without_material_or_b_solve(self):
        from make_it_stand.verification.buoyant_material import audit_saved_buoyant_case
        case=buoyant_case();case['buoyant_problem']['initial_parameters']=[1e308]*3
        with tempfile.TemporaryDirectory() as directory:
            result=load_runner().run_case(case,1,1,bundle_root=Path(directory))
            self.assertFalse(result['passed'])
            self.assertEqual(len(result['runs'][0]['stages']),1)
            self.assertEqual(result['runs'][0]['stages'][0]['function_evaluations'],1)
            audit=audit_saved_buoyant_case(result,bundle_root=directory)
            self.assertEqual(audit['qualified_runs'],0)

    def test_b_seal_failure_does_not_erase_completed_solver_ledger(self):
        from unittest.mock import patch
        from make_it_stand.verification.buoyant_material import audit_saved_buoyant_case
        # Fault only at the geometry-seal boundary, AFTER the real B optimizer.
        # A adapter already holds its own imported seal function.
        with tempfile.TemporaryDirectory() as directory, patch(
                'make_it_stand.verification.faster_material.seal_buoyant_plane_result_as_artifact',
                side_effect=ValueError('injected B sealing rejection')):
            result=load_runner().run_case(buoyant_case(),1,1,bundle_root=Path(directory))
            self.assertFalse(result['passed'])
            stages=result['runs'][0]['stages']
            self.assertEqual(len(stages),2)
            self.assertEqual(stages[1]['status'],'CONVERGED')
            self.assertGreater(stages[1]['function_evaluations'],1)
            self.assertIn('injected B sealing rejection',stages[1]['seal_error']['message'])
            self.assertEqual(audit_saved_buoyant_case(result,bundle_root=directory)['qualified_runs'],0)

    def test_native_preflight_rejection_has_zero_evaluations(self):
        from make_it_stand.verification.buoyant_material import audit_saved_buoyant_case
        case=buoyant_case();case['buoyant_problem']['waterline']['origin']=[0.,0.,-1.99]
        with tempfile.TemporaryDirectory() as directory:
            result=load_runner().run_case(case,1,1,bundle_root=Path(directory))
            self.assertEqual(result['runs'][0]['status'],'INVALID_REQUEST')
            self.assertEqual(result['runs'][0]['stages'][0]['function_evaluations'],0)
            self.assertEqual(audit_saved_buoyant_case(result,bundle_root=directory)['qualified_runs'],0)

    def test_save_collision_preserves_solver_failure_record(self):
        from make_it_stand.verification.buoyant_material import audit_saved_buoyant_case
        with tempfile.TemporaryDirectory() as directory:
            runner=load_runner()
            first=runner.run_case(buoyant_case(),1,1,bundle_root=Path(directory))
            self.assertTrue(first['passed'])
            second=runner.run_case(buoyant_case(),1,1,bundle_root=Path(directory))
            self.assertFalse(second['passed'])
            stage=second['runs'][0]['stages'][0]
            self.assertEqual(stage['seal_error']['type'],'FileExistsError')
            self.assertGreater(stage['function_evaluations'],1)
            self.assertTrue(stage['attempted_bundle'])
            self.assertEqual(audit_saved_buoyant_case(second,bundle_root=directory)['qualified_runs'],0)

    def test_forged_pass_on_failed_run_rejected(self):
        from make_it_stand.verification.buoyant_material import audit_saved_buoyant_case
        from make_it_stand.verification.spin_common import digest
        case=buoyant_case();case['buoyant_problem']['initial_parameters']=[1e308]*3
        with tempfile.TemporaryDirectory() as directory:
            result=load_runner().run_case(case,1,1,bundle_root=Path(directory))
            result['passed']=True
            result['scientific_hash']=digest({k:v for k,v in result.items() if k!='scientific_hash'})
            with self.assertRaises(ValueError):audit_saved_buoyant_case(result,bundle_root=directory)

    def test_reverse_pose_is_in_target_frame(self):
        from make_it_stand.verification.buoyant_material import audit_saved_buoyant_case
        case=buoyant_case();p=case['buoyant_problem']
        axes=np.array([[0.,0.,1.],[0.,1.,0.],[-1.,0.,0.]])
        for m in [p['outer'],p['design'],p['displacement'],*p['fixed_regions']]:
            m['vertices']=(np.asarray(m['vertices'])@axes.T).tolist()
        p['target_axes']=axes.tolist()
        for s in [p,p['stage_b']]:
            s['waterline']['origin']=(np.asarray(s['waterline']['origin'])@axes.T).tolist()
            s['waterline']['normal']=axes[:,2].tolist()
        with tempfile.TemporaryDirectory() as directory:
            result=load_runner().run_case(case,1,1,bundle_root=Path(directory))
            self.assertTrue(result['passed'],result)
            pose=result['runs'][0]['stages'][1]['pose_checks'][1]['pose']
            np.testing.assert_allclose(np.asarray(pose['rotation'])@axes[:,2],-axes[:,2],atol=1e-14)
            self.assertEqual(audit_saved_buoyant_case(result,bundle_root=directory)['qualified_runs'],1)

if __name__=="__main__":unittest.main()
