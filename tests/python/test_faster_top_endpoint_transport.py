"""Objective identity and actual Cherchi evaluations; no endpoint success claim."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

from test_faster_cherchi_table_dispatch import locked_case, ROOT
from test_spin_it_faster_table1_runner import load_runner


class TopEndpointTransportTests(unittest.TestCase):
    def test_registered_objective_rejects_wrong_mode_and_extra_weight(self):
        table=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            case=locked_case(Path(directory),"SPINNING_TOP")
            case["spinning_top_objective"]="LOWER_CENTER_OF_MASS"
            plan=table.cherchi_case_plan(case,1,1,asset_root=directory)
            self.assertEqual(plan["spinning_top_objective"],"LOWER_CENTER_OF_MASS")
            self.assertEqual(plan["objective_definition"]["expression"],"c3")
            self.assertEqual(plan["objective_definition"]["constraints"],["l1","l2","Ic13","Ic23"])
            for change in ({"center_of_mass_weight":1.},{"spinning_top_objective":"UNKNOWN"}):
                with self.subTest(change=change),self.assertRaises(ValueError):
                    table.cherchi_case_plan({**case,**change},1,1,asset_root=directory)
            yoyo=locked_case(Path(directory)/"yoyo","YO_YO")
            with self.assertRaises(ValueError):
                table.cherchi_case_plan({**yoyo,"spinning_top_objective":"LOWER_CENTER_OF_MASS"},
                    1,1,asset_root=Path(directory)/"yoyo")

    def test_actual_cherchi_callback_evaluates_height_not_inertia_ratio(self):
        table=load_runner();runner=table._cherchi_runner()
        from make_it_stand.verification.faster_material import _transform_moments
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);locked_case(root/"input","SPINNING_TOP")
            manifest=root/"input/problem.json";data=json.loads(manifest.read_text())
            # Existing frozen method geometry used by the weighted-objective
            # regression. Explicit sphere start: dispatch check, NOT random success.
            data["target"]["origin"]=[0.,0.,-2.]
            q=np.array([.25-.125**2,0.,0.,-1.,1.,1.,1.,0.,0.,0.]);q/=np.linalg.norm(q)
            data["preflight_quadric_target_normalized"]=q.tolist()
            manifest.write_text(json.dumps(data))
            build=Path(runner.spin._native_module().__file__).resolve().parent.parent
            runner.run(root/"run",build=build,
                probe=ROOT/"build/cherchi-reference-singlethread-20260917/cherchi_probe",
                problem_manifest=manifest,spinning_top_objective="LOWER_CENTER_OF_MASS",
                maximum_function_evaluations=4,maximum_restarts=1,
                icosphere_subdivisions=4,defer_finalization=True)
            inputs=json.loads((root/"run/inputs.json").read_text())
            self.assertEqual(inputs["spinning_top_objective"],"LOWER_CENTER_OF_MASS")
            payload=runner._load_terminal_audit_checkpoint(root/"run/terminal-audit-checkpoint.json")
            bridge=runner._bridge_identity_attachment(payload["run_finalization"])
            self.assertEqual(bridge["identity"]["objective_definition"],inputs["objective_definition"])
            evaluated=[r for r in payload["run_finalization"]["native_result"]["events"]
                if r["type"]=="evaluation" and r.get("objective") is not None]
            self.assertTrue(evaluated,"must exercise actual Boolean/map/Newton evaluation")
            problem,_,_=runner._load_problem_manifest(manifest)
            _,request,_=runner._prepare(problem,
                runner.spin.SpinEvaluatorConfig(derivative_method="DISCRETE_FULL_RANK_QUADRIC"),
                runner.spin.SpinSolverConfig(equation_profile="DUAL_MONOMIAL_RAW_V1"))
            axes,origin=np.asarray(request.target_axes),np.asarray(request.target_origin)
            fixed=np.asarray(_transform_moments(request.fixed_shell_moments,axes.T,-axes.T@origin))
            for event in evaluated:
                geometry=next(r for r in payload["evaluations"] if r["quadric_hash"]==event["quadric_hash"])
                s=np.asarray(geometry["moment_map"]["raw_moments"])+fixed
                self.assertAlmostEqual(event["objective"],s[3]/s[0],places=12)

    def test_seeded_saved_failure_cannot_be_relabelled_as_another_objective(self):
        table=load_runner()
        import _make_it_stand as native
        from make_it_stand.verification.cherchi_application import audit_saved_cherchi_case
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=locked_case(root/"input","SPINNING_TOP")
            case["spinning_top_objective"]="LOWER_CENTER_OF_MASS"
            result=table.run_cherchi_case(case,1,1,bundle_root=root/"runs",asset_root=root/"input",
                initial_seeds=(0,),cherchi_build=Path(native.__file__).resolve().parent.parent,
                cherchi_probe=ROOT/"build/cherchi-reference-singlethread-20260917/cherchi_probe")
            self.assertEqual(result["initialization_plan"]["spinning_top_objective"],"LOWER_CENTER_OF_MASS")
            report=audit_saved_cherchi_case(result,bundle_root=root/"runs")
            self.assertEqual(report["expected_runs"],1)
            for change in ("selector","formula"):
                bad=copy.deepcopy(result)
                if change=="selector":bad["initialization_plan"]["spinning_top_objective"]="INERTIA_AND_LOWERING"
                else:bad["initialization_plan"]["objective_definition"]["expression"]="f_yo"
                bad["scientific_hash"]=table.scientific_hash({k:v for k,v in bad.items() if k!="scientific_hash"})
                with self.subTest(change=change),self.assertRaises(ValueError):
                    audit_saved_cherchi_case(bad,bundle_root=root/"runs")


if __name__=="__main__":unittest.main()
