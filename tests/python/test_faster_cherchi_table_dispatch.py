"""Table top/yo-yo must use locked Cherchi, never the legacy fallback."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from test_spin_it_faster_table1_runner import load_runner

ROOT=Path(__file__).resolve().parents[2]

def locked_case(root, mode="YO_YO"):
    source=ROOT/"experiments/spin_it_reproduction/faster/fixtures/cherchi-centered-axis-fixed-shell-v1"
    manifest=json.loads((source/"problem.json").read_text())
    root.mkdir(parents=True,exist_ok=True)
    assets=[]
    def copy_mesh(raw,name):
        destination=root/name
        shutil.copyfile((source/raw).resolve(),destination)
        assets.append({"path":name,"sha256":hashlib.sha256(destination.read_bytes()).hexdigest()})
        return name
    manifest["outer_mesh"]=copy_mesh(manifest["outer_mesh"],"outer.obj")
    for kind in ("fixed_regions","design_regions"):
        for i,region in enumerate(manifest[kind]):
            region["mesh"]=copy_mesh(region["mesh"],f"{kind}-{i}.obj")
    manifest["target"]["mode"]=mode
    path=root/"problem.json"
    path.write_text(json.dumps(manifest,sort_keys=True)+"\n")
    return {"id":"locked-method", "mode":mode,"asset_class":"SYNTHETIC_REGRESSION",
        "cherchi_problem":{"path":"problem.json","sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                           "assets":assets},
        "cherchi_settings":{"icosphere_subdivisions":0,"maximum_radius_doublings":8,
                            "maximum_function_evaluations":1},
        "center_of_mass_weight":0.}

def strict_audit():
    spec=importlib.util.spec_from_file_location("strict_cherchi_support",ROOT/"experiments/spin_it_reproduction/strict/audit.py")
    module=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(ROOT/"experiments/spin_it_reproduction/strict"))
    try:spec.loader.exec_module(module)
    finally:sys.path.pop(0)
    return module

class CherchiTableDispatchTests(unittest.TestCase):
    def test_explicit_paper_map_reaches_native_run_and_cannot_be_relabelled(self):
        import _make_it_stand as native
        from make_it_stand.verification.cherchi_application import audit_saved_cherchi_case
        runner=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=locked_case(root/"inputs")
            case["cherchi_settings"]["paper_projective_hyperboloid"]=True
            result=runner.run_case(case,1,1,bundle_root=root/"runs",asset_root=root/"inputs",
                cherchi_build=Path(native.__file__).resolve().parent.parent,
                cherchi_probe=ROOT/"build/cherchi-reference-singlethread-20260917/cherchi_probe")
            inputs=json.loads((root/"runs/locked-method/run-000000/inputs.json").read_text())
            self.assertIs(inputs["paper_projective_hyperboloid"],True)
            self.assertEqual(result["summary"]["total_function_evaluations_all_attempts"],1)
            self.assertFalse(result["passed"],"a one-evaluation budget is not an optimization success")
            audit_saved_cherchi_case(result,bundle_root=root/"runs")
            for change in ("false", "missing"):
                altered=copy.deepcopy(result)
                settings=altered["initialization_plan"]["settings"]
                if change=="false":settings["paper_projective_hyperboloid"]=False
                else:settings.pop("paper_projective_hyperboloid")
                altered["scientific_hash"]=runner.scientific_hash(
                    {k:v for k,v in altered.items() if k!="scientific_hash"})
                with self.subTest(change=change),self.assertRaisesRegex(ValueError,"source|paper_projective"):
                    audit_saved_cherchi_case(altered,bundle_root=root/"runs")

    def test_paper_source_opt_in_rejects_non_boolean_without_changing_legacy_plan(self):
        runner=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=locked_case(root)
            old=runner.cherchi_case_plan(case,1,1,asset_root=root)
            self.assertNotIn("paper_projective_hyperboloid",old["settings"])
            for value in (True,False):
                selected=copy.deepcopy(case)
                selected["cherchi_settings"]["paper_projective_hyperboloid"]=value
                plan=runner.cherchi_case_plan(selected,1,1,asset_root=root)
                self.assertIs(plan["settings"]["paper_projective_hyperboloid"],value)
            for value in (1,"true",None):
                selected=copy.deepcopy(case)
                selected["cherchi_settings"]["paper_projective_hyperboloid"]=value
                with self.subTest(value=value),self.assertRaises(ValueError):
                    runner.cherchi_case_plan(selected,1,1,asset_root=root)

    def test_independent_native_run_files_cannot_be_reused_through_nested_root_aliases(self):
        strict=strict_audit()
        self.assertTrue(callable(getattr(strict,"_assert_distinct_cherchi_native_runs",None)))
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);nested=root/"nested";actual=nested/"case/run-000000"
            actual.mkdir(parents=True)
            (actual/"evidence.json").write_text("filesystem-identity-fixture-not-valid-scientific-evidence\n")
            (root/"case").symlink_to(nested/"case",target_is_directory=True)
            application={"cases":[{"solver_family":"CHERCHI_RAW_EQ14_V1","runs":[
                {"status":"CONVERGED","relative_output":"case/run-000000"}]}]}
            paths=[root/"a.json",root/"b.json"]
            for path in paths:path.write_text(json.dumps(application))
            with self.assertRaisesRegex(ValueError,"same native run|reused native"):
                strict._assert_distinct_cherchi_native_runs(paths,[root,nested])

    def test_fresh_process_audit_loads_bound_runtime_before_freshness_check(self):
        import _make_it_stand as native
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=locked_case(root/"inputs")
            result=load_runner().run_case(case,1,1,bundle_root=root/"runs",asset_root=root/"inputs",
                cherchi_build=Path(native.__file__).resolve().parent.parent,
                cherchi_probe=ROOT/"build/cherchi-reference-singlethread-20260917/cherchi_probe")
            application=root/"application.json";application.write_text(json.dumps({"cases":[result]}))
            command=[sys.executable,str(ROOT/"experiments/spin_it_reproduction/faster/check_saved_cherchi_application.py"),
                "--application",str(application),"--bundle-root",str(root/"runs")]
            for i in range(2):
                completed=subprocess.run(command+["--output",str(root/f"audit-{i}")],capture_output=True,text=True)
                self.assertEqual(completed.returncode,0,completed.stderr)
            code='''import sys
sys.path.insert(0,sys.argv[1])
import audit
try: audit._supporting_checks({"faster.cherchi_saved_material":sys.argv[2:]})
except ValueError as error:
    assert "no independently qualified" in str(error), str(error)
else: raise AssertionError("failed native solve granted material support")
'''
            completed=subprocess.run([sys.executable,"-c",code,str(ROOT/"experiments/spin_it_reproduction/strict"),
                str(root/"audit-0/evidence.json"),str(root/"audit-1/evidence.json")],capture_output=True,text=True)
            self.assertEqual(completed.returncode,0,completed.stderr)

    def test_certified_negative_is_rechecked_without_inventing_a_solver_attempt(self):
        from make_it_stand.verification.cherchi_application import audit_saved_cherchi_case
        runner=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=locked_case(root/"inputs")
            path=root/"inputs/problem.json"
            manifest=json.loads(path.read_text());manifest["target"]["origin"]=[.25,-.125,.125]
            path.write_text(json.dumps(manifest,sort_keys=True)+"\n")
            case["cherchi_problem"]["sha256"]=hashlib.sha256(path.read_bytes()).hexdigest()
            result=runner.run_case(case,1,2,bundle_root=root/"runs",asset_root=root/"inputs",
                cherchi_build=root/"must-not-build",cherchi_probe=root/"must-not-execute")
            self.assertEqual(result["runs"][0]["status"],"INFEASIBLE_CERTIFIED")
            self.assertEqual(result["attempts"],[])
            self.assertEqual(result["summary"]["expected_runs"],1)
            self.assertEqual(result["summary"]["failed_runs"],1)
            report=audit_saved_cherchi_case(result,bundle_root=root/"runs")
            self.assertEqual(report["qualified_material_runs"],0)
            self.assertEqual(report["total_function_evaluations"],0)
            self.assertEqual(report["runs"][0]["status"],"INFEASIBLE_CERTIFIED")
            for change in ("sha","missing","path"):
                altered=copy.deepcopy(result)
                assets=altered["initialization_plan"]["problem_identity"]["assets"]
                if change=="sha":assets[0]["sha256"]="0"*64
                elif change=="missing":assets.pop()
                else:assets[0]["path"]="../not-the-registered-asset.obj"
                altered["scientific_hash"]=runner.scientific_hash({k:v for k,v in altered.items() if k!="scientific_hash"})
                with self.subTest(change=change),self.assertRaises(ValueError):
                    audit_saved_cherchi_case(altered,bundle_root=root/"runs")

    def test_missing_locked_material_cannot_fall_back_to_legacy(self):
        runner=load_runner()
        for mode in ("SPINNING_TOP","YO_YO"):
            case={"id":"missing","mode":mode}
            planned=runner.plan_payload([case],1)["cases"][0]
            self.assertEqual(planned["status"],"BLOCKED_MISSING_CHERCHI_PROBLEM")
            self.assertEqual(planned["runs"],[])
            with self.assertRaisesRegex(ValueError,"cherchi_problem"):
                runner.run_case(case,1,1)

    def test_plan_binds_fixed_material_and_never_uses_preflight_as_start(self):
        runner=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=locked_case(root)
            plan=runner.plan_payload([case],2,asset_root=root,max_restarts=2)["cases"][0]
            self.assertEqual(plan["solver_family"],"CHERCHI_RAW_EQ14_V1")
            self.assertEqual([row["seed"] for row in plan["runs"]],[0,1])
            self.assertEqual(plan["problem"]["fixed_regions"][0]["density"],.75)
            self.assertEqual(plan["problem"]["design_regions"][0]["density"],2.)
            self.assertTrue(all(row["initialization_mode"]=="deterministic_mode_seed" for row in plan["runs"]))
            self.assertTrue(all("initial_quadric" not in row for row in plan["runs"]))
            for change in ("density","mode","missing_asset","answer"):
                with self.subTest(change=change):
                    bad=copy.deepcopy(case)
                    if change=="density":
                        mesh=root/"fixed_regions-0.obj"; data=mesh.read_bytes();mesh.write_bytes(data+b"\n")
                    elif change=="mode": bad["mode"]="SPINNING_TOP"
                    elif change=="missing_asset":bad["cherchi_problem"]["assets"].pop()
                    else:bad["initial_quadric"]=[1.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
                    try:
                        with self.assertRaises(ValueError):runner.plan_payload([bad],1,asset_root=root)
                    finally:
                        if change=="density":mesh.write_bytes(data)

    def test_actual_native_restarts_keep_failed_evaluations_and_saved_evidence(self):
        import _make_it_stand as native
        runner=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=locked_case(root/"inputs")
            result=runner.run_case(case,1,2,bundle_root=root/"runs",asset_root=root/"inputs",
                cherchi_build=Path(native.__file__).resolve().parent.parent,
                cherchi_probe=ROOT/"build/cherchi-reference-singlethread-20260917/cherchi_probe")
            self.assertEqual(result["solver_family"],"CHERCHI_RAW_EQ14_V1")
            self.assertFalse(result["passed"])
            self.assertEqual(len(result["attempts"]),2)
            self.assertEqual(result["summary"]["total_function_evaluations_all_attempts"],2)
            self.assertEqual(result["summary"]["resource_limited_runs"],1)
            for i,row in enumerate(result["attempts"]):
                self.assertEqual(row["native_restart"],i)
                self.assertEqual(row["function_evaluations"],1)
                self.assertEqual(row["solver_events"][0]["restart"],i)
                self.assertEqual(row["equation_profile"],"DUAL_MONOMIAL_RAW_V1")
                self.assertEqual(row["solver_events"][-1]["type"],"attempt_terminal")
            runner.validate_attempt_ledger(result["attempts"],expected_runs=1,maximum_restarts=2)
            forged=copy.deepcopy(result["attempts"])
            forged[1]["solver_events"][0]["restart"]=0
            with self.assertRaises(ValueError):runner.validate_attempt_ledger(forged,expected_runs=1,maximum_restarts=2)
            forged=copy.deepcopy(result["attempts"])
            forged[1]["solver_events"]=[e for e in forged[1]["solver_events"] if e["type"]!="evaluation"]
            with self.assertRaises(ValueError):runner.validate_attempt_ledger(forged,expected_runs=1,maximum_restarts=2)
            output=root/"runs/locked-method/run-000000"
            self.assertTrue((output/"evidence.json").is_file())
            inputs=json.loads((output/"inputs.json").read_text())
            self.assertEqual(inputs["initialization_mode"],"deterministic_mode_seed")
            self.assertIsNone(inputs["initial_coefficients"])
            self.assertEqual(inputs["maximum_radius_doublings"],8)
            self.assertEqual(json.loads((output/"freshness.json").read_text())["status"],"CURRENT")
            from make_it_stand.verification.cherchi_application import audit_saved_cherchi_case
            audit=audit_saved_cherchi_case(result,bundle_root=root/"runs")
            self.assertEqual(audit["total_function_evaluations"],2)
            self.assertEqual(audit["qualified_material_runs"],0)
            self.assertEqual(len(audit["runs"]),1)
            for field in ("passed","function_evaluations","seed","relative_output","asset_sha"):
                corrupt=copy.deepcopy(result)
                if field=="passed":corrupt["passed"]=True
                elif field=="function_evaluations":corrupt["attempts"][0][field]=0
                elif field=="seed":corrupt["runs"][0][field]=42
                elif field=="asset_sha":corrupt["initialization_plan"]["problem_identity"]["assets"][0]["sha256"]="0"*64
                else:corrupt["runs"][0][field]="../not-this-run"
                # Rehashing a claimed success must not bypass original run evidence.
                corrupt["scientific_hash"]=runner.scientific_hash({k:v for k,v in corrupt.items() if k!="scientific_hash"})
                with self.subTest(field=field), self.assertRaises(ValueError):
                    audit_saved_cherchi_case(corrupt,bundle_root=root/"runs")
            application=root/"application.json"
            application.write_text(json.dumps({"cases":[result]})+"\n")
            completed=subprocess.run([sys.executable,str(ROOT/"experiments/spin_it_reproduction/faster/check_saved_cherchi_application.py"),
                "--application",str(application),"--bundle-root",str(root/"runs"),"--output",str(root/"audit")],
                capture_output=True,text=True)
            self.assertEqual(completed.returncode,0,completed.stderr)
            report=json.loads((root/"audit/scientific.json").read_text())
            self.assertEqual(report["cases"][0]["qualified_material_runs"],0)
            self.assertTrue((root/"audit/evidence.json").is_file())
            strict=strict_audit()
            first=str(root/"audit/evidence.json")
            with self.assertRaisesRegex(ValueError,"not two"):
                strict._supporting_checks({"faster.cherchi_saved_material":[first,first]})
            completed=subprocess.run([sys.executable,str(ROOT/"experiments/spin_it_reproduction/faster/check_saved_cherchi_application.py"),
                "--application",str(application),"--bundle-root",str(root/"runs"),"--output",str(root/"audit-repeat")],
                capture_output=True,text=True)
            self.assertEqual(completed.returncode,0,completed.stderr)
            with self.assertRaisesRegex(ValueError,"no independently qualified"):
                strict._supporting_checks({"faster.cherchi_saved_material":[first,str(root/"audit-repeat/evidence.json")]})

if __name__=="__main__":unittest.main()
