"""Actual static-plane application routing, not mock solver dispatch checks."""
import copy
import hashlib
import json
import importlib.util
import os
import subprocess
import sys
from pathlib import Path
import tempfile
import unittest

import numpy as np
from test_spin_it_faster_table1_runner import load_runner
from test_spin_it_faster import cube


def method_case(root, objective="LOWER_CENTER_OF_MASS", weight=0., **settings):
    entries = {}
    for name, low, high in (("outer", -2., 2.), ("design_domain", -1., 1.)):
        vertices, faces = cube(low, high)
        data = "".join("v "+" ".join(map(str, point))+"\n" for point in vertices)
        data += "".join("f "+" ".join(str(i+1) for i in face)+"\n" for face in faces)
        path = root/f"{name}.obj"
        path.write_text(data)
        entries[name] = {"path": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    return {"id": "static-plane-method", "mode": "STATIC_BALANCE", **entries,
            "asset_class": "SYNTHETIC_REGRESSION", "material_density": 1.7,
            "target_origin": [0.,0.,0.], "target_axes": np.eye(3).tolist(),
            "minimum_wall_thickness": 0., "static_objective": objective,
            "static_pareto_weight": weight, "static_plane_settings": settings}


class StaticTableDispatchTests(unittest.TestCase):
    def test_actual_iteration_limits_are_preserved_in_resource_statistics(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case = method_case(root, maximum_iterations=0)
            result = runner.run_case(case, 2, 2, asset_root=root, initial_seeds=(0, 1))
            self.assertEqual(len(result['attempts']), 4)
            self.assertEqual([row['status'] for row in result['attempts']],
                             ['ITERATION_LIMIT'] * 4,
                             [row.get('solver_result', {}).get('diagnostic')
                              for row in result['attempts']])
            self.assertFalse(result['passed'])
            summary = result['summary']
            self.assertEqual(summary['resource_limited_runs'], 2)
            self.assertEqual(summary['resource_limited_attempts'], 4)
            self.assertEqual(summary['runs_with_resource_limited_attempts'], 2)
            self.assertEqual(summary['failed_runs'], 2)
            self.assertEqual(summary['strict_terminal_numerics']['numerically_qualified_runs'], 0)

    def test_independent_normalization_accepts_roundoff_not_changed_input_scale(self):
        from make_it_stand.verification.static_material import audit_static_case_material
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); case = method_case(root, maximum_iterations=0)
            case["target_origin"] = [.7262434559120877,.5927419911039958,.7841737089762773]
            case["target_axes"] = [[-.04100874821883127,-.7168823741180277,-.6959870287932273],
                [.19245148875083423,-.6891967610366285,.6985486733608532],
                [-.9804492373637259,-.10529713317807243,.16622817660483225]]
            actual = runner.run_case(case,1,1,asset_root=root)
            self.assertTrue(audit_static_case_material(actual)["attempts"][0]["moments_match"])
            actual["static_problem"]["target_length"] *= 1.001
            actual["static_problem_sha256"] = runner.scientific_hash(actual["static_problem"])
            with self.assertRaisesRegex(ValueError, "normalization"):
                audit_static_case_material(actual)

    def test_real_failed_initial_geometry_is_counted_not_discarded(self):
        from dataclasses import replace
        from make_it_stand.spin import solve_spin_it_faster_static_plane_family
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); case = method_case(root)
            request = runner.build_request(case, asset_root=root)
            # Native API permits supplied moments; this deliberately empty
            # mass is an exception-path test, not a valid material fixture.
            request = replace(request, fixed_shell_moments=(0.,)*10)
            result = solve_spin_it_faster_static_plane_family(request,
                initial_coefficients=(-10.,0.,0.,-1.))
            self.assertEqual(result["status"], "GEOMETRY_EVALUATION_FAILED")
            self.assertEqual(len(result["trace"]), 1)
            event = result["trace"][0]
            self.assertEqual(event["evaluation"], 0)
            self.assertFalse(event["evaluation_succeeded"])
            self.assertFalse(event["accepted"])
            self.assertEqual(event["diagnostic"], result["diagnostic"])
            from make_it_stand.verification.static_material import audit_static_case_material
            application = runner.run_case(case, 1, 1, asset_root=root)
            row = application["attempts"][0]
            row.update(status=result["status"], solver_result=runner._planar_json_value(result),
                       function_evaluations=1)
            runner.validate_static_plane_ledger(row)
            checked = audit_static_case_material(application)["attempts"][0]
            self.assertFalse(checked["moments_match"])
            self.assertFalse(checked["terminal_numerics_qualified"])
            self.assertEqual(checked["geometry_status"], "UNAVAILABLE_INITIAL_EVALUATION_FAILED")

    def test_saved_checker_cli_uses_existing_support_registry_without_promotion(self):
        from test_spin_it_strict_audit import load_audit
        root = Path(__file__).resolve().parents[2]
        cli = root/"experiments/spin_it_reproduction/faster/check_saved_static_material.py"
        self.assertTrue(cli.is_file())
        audit = load_audit()
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory); paths = []
            for index in (1,2):
                application = folder/f"application-{index}.json"
                case = load_runner().run_case(method_case(folder), 1, 2, asset_root=folder)
                application.write_text(json.dumps({"cases": [case]}))
                output = folder/f"check-{index}"
                completed = subprocess.run([sys.executable, str(cli), "--application", str(application),
                    "--output", str(output)], env={**os.environ, "PYTHONPATH": str(root/"python")},
                    capture_output=True, text=True, timeout=60)
                self.assertEqual(completed.returncode, 0, completed.stdout+completed.stderr)
                paths.append(str(output/"evidence.json"))
            checked = audit._supporting_checks({"faster.static_saved_material": paths})
            self.assertEqual(checked["faster.static_saved_material"]["verification_effect"],
                             "SUPPORTING_ONLY_NO_CATALOGUE_REQUIREMENT_VERIFIED")
            with self.assertRaises(ValueError):
                audit._supporting_checks({"faster.static_saved_material": [paths[0]]*2})
            application.write_text("{}")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit._supporting_checks({"faster.static_saved_material": paths})

    def test_static_pareto_routes_all_registered_weights_and_independent_endpoints(self):
        from test_spin_it_faster_paper_figures import load_runner as load_figures
        figures = load_figures()
        with tempfile.TemporaryDirectory() as directory:
            result = figures.run_pareto_sweep(method_case(Path(directory)), [0., .5, 1.],
                                              asset_root=Path(directory))
            self.assertEqual(result["solver_family"], "STATIC_PLANE_FAMILY_V1")
            self.assertFalse(result["warm_started"])
            self.assertEqual([row["weight"] for row in result["rows"]], [0., .5, 1.])
            self.assertEqual([row["application"]["static_problem"]["static_objective"]
                              for row in result["rows"]],
                             ["LOWER_CENTER_OF_MASS", "PARETO_TRADEOFF", "MAXIMUM_MASS"])
            total = sum(a["function_evaluations"] for row in result["rows"]
                        for a in row["application"]["attempts"])
            self.assertEqual(result["cost_ledger"]["function_evaluations_all_solves"], total)
            self.assertTrue(all(row["application"]["seed_sequence"] == [0]
                                for row in result["rows"]))
            json.dumps(result, allow_nan=False)

    def test_terminal_plane_cannot_be_detached_from_accepted_state(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            for objective in ("MAXIMUM_MASS", "LOWER_CENTER_OF_MASS"):
                result = runner.run_case(method_case(Path(directory), objective), 1, 2,
                                         asset_root=Path(directory))
                row = copy.deepcopy(result["attempts"][-1])
                row["solver_result"]["plane_coefficients"][1] += .125
                self.assertFalse(runner.numerical_attempt_acceptance(row)["passed"])

    def test_saved_numerical_evidence_rejects_unregistered_or_broadcast_values(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            result = runner.run_case(method_case(Path(directory)), 1, 2, asset_root=Path(directory))
            for field,value in (("static_objective", "UNREGISTERED_OBJECTIVE"),
                                ("static_pareto_weight", 2.), ("material_density", True)):
                row = copy.deepcopy(result["attempts"][-1]); row[field] = value
                self.assertFalse(runner.numerical_attempt_acceptance(row)["passed"], field)
            row = copy.deepcopy(result["attempts"][-1])
            row["solver_result"]["trace"][-1]["normalized_constraints"] = 0.
            self.assertFalse(runner.numerical_attempt_acceptance(row)["passed"])

    def test_saved_geometry_is_integrated_independently_and_tampering_rejected(self):
        spec = importlib.util.find_spec("make_it_stand.verification.static_material")
        self.assertIsNotNone(spec, "independent saved static geometry checker is missing")
        from make_it_stand.verification.static_material import audit_saved_static_application
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = runner.run_case(method_case(root), 1, 2, asset_root=root)
            path = root/"application.json"
            path.write_text(json.dumps({"cases": [result]}))
            checked = audit_saved_static_application(path)
            self.assertTrue(checked["cases"][0]["attempts"][-1]["moments_match"])
            self.assertFalse(checked["material_partition_verified"])
            for field in ("total_raw_moments", "variable_material_vertices"):
                poisoned = copy.deepcopy(result)
                target = poisoned["attempts"][-1]["solver_result"][field]
                if field == "total_raw_moments": target[0] += 1.
                else: target[0][0] += .25
                # Hash agreement alone must not make the geometry valid.
                poisoned["scientific_hash"] = runner.scientific_hash({
                    k:v for k,v in poisoned.items() if k != "scientific_hash"})
                path.write_text(json.dumps({"cases": [poisoned]}))
                with self.assertRaises(ValueError): audit_saved_static_application(path)

    def test_plan_freezes_four_coefficient_starts_without_answer(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); case = method_case(root)
            plan = runner.plan_payload([case], 2, max_restarts=2, asset_root=root)
            item = plan["cases"][0]
            self.assertEqual(item.get("solver_family"), "STATIC_PLANE_FAMILY_V1")
            self.assertNotIn("initial_quadrics", item)
            self.assertEqual(len(item["runs"]), 2)
            for row in item["runs"]:
                self.assertEqual(len(row["initial_planes"]), 2)
                for initial in row["initial_planes"]:
                    self.assertEqual(len(initial), 4)
                    self.assertLess(initial[3], 0.)
            poisoned = copy.deepcopy(case)
            poisoned["static_plane_settings"]["initial_coefficients"] = [0.,0.,0.,-1.]
            with self.assertRaises(ValueError):
                runner.plan_payload([poisoned], 1, asset_root=root)

    def test_all_static_objectives_use_actual_plane_solver_and_keep_geometry(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for objective,weight in (("MAXIMUM_MASS",0.),("LOWER_CENTER_OF_MASS",0.),
                                     ("PARETO_TRADEOFF",0.),("PARETO_TRADEOFF",1.)):
                with self.subTest(objective=objective, weight=weight):
                    case = method_case(root, objective, weight)
                    plan = runner.plan_payload([case], 1, max_restarts=2, asset_root=root)
                    actual = runner.run_case(case, 1, 2, asset_root=root)
                    self.assertEqual(actual.get("solver_family"), "STATIC_PLANE_FAMILY_V1")
                    self.assertEqual(actual["initialization_plan"], plan["cases"][0]["runs"])
                    self.assertEqual(actual["static_problem"]["static_pareto_weight"], weight)
                    self.assertTrue(actual["passed"], str(actual["attempts"]))
                    row = actual["attempts"][-1]
                    self.assertNotIn("raw_equation_residual", row)
                    self.assertNotIn("initial_quadric", row)
                    self.assertGreater(row["solver_result"]["accepted_steps"], 0)
                    self.assertGreater(len(row["solver_result"]["variable_material_faces"]), 0)
                    self.assertEqual(row["function_evaluations"], len(row["solver_result"]["trace"]))
                    self.assertEqual(actual["formal_paper_admission"], "NOT_GRANTED_BY_METHOD_RUNNER")
                    json.dumps(actual, allow_nan=False)

    def test_zero_iteration_budget_keeps_failed_restarts_and_rejects_forged_pass(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            actual = runner.run_case(method_case(root, maximum_iterations=0), 2, 2, asset_root=root)
            self.assertFalse(actual["passed"])
            self.assertEqual(len(actual["attempts"]), 4)
            self.assertEqual(actual["summary"]["total_function_evaluations_all_attempts"], 4)
            row = copy.deepcopy(actual["attempts"][0]); row["status"] = "CONVERGED"
            row["numerical_acceptance"] = {"passed": True}
            self.assertFalse(runner.numerical_attempt_acceptance(row)["passed"])
            damaged = copy.deepcopy(actual["attempts"])
            damaged[0]["solver_result"]["trace"] = []
            with self.assertRaises(ValueError):
                runner.validate_attempt_ledger(damaged, expected_runs=2, maximum_restarts=2)


if __name__ == "__main__": unittest.main()
