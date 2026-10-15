"""Actual public dogleg call accounting, separate from material approval."""
import json
import math
import os
import subprocess
import sys
import unittest

from make_it_stand.spin import solve_spin_it_faster_planar_fidget


def request(**overrides):
    outer = ((-3., -3.), (3., -3.), (3., 3.), (-3., 3.))
    patch = ((math.pi/4-.1, -.1), (math.pi/4+.1, -.1),
             (math.pi/4+.1, .1), (math.pi/4-.1, .1))
    values = dict(design_region_loops=(outer, tuple(reversed(patch))),
                  fixed_material_loops=(patch,), fixed_areal_density=100.,
                  design_areal_density=1., initial_center=(.75, .2))
    values.update(overrides)
    return values


class PlanarEvaluationLedgerTests(unittest.TestCase):
    def check_ledger(self, result):
        ledger = result.get("evaluation_ledger")
        self.assertIsInstance(ledger, dict, "actual dogleg call ledger is missing")
        self.assertEqual(ledger["schema"], "planar-dogleg-evaluation-ledger-v1")
        self.assertEqual(ledger["scope"],
                         "FUNCTION_EVALUATION_LEDGER_NOT_MATERIAL_OR_PAPER_APPROVAL")
        events = ledger["solver_evaluations"]
        self.assertEqual(len(events), result["function_evaluations"])
        self.assertEqual([event["evaluation_index"] for event in events],
                         list(range(1, len(events)+1)))
        self.assertEqual(events[0]["center"], [.75, .2])
        self.assertEqual(events[0]["purpose"], "initial")
        self.assertEqual(sum(event["accepted_trial"] for event in events),
                         result["accepted_steps"])
        for event in events:
            self.assertEqual(event["outcome"], "evaluated")
            self.assertTrue(event["has_residual"])
            self.assertTrue(all(math.isfinite(x) for x in event["center"]))
            self.assertTrue(math.isfinite(event["normalized_residual"]))
        auxiliary = ledger["auxiliary_evaluations"]
        self.assertEqual([event["purpose"] for event in auxiliary],
                         ["binding_initial_assessment", "binding_final_assessment"])
        self.assertEqual(ledger["auxiliary_geometry_evaluations"], 2)
        self.assertEqual(ledger["solver_function_evaluations"], len(events))
        self.assertEqual(ledger["total_residual_evaluations"], len(events)+2)
        self.assertEqual(ledger["terminal_phase"], "complete")
        json.dumps(ledger, allow_nan=False)
        return events

    def test_nontrivial_solve_keeps_finite_difference_and_accepted_calls(self):
        result = solve_spin_it_faster_planar_fidget(**request())
        self.assertEqual(result["status"], "CONVERGED")
        self.assertEqual(result["accepted_steps"], 5)
        self.assertEqual(result["function_evaluations"], 26)
        self.assertLessEqual(result["normalized_residual"], 1e-10)
        self.assertAlmostEqual(float(result["center"][0]), 1., places=9)
        self.assertAlmostEqual(float(result["center"][1]), 0., places=9)
        events = self.check_ledger(result)
        finite_differences = [event for event in events
                              if event["purpose"].startswith("finite_difference_")]
        self.assertEqual(len(finite_differences), 20)
        self.assertEqual({event["difference_axis"] for event in finite_differences}, {0, 1})

    def test_budget_failure_keeps_initial_call_and_two_auxiliary_assessments(self):
        result = solve_spin_it_faster_planar_fidget(
            **request(maximum_function_evaluations=1))
        self.assertEqual(result["status"], "MAXIMUM_FUNCTION_EVALUATIONS")
        self.assertEqual(result["function_evaluations"], 1)
        self.assertEqual(result["accepted_steps"], 0)
        self.check_ledger(result)

    def test_failed_initial_assessment_retains_exception_and_no_fabricated_residual(self):
        with self.assertRaises(ValueError) as caught:
            solve_spin_it_faster_planar_fidget(**request(initial_center=(1e200, 1e200)))
        ledger = getattr(caught.exception, "evaluation_ledger", None)
        self.assertIsInstance(ledger, dict, "failed residual call lost its ledger")
        self.assertEqual(ledger["terminal_phase"], "binding_initial_assessment")
        self.assertEqual(ledger["solver_function_evaluations"], 0)
        self.assertEqual(ledger["auxiliary_geometry_evaluations"], 1)
        self.assertEqual(ledger["total_residual_evaluations"], 1)
        event = ledger["auxiliary_evaluations"][0]
        self.assertEqual(event["outcome"], "failed")
        self.assertFalse(event["has_residual"])
        self.assertIsNone(event["normalized_residual"])
        self.assertIsNone(event["normalized_residual_vector"])
        self.assertEqual(event["diagnostic"], str(caught.exception))
        self.assertEqual(event["center"], [1e200, 1e200])
        json.dumps(ledger, allow_nan=False)

    def test_post_solve_export_error_keeps_actual_solver_calls(self):
        # Exercise the native boundary after all genuine evaluations, bypassing
        # the wrapper's earlier density validation; no solver mocks are used.
        from make_it_stand.spin import _native_module
        config = request()
        with self.assertRaises(ValueError) as caught:
            _native_module().solve_planar_fidget_spinner(
                config["design_region_loops"], config["fixed_material_loops"],
                100., 1., (.75, .2), 100, 60, 1e-10, .25, 1.,
                (0., 2.), 1., [])
        ledger = getattr(caught.exception, "evaluation_ledger", None)
        self.assertIsInstance(ledger, dict, "export failure discarded completed solver calls")
        self.assertEqual(ledger["terminal_phase"], "boundary_and_extrusion")
        self.assertEqual(ledger["solver_function_evaluations"], 26)
        self.assertEqual(ledger["total_residual_evaluations"], 28)
        self.assertEqual(sum(row["accepted_trial"] for row in ledger["solver_evaluations"]), 5)

    def test_two_independent_processes_return_identical_scientific_ledgers(self):
        code = (
            "import json; from test_planar_evaluation_ledger import request; "
            "from make_it_stand.spin import solve_spin_it_faster_planar_fidget; "
            "r=solve_spin_it_faster_planar_fidget(**request()); "
            "print(json.dumps({'ledger':r.get('evaluation_ledger'),"
            "'status':r['status'],'center':r['center'].tolist(),"
            "'residual':r['normalized_residual']},sort_keys=True,allow_nan=False))"
        )
        environment = dict(os.environ)
        test_root = os.path.dirname(__file__)
        environment["PYTHONPATH"] = test_root + os.pathsep + environment.get("PYTHONPATH", "")
        outputs = [subprocess.run([sys.executable, "-c", code], env=environment,
                                  check=True, capture_output=True, text=True).stdout
                   for _ in range(2)]
        self.assertIsInstance(json.loads(outputs[0])["ledger"], dict)
        self.assertEqual(outputs[0], outputs[1])


if __name__ == "__main__":
    unittest.main()
