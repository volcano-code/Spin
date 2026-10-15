"""Finite top weights must use saved actual Cherchi runs, not warm starts."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from test_faster_cherchi_table_dispatch import locked_case, ROOT
from test_spin_it_faster_paper_figures import load_runner


class CherchiParetoDispatchTests(unittest.TestCase):
    def test_finite_sweep_cannot_relabel_pure_endpoint_as_zero_weight_inertia(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=locked_case(root/"inputs","SPINNING_TOP")
            case["spinning_top_objective"]="LOWER_CENTER_OF_MASS"
            with self.assertRaisesRegex(ValueError,"finite.*objective|endpoint"):
                load_runner().run_cherchi_top_pareto_sweep(case,[0.],asset_root=root/"inputs")
            self.assertFalse((root/"runs").exists())

    def test_method_envelope_cannot_claim_default_41_weight_paper_scope(self):
        runner=load_runner()
        result={"weights":[0.,.5],"endpoint_status":{"pure_inertia":"REQUESTED",
                "pure_lowering":"BLOCKED_MISSING_OBJECTIVE"}}
        payload=runner.build_payload()
        runner.attach_top_pareto_method_result(payload,result,Path("/registered/bundles"))
        self.assertEqual(payload["figure12"]["weights"],[0.,.5])
        self.assertEqual(payload["figure12"]["status"],"PARTIAL_METHOD_ONLY")
        self.assertEqual(payload["figure12"]["checker_input"],
                         {"relative_path":"pareto-result.json","base":"--pareto-bundle-root"})
        self.assertEqual(payload["figure12"]["endpoint_status"],result["endpoint_status"])

    def test_independent_metrics_use_target_axis_origin_and_input_scale(self):
        from make_it_stand.verification.cherchi_application import _top_material_metrics
        facts={"physical":{"mass":2.,"center_of_mass":[3.,4.,5.],
                "central_inertia":[[2.,0.,0.],[0.,3.,.25],[0.,.25,4.]]},
               "normalization":{"L":2.}}
        result=_top_material_metrics(facts,{"axis":[2.,0.,0.],"origin":[1.,4.,5.]},.5)
        self.assertAlmostEqual(result["f_yo"],6.28125)
        self.assertEqual(result["center_height_input_normalized"],1.)
        self.assertAlmostEqual(result["objective"],6.78125)
        self.assertEqual(result["physical_mass"],2.)

    def test_display_output_directory_is_not_scientific_identity(self):
        runner=load_runner()
        result={"weights":[0.,.5],"endpoint_status":{"pure_inertia":"REQUESTED",
                "pure_lowering":"BLOCKED_MISSING_OBJECTIVE"}}
        outputs=[]
        for directory in ("/repeat-one/bundles","/repeat-two/bundles"):
            payload=runner.build_payload()
            runner.attach_top_pareto_method_result(payload,result,Path(directory))
            outputs.append(runner.finalize_payload(payload)["scientific_hash"])
        self.assertEqual(*outputs)

    def test_no_locked_geometry_is_not_a_legacy_fallback(self):
        with self.assertRaisesRegex(ValueError, "locked|cherchi_problem"):
            load_runner().run_pareto_sweep({"id":"missing", "mode":"SPINNING_TOP"}, [0., 1.])

    def test_invalid_weights_rejected_before_any_run(self):
        runner=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=locked_case(root/"inputs", "SPINNING_TOP")
            for weights in ([], [0.,0.], [0.,float("nan")], [0.,-1.], [True]):
                with self.subTest(weights=weights), self.assertRaises(ValueError):
                    runner.run_pareto_sweep(case,weights,asset_root=root/"inputs",bundle_root=root/"runs")
                self.assertFalse((root/"runs").exists())

    def test_real_native_geometry_failures_are_retained_and_rechecked(self):
        import _make_it_stand as native
        from make_it_stand.verification.cherchi_application import audit_saved_cherchi_application
        from test_spin_it_faster_table1_runner import load_runner as table_runner
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);case=locked_case(root/"inputs", "SPINNING_TOP")
            result=load_runner().run_pareto_sweep(case,[0.,.5],asset_root=root/"inputs",
                bundle_root=root/"runs",cherchi_build=Path(native.__file__).resolve().parent.parent,
                cherchi_probe=ROOT/"build/cherchi-reference-singlethread-20260917/cherchi_probe",
                initial_seeds=(0,1),max_restarts=1)
            self.assertEqual(result["solver_family"],"CHERCHI_RAW_EQ14_V1")
            self.assertFalse(result["warm_started"])
            self.assertEqual(result["endpoint_status"]["pure_lowering"],"BLOCKED_MISSING_GEOMETRY_ENDPOINT_VALIDATION")
            self.assertEqual(len(result["cases"]),2)
            self.assertEqual([r["seed_sequence"] for r in result["cases"]],[[0,1],[0,1]])
            self.assertEqual(result["cost_ledger"]["function_evaluations_all_solves"],4)
            self.assertTrue(all(not c["passed"] for c in result["cases"]))
            path=root/"sweep.json";path.write_text(json.dumps(result))
            report=audit_saved_cherchi_application(path,bundle_root=root/"runs")
            self.assertEqual(sum(c["expected_runs"] for c in report["cases"]),4)
            self.assertEqual(sum(c["qualified_material_runs"] for c in report["cases"]),0)
            for kind in ("weight","missing","cost","endpoint"):
                altered=copy.deepcopy(result)
                if kind=="weight":altered["weights"][1]=.75
                elif kind=="missing":altered["cases"].pop()
                elif kind=="cost":altered["cost_ledger"]["function_evaluations_all_solves"]-=1
                else:altered["endpoint_status"]["pure_lowering"]="PASSED"
                altered["scientific_hash"]=table_runner().scientific_hash({k:v for k,v in altered.items() if k!="scientific_hash"})
                path.write_text(json.dumps(altered))
                with self.subTest(kind=kind), self.assertRaises(ValueError):
                    audit_saved_cherchi_application(path,bundle_root=root/"runs")


if __name__=="__main__":unittest.main()
