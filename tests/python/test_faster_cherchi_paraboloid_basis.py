"""Reduced source API and actual reference-Boolean paraboloid derivatives."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import test_faster_cherchi_basis as full_basis_test
from make_it_stand import spin

ROOT, PROBE = full_basis_test.ROOT, full_basis_test.PROBE
FIXTURE = ROOT/"experiments/spin_it_reproduction/faster/fixtures/cherchi-paraboloid-basis-method-v1/problem.json"


class CherchiParaboloidBasisTests(unittest.TestCase):
    def test_pure_height_active_tangent_preserves_required_singular_structure(self):
        runner = full_basis_test.CherchiBasisTests().module()
        problem = json.loads((FIXTURE.parent.parent/"cherchi-pure-height-basis-method-v1/problem.json").read_text())
        protocol = runner._protocol(problem)
        basis = np.asarray(protocol['parameter_basis'])
        q = np.asarray(protocol['normalized_coefficients'])
        self.assertEqual(basis.shape, (10,5))
        np.testing.assert_allclose(basis.T@basis,np.eye(5),atol=1e-13)
        np.testing.assert_allclose(q@basis,0,atol=1e-13)
        np.testing.assert_array_equal(basis[[4,5,6,9]],0)
        projector=np.diag([1,1,1,1,0,0,0,1,1,0])-np.outer(q,q)
        np.testing.assert_allclose(basis@basis.T,projector,atol=1e-13)
        self.assertEqual(len(protocol['samples']),31)
        for sample in protocol['samples']:
            np.testing.assert_array_equal(np.asarray(sample['coefficients'])[[4,5,6,9]],0)
        bad=copy.deepcopy(problem);bad['coefficients'][4]=1e-8
        with self.assertRaisesRegex(ValueError,'active|support'):
            runner._protocol(bad)

    def test_pure_height_actual_material_derivative_has_five_not_ten_columns(self):
        runner=full_basis_test.CherchiBasisTests().module()
        fixture=FIXTURE.parent.parent/'cherchi-pure-height-basis-method-v1/problem.json'
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'output'
            record=runner.run(output,fixture=fixture,probe=PROBE)
            self.assertEqual(record['status'],'PASSED',json.dumps(record))
            self.assertEqual(record['evaluation_count'],31)
            self.assertFalse(record['optimizer_executed'])
            self.assertEqual(len(record['checks']),5)
            rows=[json.loads(x) for x in (output/'evaluations.jsonl').read_text().splitlines()]
            self.assertEqual(np.asarray(rows[0]['jacobian']).shape,(10,5))
            self.assertTrue(all(x['independent_moment_pass'] for x in rows))
            for defect in ('padded','wrong_column','missing'):
                bad=copy.deepcopy(rows)
                if defect=='padded':bad[0]['jacobian']=np.pad(bad[0]['jacobian'],((0,0),(0,5))).tolist()
                elif defect=='wrong_column':bad[0]['jacobian'][0][4]+=.5
                else:bad.pop()
                self.assertNotEqual(runner.capture_comparison(record['protocol'],bad)['status'],'PASSED')

    def test_lpi_line_endpoint_order_does_not_change_construction_identity(self):
        runner = full_basis_test.CherchiBasisTests().module()
        row = {"point_type": "LPI", "operands": [
            {"source_vertices": [{"input_object": obj, "input_vertex": index}]}
            for obj,index in ((1,164),(1,181),(0,4),(0,5),(0,6))],
            "surface_planes": [{"source_triangles": [{"input_object": 0,"input_triangle": 2}]}]}
        swapped = copy.deepcopy(row)
        swapped["operands"][:2] = row["operands"][1::-1]
        first = runner._vertex_identities({"records": [row]})
        self.assertEqual(first, runner._vertex_identities({"records": [swapped]}))
        wrong = copy.deepcopy(row)
        wrong["operands"][0], wrong["operands"][2] = wrong["operands"][2], wrong["operands"][0]
        self.assertNotEqual(first, runner._vertex_identities({"records": [wrong]}))

    def test_stationary_exact_plane_diagonal_flip_is_not_a_material_event(self):
        runner = full_basis_test.CherchiBasisTests().module()
        signature = getattr(runner, "stationary_boundary_signature", None)
        self.assertTrue(callable(signature), "exact stationary patch identity missing")
        # Rational square, with two triangulations and a genuinely bent variant.
        points = [["0","0","1/3"],["1","0","1/3"],["1","1","1/3"],["0","1","1/3"]]
        a, b = [[0,1,2],[0,2,3]], [[0,1,3],[1,2,3]]
        identities = [str(i) for i in range(4)]
        first = signature(identities, points, a, points, a)
        self.assertEqual(first, signature(identities, points, b, points, a))
        bent = copy.deepcopy(points); bent[3][2] = "1000000001/3000000000"
        self.assertNotEqual(first, signature(identities, bent, b, points, a))
        # Winding and duplicate boundary coverage must not be normalized away.
        self.assertNotEqual(first, signature(identities, points, [f[::-1] for f in b], points, a))
        with self.assertRaises(ValueError):
            signature(identities, points, a+a, points, a)

    def test_reduced_source_motion_retains_eight_directions_without_full_rank_claim(self):
        native = spin._native_module()
        evaluate = getattr(native, "evaluate_spin_it_faster_stratum_source_motion", None)
        self.assertTrue(callable(evaluate), "reduced native source differential API missing")
        problem = json.loads(FIXTURE.read_text())
        settings = native.SpinTriangulationSettings()
        settings.paraboloid_grid_resolution = 16
        settings.truncation_radius = 2.
        source = evaluate(problem["coefficients"], problem["design_vertices"], settings)
        self.assertTrue(source["defined"], source["diagnostic"])
        basis = np.asarray(source["coefficient_basis"])
        velocity = np.asarray(source["vertex_velocity"])
        self.assertEqual(basis.shape, (10, 8))
        self.assertEqual(velocity.shape, (3*578, 8))
        np.testing.assert_allclose(basis.T@basis, np.eye(8), atol=1e-14)
        q = np.asarray(problem["coefficients"])/np.linalg.norm(problem["coefficients"])
        np.testing.assert_allclose(q@basis, np.zeros(8), atol=1e-14)
        # Null z curvature must not become a fictitious ninth tangent direction.
        np.testing.assert_allclose(basis[6], np.zeros(8), atol=1e-14)
        self.assertTrue(source["side_selection_defined"])
        full = native.evaluate_spin_it_faster_source_motion(
            problem["coefficients"], problem["design_vertices"], settings)
        self.assertFalse(full["defined"])
        with self.assertRaises(ValueError):
            evaluate([float("nan")]*10, problem["design_vertices"], settings)

    def test_actual_paraboloid_boolean_checks_all_eight_directions(self):
        runner = full_basis_test.CherchiBasisTests().module()
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/"output"
            record = runner.run(output, fixture=FIXTURE, probe=PROBE)
            self.assertEqual(record["status"], "PASSED", json.dumps(record))
            self.assertFalse(record["optimizer_executed"])
            self.assertEqual(record["evaluation_count"], 49)
            self.assertEqual(len(record["checks"]), 8)
            self.assertEqual(np.asarray(record["parameter_basis"]).shape, (10, 8))
            evaluations = [json.loads(line) for line in (output/"evaluations.jsonl").read_text().splitlines()]
            self.assertEqual(np.asarray(evaluations[0]["jacobian"]).shape, (10, 8))
            self.assertIn("LPI", evaluations[0]["provenance_types"])
            self.assertTrue(all(row["independent_moment_pass"] for row in evaluations))
            for item in record["protocol"]["samples"]:
                q = np.asarray(item["coefficients"])
                matrix = np.array([[q[4],q[9],q[8]],[q[9],q[5],q[7]],[q[8],q[7],q[6]]])
                self.assertLess(abs(np.linalg.eigvalsh(matrix)).min(), 1e-14)
            for defect in ("padded", "bad_column", "branch", "missing"):
                bad = copy.deepcopy(evaluations)
                if defect == "padded":
                    bad[0]["jacobian"] = np.pad(bad[0]["jacobian"], ((0,0),(0,2))).tolist()
                elif defect == "bad_column": bad[0]["jacobian"][0][7] += .5
                elif defect == "branch": bad[1]["branch"]["truncation_radius_hex"] = "changed"
                else: bad.pop()
                self.assertNotEqual(runner.capture_comparison(record["protocol"], bad)["status"], "PASSED")


if __name__ == "__main__":
    unittest.main()
