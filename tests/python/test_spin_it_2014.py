import json
import tempfile
import unittest

import make_it_stand.spin2014 as spin2014_module
from make_it_stand import (
    SpinIt2014Request,
    solve_spin_it_2014,
    solve_spin_it_2014_cage,
)


class SpinIt2014PythonTests(unittest.TestCase):
    def test_audited_profile_exports_independently_checkable_kkt(self) -> None:
        self.assertIn("slqp_algorithm", SpinIt2014Request.__dataclass_fields__)
        common = dict(dimensions=(4,4,4), minimum=(-1.5,-3.,-1.5), maximum=(1.5,3.,1.5),
            material_mask=(1,)*64, carveable_mask=tuple(int(x in (1,2) and y in (1,2) and z in (1,2))
                for z in range(4) for y in range(4) for x in range(4)),
            initial_level=2, maximum_level=2, initial_fill=.2, laplacian_weight=.03,
            objective="TOP_INERTIA_ONLY", formulation="PAPER_GIVENS", initial_phi=1.7,
            slqp_algorithm="AUDITED_KKT_RESTORATION_V1")
        a=solve_spin_it_2014(SpinIt2014Request(**common)).diagnostics
        b=solve_spin_it_2014(SpinIt2014Request(**common)).diagnostics
        self.assertEqual(a,b)
        self.assertEqual(a["slqp_algorithm"],common["slqp_algorithm"])
        self.assertEqual(a["optimality_scope"],"FIRST_ORDER_KKT_ONLY")
        self.assertEqual(a["slqp_audit_controls"]["restoration_trust_bound"],1.)
        trace=a["adaptive_passes"][0]["slqp_trace"]
        self.assertTrue(trace[-1]["final_state"])
        for t in trace:
            k=t["kkt"]; x=t["state"]; g=t["gradient"]; j=t["constraint_jacobian"]
            self.assertTrue(k["valid"])
            residual=[g[i]+sum(row[i]*lam for row,lam in zip(j,k["lambda"]))-k["alpha"][i]+k["eta"][i] for i in range(len(x))]
            values={"stationarity":max(map(abs,residual),default=0),
                "equality_primal":max(map(abs,t["constraints"]),default=0),
                "bound_primal":max([0]+[-v for v in x[:-1]]+[v-1 for v in x[:-1]]),
                "dual_sign":max([0]+[-v for v in k["alpha"]+k["eta"]]),
                "complementarity":max([0]+[abs(k["alpha"][i]*v) for i,v in enumerate(x[:-1])]+[abs(k["eta"][i]*(1-v)) for i,v in enumerate(x[:-1])])}
            for key,value in values.items(): self.assertAlmostEqual(value,k[key],delta=1e-10)
            self.assertAlmostEqual(max(values.values()),k["full"],delta=1e-10)
            self.assertEqual(k["alpha"][-1],0); self.assertEqual(k["eta"][-1],0)
        self.assertLessEqual(trace[-1]["kkt"]["full"],1e-7)
        changed=solve_spin_it_2014(SpinIt2014Request(**(common|{"slqp_restoration_trust_bound":.5}))).diagnostics
        self.assertNotEqual(changed["input_hash"],a["input_hash"])
        for update in ({"slqp_algorithm":"unknown"},{"formulation":"LEGACY_ANGLE_ELIMINATED","initial_phi":0},
            {"slqp_restoration_trust_bound":float("nan")},{"slqp_lp_iteration_limit":0}):
            with self.assertRaises(ValueError): SpinIt2014Request(**(common|update))

    def test_paper_givens_native_joint_angle_and_rejections(self) -> None:
        self.assertIn("formulation", SpinIt2014Request.__dataclass_fields__,
                      "explicit native Givens path is unavailable")
        common = dict(dimensions=(2, 2, 2), minimum=(-1., -2., 0.),
                      maximum=(1., 2., 2.), material_mask=(1,)*8,
                      carveable_mask=(1,)*8, initial_level=1, maximum_level=1,
                      initial_fill=0.2, objective="TOP_INERTIA_ONLY",
                      formulation="PAPER_GIVENS", initial_phi=1.7,
                      slqp_maximum_iterations=120)
        first = solve_spin_it_2014(SpinIt2014Request(**common))
        second = solve_spin_it_2014(SpinIt2014Request(**common))
        d = first.diagnostics
        self.assertEqual(d["formulation"], "PAPER_GIVENS")
        self.assertEqual(d["initial_phi"], 1.7)
        self.assertTrue(abs(d["final_phi"]-1.7)>1e-6)
        trace = d["adaptive_passes"][0]["slqp_trace"]
        self.assertTrue(any(t["accepted_step"] > 0 for t in trace))
        self.assertEqual(trace[0]["phi"], 1.7)
        self.assertEqual(d["output_hash"], second.diagnostics["output_hash"])
        self.assertEqual(first.formulation, "PAPER_GIVENS")
        self.assertEqual(first.final_phi, d["final_phi"])
        with self.assertRaises(ValueError):
            solve_spin_it_2014_cage(first, SpinIt2014Request(**common))
        native = spin2014_module._native_module()
        with self.assertRaisesRegex(ValueError, "Givens"):
            native.solve_spin_it_2014_cage(first._native, common["dimensions"],
                common["minimum"], common["maximum"], [1]*8, [1]*8)
        for update in ({"initial_phi": float("nan")}, {"formulation":"unknown"},
                       {"formulation":"LEGACY_ANGLE_ELIMINATED"},
                       {"mode":"STATIC_BALANCE","objective":"STATIC_LOWERING"}):
            with self.assertRaises(ValueError):
                SpinIt2014Request(**(common | update))

    def test_objective_selector_must_match_mode(self) -> None:
        self.assertEqual(
            getattr(spin2014_module, "SPIN_2014_OBJECTIVES", None),
            frozenset({
                "YOYO",
                "TOP_INERTIA_ONLY",
                "TOP_LOWERING_ONLY",
                "TOP_MASS_AND_LOWERING",
                "STATIC_LOWERING",
            }),
        )
        common = dict(
            dimensions=(1, 1, 1),
            minimum=(0.0, 0.0, 0.0),
            maximum=(1.0, 1.0, 1.0),
            material_mask=(1,),
            carveable_mask=(1,),
        )
        with self.assertRaises(ValueError):
            SpinIt2014Request(
                **common,
                mode="STATIC_BALANCE",
                objective="TOP_INERTIA_ONLY",
            )
        with self.assertRaises(ValueError):
            SpinIt2014Request(
                **common,
                mode="YO_YO",
                objective="TOP_LOWERING_ONLY",
            )

    def test_objective_weights_must_be_finite_and_nonnegative(self) -> None:
        common = dict(
            dimensions=(1, 1, 1),
            minimum=(0.0, 0.0, 0.0),
            maximum=(1.0, 1.0, 1.0),
            material_mask=(1,),
            carveable_mask=(1,),
        )
        for field, value in (
            ("inertia_weight", -1.0),
            ("inertia_weight", float("nan")),
            ("center_of_mass_weight", float("inf")),
        ):
            with self.subTest(field=field, value=value):
                with self.assertRaises(ValueError):
                    SpinIt2014Request(**common, **{field: value})

    def test_cage_api_rejects_non_artifact_before_native_execution(self) -> None:
        request = SpinIt2014Request(
            dimensions=(1, 1, 1), minimum=(0.0, 0.0, 0.0), maximum=(1.0, 1.0, 1.0),
            material_mask=(1,), carveable_mask=(1,),
        )
        with self.assertRaises(TypeError):
            solve_spin_it_2014_cage(object(), request)  # type: ignore[arg-type]

    def test_binary_constraint_tolerance_must_be_positive(self) -> None:
        with self.assertRaises(ValueError):
            SpinIt2014Request(
                dimensions=(1, 1, 1), minimum=(0.0, 0.0, 0.0), maximum=(1.0, 1.0, 1.0),
                material_mask=(1,), carveable_mask=(1,),
                binary_constraint_tolerance=0.0,
            )

    def test_slqp_iteration_limit_must_be_positive(self) -> None:
        with self.assertRaises(ValueError):
            SpinIt2014Request(
                dimensions=(1, 1, 1), minimum=(0.0, 0.0, 0.0), maximum=(1.0, 1.0, 1.0),
                material_mask=(1,), carveable_mask=(1,),
                slqp_maximum_iterations=0,
            )

    def test_slqp_gradient_tolerance_must_be_positive(self) -> None:
        with self.assertRaises(ValueError):
            SpinIt2014Request(
                dimensions=(1, 1, 1), minimum=(0.0, 0.0, 0.0), maximum=(1.0, 1.0, 1.0),
                material_mask=(1,), carveable_mask=(1,),
                slqp_projected_gradient_tolerance=0.0,
            )

    def test_sealed_deterministic_static_artifact_and_viewer_bundle(self) -> None:
        dimensions = (4, 4, 4)
        material = [1] * 64
        carveable = [0] * 64
        for z in range(1, 3):
            for y in range(1, 3):
                for x in range(1, 3):
                    carveable[x + 4 * (y + 4 * z)] = 1
        request = SpinIt2014Request(
            dimensions=dimensions,
            minimum=(-1.0, -1.0, 0.0),
            maximum=(1.0, 1.0, 2.0),
            material_mask=tuple(material),
            carveable_mask=tuple(carveable),
            mode="STATIC_BALANCE",
            objective="STATIC_LOWERING",
            inertia_weight=2.5,
            initial_level=1,
            maximum_level=2,
            maximum_adaptive_passes=3,
        )
        first = solve_spin_it_2014(request)
        second = solve_spin_it_2014(request)
        self.assertTrue(first.succeeded)
        self.assertTrue(first.geometry_sealed)
        self.assertEqual(first.diagnostics["objective"], "STATIC_LOWERING")
        self.assertEqual(first.diagnostics["objective_resolution"], "EXPLICIT")
        self.assertEqual(first.diagnostics["inertia_weight"], 2.5)
        self.assertEqual(first.diagnostics["center_of_mass_weight"], 1.0)
        self.assertEqual(
            first.diagnostics["objective_semantics"],
            "ANGLE_ELIMINATED_TRANSVERSE_INERTIA",
        )
        self.assertFalse(first.diagnostics["original_phi_trajectory"])
        self.assertEqual(first.diagnostics["output_hash"], second.diagnostics["output_hash"])
        self.assertFalse(first.outer_vertices.flags.writeable)
        with tempfile.TemporaryDirectory() as directory:
            paths = first.save_viewer_bundle(directory)
            self.assertEqual(len(paths), 5)
            with open(paths[-1], encoding="utf-8") as stream:
                self.assertEqual(json.load(stream)["schema"], "spin-it-2014-artifact-v1")


if __name__ == "__main__":
    unittest.main()
