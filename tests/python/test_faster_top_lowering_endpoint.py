"""Pure-height endpoint: algebraic Newton validation, not a geometry case."""
from dataclasses import replace
import math
import hashlib
import unittest
import numpy as np
from make_it_stand.spin import SpinItFasterRequest, SpinSolverConfig, solve_spin_it_faster_with_external_moment_evaluator
from test_spin_it_faster import cube


class TopLoweringEndpointTests(unittest.TestCase):
    def test_directional_native_bridge_runs_real_updates_and_preserves_failed_trial(self):
        # Nonphysical analytic map: isolates bridge direction, state and ledger
        # transmission. Missing/padded directional data must not pass this test.
        u=np.eye(10)[:,[0,1,2,3,7,8]]
        start=np.array([.125,.3,.7,-1,0,0,0,.2,.4,0.])
        a=np.diag([.1,.2,.2,.3,0,0,0,.1,.1,0])
        offset=np.array([2,0,0,.8,3,4,5,0,0,0.])
        calls=[]
        def evaluate(p):
            p=np.asarray(p);calls.append(p.copy())
            if len(calls)==2:
                raise ArithmeticError("directional trial rejected")
            v=u.T@p;v[0]+=1 if v[0]>=0 else -1;v/=np.linalg.norm(v)
            b=u@(np.eye(6)-2*np.outer(v,v))[:,1:]
            return dict(raw_moments=offset+a@p,native_tangent_basis=b,
                native_directional_jacobian=a@b,derivative_defined=True,
                geometry_hash=hashlib.sha256(p.tobytes()).hexdigest(),topology_hash="analytic",
                evaluation_map_identity="affine-test",boolean_backend_identity="NONE",
                derivative_identity="analytic-directional")
        request=replace(self.request(spinning_top_objective="LOWER_CENTER_OF_MASS"),
            settings=SpinSolverConfig(equation_profile="DUAL_MONOMIAL_RAW_V1",maximum_restarts=1))
        result=solve_spin_it_faster_with_external_moment_evaluator(request,evaluate,
            initial_quadrics=[start],support_basis=u)
        self.assertEqual(result["status"],"CONVERGED")
        self.assertLessEqual(result["raw_equation_residual"],1e-12)
        self.assertGreater(len(result["trace"]),2)
        attachment=result["supported_space"]
        self.assertEqual(attachment["version"],1)
        np.testing.assert_array_equal(attachment["support_basis"],u)
        records=attachment["evaluations"]
        self.assertEqual(len(records),len(calls))
        self.assertEqual(len(calls),result["function_evaluations"])
        self.assertIn("ArithmeticError",records[1]["error"])
        self.assertIn("directional trial rejected",records[1]["error"])
        failure=result["supported_space"]["python_callback_diagnostics"][0]
        self.assertEqual(failure["evaluation_index"],2)
        self.assertEqual(failure["exception_type"],"builtins.ArithmeticError")
        self.assertEqual(failure["message"],"directional trial rejected")
        self.assertIn("test_directional_native_bridge",failure["traceback"])
        self.assertIsNone(records[1]["evaluation"])
        self.assertIsNone(records[1]["linearization"])
        for record,p in zip(records,calls):
            np.testing.assert_array_equal(record["quadric_coefficients"],p)
            if record["evaluation"] is not None:
                self.assertEqual(np.shape(record["evaluation"]["native_directional_jacobian"]),(10,5))
                np.testing.assert_allclose(record["evaluation"]["raw_moments"],offset+a@p,rtol=0,atol=1e-14)
                self.assertEqual(np.shape(record["linearization"]["full_residual"]),(14,))
        self.assertFalse(result["final_geometry_verified"])
        original_hash=attachment["scientific_hash"]
        calls.clear()
        again=solve_spin_it_faster_with_external_moment_evaluator(request,evaluate,
            initial_quadrics=[start],support_basis=u)
        self.assertEqual(again["supported_space"]["scientific_hash"],original_hash)

    def test_directional_bridge_rejects_padded_derivative_and_invalid_basis(self):
        from make_it_stand.spin import _native_request
        native,request=_native_request(self.request(spinning_top_objective="LOWER_CENTER_OF_MASS"))
        u=np.eye(10)[:,[0,1,2,3,7,8]]
        start=np.array([[.125,.3,.7,-1,0,0,0,.2,.4,0.]])
        def bad(_):
            return dict(raw_moments=np.ones(10),moment_jacobian=np.zeros((10,10)),
                derivative_defined=True,geometry_hash="analytic",topology_hash="analytic",
                evaluation_map_identity="bad-map",boolean_backend_identity="NONE",derivative_identity="bad")
        result=native.solve_spin_it_faster_external_moment_evaluator(request,bad,start,None,u)
        self.assertEqual(result["status"],"GEOMETRY_EVALUATION_FAILED")
        self.assertEqual(result["function_evaluations"],1)
        self.assertIn("ambient",result["supported_space"]["evaluations"][0]["error"])
        for malformed in (u.T,np.zeros((10,6)),np.full((10,6),np.nan)):
            with self.subTest(shape=malformed.shape),self.assertRaises(ValueError):
                native.solve_spin_it_faster_external_moment_evaluator(request,bad,start,None,malformed)

    def test_mapping_extraction_exception_keeps_diagnostics_without_caller_identity(self):
        class BrokenMapping(dict):
            def __getitem__(self,key):
                if key=="raw_moments":raise ArithmeticError("mapping extraction failed")
                return super().__getitem__(key)
        def evaluate(_):
            return BrokenMapping(raw_moments=[],geometry_hash="analytic",topology_hash="analytic",
                evaluation_map_identity="bad-map",boolean_backend_identity="NONE",derivative_identity="bad")
        req=self.request(spinning_top_objective="LOWER_CENTER_OF_MASS")
        options=dict(initial_quadrics=[[.125,.3,.7,-1,0,0,0,.2,.4,0]],support_basis=np.eye(10)[:,[0,1,2,3,7,8]])
        first=solve_spin_it_faster_with_external_moment_evaluator(req,evaluate,**options)
        second=solve_spin_it_faster_with_external_moment_evaluator(req,evaluate,**options)
        self.assertEqual(first["supported_space"]["scientific_hash"],second["supported_space"]["scientific_hash"])
        self.assertEqual(first["status"],"GEOMETRY_EVALUATION_FAILED")
        diagnostic=first["supported_space"]["python_callback_diagnostics"][0]
        self.assertEqual(diagnostic["exception_type"],"builtins.ArithmeticError")
        self.assertIn("__getitem__",diagnostic["traceback"])

    def test_supported_entry_requires_explicit_starts_without_calling_geometry(self):
        from make_it_stand.spin import _native_request
        native,request=_native_request(self.request(spinning_top_objective="LOWER_CENTER_OF_MASS"))
        def forbidden(_):raise AssertionError("geometry must not run")
        with self.assertRaisesRegex(ValueError,"explicit.*initial"):
            native.solve_spin_it_faster_external_moment_evaluator(request,forbidden,
                np.empty((0,10)),None,np.eye(10)[:,[0,1,2,3,7,8]])

    def request(self, **changes):
        vertices,faces=cube(-1.,1.)
        return SpinItFasterRequest(outer_vertices=vertices,outer_faces=faces,
            design_vertices=vertices,design_faces=faces,mode="SPINNING_TOP",
            settings=SpinSolverConfig(equation_profile="DUAL_MONOMIAL_RAW_V1",
                maximum_restarts=1,maximum_function_evaluations_per_restart=1),**changes)

    def test_explicit_endpoint_rejects_wrong_mode_and_extra_weight(self):
        request=self.request(spinning_top_objective="LOWER_CENTER_OF_MASS")
        for changes in ({"mode":"YO_YO"},{"center_of_mass_weight":1.},
                        {"spinning_top_objective":"UNKNOWN"}):
            with self.subTest(changes=changes),self.assertRaises(ValueError):replace(request,**changes)

    def test_actual_raw_eq14_uses_pure_height_with_four_top_constraints(self):
        # This is explicitly KNOWN_SOLUTION_VALIDATION of algebraic dispatch.
        # r has mass8, c3=1/8, and zero cross-axis moments. grad(c3)=(-1/64,0,0,1/8,...).
        r=np.array([8.,0.,0.,1.,3.,3.,1.,0.,0.,0.])
        q=np.array([1.,0.,0.,-8.,0.,0.,0.,0.,0.,0.])/math.sqrt(65.)
        def evaluator(p):
            raw=r+.01*(np.asarray(p)-q)
            identity=hashlib.sha256(np.asarray(p,dtype="<f8").tobytes()+raw.astype("<f8").tobytes()).hexdigest()
            return {"raw_moments":raw,"moment_jacobian":.01*np.eye(10),"derivative_defined":True,
                    "geometry_hash":"ALGEBRAIC_NO_MESH_"+identity,"topology_hash":"ALGEBRAIC_NO_TOPOLOGY",
                    "evaluation_map_identity":"ALGEBRAIC_LINEAR_MAP_"+identity,
                    "boolean_backend_identity":"NONE_ALGEBRAIC_TEST",
                    "derivative_identity":"ALGEBRAIC_0.01_IDENTITY"}
        pure=self.request(spinning_top_objective="LOWER_CENTER_OF_MASS")
        result=solve_spin_it_faster_with_external_moment_evaluator(pure,evaluator,initial_quadrics=[q])
        self.assertEqual(result["status"],"CONVERGED")
        self.assertLessEqual(result["raw_equation_residual"],1e-12)
        self.assertAlmostEqual(result["multiplier_mu"],math.sqrt(65.)/64.,places=13)
        self.assertEqual(len(result["multipliers"]),4)
        self.assertFalse(result["final_geometry_verified"])
        self.assertAlmostEqual(result["trace"][0]["objective"],.125,places=14)
        legacy=solve_spin_it_faster_with_external_moment_evaluator(self.request(),evaluator,initial_quadrics=[q])
        self.assertNotEqual(legacy["status"],"CONVERGED")

    def test_native_external_entry_rejects_invalid_endpoint_before_callback(self):
        from make_it_stand.spin import _native_request
        native,request=_native_request(self.request(spinning_top_objective="LOWER_CENTER_OF_MASS"))
        def should_not_evaluate(_):
            raise AssertionError("invalid endpoint reached evaluator")
        request.mode=native.SpinMode.YO_YO
        with self.assertRaises(ValueError):
            native.solve_spin_it_faster_external_moment_evaluator(request,should_not_evaluate,np.empty((0,10)),None)
        request.mode=native.SpinMode.SPINNING_TOP
        request.center_of_mass_weight=1.
        with self.assertRaises(ValueError):
            native.solve_spin_it_faster_external_moment_evaluator(request,should_not_evaluate,np.empty((0,10)),None)


if __name__=="__main__":unittest.main()
