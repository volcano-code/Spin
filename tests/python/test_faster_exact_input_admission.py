"""Exact input admission must not confuse small area with zero area."""
import unittest
from dataclasses import replace
import numpy as np
import _make_it_stand as native
from make_it_stand.spin import (
    SpinItFasterRequest, SpinEvaluatorConfig, SpinSolverConfig,
    evaluate_spin_it_faster_fullrank_quadric,
    solve_spin_it_faster_fullrank_quadric,
)


def box(a, b, split=False):
    v = np.array([(a,a,a),(b,a,a),(b,b,a),(a,b,a),
                  (a,a,b),(b,a,b),(b,b,b),(a,b,b)], dtype=float)
    f = [(0,2,1),(0,3,2),(4,5,6),(4,6,7),(0,1,5),(0,5,4),
         (3,7,6),(3,6,2),(0,4,7),(0,7,3),(1,2,6),(1,6,5)]
    if split:
        # Exact dyadic subdivision of one edge. Same cube, not a repaired mesh.
        v = np.vstack([v, v[0] + 2.**-45*(v[1]-v[0])])
        f = [(0,2,8),(8,2,1),(0,8,5),(8,1,5)] + [x for i,x in enumerate(f) if i not in (0,4)]
    return v, np.array(f, dtype=np.int32)


class ExactInputAdmission(unittest.TestCase):
    def setUp(self):
        self.config = SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
            icosphere_subdivisions=1, exact_binary64_source_boundary=True,
            exact_binary64_material_boundary=True)
        self.q = (-.25, 0., 0., 0., 1., 1., 1., 0., 0., 0.)

    def request(self, split):
        outer, of = box(-2., 2., split)
        design, df = box(-1., 1., split)
        return SpinItFasterRequest(outer.tolist(), of.tolist(), design.tolist(), df.tolist())

    def evaluate(self, request, config=None):
        return evaluate_spin_it_faster_fullrank_quadric(request, self.q,
            evaluator_settings=config or self.config)

    def test_exact_input_policy_accepts_same_cube_with_tiny_nonzero_faces(self):
        # Break caught: using the area-tolerance gate before the exact policy.
        request = self.request(True)
        audit = native.audit_exact_triangle_boundary(np.asarray(request.outer_vertices),
            np.asarray(request.outer_faces, dtype=np.int32), 'OutwardPositive')
        self.assertTrue(audit['boundary_valid'])
        self.assertGreater(native.audit_triangle_boundary(np.asarray(request.outer_vertices),
            np.asarray(request.outer_faces, dtype=np.int32))['degenerate_face_count'], 0)
        expected = self.evaluate(self.request(False))
        actual = self.evaluate(request)
        self.assertTrue(actual.derivative_defined, actual.derivative_diagnostic)
        np.testing.assert_allclose(actual.total_raw_moments, expected.total_raw_moments,
                                   rtol=1e-12, atol=1e-14)

    def test_legacy_policy_still_rejects_tiny_faces(self):
        with self.assertRaises(ValueError):
            self.evaluate(self.request(True), replace(self.config, exact_binary64_source_boundary=False))

    def test_stationary_skinny_design_facet_preserves_moments_and_derivative(self):
        # Break caught: applying moving-facet conditioning to an exact fixed
        # facet. This dyadic subdivision changes neither the cube nor its motion.
        request = self.request(False)
        v = np.array(request.design_vertices)
        f = np.array(request.design_faces)
        a, b, c = f[0]
        point = .5*(v[a]+v[b]) + 2.**-45*(v[c]-v[a])
        skinny = replace(request, design_vertices=np.vstack([v,point]).tolist(),
            design_faces=np.vstack([[(a,b,8),(b,c,8),(c,a,8)], f[1:]]).tolist())
        expected = self.evaluate(request)
        actual = self.evaluate(skinny)
        self.assertTrue(actual.derivative_defined, actual.derivative_diagnostic)
        np.testing.assert_allclose(actual.total_raw_moments, expected.total_raw_moments,
                                   rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(actual.coefficient_jacobian, expected.coefficient_jacobian,
                                   rtol=1e-10, atol=1e-12)

    def test_solver_enters_real_evaluation_under_exact_input_policy(self):
        request = replace(self.request(True), mode='YO_YO', settings=SpinSolverConfig(
            equation_profile='DUAL_MONOMIAL_RAW_V1', maximum_restarts=1,
            maximum_function_evaluations_per_restart=1))
        result = solve_spin_it_faster_fullrank_quadric(request, evaluator_settings=self.config)
        self.assertEqual(result.function_evaluations, 1)
        self.assertTrue(result.events)
        # This checks entry and real evaluation, NOT convergence or a seal.

    def test_exact_policy_still_rejects_actual_zero_area(self):
        request = self.request(True)
        v = np.array(request.outer_vertices); v[8] = v[0]
        with self.assertRaises(ValueError):
            self.evaluate(replace(request, outer_vertices=v))

    def test_exact_policy_still_rejects_open_boundary(self):
        request = self.request(True)
        with self.assertRaises(ValueError):
            self.evaluate(replace(request, outer_faces=request.outer_faces[:-1]))

    def test_exact_policy_still_rejects_inward_boundary(self):
        request = self.request(True)
        with self.assertRaises(ValueError):
            self.evaluate(replace(request, outer_faces=[(c,b,a) for a,b,c in request.outer_faces]))

    def test_exact_policy_still_rejects_intersecting_closed_shells(self):
        request = self.request(False)
        v, f = box(-2., 2.)
        v2 = v + np.array([.3, .4, .5])
        with self.assertRaises(ValueError):
            self.evaluate(replace(request, outer_vertices=np.vstack([v,v2]).tolist(),
                outer_faces=np.vstack([f,f+8]).tolist()))


if __name__ == '__main__':
    unittest.main()
