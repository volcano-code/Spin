"""Same-call authority transport, not an optimizer or manufacturing test."""
import json
import unittest
from fractions import Fraction

import numpy as np

from make_it_stand.spin import (SpinEvaluatorConfig, SpinItFasterRequest, _native_request,
    evaluate_spin_it_faster_quadric, evaluate_spin_it_faster_fullrank_quadric)
from make_it_stand.verification.mass_moments import integrate_exact_rational_raw_moments
from test_spin2014_quadric_interface import box


class RationalSourceProduction(unittest.TestCase):
    def test_saved_candidate_consumes_same_call_authority_and_existing_bundle_io(self):
        # Geometry adapter test only: no fabricated converged optimization,
        # KKT, manufacturing or repeat claim. The candidate is an ACTUAL
        # production evaluation and uses the existing exact bundle mesh I/O.
        from dataclasses import replace
        from types import SimpleNamespace
        import os
        import tempfile
        from make_it_stand.verification.faster_material import (
            _exact_spin2014_saved_retained, _canonical_exact_boolean_boundary)
        from make_it_stand.verification.spin_bundle import _Writer, _Reader, SCHEMA
        outer, faces = box(-2., 2.)
        design, df = box(-1., 1.)
        request = SpinItFasterRequest(outer, faces, design, df, material_density=2.7)
        config = SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
            icosphere_subdivisions=1, rational_source_coordinates=True)
        q = np.array([-.25, 0., 0., 0., 1., 1., 1., 0., 0., 0.])
        q /= np.linalg.norm(q)
        value = evaluate_spin_it_faster_quadric(request, q, evaluator_settings=config,
                                                coefficients_are_normalized=True)
        endpoint = SimpleNamespace(quadric_coefficients=tuple(q),
            final_geometry_hash=value.geometry_hash, material_vertices=value.material_vertices,
            material_faces=value.material_faces, rational_authority=value.rational_authority)
        candidate = SimpleNamespace(measured_request=request, solve=endpoint)
        retained, binding = _exact_spin2014_saved_retained(candidate, config)
        self.assertEqual(binding['schema'], 'spin2014-saved-rational-authority-v1')
        self.assertEqual(binding['source_construction_identity'],
                         value.rational_authority['source_construction_identity'])
        self.assertTrue(binding['independent_exact_moments_identical'])
        native, _ = _native_request(request)
        complement = native.libigl_exact_boolean_construction_from_rationals(
            value.rational_authority['design_vertices'], df,
            retained.rational_vertices, retained.arrays()[1],
            'difference')
        self.assertTrue(complement['succeeded'])
        cavity = _canonical_exact_boolean_boundary(complement['vertices'], complement['faces'])
        with tempfile.TemporaryDirectory() as directory:
            fd = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
            try:
                writer = _Writer(fd)
                records = [writer.mesh(mesh) for mesh in (retained, cavity)]
                reader = _Reader(fd)
                restored = [reader.mesh(record, schema=SCHEMA) for record in records]
                self.assertEqual([mesh.identity for mesh in restored],
                                 [mesh.identity for mesh in (retained, cavity)])
                self.assertEqual(restored[0].rational_vertices, retained.rational_vertices)
            finally:
                os.close(fd)
        for field in ('source_construction_identity', 'variable_moments', 'material_vertices'):
            bad = value.rational_authority
            if field == 'source_construction_identity':
                bad[field] = 'stale-source'
            elif field == 'variable_moments':
                bad[field][0] = '1/1'
            else:
                # Finite doubles masquerading as authority are not equivalent
                # to the certified rational approximation, even if close.
                bad[field] = [[str(Fraction(float(Fraction(x)))) for x in row] for row in bad[field]]
            damaged = SimpleNamespace(**vars(endpoint))
            damaged.rational_authority = bad
            with self.assertRaises(ValueError):
                _exact_spin2014_saved_retained(
                    SimpleNamespace(measured_request=request, solve=damaged), config)
        with self.assertRaises(ValueError):
            _exact_spin2014_saved_retained(
                SimpleNamespace(measured_request=replace(request, material_density=3.),
                                solve=endpoint), config)

    def test_authority_target_transform_uses_the_actual_affine_inverse(self):
        # Break caught: treating the binary64 rotation as EXACTLY orthogonal
        # displaces rational authority and can create a false containment gap.
        from dataclasses import replace
        from make_it_stand.verification.faster_material import (
            _target_rationals_to_normalized_model, verify_rational_target_affine_binding)
        outer, faces = box(-2., 2.)
        design, df = box(-1., 1.)
        request = SpinItFasterRequest(outer, faces, design, df,
            target_axes=((.8, -.6, 0.), (.6, .8, 0.), (0., 0., 1.)),
            target_origin=(.137, -.219, .05))
        c, s = Fraction(.8), Fraction(.6)
        origin = tuple(Fraction(x) for x in request.target_origin)
        model = (Fraction(1, 7), Fraction(-2, 9), Fraction(1, 13))
        x, y, z = tuple(a-b for a, b in zip(model, origin))
        target = (c*x+s*y, -s*x+c*y, z)
        encode = lambda v: f'{v.numerator}/{v.denominator}'
        restored, binding = _target_rationals_to_normalized_model(
            (tuple(encode(v) for v in target),), request,
            exact_inverse=True, return_binding=True)
        self.assertEqual(restored, (tuple(encode(v) for v in model),))
        verify_rational_target_affine_binding(binding, request)
        bad = json.loads(json.dumps(binding))
        bad['target_to_model_matrix'][0][0] = encode(c)
        with self.assertRaises(ValueError):
            verify_rational_target_affine_binding(bad, request)
        with self.assertRaises(ValueError):
            verify_rational_target_affine_binding(binding,
                replace(request, target_origin=(.138, -.219, .05)))
        with self.assertRaises(ValueError):
            verify_rational_target_affine_binding(binding,
                replace(request, original_units_per_normalized_unit=2.))

    def test_public_profile_preserves_authority_without_changing_legacy_result_shape(self):
        # Break caught: opt-in data is silently lost by the public facade, or
        # adding its default changes historical result/setting identities.
        from dataclasses import asdict, replace
        outer, faces = box(-2., 2.)
        design, df = box(-1., 1.)
        request = SpinItFasterRequest(outer, faces, design, df, material_density=2.7)
        legacy_config = SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
                                             icosphere_subdivisions=1)
        config = replace(legacy_config, rational_source_coordinates=True)
        self.assertNotIn('rational_source_coordinates', legacy_config.identity_dict())
        self.assertTrue(config.identity_dict()['rational_source_coordinates'])
        for invalid in (1, 'true', None):
            with self.assertRaises(ValueError):
                replace(config, rational_source_coordinates=invalid)
        for mixed in ('algebraic_finite_map', 'paper_projective_hyperboloid', 'fixed_normal_paraboloid'):
            with self.assertRaises(ValueError):
                replace(config, **{mixed: True})
        q = (-.25, 0., 0., 0., 1., 1., 1., 0., 0., 0.)
        value = evaluate_spin_it_faster_quadric(request, q, evaluator_settings=config)
        full = evaluate_spin_it_faster_fullrank_quadric(request, q, evaluator_settings=config)
        self.assertEqual(value.rational_authority, full.rational_authority)
        # The five-point protocol must use the actual normalization performed
        # by the producer, not NumPy's potentially different reduction order.
        # norm(-1/4,0,0,0,1,1,1,0,0,0)=7/4 exactly.
        self.assertEqual(full.normalized_coefficients,
                         (-1/7,0.,0.,0.,4/7,4/7,4/7,0.,0.,0.))
        self.assertEqual(len(full.branch_identity),64)
        authority = value.rational_authority
        material = integrate_exact_rational_raw_moments(authority['material_vertices'],
                                                       [list(f) for f in value.material_faces])
        self.assertEqual(tuple(Fraction(x) for x in authority['variable_moments']),
                         tuple(Fraction(2.7)*x for x in material))
        payload = asdict(value)
        self.assertEqual(json.loads(payload['rational_authority_json']), authority)
        authority['source_vertices'][0][0] = '999/1'
        self.assertNotEqual(authority, value.rational_authority)
        legacy = evaluate_spin_it_faster_quadric(request, q, evaluator_settings=legacy_config)
        self.assertNotIn('rational_authority_json', asdict(legacy))
        # Serialization-only NOT_RUN fixture carrying an actual production
        # packet. This tests the existing complete-run parser without inventing
        # a converged optimizer, a KKT check, or upgrading a historical solve.
        from make_it_stand.spin import (Spin2014QuadricMaterialRun,
            Spin2014QuadricObjectiveConfig, _fullrank_solve_result_from_native)
        from make_it_stand.verification.spin_bundle import _rotational_source_run
        from make_it_stand.verification.faster_material import _json_value
        solver = _fullrank_solve_result_from_native(dict(
            status='NOT_RUN', diagnostic='serialization-only production packet',
            final_geometry_verified=False, verification_diagnostic='NOT_VERIFIED',
            input_hash='', output_hash='', function_evaluations=0, restart_count=0,
            normalized_equation_residual=0., physical_relative_residual=0.,
            independent_raw_residual=None, multiplier_mu=0., multipliers=[],
            quadric_coefficients=None, equation_profile='DUAL_MONOMIAL_RAW_V1',
            event_hash='', trace_hash='', trace=[], events=[],
            restart_function_evaluations=[], evaluation_timings=[],
            final_geometry_hash=value.geometry_hash, final_derivative_defined=False,
            final_derivative_diagnostic='NOT_VERIFIED', final_variable_raw_moments=None,
            final_total_raw_moments=None, independent_total_raw_moments=None,
            material_vertices=value.material_vertices, material_faces=value.material_faces,
            rational_authority_json=asdict(value)['rational_authority_json']))
        run = Spin2014QuadricMaterialRun(objective_profile='SPIN2014_OBJECTIVES_QUADRIC_V1',
            objective_identity='serialization-only',
            objective_settings=Spin2014QuadricObjectiveConfig(), measured_request=request,
            fixed_material_density=1., solve=solver, material_export_complete=False,
            material_export_diagnostic='NOT_RUN', cavity_vertices=None,
            cavity_faces=None, cavity_raw_moments=None)
        saved = _json_value(asdict(run))
        restored = _rotational_source_run(saved)
        self.assertEqual(_json_value(asdict(restored)), saved)
        self.assertEqual(restored.solve.rational_authority, value.rational_authority)
        self.assertFalse(restored.solve.final_geometry_verified)
        damaged = json.loads(json.dumps(saved))
        damaged['solve']['final_geometry_hash'] = 'wrong-saved-geometry'
        with self.assertRaises(RuntimeError):
            _rotational_source_run(damaged)

    def test_native_opt_in_and_same_operation_authority_transport(self):
        outer, faces = box(-2., 2.)
        design, design_faces = box(-1., 1.)
        request = SpinItFasterRequest(outer, faces, design, design_faces,
                                      material_density=2.7)
        native, nr = _native_request(request)
        settings = native.SpinEvaluatorSettings()
        tri = native.SpinTriangulationSettings()
        self.assertFalse(tri.rational_source_coordinates)
        tri.rational_source_coordinates = True
        tri.icosphere_subdivisions = 1
        settings.triangulation = tri
        settings.derivative_method = native.SpinDerivativeMethod.DISCRETE_FULL_RANK_QUADRIC
        q = np.array([-.25, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype=float)
        packet = native.evaluate_spin_it_faster_fullrank_quadric(nr, q, settings)
        self.assertTrue(packet['derivative_defined'], packet['derivative_diagnostic'])
        authority = packet['rational_authority']
        self.assertEqual(authority['schema'], 'rational-quadric-production-authority-v1')
        self.assertEqual(authority['coordinate_space'], 'TARGET_FRAME')
        self.assertFalse(authority['analytic_coordinates_exact'])
        self.assertEqual(authority['geometry_hash'], packet['geometry_hash'])
        self.assertTrue(authority['source_construction_identity'])
        self.assertEqual(json.loads(json.dumps(authority)), authority)
        self.assertEqual(len(authority['material_vertices']), len(packet['material_vertices']))
        self.assertEqual(len(authority['material_birth_faces']), len(packet['material_faces']))
        self.assertEqual(len(authority['source_coordinate_lower']), len(authority['source_vertices']))
        self.assertEqual(len(authority['source_coordinate_upper']), len(authority['source_vertices']))
        for point, lower, upper in zip(authority['source_vertices'],
                                      authority['source_coordinate_lower'],
                                      authority['source_coordinate_upper']):
            for x, lo, hi in zip(point, lower, upper):
                self.assertLessEqual(Fraction(lo), Fraction(x))
                self.assertLessEqual(Fraction(x), Fraction(hi))
                self.assertLessEqual(Fraction(hi)-Fraction(lo), Fraction(1, 2**384))
        material_faces = packet['material_faces'].tolist()
        raw = integrate_exact_rational_raw_moments(authority['material_vertices'], material_faces)
        self.assertEqual(tuple(Fraction(x) for x in authority['variable_moments']),
                         tuple(Fraction(request.material_density)*x for x in raw))
        regular = native.evaluate_spin_it_faster_quadric(nr, q, settings)
        self.assertEqual(regular['geometry_hash'], packet['geometry_hash'])
        self.assertEqual(regular['rational_authority'], authority)
        legacy = native.evaluate_spin_it_faster_quadric(nr, q)
        self.assertNotIn('rational_authority', legacy)
        # The complement consumes the saved authority, never the display OBJ.
        complement = native.libigl_exact_boolean_construction_from_rationals(
            authority['design_vertices'], np.asarray(design_faces, dtype=np.int32),
            authority['material_vertices'], packet['material_faces'], 'difference')
        self.assertTrue(complement['succeeded'])
        cavity_raw = integrate_exact_rational_raw_moments(complement['vertices'],
                                                         complement['faces'].tolist())
        design_raw = integrate_exact_rational_raw_moments(authority['design_vertices'],
                                                         [list(f) for f in design_faces])
        self.assertEqual(tuple(a+b for a, b in zip(raw, cavity_raw)), design_raw)
        whole = native.make_exact_material_region_from_rationals(authority['design_vertices'],
                                                                 design_faces)
        material = native.make_exact_material_region_from_rationals(authority['material_vertices'],
                                                                    packet['material_faces'])
        cavity = native.make_exact_material_region_from_rationals(complement['vertices'],
                                                                  complement['faces'])
        self.assertTrue(whole.facts['valid'])
        self.assertTrue(material.facts['valid'])
        self.assertTrue(cavity.facts['valid'])
        partition = native.check_exact_material_partition(whole, [material, cavity])
        self.assertTrue(partition['completed'], partition)
        self.assertTrue(partition['complete_partition'], partition)


if __name__ == '__main__':
    unittest.main()
