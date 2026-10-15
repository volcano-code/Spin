"""Actual native Faster geometry to independent complete material sealing."""
from dataclasses import FrozenInstanceError, replace
import importlib
import json
import os
from pathlib import Path
import tempfile
import unittest
from test_exact_boolean_coordinate_reader import CoordinateReaderTests
from unittest.mock import patch

import numpy as np

from test_exact_material_region import box, combine
from make_it_stand import spin
from make_it_stand.verification.spin_problem import SpinArtifactProblem
from make_it_stand.verification.mass_moments import TriangleMesh, compute_mass_properties
from make_it_stand.verification.spin_artifact import SealedSpinArtifact
from make_it_stand.verification.spin_common import FRAME, ORDER
from make_it_stand.verification import faster_material as faster_material_module
from make_it_stand.verification.spin_problem import MeshSnapshot, solver_candidate_snapshot
from make_it_stand.verification.rational_text import (
    fraction_from_rational_text as _fraction_from_rational_text,
    fraction_to_rational_text as _fraction_to_rational_text)


COEFFICIENTS = (-.01, .01, .005, -.003, 1.1, .8, 1.3, .12, -.08, .17)


def problem(density=1.7, fixed_density=1.7, scale=1., offset=(0., 0., 0.), axis=(0., 0., 1.), target_origin=None):
    def transformed(mesh):
        v, f = mesh
        return v * scale + np.array(offset), f
    outer = transformed(box((-2.,)*3, (2.,)*3))
    design = transformed(box((-1.,)*3, (1.,)*3))
    shell = combine(outer, (design[0], design[1][:, ::-1]))
    return SpinArtifactProblem(original=outer, fixed_regions=[('shell', shell, fixed_density)],
        design_regions=[('design', design, density)], length_unit='m', mass_unit='kg',
        target={'mode': 'YO_YO', 'axis': list(axis), 'origin': list(offset if target_origin is None else target_origin), 'support_points': []},
        minimum_outer_wall=.25*scale)


class DirectionalMaterialBridgeTests(unittest.TestCase):
    """Geometry/transport tests, not a Cherchi optimization success fixture."""

    def test_directional_science_digest_excludes_only_named_runtime_diagnostics(self):
        def recorded(seconds, traceback, mass=1.):
            return faster_material_module._result('RESOURCE_LIMIT_REACHED', {
                'directional_mapping': {'schema': 'faster-material-supported-directions-v1'},
                'native_run': {'evaluation_timings': [{'wall_seconds': seconds}],
                    'raw_moments': [mass], 'supported_space': {
                        'evaluations': [{'error': 'ArithmeticError: rejected'}],
                        'python_callback_diagnostics': [{'traceback': traceback}]}}})
        first, second = recorded(.1, 'caller A'), recorded(.2, 'caller B')
        self.assertEqual(first.scientific_digest, second.scientific_digest)
        self.assertNotEqual(first.identity, second.identity)
        self.assertEqual(second.identity['native_run']['evaluation_timings'][0]['wall_seconds'], .2)
        self.assertNotEqual(first.scientific_digest, recorded(.1, 'caller A', mass=2.).scientific_digest)

    def fixture(self):
        p = problem(density=2., fixed_density=1.)
        p = SpinArtifactProblem(original=p.original.arrays(),
            fixed_regions=[(r.region_id, r.mesh.arrays(), r.density) for r in p.fixed_regions],
            design_regions=[(r.region_id, r.mesh.arrays(), r.density) for r in p.design_regions],
            length_unit='m', mass_unit='kg', minimum_outer_wall=p.minimum_outer_wall,
            target={'mode': 'SPINNING_TOP', 'axis': [0., 0., 1.],
                    'origin': [0., 0., 0.], 'support_points': []})
        evaluator = spin.SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC')
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1',
            maximum_restarts=1, maximum_function_evaluations_per_restart=1)
        _, request, _ = faster_material_module._prepare(p, evaluator, solver)
        u = np.eye(10)[:, [0, 1, 2, 3, 7, 8]]
        q = np.array([.125, .3, .7, -1., 0., 0., 0., .2, .4, 0.])
        q /= np.linalg.norm(q)
        return p, request, evaluator, solver, u, q

    @staticmethod
    def geometry(q, u):
        # An actual cube of side 1/4 translated by q3/20. Literal polynomial
        # volume moments and their derivatives, not the production integrator.
        q = np.asarray(q); t = q[3] / 20.; m = 1. / 64.
        vertices, faces = box((-.125, -.125, -.125+t), (.125, .125, .125+t))
        moments = np.array([m, 0., 0., m*t, m/192., m/192., m*(1/192.+t*t), 0., 0., 0.])
        h = u.T @ q; h[0] += 1 if h[0] >= 0 else -1; h /= np.linalg.norm(h)
        b = u @ (np.eye(6)-2*np.outer(h, h))[:, 1:]
        derivative = np.zeros((10, 10)); derivative[3, 3] = m/20.; derivative[6, 3] = m*t/10.
        return dict(variable_vertices=vertices, variable_vertices_model=vertices.copy(),
            variable_faces=faces, variable_raw_moments=moments,
            native_tangent_basis=b, native_directional_jacobian=derivative @ b,
            quadric_coefficients=q.copy(), derivative_defined=True,
            geometry_hash='analytic-translation-interface-fixture', topology_hash='cube',
            evaluation_map_identity='literal-translated-cube-not-quadric-boolean',
            boolean_backend_identity='cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:interface-fixture-only',
            derivative_identity='analytic-cube-directional')

    def test_directional_mesh_map_preserves_real_columns_and_adds_shell_once(self):
        p, request, _, _, u, q = self.fixture()
        source = self.geometry(q, u)
        mapped = faster_material_module._external_geometry_moment_map(
            p, request, source, support_basis=u, expected_coefficients=q)
        # Normalized outer side1 minus design side1/2, shell density1/2.
        fixed = np.array([7/16., 0, 0, 0, 31/768., 31/768., 31/768., 0, 0, 0])
        np.testing.assert_allclose(mapped['raw_moments'], fixed+source['variable_raw_moments'], atol=1e-14, rtol=0)
        self.assertNotIn('moment_jacobian', mapped)
        np.testing.assert_array_equal(mapped['native_directional_jacobian'], source['native_directional_jacobian'])
        expected = mapped['native_tangent_basis'].copy()
        source['native_tangent_basis'][:] = 0
        np.testing.assert_array_equal(mapped['native_tangent_basis'], expected)

    def test_directional_map_rejects_padding_stale_coefficients_and_wrong_mass(self):
        p, request, _, _, u, q = self.fixture()
        cases = [('ambient', {'moment_jacobian': np.zeros((10, 10))}),
                 ('dimensions', {'native_directional_jacobian': np.zeros((10, 10))}),
                 ('tangent', {'native_tangent_basis': np.ones((10, 5))}),
                 ('coefficients', {'quadric_coefficients': -q}),
                 ('moments', {'variable_raw_moments': self.geometry(q, u)['variable_raw_moments']*2})]
        for reason, change in cases:
            with self.subTest(reason=reason), self.assertRaisesRegex(ValueError, reason):
                faster_material_module._external_geometry_moment_map(p, request,
                    {**self.geometry(q, u), **change}, support_basis=u, expected_coefficients=q)

    def test_complement_geometry_needs_no_fabricated_directional_jacobian(self):
        p, request, _, _, u, q = self.fixture()
        source = self.geometry(q, u)
        del source['native_tangent_basis'], source['native_directional_jacobian']
        source['derivative_defined'] = False
        result = faster_material_module._external_geometry_moment_map(
            p, request, source, support_basis=u, require_derivative=False)
        self.assertNotIn('moment_jacobian', result)
        self.assertNotIn('native_directional_jacobian', result)
        np.testing.assert_allclose(result['variable_raw_moments'], source['variable_raw_moments'], atol=1e-15)

    def test_supported_material_request_requires_explicit_starts(self):
        p, _, evaluator, solver, u, _ = self.fixture()
        def forbidden(*args):
            raise AssertionError('geometry must not run before explicit starts')
        with self.assertRaisesRegex(ValueError, 'explicit.*initial'):
            faster_material_module.solve_cherchi_geometry_and_seal_faster(p,
                geometry_evaluator=forbidden, complement_evaluator=forbidden,
                initial_quadrics=(), evaluator_settings=evaluator, solver_settings=solver,
                spinning_top_objective='LOWER_CENTER_OF_MASS', support_basis=u)

    def test_geometry_owned_request_reaches_native_directional_ledger_without_padding(self):
        p, _, evaluator, solver, u, q = self.fixture()
        frozen = u.copy()
        def geometry(coefficients, context):
            # Mutation of the caller-owned basis must not affect this run.
            u[:] = 0
            return self.geometry(coefficients, frozen)
        def forbidden(*args):
            raise AssertionError('budget-limited run must not construct terminal material')
        result = faster_material_module.solve_cherchi_geometry_and_seal_faster(p,
            geometry_evaluator=geometry, complement_evaluator=forbidden,
            initial_quadrics=[q], evaluator_settings=evaluator, solver_settings=solver,
            spinning_top_objective='LOWER_CENTER_OF_MASS', support_basis=u)
        self.assertNotEqual(result.status, 'CONVERGED')
        identity = result.identity
        np.testing.assert_array_equal(identity['directional_mapping']['support_basis'], frozen)
        run = identity['native_run']; records = run['supported_space']['evaluations']
        self.assertEqual(run['function_evaluations'], 1)
        self.assertEqual(len(records), 1)
        self.assertFalse(records[0]['error'])
        self.assertEqual(np.asarray(records[0]['evaluation']['native_directional_jacobian']).shape, (10, 5))
        self.assertAlmostEqual(records[0]['evaluation']['raw_moments'][0], 29/64., places=13)


class FasterSettingsTests(unittest.TestCase):
    def test_exact_boundary_can_be_valid_when_binary64_derivative_is_degenerate(self):
        # Terminal complement evidence contains genuine rational slivers whose
        # correctly rounded binary64 coordinates are collinear. Keep exact
        # boundary validity separate from display/derivative admissibility.
        from fractions import Fraction
        from make_it_stand.verification.cherchi_motion_differential import (
            CherchiDifferentialUnavailable, normalize_closed_component_orientation)

        rows = ((Fraction(0), Fraction(0), Fraction(1)),
                (Fraction(1, 2), Fraction(0), Fraction(1) + Fraction(1, 2**60)),
                (Fraction(1), Fraction(0), Fraction(1)),
                (Fraction(0), Fraction(1), Fraction(0)))
        rationals = tuple(tuple(_fraction_to_rational_text(x) for x in row) for row in rows)
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)), dtype=np.int32)
        exact = MeshSnapshot.from_rationals(rationals, faces)
        self.assertTrue(exact.exact_region().facts['valid'])
        self.assertTrue(exact.exact_region().facts['completed'])
        display, triangles = exact.arrays()
        np.testing.assert_array_equal(display[:3, 2], np.ones(3))
        with self.assertRaisesRegex(ValueError, 'zero-area'):
            MeshSnapshot((display, triangles))
        with self.assertRaisesRegex(CherchiDifferentialUnavailable, 'zero-area face'):
            normalize_closed_component_orientation(display, triangles)
        # Analytic tetrahedron integration is independent of the checker used
        # by the bridge. Density is applied once, after geometric integration.
        density = Fraction(5, 2)
        volume = Fraction(1, 2**60) / 6
        sums = [sum(row[i] for row in rows) for i in range(3)]
        first = [volume*s/4 for s in sums]
        second = [[volume*(sums[i]*sums[j] + sum(row[i]*row[j] for row in rows))/20
                   for j in range(3)] for i in range(3)]
        expected = [float(density*x) for x in (volume, *first,
                    second[0][0], second[1][1], second[2][2],
                    second[1][2], second[0][2], second[0][1])]
        actual = faster_material_module._raw_moments(exact, float(density))
        np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0)

    def test_huge_rational_decimal_roundtrip_bypasses_python_integer_digit_guard(self):
        numerator = '1' + '0' * 5000 + '7'
        denominator = '3' + '0' * 5000 + '1'
        token = numerator + '/' + denominator
        value = _fraction_from_rational_text(token)
        canonical = _fraction_to_rational_text(value)
        self.assertEqual(_fraction_from_rational_text(canonical), value)
        self.assertEqual(len(canonical.split('/')[0]), 5002)
        self.assertEqual(len(canonical.split('/')[1]), 5002)

    def test_mesh_snapshot_accepts_large_exact_rational_coordinates(self):
        numerator = '1' + '0' * 5000 + '7'
        denominator = '3' + '0' * 5000 + '1'
        token = numerator + '/' + denominator
        value = _fraction_from_rational_text(token)
        display = float(value)
        vertices = ((0, 0, 0), (display, 0, 0), (0, 1, 0), (0, 0, 1))
        faces = ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))
        exact = (('0/1', '0/1', '0/1'), (token, '0/1', '0/1'),
                 ('0/1', '1/1', '0/1'), ('0/1', '0/1', '1/1'))
        class ValidRegion:
            facts = {'completed': True, 'valid': True, 'error': ''}

        class ExactNative:
            def make_exact_material_region_from_rationals(self, rows, triangles):
                self.rows = rows
                return ValidRegion()

        with patch('make_it_stand.verification.spin_problem.native', return_value=ExactNative()):
            snapshot = MeshSnapshot.from_rationals(exact, faces)
        self.assertEqual(snapshot.rational_vertices[1][0], _fraction_to_rational_text(value))
        self.assertEqual(snapshot.arrays()[0][1, 0], display)

    def test_mesh_snapshot_reuses_the_exact_region_that_authenticated_its_bytes(self):
        spin_problem_module = importlib.import_module(
            'make_it_stand.verification.spin_problem')

        class AuthenticatedRegion:
            facts = {'completed': True, 'valid': True, 'error': '',
                     'region_digest': 'authenticated-once'}

        class CountingNative:
            def __init__(self):
                self.builds = 0

            def make_exact_material_region(self, vertices, faces):
                self.builds += 1
                return AuthenticatedRegion()

            def empty_exact_material_region(self):
                self.builds += 1
                return AuthenticatedRegion()

        fake_native = CountingNative()
        mesh = box((-1., -1., -1.), (1., 1., 1.))
        with patch.object(spin_problem_module, 'native', return_value=fake_native):
            snapshot = MeshSnapshot(mesh)
            first_region = snapshot.exact_region()
            second_region = snapshot.exact_region()
            copied_snapshot = MeshSnapshot(snapshot)
            independently_authenticated = MeshSnapshot(mesh)

        self.assertIs(first_region, second_region)
        self.assertIs(first_region, copied_snapshot.exact_region())
        self.assertEqual(snapshot, independently_authenticated)
        self.assertEqual(fake_native.builds, 2)

        with patch.object(spin_problem_module, 'native', return_value=fake_native):
            empty_snapshot = MeshSnapshot(None)
            empty_region = empty_snapshot.exact_region()
            empty_copy = MeshSnapshot(empty_snapshot)

        self.assertIs(empty_region, empty_copy.exact_region())
        self.assertEqual(fake_native.builds, 3)

    def test_exact_boolean_coordinate_sidecar_builds_authoritative_rational_snapshot(self):
        from fractions import Fraction

        vertices = np.asarray(((0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
                               (0., 0., float(Fraction(1, 3)))), dtype=np.float64)
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)), dtype=np.int32)
        with tempfile.TemporaryDirectory() as directory:
            sidecar = Path(directory) / 'non_dyadic.exact.tsv'
            header = ('vertex_index\tstatus\texact_x\texact_y\texact_z\t'
                      'memory_x\tmemory_y\tmemory_z\tobj_x\tobj_y\tobj_z\n')
            rows = []
            exact_rows = ((Fraction(0), Fraction(0), Fraction(0)),
                          (Fraction(1), Fraction(0), Fraction(0)),
                          (Fraction(0), Fraction(1), Fraction(0)),
                          (Fraction(0), Fraction(0), Fraction(1, 3)))
            for index, (exact_row, display_row) in enumerate(zip(exact_rows, vertices)):
                exact_tokens = [f'{value.numerator}/{value.denominator}' for value in exact_row]
                display_tokens = [format(float(value), '.17g') for value in display_row]
                rows.append('\t'.join([str(index), 'MATCHED', *exact_tokens,
                                       *display_tokens, *display_tokens]) + '\n')
            sidecar.write_text(header + ''.join(rows), encoding='ascii')

            snapshot = faster_material_module.mesh_snapshot_from_exact_boolean_sidecar(
                vertices, faces, sidecar,
            )

        self.assertEqual(snapshot.rational_vertices[3][2], '1/3')
        self.assertNotEqual(Fraction(snapshot.arrays()[0][3, 2]), Fraction(1, 3))
        self.assertEqual(snapshot.arrays()[1].tolist(), faces.tolist())

    def test_exact_boolean_sidecar_preserves_two_thirds_intersection_coordinate(self):
        from fractions import Fraction

        display = np.asarray(((0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
                              (0., 0., float(Fraction(2, 3)))), dtype=np.float64)
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)), dtype=np.int32)
        header = ('vertex_index\tstatus\texact_x\texact_y\texact_z\t'
                  'memory_x\tmemory_y\tmemory_z\tobj_x\tobj_y\tobj_z\n')
        rows = []
        exact = ((Fraction(0), Fraction(0), Fraction(0)),
                 (Fraction(1), Fraction(0), Fraction(0)),
                 (Fraction(0), Fraction(1), Fraction(0)),
                 (Fraction(0), Fraction(0), Fraction(2, 3)))
        for index, (rational, point) in enumerate(zip(exact, display)):
            exact_tokens = [f'{value.numerator}/{value.denominator}' for value in rational]
            floats = [format(float(value), '.17g') for value in point]
            rows.append('\t'.join([str(index), 'APPROXIMATION_DISAGREES',
                                    *exact_tokens, *floats, *floats]) + '\n')
        class ValidRegion:
            facts = {'completed': True, 'valid': True, 'error': ''}
        class ExactNative:
            def make_exact_material_region_from_rationals(self, vertices, triangles):
                return ValidRegion()
        with tempfile.TemporaryDirectory() as directory:
            sidecar = Path(directory) / 'two_thirds.exact.tsv'
            sidecar.write_text(header + ''.join(rows), encoding='ascii')
            with patch('make_it_stand.verification.spin_problem.native', return_value=ExactNative()):
                snapshot = faster_material_module.mesh_snapshot_from_exact_boolean_sidecar(
                    display, faces, sidecar)
        self.assertEqual(snapshot.rational_vertices[3][2], '2/3')
        self.assertEqual(snapshot.arrays()[0][3, 2], float(Fraction(2, 3)))

    def test_exact_boolean_coordinate_sidecar_rejects_unbound_or_reordered_rows(self):
        vertices = np.asarray(((0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
                               (0., 0., 1.)), dtype=np.float64)
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)), dtype=np.int32)
        header = ('vertex_index\tstatus\texact_x\texact_y\texact_z\t'
                  'memory_x\tmemory_y\tmemory_z\tobj_x\tobj_y\tobj_z\n')
        rows = []
        for index, row in enumerate(vertices):
            tokens = [format(float(value), '.17g') for value in row]
            rows.append('\t'.join([str(index), 'MATCHED',
                                   *(f'{int(value)}/1' for value in row),
                                   *tokens, *tokens]) + '\n')
        with tempfile.TemporaryDirectory() as directory:
            sidecar = Path(directory) / 'tampered.exact.tsv'
            unbound = list(rows)
            unbound[3] = unbound[3].replace('\t1/1\t', '\t2/1\t', 1)
            sidecar.write_text(header + ''.join(unbound), encoding='ascii')
            with self.assertRaisesRegex(ValueError, 'does not bind'):
                faster_material_module.mesh_snapshot_from_exact_boolean_sidecar(
                    vertices, faces, sidecar,
                )

            reordered = list(rows)
            reordered[0], reordered[1] = reordered[1], reordered[0]
            sidecar.write_text(header + ''.join(reordered), encoding='ascii')
            with self.assertRaisesRegex(ValueError, 'reordered'):
                faster_material_module.mesh_snapshot_from_exact_boolean_sidecar(
                    vertices, faces, sidecar,
                )

    def test_exact_boolean_sidecar_audits_bounded_float_approximation(self):
        from fractions import Fraction

        vertices = np.asarray(((0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
                               (0., 0., np.nextafter(float(Fraction(1, 3)), np.inf))), dtype=np.float64)
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)), dtype=np.int32)
        header = ('vertex_index\tstatus\texact_x\texact_y\texact_z\t'
                  'memory_x\tmemory_y\tmemory_z\tobj_x\tobj_y\tobj_z\n')
        rows = []
        for index, row in enumerate(vertices):
            tokens = [format(float(value), '.17g') for value in row]
            exact = [f'{int(value)}/1' for value in row]
            if index == 3:
                exact[2] = '1/3'
            status = 'APPROXIMATION_DISAGREES' if index == 3 else 'MATCHED'
            rows.append('\t'.join([str(index), status, *exact, *tokens, *tokens]) + '\n')
        with tempfile.TemporaryDirectory() as directory:
            sidecar = Path(directory) / 'bounded_approximation.exact.tsv'
            sidecar.write_text(header + ''.join(rows), encoding='ascii')
            result = faster_material_module.mesh_snapshot_from_exact_boolean_sidecar(
                vertices, faces, sidecar, return_audit=True,
            )

        snapshot, audit = result
        self.assertEqual(snapshot.rational_vertices[3][2], '1/3')
        self.assertEqual(audit['status'], 'BOUNDED_FLOAT_DERIVATION')
        self.assertGreater(audit['maximum_absolute_error'], 0.)
        self.assertLessEqual(audit['maximum_absolute_error'], audit['absolute_error_bound'])
        faster_material_module.canonical(audit)

        rows[3] = rows[3].replace('1/3', '1003/3000', 1)
        with tempfile.TemporaryDirectory() as directory:
            sidecar = Path(directory) / 'out_of_bound.exact.tsv'
            sidecar.write_text(header + ''.join(rows), encoding='ascii')
            with self.assertRaisesRegex(ValueError, 'bounded derivation error'):
                faster_material_module.mesh_snapshot_from_exact_boolean_sidecar(
                    vertices, faces, sidecar, return_audit=True,
                )

    def test_nonfinite_native_diagnostic_values_remain_hashable_and_explicit(self):
        converted = faster_material_module._json_value({
            'nan': float('nan'),
            'positive_infinity': float('inf'),
            'negative_infinity': float('-inf'),
            'finite': 1.25,
        })

        self.assertEqual(converted['nan'], {'$nonfinite_float': 'nan'})
        self.assertEqual(converted['positive_infinity'], {'$nonfinite_float': '+inf'})
        self.assertEqual(converted['negative_infinity'], {'$nonfinite_float': '-inf'})
        self.assertEqual(converted['finite'], 1.25)

    def test_candidate_snapshot_keeps_full_nef_check_for_nested_material_shells(self):
        outer = box((-2., -2., -2.), (2., 2., 2.))
        inner = box((-1., -1., -1.), (1., 1., 1.))
        nested = combine(outer, (inner[0], inner[1][:, ::-1]))

        snapshot, scope = solver_candidate_snapshot(nested)

        self.assertIsInstance(snapshot, MeshSnapshot)
        self.assertEqual(scope, 'EXACT_NEF_SIGNED_REGION')
        facts = snapshot.exact_region().facts
        self.assertTrue(facts['valid'], facts['error'])
        self.assertEqual(facts['bounded_void_components'], 1)

    def test_external_cherchi_seal_preserves_model_frame_source_vertices_under_rotation(self):
        p = problem(density=2.0, fixed_density=1.0, scale=2.0,
                    offset=(1.0, -3.0, 2.0), axis=(1.0, 2.0, 3.0),
                    target_origin=(0.0, -2.0, 2.0))
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
        _, locked_request, _ = faster_material_module._prepare(p, evaluator, solver)
        initial = (1.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0)

        def callback(coefficients, context):
            positive = coefficients[0] > 0.0
            target_vertices = context['design_vertices'] if positive else np.empty((0, 3))
            model_vertices = np.asarray(locked_request.design_vertices) if positive else np.empty((0, 3))
            faces = context['design_faces'] if positive else np.empty((0, 3), dtype=np.int32)
            moments = spin.integrate_raw_mass_moments(
                np.asarray(target_vertices).tolist(), np.asarray(faces).tolist(),
                context['material_density']) if positive else (0.0,) * 10
            return {
                'variable_vertices': target_vertices,
                'variable_vertices_model': model_vertices,
                'variable_faces': faces,
                'variable_raw_moments': moments,
                'quadric_coefficients': coefficients,
                'moment_jacobian': np.eye(10),
                'derivative_defined': True,
                'geometry_hash': 'rotated-geometry-' + str(int(positive)),
                'topology_hash': 'rotated-topology-' + str(int(positive)),
                'evaluation_map_identity': 'rotated-map-' + str(int(positive)),
                'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:test',
                'derivative_identity': 'source-identity-fixture-v1',
            }

        def fake_native_solver(request, mapped_evaluator, *, initial_quadrics,
                               event_sink=None):
            mapped = mapped_evaluator(initial)
            return {'status': 'CONVERGED', 'quadric_coefficients': initial,
                    'trace': [{'material_geometry_hash': mapped['geometry_hash']}],
                    'function_evaluations': 1, 'events': [], 'event_hash': 'fixture',
                    'trace_hash': 'fixture', 'final_geometry_verified': False,
                    'scope': 'numerical_newton_only_external_geometry_not_sealed'}

        def complement_callback(retained, context):
            self.assertGreater(len(retained['variable_faces']), 0)
            return {
                'variable_vertices': np.empty((0, 3)),
                'variable_vertices_model': np.empty((0, 3)),
                'variable_faces': np.empty((0, 3), dtype=np.int32),
                'variable_raw_moments': (0.0,) * 10,
                'moment_jacobian': np.zeros((10, 10)),
                'derivative_defined': False,
                'geometry_hash': 'complement-empty',
                'topology_hash': 'complement-empty-topology',
                'evaluation_map_identity': 'complement-empty-map',
                'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:test',
                'derivative_identity': 'not-applicable-complement',
            }

        timed_stages = []

        def stage_observer(name, operation):
            timed_stages.append(name)
            return operation()

        with patch.object(faster_material_module,
                          'solve_spin_it_faster_with_external_moment_evaluator',
                          fake_native_solver):
            result = faster_material_module.solve_cherchi_geometry_and_seal_faster(
                p, geometry_evaluator=callback, initial_quadrics=(initial,),
                complement_evaluator=complement_callback,
                evaluator_settings=evaluator, solver_settings=solver,
                stage_observer=stage_observer)
        self.assertEqual(result.status, 'CONVERGED', result.seal_error)
        self.assertIsNotNone(result.artifact)
        self.assertTrue(result.artifact.revalidate())
        self.assertEqual(timed_stages, [
            'native_newton_solver', 'external_geometry_moment_map',
            'map_variable_mesh_snapshot',
            'map_variable_independent_integration',
            'map_fixed_material_integration_0',
            'map_target_model_frame_identity',
            'terminal_complement_boolean',
        ])

    def test_external_cherchi_bridge_solves_callback_and_seals_actual_material(self):
        p = problem(density=2.0, fixed_density=1.0, scale=2.0,
                    offset=(1.0, -3.0, 2.0), axis=(0.0, 0.0, 1.0),
                    target_origin=(1.0, -3.0, 2.0))
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
        initial = (1.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0)
        calls = []

        def callback(coefficients, context):
            calls.append(tuple(coefficients))
            positive = coefficients[0] > 0.0
            vertices = context['design_vertices'] if positive else np.empty((0, 3))
            faces = context['design_faces'] if positive else np.empty((0, 3), dtype=np.int32)
            moments = spin.integrate_raw_mass_moments(
                np.asarray(vertices).tolist(), np.asarray(faces).tolist(),
                context['material_density']) if positive else (0.0,) * 10
            from fractions import Fraction
            rational_vertices = tuple(tuple(
                f'{Fraction.from_float(float(value)).numerator}/'
                f'{Fraction.from_float(float(value)).denominator}' for value in row)
                for row in vertices)
            return {
                'variable_vertices': vertices,
                'variable_rational_vertices': rational_vertices,
                'variable_faces': faces,
                'variable_raw_moments': moments,
                'quadric_coefficients': coefficients,
                'moment_jacobian': np.eye(10),
                'derivative_defined': True,
                'geometry_hash': 'geometry-' + str(int(positive)),
                'topology_hash': 'topology-' + str(int(positive)),
                'evaluation_map_identity': 'map-' + str(int(positive)),
                'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:test-fixture',
                'derivative_identity': 'test-fixture-jacobian-v1',
            }

        def fake_native_solver(request, mapped_evaluator, *, initial_quadrics,
                               event_sink=None):
            self.assertEqual(tuple(initial_quadrics[0]), initial)
            mapped = mapped_evaluator(initial)
            self.assertGreater(mapped['raw_moments'][0], mapped['variable_raw_moments'][0])
            return {'status': 'CONVERGED', 'quadric_coefficients': initial,
                    'trace': [{'material_geometry_hash': mapped['geometry_hash']}],
                    'function_evaluations': 1, 'events': [], 'event_hash': 'fixture',
                    'trace_hash': 'fixture', 'final_geometry_verified': False,
                    'scope': 'numerical_newton_only_external_geometry_not_sealed'}

        def complement_callback(retained, context):
            self.assertGreater(len(retained['variable_faces']), 0)
            return {
                'variable_vertices': np.empty((0, 3)),
                'variable_vertices_model': np.empty((0, 3)),
                'variable_faces': np.empty((0, 3), dtype=np.int32),
                'variable_raw_moments': (0.0,) * 10,
                'moment_jacobian': np.zeros((10, 10)),
                'derivative_defined': False,
                'geometry_hash': 'rotated-complement-empty',
                'topology_hash': 'rotated-complement-empty-topology',
                'evaluation_map_identity': 'rotated-complement-empty-map',
                'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:test',
                'derivative_identity': 'not-applicable-complement',
            }

        with patch.object(faster_material_module,
                          'solve_spin_it_faster_with_external_moment_evaluator',
                          fake_native_solver):
            result = faster_material_module.solve_cherchi_geometry_and_seal_faster(
                p, geometry_evaluator=callback, initial_quadrics=(initial,),
                complement_evaluator=complement_callback,
                evaluator_settings=evaluator, solver_settings=solver)

        self.assertEqual(result.status, 'CONVERGED', result.seal_error)
        self.assertIsNotNone(result.artifact)
        self.assertIsNotNone(result.artifact.design_results[0][1].rational_vertices)
        self.assertEqual(len(calls), 1)  # the complement has its own Boolean path
        self.assertTrue(result.artifact.revalidate())
        self.assertEqual(result.identity['operation'], 'SOLVE_CHERCHI_EXTERNAL_GEOMETRY')

    def test_external_cherchi_bridge_accepts_native_seeded_empty_initial_sequence(self):
        p = problem(density=1.0, fixed_density=0.5, scale=2.0)
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(
            equation_profile='DUAL_MONOMIAL_RAW_V1', deterministic_seed=17)
        observed = {}

        def fake_native_solver(request, mapped_evaluator, *, initial_quadrics,
                               event_sink=None):
            observed['seed'] = request.settings.deterministic_seed
            observed['initial_shape'] = np.asarray(initial_quadrics).shape
            return {'status': 'RESOURCE_LIMIT_REACHED', 'quadric_coefficients': None,
                    'trace': [], 'function_evaluations': 0, 'events': [],
                    'event_hash': 'empty', 'trace_hash': 'empty',
                    'final_geometry_verified': False,
                    'scope': 'numerical_newton_only_external_geometry_not_sealed'}

        def unused_geometry(*_args):
            self.fail('geometry callback must not run when the fake solver has no trials')

        with patch.object(faster_material_module,
                          'solve_spin_it_faster_with_external_moment_evaluator',
                          fake_native_solver):
            result = faster_material_module.solve_cherchi_geometry_and_seal_faster(
                p, geometry_evaluator=unused_geometry,
                complement_evaluator=unused_geometry, initial_quadrics=(),
                evaluator_settings=evaluator, solver_settings=solver)

        self.assertEqual(result.status, 'RESOURCE_LIMIT_REACHED')
        self.assertEqual(observed, {'seed': 17, 'initial_shape': (0,)})

    def test_external_cherchi_bridge_records_complement_failure_without_raising(self):
        p = problem(density=1.0, fixed_density=0.5)
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
        initial = (1.0,) + (0.0,) * 9

        def geometry(coefficients, context):
            vertices = np.asarray(context['design_vertices'])
            faces = np.asarray(context['design_faces'], dtype=np.int32)
            moments = spin.integrate_raw_mass_moments(
                vertices.tolist(), faces.tolist(), context['material_density'])
            return {
                'variable_vertices': vertices,
                'variable_vertices_model': vertices,
                'variable_faces': faces,
                'variable_raw_moments': moments,
                'quadric_coefficients': coefficients,
                'moment_jacobian': np.eye(10),
                'derivative_defined': True,
                'geometry_hash': 'converged-retained-geometry',
                'topology_hash': 'converged-retained-topology',
                'evaluation_map_identity': 'converged-retained-map',
                'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:test-fixture',
                'derivative_identity': 'test-fixture-derivative-v1',
            }

        def native_solver(request, mapped_evaluator, *, initial_quadrics,
                          event_sink=None):
            mapped = mapped_evaluator(initial)
            return {'status': 'CONVERGED', 'quadric_coefficients': initial,
                    'trace': [{'material_geometry_hash': mapped['geometry_hash']}],
                    'function_evaluations': 1, 'events': [], 'event_hash': 'fixture',
                    'trace_hash': 'fixture', 'final_geometry_verified': False,
                    'scope': 'numerical_newton_only_external_geometry_not_sealed'}

        def failed_complement(*_args):
            raise ValueError('Cherchi emitted an invalid nested complement boundary')

        with patch.object(faster_material_module,
                          'solve_spin_it_faster_with_external_moment_evaluator',
                          native_solver):
            result = faster_material_module.solve_cherchi_geometry_and_seal_faster(
                p, geometry_evaluator=geometry, complement_evaluator=failed_complement,
                initial_quadrics=(initial,), evaluator_settings=evaluator,
                solver_settings=solver)

        self.assertEqual(result.status, 'FINAL_COMPLEMENT_AUDIT_FAILED')
        self.assertIsNone(result.artifact)
        self.assertEqual(result.seal_error,
                         'ValueError: Cherchi emitted an invalid nested complement boundary')

    def test_external_cherchi_bridge_binds_complement_to_same_quadric_source_and_side(self):
        p = problem(density=1.0, fixed_density=0.5)
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
        initial = (1.0,) + (0.0,) * 9

        def geometry(coefficients, context):
            vertices = np.asarray(context['design_vertices'])
            faces = np.asarray(context['design_faces'], dtype=np.int32)
            moments = spin.integrate_raw_mass_moments(
                vertices.tolist(), faces.tolist(), context['material_density'])
            return {
                'variable_vertices': vertices,
                'variable_vertices_model': vertices,
                'variable_faces': faces,
                'variable_raw_moments': moments,
                'quadric_coefficients': coefficients,
                'moment_jacobian': np.eye(10),
                'derivative_defined': True,
                'geometry_hash': 'retained-geometry',
                'topology_hash': 'retained-topology',
                'evaluation_map_identity': 'retained-map',
                'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:test-fixture',
                'derivative_identity': 'test-fixture-derivative-v1',
                'boolean_operation': 'intersection',
                'quadric_source_identity': 'same-source-sha256',
            }

        def native_solver(request, mapped_evaluator, *, initial_quadrics,
                          event_sink=None):
            mapped = mapped_evaluator(initial)
            return {'status': 'CONVERGED', 'quadric_coefficients': initial,
                    'trace': [{'material_geometry_hash': mapped['geometry_hash']}],
                    'function_evaluations': 1, 'events': [], 'event_hash': 'fixture',
                    'trace_hash': 'fixture', 'final_geometry_verified': False,
                    'scope': 'numerical_newton_only_external_geometry_not_sealed'}

        def wrong_complement(retained, context):
            return {
                'variable_vertices': np.empty((0, 3)),
                'variable_vertices_model': np.empty((0, 3)),
                'variable_faces': np.empty((0, 3), dtype=np.int32),
                'variable_raw_moments': (0.0,) * 10,
                'moment_jacobian': np.zeros((10, 10)),
                'derivative_defined': False,
                'geometry_hash': 'complement-geometry',
                'topology_hash': 'complement-topology',
                'evaluation_map_identity': 'complement-map',
                'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:test-fixture',
                'derivative_identity': 'not-applicable-complement',
                'boolean_operation': 'intersection',
                'quadric_source_identity': 'different-source-sha256',
            }

        with patch.object(faster_material_module,
                          'solve_spin_it_faster_with_external_moment_evaluator',
                          native_solver):
            result = faster_material_module.solve_cherchi_geometry_and_seal_faster(
                p, geometry_evaluator=geometry, complement_evaluator=wrong_complement,
                initial_quadrics=(initial,), evaluator_settings=evaluator,
                solver_settings=solver)

        self.assertEqual(result.status, 'FINAL_COMPLEMENT_AUDIT_FAILED')
        self.assertIn('did not use the retained material quadric source', result.seal_error)

    def test_external_cherchi_bridge_rejects_same_identity_with_different_source_mesh(self):
        """A copied source ID must not substitute for the actual terminal mesh."""
        p = problem(density=1.0, fixed_density=0.5)
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
        initial = (1.0,) + (0.0,) * 9

        def geometry(coefficients, context):
            vertices = np.asarray(context['design_vertices'])
            faces = np.asarray(context['design_faces'], dtype=np.int32)
            moments = spin.integrate_raw_mass_moments(
                vertices.tolist(), faces.tolist(), context['material_density'])
            return {
                'variable_vertices': vertices,
                'variable_vertices_model': vertices,
                'variable_faces': faces,
                'variable_raw_moments': moments,
                'quadric_coefficients': coefficients,
                'moment_jacobian': np.eye(10),
                'derivative_defined': True,
                'geometry_hash': 'retained-geometry',
                'topology_hash': 'retained-topology',
                'evaluation_map_identity': 'retained-map',
                'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:test-fixture',
                'derivative_identity': 'test-fixture-derivative-v1',
                'boolean_operation': 'intersection',
                'quadric_source_identity': 'copied-source-id',
                'quadric_source_vertices': vertices,
                'quadric_source_faces': faces,
            }

        def native_solver(request, mapped_evaluator, *, initial_quadrics,
                          event_sink=None):
            mapped = mapped_evaluator(initial)
            return {'status': 'CONVERGED', 'quadric_coefficients': initial,
                    'trace': [{'material_geometry_hash': mapped['geometry_hash']}],
                    'function_evaluations': 1, 'events': [], 'event_hash': 'fixture',
                    'trace_hash': 'fixture', 'final_geometry_verified': False,
                    'scope': 'numerical_newton_only_external_geometry_not_sealed'}

        def forged_complement(retained, context):
            vertices = np.asarray(retained['quadric_source_vertices'])
            faces = np.asarray(retained['quadric_source_faces'], dtype=np.int32)
            changed_vertices = vertices.copy()
            changed_vertices[0, 0] += 0.125
            design_vertices = np.asarray(context['design_vertices'])
            design_faces = np.asarray(context['design_faces'], dtype=np.int32)
            moments = spin.integrate_raw_mass_moments(
                design_vertices.tolist(), design_faces.tolist(),
                context['material_density'])
            return {
                'variable_vertices': design_vertices,
                'variable_vertices_model': design_vertices,
                'variable_faces': design_faces,
                'variable_raw_moments': moments,
                'moment_jacobian': np.zeros((10, 10)),
                'derivative_defined': False,
                'geometry_hash': 'complement-geometry',
                'topology_hash': 'complement-topology',
                'evaluation_map_identity': 'complement-map',
                'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:test-fixture',
                'derivative_identity': 'not-applicable-complement',
                'boolean_operation': 'subtraction',
                'quadric_source_identity': 'copied-source-id',
                'quadric_source_vertices': changed_vertices,
                'quadric_source_faces': faces,
            }

        with patch.object(faster_material_module,
                          'solve_spin_it_faster_with_external_moment_evaluator',
                          native_solver):
            result = faster_material_module.solve_cherchi_geometry_and_seal_faster(
                p, geometry_evaluator=geometry, complement_evaluator=forged_complement,
                initial_quadrics=(initial,), evaluator_settings=evaluator,
                solver_settings=solver)

        self.assertEqual(result.status, 'FINAL_COMPLEMENT_AUDIT_FAILED')
        self.assertIn('source mesh differs', result.seal_error)

    def test_external_geometry_moment_map_adds_transformed_fixed_material_once(self):
        p = problem(density=2.3, fixed_density=1.4, scale=2.0,
                    offset=(3.0, -2.0, 1.0), axis=(1.0, 2.0, 3.0),
                    target_origin=(2.0, -1.0, 0.5))
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
        _, request, _ = faster_material_module._prepare(p, evaluator, solver)
        axes = np.asarray(request.target_axes)
        origin = np.asarray(request.target_origin)
        model_vertices = np.asarray(request.design_vertices)
        target_vertices = (model_vertices - origin) @ axes
        variable = spin.integrate_raw_mass_moments(
            target_vertices, request.design_faces, request.material_density)
        jacobian = np.arange(100, dtype=float).reshape(10, 10)
        raw = faster_material_module._external_geometry_moment_map(
            p, request,
            {'variable_vertices': target_vertices,
             'variable_vertices_model': model_vertices,
             'variable_faces': request.design_faces,
             'variable_raw_moments': variable,
             'moment_jacobian': jacobian,
             'derivative_defined': True,
             'geometry_hash': 'geometry-fixture',
             'topology_hash': 'topology-fixture',
             'evaluation_map_identity': 'map-fixture',
             'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:fixture',
             'derivative_identity': 'discrete-fixture'})

        fixed_model = []
        normalization = p.normalization
        for region in p.fixed_regions:
            vertices, faces = region.mesh.arrays()
            model = (vertices - np.asarray(normalization['origin'])) / normalization['L']
            target = (model - origin) @ axes
            fixed_model.append(spin.integrate_raw_mass_moments(
                target.tolist(), faces.tolist(), region.density / normalization['R']))
        expected_fixed = tuple(sum(row[i] for row in fixed_model) for i in range(10))
        np.testing.assert_allclose(raw['raw_moments'], np.asarray(expected_fixed) + variable,
                                   rtol=1e-12, atol=1e-12)
        np.testing.assert_array_equal(raw['moment_jacobian'], jacobian)
        np.testing.assert_allclose(raw['variable_raw_moments'], variable,
                                   rtol=1e-12, atol=1e-12)
        self.assertEqual(raw['boolean_backend_identity'],
                         'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:fixture')
        self.assertEqual(raw['variable_geometry_validation'],
                         'EXACT_SINGLE_OUTWARD_BOUNDARY')

    def test_external_geometry_moment_map_rejects_geometry_from_another_quadric(self):
        p = problem()
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
        _, request, _ = faster_material_module._prepare(p, evaluator, solver)
        vertices = np.asarray(request.design_vertices)
        faces = np.asarray(request.design_faces, dtype=np.int32)
        moments = spin.integrate_raw_mass_moments(
            vertices.tolist(), faces.tolist(), request.material_density)
        raw = {
            'variable_vertices': vertices,
            'variable_vertices_model': vertices,
            'variable_faces': faces,
            'variable_raw_moments': moments,
            'moment_jacobian': np.eye(10),
            'derivative_defined': True,
            'geometry_hash': 'stale-geometry',
            'topology_hash': 'stale-topology',
            'evaluation_map_identity': 'stale-map',
            'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:fixture',
            'derivative_identity': 'stale-derivative',
            'quadric_coefficients': np.zeros(10),
        }

        with self.assertRaisesRegex(ValueError, 'quadric coefficients do not match'):
            faster_material_module._external_geometry_moment_map(
                p, request, raw, expected_coefficients=np.ones(10))

    def test_external_cherchi_bridge_binds_callback_geometry_to_requested_quadric(self):
        p = problem()
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
        initial = (1.0,) + (0.0,) * 9

        def stale_geometry(_coefficients, context):
            vertices = np.asarray(context['design_vertices'])
            faces = np.asarray(context['design_faces'], dtype=np.int32)
            moments = spin.integrate_raw_mass_moments(
                vertices.tolist(), faces.tolist(), context['material_density'])
            return {
                'variable_vertices': vertices,
                'variable_vertices_model': vertices,
                'variable_faces': faces,
                'variable_raw_moments': moments,
                'moment_jacobian': np.eye(10),
                'derivative_defined': True,
                'geometry_hash': 'stale-geometry',
                'topology_hash': 'stale-topology',
                'evaluation_map_identity': 'stale-map',
                'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:fixture',
                'derivative_identity': 'stale-derivative',
                'quadric_coefficients': np.zeros(10),
            }

        def fake_native_solver(request, mapped_evaluator, *, initial_quadrics,
                               event_sink=None):
            return mapped_evaluator(initial_quadrics[0])

        with patch.object(faster_material_module,
                          'solve_spin_it_faster_with_external_moment_evaluator',
                          fake_native_solver):
            with self.assertRaisesRegex(ValueError, 'quadric coefficients do not match'):
                faster_material_module.solve_cherchi_geometry_and_seal_faster(
                    p, geometry_evaluator=stale_geometry,
                    complement_evaluator=stale_geometry, initial_quadrics=(initial,),
                    evaluator_settings=evaluator, solver_settings=solver)

    def test_external_geometry_moment_map_requires_consistent_source_frame_mesh(self):
        p = problem(density=2.0, fixed_density=1.0, scale=2.0,
                    offset=(1.0, -3.0, 2.0), axis=(1.0, 2.0, 3.0),
                    target_origin=(0.0, -2.0, 2.0))
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
        _, request, _ = faster_material_module._prepare(p, evaluator, solver)
        model_vertices = np.asarray(request.design_vertices)
        axes = np.asarray(request.target_axes)
        origin = np.asarray(request.target_origin)
        target_vertices = (model_vertices - origin) @ axes
        moments = spin.integrate_raw_mass_moments(
            target_vertices.tolist(), request.design_faces, request.material_density)
        base = {
            'variable_vertices': target_vertices,
            'variable_faces': request.design_faces,
            'variable_raw_moments': moments,
            'moment_jacobian': np.eye(10),
            'derivative_defined': True,
            'geometry_hash': 'frame-geometry',
            'topology_hash': 'frame-topology',
            'evaluation_map_identity': 'frame-map',
            'boolean_backend_identity': 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:fixture',
            'derivative_identity': 'frame-derivative',
        }
        with self.assertRaisesRegex(ValueError, 'non-identity target frames require'):
            faster_material_module._external_geometry_moment_map(p, request, base)
        inconsistent = dict(base, variable_vertices_model=model_vertices + (1e-5, 0.0, 0.0))
        with self.assertRaisesRegex(ValueError, 'coordinates do not match'):
            faster_material_module._external_geometry_moment_map(p, request, inconsistent)

    def test_external_geometry_moment_map_rejects_callback_that_double_counts_fixed_shell(self):
        p = problem(density=2.3, fixed_density=1.4)
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
        _, request, _ = faster_material_module._prepare(p, evaluator, solver)
        vertices = np.asarray(request.design_vertices)
        variable = spin.integrate_raw_mass_moments(
            vertices, request.design_faces, request.material_density)
        bad = dict(variable_raw_moments=np.asarray(variable) + 1.0,
                   variable_vertices=vertices, variable_faces=request.design_faces,
                   moment_jacobian=np.zeros((10, 10)), derivative_defined=True,
                   geometry_hash='g', topology_hash='t',
                   evaluation_map_identity='m',
                   boolean_backend_identity='cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:bad',
                   derivative_identity='d')
        with self.assertRaisesRegex(ValueError, 'independently integrated variable material'):
            faster_material_module._external_geometry_moment_map(p, request, bad)
        bad_indices = dict(bad, variable_raw_moments=variable,
                           variable_faces=np.asarray(request.design_faces, dtype=float))
        with self.assertRaisesRegex(TypeError, 'faces must use an integer dtype'):
            faster_material_module._external_geometry_moment_map(p, request, bad_indices)

    def test_native_external_solver_receives_total_moments_and_variable_jacobian_split(self):
        p = problem(density=2.0, fixed_density=1.3,
                    axis=(1.0, 2.0, 3.0), target_origin=(0.0, 0.0, 0.0))
        settings = spin.SpinSolverConfig(
            equation_profile='DUAL_MONOMIAL_RAW_V1', capture_evaluation_payload=True,
            maximum_restarts=1, maximum_function_evaluations_per_restart=1)
        evaluator = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=1)
        _, request, _ = faster_material_module._prepare(p, evaluator, settings)
        axes = np.asarray(request.target_axes)
        origin = np.asarray(request.target_origin)
        normalization = p.normalization
        fixed_model = []
        for region in p.fixed_regions:
            vertices, faces = region.mesh.arrays()
            model = (vertices - np.asarray(normalization['origin'])) / normalization['L']
            target = (model - origin) @ axes
            fixed_model.append(spin.integrate_raw_mass_moments(
                target.tolist(), faces.tolist(), region.density / normalization['R']))
        fixed = np.asarray(tuple(sum(row[i] for row in fixed_model)
                                 for i in range(10)))
        fixed_target = np.asarray(faster_material_module._transform_moments(
            request.fixed_shell_moments, axes.T, -axes.T @ origin))
        np.testing.assert_allclose(fixed, fixed_target, rtol=1e-12, atol=1e-12)
        variable = np.asarray((0.0,) * 10)
        jacobian = np.arange(100, dtype=float).reshape(10, 10)
        initial = (1.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0)

        def callback(coefficients):
            return {
                'raw_moments': fixed_target + variable,
                'variable_raw_moments': variable,
                'moment_jacobian': jacobian,
                'derivative_defined': True,
                'geometry_hash': 'fixed-plus-empty-variable',
                'topology_hash': 'fixed-topology',
                'evaluation_map_identity': 'map-fixed-plus-empty',
                'boolean_backend_identity': 'cherchi-test-contract-v1',
                'derivative_identity': 'zero-variable-derivative-v1',
            }

        result = spin.solve_spin_it_faster_with_external_moment_evaluator(
            request, callback,
            initial_quadrics=(initial,))
        payload = next(event['evaluation_payload'] for event in result['events']
                       if event.get('evaluation_payload') is not None)
        np.testing.assert_allclose(payload['solver_raw_moments'], fixed_target,
                                   rtol=1e-12, atol=1e-12)
        np.testing.assert_array_equal(payload['variable_raw_moments'], variable)
        dual = np.asarray(initial, dtype=np.float64)
        dual[7:] *= 2.0
        dual /= np.linalg.norm(dual)
        half = dual.copy()
        half[7:] *= 0.5
        geometry_direction = half / np.linalg.norm(half)
        coefficient_scale = np.eye(10)
        coefficient_scale[7:, 7:] *= 0.5
        geometry_jacobian = ((np.eye(10) - np.outer(geometry_direction,
                                                    geometry_direction)) @
                             coefficient_scale / np.linalg.norm(half))
        np.testing.assert_allclose(payload['solver_moment_jacobian'],
                                   jacobian @ geometry_jacobian,
                                   rtol=1e-12, atol=1e-12)

    def test_bounded_normal_resolution_policy_is_validated_and_sent_to_native(self):
        config = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC',
            icosphere_subdivisions=3,
            maximum_icosphere_subdivisions=5,
        )
        vertices, faces = box((-1.,)*3, (1.,)*3)
        native, _ = spin._native_request(spin.SpinItFasterRequest(
            vertices.tolist(), faces.tolist(), vertices.tolist(), faces.tolist()
        ))
        native_settings = spin._native_evaluator_settings(native, config)
        self.assertEqual(
            native_settings.triangulation.maximum_icosphere_subdivisions, 5
        )
        for maximum in (-2, 2, 8):
            with self.subTest(maximum=maximum):
                with self.assertRaises(ValueError):
                    replace(config, maximum_icosphere_subdivisions=maximum)

    def test_solver_exposes_immutable_actual_export_and_derivative_metadata(self):
        self.assertTrue(hasattr(spin.SpinItFasterResult, 'derivative_metadata'))
        v, f = box((-2.,)*3, (2.,)*3)
        d, g = box((-1.,)*3, (1.,)*3)
        request = spin.SpinItFasterRequest(v.tolist(), f.tolist(), d.tolist(), g.tolist(), mode='YO_YO',
            settings=spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1',
                maximum_restarts=1, maximum_function_evaluations_per_restart=1))
        run = spin.solve_spin_it_faster(request, initial_quadrics=(COEFFICIENTS,),
            evaluator_settings=spin.SpinEvaluatorConfig(derivative_method='DISCRETE_ELLIPSOID', icosphere_subdivisions=1))
        self.assertEqual(run.derivative_metadata['method'], 'discrete-fullrank-quadric-boolean-transport-v2')
        copied = run.derivative_metadata
        copied['method'] = 'INVALID'
        self.assertNotEqual(run.derivative_metadata['method'], 'INVALID')
        self.assertIsNone(run.final_raw_moments)
        self.assertIsNone(run.solution_coefficients)

    def test_native_invalid_evaluator_controls_reject_before_solver(self):
        v, f = box((-2.,)*3, (2.,)*3)
        d, g = box((-1.,)*3, (1.,)*3)
        n, request = spin._native_request(spin.SpinItFasterRequest(v.tolist(), f.tolist(), d.tolist(), g.tolist(),
            settings=spin.SpinSolverConfig(maximum_restarts=1)))
        for name, value in (('derivative_method', n.SpinDerivativeMethod(99)),
                            ('minimum_interface_gradient', 0.)):
            config = n.SpinEvaluatorSettings()
            setattr(config, name, value)
            with self.assertRaises(ValueError):
                n.solve_spin_it_faster(request, n.SpinItFasterEvidence.PAPER_FORMULA_REPRODUCTION,
                                      np.array([COEFFICIENTS]), config)
        config = n.SpinEvaluatorSettings()
        triangle = n.SpinTriangulationSettings()
        triangle.icosphere_subdivisions = 8
        config.triangulation = triangle
        with self.assertRaises(ValueError):
            n.solve_spin_it_faster(request, n.SpinItFasterEvidence.PAPER_FORMULA_REPRODUCTION,
                                  np.array([COEFFICIENTS]), config)

    def test_public_settings_and_legacy_native_dictionary_compatibility(self):
        import make_it_stand
        self.assertTrue(hasattr(make_it_stand, 'SpinEvaluatorConfig'))
        v, f = box((-2.,)*3, (2.,)*3)
        d, g = box((-1.,)*3, (1.,)*3)
        request = spin.SpinItFasterRequest(v.tolist(), f.tolist(), d.tolist(), g.tolist())
        n, r = spin._native_request(request)
        value = n.evaluate_spin_it_faster_quadric(r, np.array(COEFFICIENTS))
        self.assertEqual(set(value), {'variable_raw_moments', 'total_raw_moments', 'geometry_hash',
                                    'topology_hash', 'interface_type', 'interface_subtype'})

    def test_opt_in_exposes_actual_discrete_geometry_and_settings_identity(self):
        self.assertTrue(hasattr(spin, 'SpinEvaluatorConfig'), 'immutable evaluator settings missing')
        config = spin.SpinEvaluatorConfig(derivative_method='DISCRETE_ELLIPSOID', icosphere_subdivisions=1)
        p = problem()
        v, f = p.original.arrays()
        d, g = p.design_regions[0].mesh.arrays()
        request = spin.SpinItFasterRequest((v/4).tolist(), f.tolist(), (d/4).tolist(), g.tolist(),
            mode='YO_YO', settings=spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1'))
        e = spin.evaluate_spin_it_faster_quadric(request, COEFFICIENTS, evaluator_settings=config)
        self.assertTrue(e.derivative_defined, e.derivative_diagnostic)
        self.assertEqual(e.derivative_method, 'discrete-fullrank-quadric-boolean-transport-v2')
        self.assertRegex(e.boolean_backend_identity,
            r'^libigl::copyleft::cgal::mesh_boolean;libigl_source=revision=[0-9a-f]+,tree-sha256=[0-9a-f]{64}$')
        self.assertEqual(e.interface_type, 'ELLIPSOID')
        self.assertIsInstance(e.quadric_solid_is_positive_region, bool)
        self.assertGreater(len(e.interface_vertices), 0)
        self.assertGreater(len(e.interface_faces), 0)
        self.assertGreater(len(e.material_faces), 0)
        self.assertEqual(len(e.vertex_provenance), len(e.material_vertices))
        self.assertRegex(e.derivative_scientific_hash, r'^[0-9a-f]{64}$')
        changed = spin.evaluate_spin_it_faster_quadric(request, COEFFICIENTS,
            evaluator_settings=replace(config, minimum_interface_gradient=2e-12))
        self.assertNotEqual(e.request_hash, changed.request_hash)
        self.assertNotEqual(e.derivative_scientific_hash, changed.derivative_scientific_hash)
        geometry_changed = spin.evaluate_spin_it_faster_quadric(request, COEFFICIENTS,
            evaluator_settings=replace(config, icosphere_subdivisions=2))
        self.assertEqual(e.base_request_hash, geometry_changed.base_request_hash)
        self.assertNotEqual(e.request_hash, geometry_changed.request_hash)
        self.assertNotEqual(e.geometry_hash, geometry_changed.geometry_hash)
        adaptive_level_3 = spin.evaluate_spin_it_faster_quadric(
            request, COEFFICIENTS,
            evaluator_settings=replace(
                config, maximum_icosphere_subdivisions=3
            ),
        )
        adaptive_level_4 = spin.evaluate_spin_it_faster_quadric(
            request, COEFFICIENTS,
            evaluator_settings=replace(
                config, maximum_icosphere_subdivisions=4
            ),
        )
        self.assertNotEqual(adaptive_level_3.request_hash, adaptive_level_4.request_hash)
        hyperbolic_resolution_changed = spin.evaluate_spin_it_faster_quadric(
            request, COEFFICIENTS,
            evaluator_settings=replace(config, hyperbolic_cylinder_grid_resolution=25))
        self.assertNotEqual(e.request_hash, hyperbolic_resolution_changed.request_hash)
        parabolic_resolution_changed = spin.evaluate_spin_it_faster_quadric(
            request, COEFFICIENTS,
            evaluator_settings=replace(config, parabolic_cylinder_grid_resolution=25))
        self.assertNotEqual(e.request_hash, parabolic_resolution_changed.request_hash)
        with self.assertRaises(FrozenInstanceError):
            config.icosphere_subdivisions = 2
        for kw in ({'derivative_method': 'UNKNOWN'}, {'icosphere_subdivisions': True},
                   {'icosphere_subdivisions': 8}, {'minimum_interface_gradient': 0.},
                   {'maximum_icosphere_subdivisions': -2},
                   {'maximum_icosphere_subdivisions': 0},
                   {'maximum_icosphere_subdivisions': 8},
                   {'retain_fixed_normal_provenance': 1},
                   {'hyperbolic_cylinder_grid_resolution': 1},
                   {'parabolic_cylinder_grid_resolution': 0}):
            with self.assertRaises(ValueError):
                replace(config, **kw)
        plane = spin.evaluate_spin_it_faster_quadric(request, (.01,1.,0.,0.,0.,0.,0.,0.,0.,0.),
            evaluator_settings=config)
        self.assertFalse(plane.derivative_defined)
        self.assertEqual(plane.derivative_method, 'discrete-fullrank-quadric-boolean-transport-v2')
        reduced_config = replace(config, derivative_method='DISCRETE_REDUCED_STRATUM')
        reduced_plane = spin.evaluate_spin_it_faster_quadric(
            request, (.14,.31,-.47,.826,0.,0.,0.,0.,0.,0.),
            evaluator_settings=reduced_config)
        self.assertFalse(reduced_plane.derivative_defined)
        self.assertTrue(reduced_plane.stratum_derivative_defined,
                        reduced_plane.stratum_derivative_diagnostic)
        self.assertEqual(reduced_plane.derivative_method,
                         'discrete-rank-deficient-stratum-boolean-transport-v1')
        self.assertEqual(len(reduced_plane.stratum_moment_jacobian), 10)
        self.assertEqual(len(reduced_plane.stratum_moment_jacobian[0]), 3)
        self.assertEqual(len(reduced_plane.stratum_coefficient_basis), 10)
        self.assertEqual(len(reduced_plane.stratum_coefficient_basis[0]), 3)
        stratum_offset = (0.013, -0.017, 0.009)
        stratum_request = replace(
            request,
            design_vertices=tuple(
                tuple(value + stratum_offset[axis] for axis, value in enumerate(vertex))
                for vertex in request.design_vertices))
        for coefficients, tangent_dimension, triangulation in (
            ((-.04,0.,0.,0.,1.,2.,0.,0.,0.,0.), 7,
             {'cylinder_angular_resolution': 24, 'truncation_radius': 8.}),
            ((-.2,0.,0.,-1.,1.,2.,0.,0.,0.,0.), 8,
             {'paraboloid_grid_resolution': 16, 'truncation_radius': 2.}),
        ):
            stratum = spin.evaluate_spin_it_faster_quadric(
                stratum_request, coefficients,
                evaluator_settings=replace(reduced_config, **triangulation))
            self.assertFalse(stratum.derivative_defined)
            self.assertTrue(stratum.stratum_derivative_defined,
                            stratum.stratum_derivative_diagnostic)
            self.assertEqual(len(stratum.stratum_moment_jacobian), 10)
            self.assertEqual(len(stratum.stratum_moment_jacobian[0]), tangent_dimension)
            self.assertEqual(len(stratum.stratum_coefficient_basis), 10)
            self.assertEqual(len(stratum.stratum_coefficient_basis[0]), tangent_dimension)
        with self.assertRaisesRegex(ValueError, 'evaluation-only'):
            spin.solve_spin_it_faster(
                request, initial_quadrics=((.14,.31,-.47,.826,0.,0.,0.,0.,0.,0.),),
                evaluator_settings=reduced_config)


class FasterMaterialBridgeTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('make_it_stand.verification.faster_material'),
                             'locked Faster material bridge missing')
        self.bridge = importlib.import_module('make_it_stand.verification.faster_material')
        self.config = spin.SpinEvaluatorConfig(derivative_method='DISCRETE_ELLIPSOID', icosphere_subdivisions=1)
        self.solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1',
            maximum_iterations=1, maximum_restarts=2, maximum_function_evaluations_per_restart=1)

    def evaluate(self, p, coefficients=COEFFICIENTS, **kw):
        return self.bridge.evaluate_and_seal_faster_candidate(p, coefficients,
            evaluator_settings=self.config, solver_settings=self.solver, **kw)

    def _sealed_source_material(self):
        """Create a sealed geometry-only A material fixture for composition tests."""
        from make_it_stand.verification.spin_artifact import SealedSpinArtifact
        left = box((-2., -2., -2.), (0., 2., 2.))
        source = SpinArtifactProblem(original=left, fixed_regions=[],
            design_regions=[('editable', left, 1.0)], length_unit='m', mass_unit='kg',
            target={'mode':'YO_YO','axis':[0.,0.,1.],'origin':[0.,0.,0.],
                    'support_points':[]}, minimum_outer_wall=.1)
        v, f = source.design_regions[0].mesh.arrays()
        frame = source.normalization
        normalized = spin.integrate_raw_mass_moments(
            ((v-np.asarray(frame['origin']))/frame['L']).tolist(), f.tolist(),
            density=1.0/frame['R'])
        return SealedSpinArtifact.create(source, outer=source.original,
            design_results=[('editable', left, None)], final_material=left,
            newly_carved=None,
            solver_claim={'frame':'INPUT_NORMALIZED_V1',
                'order':['M','LX','LY','LZ','QXX','QYY','QZZ','QYZ','QXZ','QXY'],
                'moments':list(normalized)})

    def _floating_b_target(self, **overrides):
        outer = box((-2., -2., -2.), (2., 2., 2.))
        design = box((0., -2., -2.), (2., 2., 2.))
        values = dict(mode='YO_YO', fixed_shell_moments=(0.,)*10,
            material_density=1.0, original_units_per_normalized_unit=1.0)
        values.update(overrides)
        return spin.SpinItFasterRequest(outer[0].tolist(), outer[1].tolist(),
            design[0].tolist(), design[1].tolist(), **values)

    def test_floating_island_material_a_is_exactly_frozen_as_b_fixed_input(self):
        source = self._sealed_source_material()

        composed = self.bridge.compose_sealed_faster_material_as_fixed_region(
            self._floating_b_target(), source, region_prefix='island_a',
            source_to_request_rotation=np.eye(3),
            source_to_request_translation=(0.,0.,0.),
            source_length_to_request_unit=1.0, reference_density=1.0,
            source_length_unit='m', target_length_unit='m',
            source_mass_unit='kg', target_mass_unit='kg')

        self.assertEqual(composed.request.fixed_shell_moments[0], 32.0)
        self.assertEqual(composed.request.fixed_shell_moments[1], -32.0)
        self.assertNotEqual(composed.request.fixed_shell_moments, (0.,)*10)
        binding = composed.binding
        self.assertEqual(binding['source_artifact_scientific_digest'], source.scientific_digest)
        self.assertEqual(binding['source_final_material'], source.final_material.identity)
        self.assertEqual(binding['request_digest'], composed.request.request_digest)
        self.assertEqual(binding['status'], 'FROZEN_AS_FIXED_GEOMETRY')
        self.assertEqual(binding['rigid_transform']['rotation'], np.eye(3).tolist())
        self.assertEqual(binding['rigid_transform']['translation'], [0.,0.,0.])
        half_density = self.bridge.compose_sealed_faster_material_as_fixed_region(
            self._floating_b_target(), source, region_prefix='island_a',
            source_to_request_rotation=np.eye(3), source_to_request_translation=(0.,0.,0.),
            source_length_to_request_unit=1., reference_density=2.,
            source_length_unit='m', target_length_unit='m',
            source_mass_unit='kg', target_mass_unit='kg')
        self.assertEqual(half_density.request.fixed_shell_moments[0], 16.)
        self.assertEqual(half_density.request.fixed_shell_moments[1], -16.)

    def test_floating_island_rejects_a_with_mismatched_scientific_summary(self):
        source = self._sealed_source_material()
        object.__setattr__(source, 'scientific_digest', '0' * 64)

        with self.assertRaisesRegex(ValueError, 'recomputed scientific facts mismatch'):
            self.bridge.compose_sealed_faster_material_as_fixed_region(
                self._floating_b_target(), source, region_prefix='island_a',
                source_to_request_rotation=np.eye(3),
                source_to_request_translation=(0., 0., 0.),
                source_length_to_request_unit=1.0, reference_density=1.0,
                source_length_unit='m', target_length_unit='m',
                source_mass_unit='kg', target_mass_unit='kg')

    def test_floating_island_buoyant_b_rejects_design_touching_outer_surface(self):
        source = self._sealed_source_material()
        target = self._floating_b_target(mode='STATIC_BALANCE', minimum_wall_thickness=0.1)
        displacement = box((-2., -2., -2.), (2., 2., 2.))

        result = self.bridge.solve_floating_island_b_buoyant_from_sealed_a(
            target, source, region_prefix='floating_a',
            source_to_request_rotation=np.eye(3),
            source_to_request_translation=(0., 0., 0.),
            source_length_to_request_unit=1., reference_density=2.,
            source_length_unit='m', target_length_unit='m',
            source_mass_unit='kg', target_mass_unit='kg',
            displacement_vertices=displacement[0].tolist(),
            displacement_faces=displacement[1].tolist(),
            waterline_origin=(0., 0., 0.), waterline_normal=(0., 0., 1.),
            liquid_density=1., initial_parameters=(0., 0., 0.),
            maximum_iterations=60, constraint_tolerance=1e-10)

        self.assertEqual(result.status, 'INVALID_REQUEST')
        self.assertIsNone(result.artifact)
        self.assertEqual(result.identity['operation'], 'FLOATING_ISLAND_BUOYANT_B_FROM_SEALED_A')
        self.assertEqual(result.identity['source_a_scientific_digest'], source.scientific_digest)
        self.assertEqual(result.identity['composition_binding']['status'], 'FROZEN_AS_FIXED_GEOMETRY')
        self.assertIn('sufficiently clear of the outer surface', result.identity['solver_diagnostic'])

    def test_floating_island_buoyant_b_solves_and_seals_with_legal_partition(self):
        a_mesh = box((-1., -1., -1.5), (1., 1., -0.5))
        problem_a = SpinArtifactProblem(
            original=a_mesh, fixed_regions=[],
            design_regions=[('a_core', a_mesh, 1.0)],
            length_unit='m', mass_unit='kg',
            target={'mode': 'YO_YO', 'axis': [0., 0., 1.],
                    'origin': [0., 0., 0.], 'support_points': []},
            minimum_outer_wall=0.1)
        a_vertices, a_faces = problem_a.design_regions[0].mesh.arrays()
        frame = problem_a.normalization
        a_claim = spin.integrate_raw_mass_moments(
            ((a_vertices - np.asarray(frame['origin'])) / frame['L']).tolist(),
            a_faces.tolist(), density=1.0 / frame['R'])
        source = SealedSpinArtifact.create(
            problem_a, outer=problem_a.original,
            design_results=[('a_core', problem_a.design_regions[0].mesh, None)],
            final_material=a_mesh, newly_carved=None,
            solver_claim={'frame': FRAME, 'order': list(ORDER), 'moments': list(a_claim)})

        outer = box((-2., -2., -2.), (2., 2., 2.))
        design = box((-1., -1., -0.5), (1., 1., 1.5))
        fixed = [
            ('outer_x_low', box((-2., -2., -2.), (-1., 2., 2.)), 0.2),
            ('outer_x_high', box((1., -2., -2.), (2., 2., 2.)), 0.2),
            ('outer_y_low', box((-1., -2., -2.), (1., -1., 2.)), 0.2),
            ('outer_y_high', box((-1., 1., -2.), (1., 2., 2.)), 0.2),
            ('central_z_low', box((-1., -1., -2.), (1., 1., -1.5)), 0.2),
            ('central_z_high', box((-1., -1., 1.5), (1., 1., 2.)), 0.2),
        ]
        fixed_moments = np.zeros(10)
        for _region_id, mesh, density in fixed:
            fixed_moments += np.asarray(spin.integrate_raw_mass_moments(
                mesh[0].tolist(), mesh[1].tolist(), density=density))
        request_b = spin.SpinItFasterRequest(
            outer[0].tolist(), outer[1].tolist(), design[0].tolist(), design[1].tolist(),
            fixed_shell_moments=tuple(fixed_moments),
            material_density=1.0, minimum_wall_thickness=0.1)
        result = self.bridge.solve_floating_island_b_buoyant_from_sealed_a(
            request_b, source, region_prefix='floating_a',
            source_to_request_rotation=np.eye(3),
            source_to_request_translation=(0., 0., 0.),
            source_length_to_request_unit=1., reference_density=5.,
            source_length_unit='m', target_length_unit='m',
            source_mass_unit='kg', target_mass_unit='kg',
            existing_fixed_regions=[(region_id, mesh[0].tolist(), mesh[1].tolist(), density)
                                    for region_id, mesh, density in fixed],
            displacement_vertices=outer[0].tolist(), displacement_faces=outer[1].tolist(),
            waterline_origin=(0., 0., 1.2), waterline_normal=(0., 0., 1.),
            liquid_density=0.296875, initial_parameters=(0., 0., 0.),
            maximum_iterations=60, constraint_tolerance=1e-10)

        self.assertEqual(result.status, 'CONVERGED', result.identity)
        self.assertIsNotNone(result.artifact)
        self.assertEqual(result.identity['operation'], 'FLOATING_ISLAND_BUOYANT_B_FROM_SEALED_A')
        self.assertEqual(result.identity['source_a_scientific_digest'], source.scientific_digest)
        self.assertEqual(result.identity['composition_binding']['status'], 'FROZEN_AS_FIXED_GEOMETRY')
        self.assertTrue(result.identity['b_seal']['plane_partition']['complete_partition'])
        self.assertTrue(result.identity['b_seal']['prescribed_waterline_stable'])
        self.assertLessEqual(result.identity['b_seal']['independent_moment_error'], 1e-10)
        self.assertGreater(result.identity['solver_accepted_steps'], 0)
        np.testing.assert_allclose(result.identity['plane_coefficients'],
                                   (0.5, 0., 0., -1.), atol=1e-9)
        self.assertAlmostEqual(result.identity['fixed_material_measurement']['mass'], 11.2)
        self.assertAlmostEqual(result.identity['hydrostatic_measurement']['displaced_volume'], 51.2)
        self.assertTrue(result.identity['solver_trace'])

    def test_two_stage_buoyancy_optimizes_a_then_consumes_its_sealed_material_in_b(self):
        outer = box((-2., -2., -2.), (2., 2., 2.))
        a_design = box((-1.9, -1.9, -1.9), (1.9, 1.9, 1.9))
        a_shell = combine(outer, (a_design[0], a_design[1][:, ::-1]))
        problem_a = SpinArtifactProblem(
            original=outer,
            fixed_regions=[('a_shell', a_shell, 0.1)],
            design_regions=[('a_design', a_design, 1.0)],
            length_unit='m', mass_unit='kg',
            target={'mode': 'BUOYANT_STABILITY', 'axis': [0., 0., 1.],
                    'origin': [0., 0., 0.], 'support_points': []},
            minimum_outer_wall=0.05)
        request_a = spin.SpinItFasterRequest(
            outer[0].tolist(), outer[1].tolist(), a_design[0].tolist(), a_design[1].tolist(),
            material_density=1.0, minimum_wall_thickness=0.05)
        request_a = spin.with_measured_buoyancy_geometry(
            request_a, outer[0].tolist(), outer[1].tolist(),
            a_shell[0].tolist(), a_shell[1].tolist(),
            waterline_origin=(0., 0., 1.2), waterline_normal=(0., 0., 1.),
            liquid_density=(0.9128 + 28.88) / 51.2, fixed_shell_density=0.1)
        solved_a = spin.solve_spin_it_faster_buoyant_plane(
            request_a, initial_parameters=(0., 0., 0.), maximum_iterations=60,
            constraint_tolerance=1e-10)
        self.assertEqual(solved_a['status'], 'CONVERGED', solved_a['diagnostic'])
        self.assertGreater(solved_a['accepted_steps'], 0)
        sealed_a = self.bridge.seal_buoyant_plane_result_as_artifact(
            problem_a, request_a, solved_a,
            displacement_vertices=outer[0].tolist(), displacement_faces=outer[1].tolist(),
            waterline_origin=(0., 0., 1.2), waterline_normal=(0., 0., 1.))
        self.assertEqual(sealed_a.status, 'CONVERGED')

        b_design = sealed_a.artifact.design_results[0][2]
        b_vertices, b_faces = b_design.arrays()
        request_b = spin.SpinItFasterRequest(
            outer[0].tolist(), outer[1].tolist(), b_vertices.tolist(), b_faces.tolist(),
            material_density=1.0, minimum_wall_thickness=0.05)
        b_target_mass = 0.9128 + 28.88 + 12.996
        result_b = self.bridge.solve_floating_island_b_buoyant_from_sealed_a(
            request_b, sealed_a.artifact, region_prefix='actual_a',
            source_to_request_rotation=np.eye(3),
            source_to_request_translation=(0., 0., 0.),
            source_length_to_request_unit=1., reference_density=1.,
            source_length_unit='m', target_length_unit='m',
            source_mass_unit='kg', target_mass_unit='kg',
            displacement_vertices=outer[0].tolist(), displacement_faces=outer[1].tolist(),
            waterline_origin=(0., 0., 1.2), waterline_normal=(0., 0., 1.),
            liquid_density=b_target_mass / 51.2,
            initial_parameters=(0.8, 0., 0.), maximum_iterations=60,
            constraint_tolerance=1e-10)

        self.assertEqual(result_b.status, 'CONVERGED', result_b.identity)
        self.assertIsNotNone(result_b.artifact)
        self.assertEqual(result_b.identity['source_a_scientific_digest'],
                         sealed_a.artifact.scientific_digest)
        self.assertTrue(result_b.identity['composition_binding']['exact_geometry_check'][
            'complete_partition'])
        self.assertTrue(result_b.identity['b_seal']['plane_partition']['complete_partition'])
        self.assertTrue(result_b.identity['b_seal']['prescribed_waterline_stable'])
        self.assertGreater(result_b.identity['solver_accepted_steps'], 0)
        self.assertTrue(sealed_a.artifact.revalidate())
        self.assertTrue(result_b.artifact.revalidate())

        # The final application result must survive the same persistence
        # boundary used by downstream tools; reload and independently recompute
        # its geometry and scientific facts rather than trusting the in-memory
        # solver result or its caller-supplied identity.
        from make_it_stand.verification.spin_bundle import save_spin_artifact, load_spin_artifact
        with tempfile.TemporaryDirectory() as directory:
            bundle = Path(directory) / 'floating-island-b'
            save_spin_artifact(result_b.artifact, bundle)
            restored_b = load_spin_artifact(
                bundle, expected_problem_digest=result_b.artifact.problem.input_digest)
            self.assertEqual(restored_b.scientific_digest,
                             result_b.artifact.scientific_digest)
            self.assertEqual(restored_b.evidence_digest,
                             result_b.artifact.evidence_digest)
            self.assertTrue(restored_b.revalidate())

        a_pose_stability = self.bridge.recompute_sealed_buoyant_stability_at_pose(
            sealed_a.artifact, displacement_vertices=outer[0].tolist(),
            displacement_faces=outer[1].tolist(), rotation=np.eye(3),
            translation=(0., 0., 0.), reference_origin=(0., 0., 0.),
            waterline_normal=(0., 0., 1.),
            liquid_density=(0.9128 + 28.88) / 51.2)
        self.assertAlmostEqual(a_pose_stability['target_displaced_volume'], 51.2, places=9)
        self.assertTrue(a_pose_stability['stable_at_pose'])

        # A prescribed-waterline solver result is not a free-waterline test.
        # Recompute the displacement geometry after a small rigid-body pitch
        # using the actual B displacement hull and the same displaced volume.
        angle = np.deg2rad(2.0)
        pitch = np.asarray(((1., 0., 0.),
                            (0., np.cos(angle), -np.sin(angle)),
                            (0., np.sin(angle), np.cos(angle))))
        tilted_outer = outer[0] @ pitch.T
        equilibrium = spin.find_hydrostatic_equilibrium(
            tilted_outer.tolist(), outer[1].tolist(),
            reference_origin=(0., 0., 0.), waterline_normal=(0., 0., 1.),
            target_displaced_volume=51.2, relative_volume_tolerance=1e-10,
            maximum_iterations=80)
        self.assertLessEqual(abs(equilibrium.volume_residual), 64.0e-10)
        self.assertAlmostEqual(equilibrium.geometry.displaced_volume, 51.2, places=8)
        self.assertGreater(equilibrium.geometry.waterline_area, 0.0)
        stability = self.bridge.recompute_sealed_buoyant_stability_at_pose(
            result_b.artifact, displacement_vertices=outer[0].tolist(),
            displacement_faces=outer[1].tolist(), rotation=pitch,
            translation=(0., 0., 0.), reference_origin=(0., 0., 0.),
            waterline_normal=(0., 0., 1.), liquid_density=b_target_mass / 51.2)
        self.assertEqual(stability['artifact_scientific_digest'],
                         result_b.artifact.scientific_digest)
        self.assertEqual(stability['artifact_evidence_digest'],
                         result_b.artifact.evidence_digest)
        self.assertAlmostEqual(stability['target_displaced_volume'], 51.2, places=9)
        self.assertLessEqual(stability['volume_residual'], 64.0e-10)
        self.assertGreater(stability['minimum_metacentric_margin'], 0.)
        self.assertTrue(stability['stable_at_pose'])
        expected_matrix = (
            np.asarray(stability['waterplane_central_second_moment']) /
            stability['recomputed_displaced_volume'] -
            np.dot(np.asarray(stability['center_of_gravity']) -
                   np.asarray(stability['center_of_buoyancy']), (0., 0., 1.)) * np.eye(2))
        np.testing.assert_allclose(stability['metacentric_stability_matrix'],
                                   expected_matrix, rtol=0., atol=1e-12)
        self.assertEqual(stability['scope'],
                         'SEALED_OUTER_HULL_SMALL_ANGLE_HYDROSTATIC_LINEARIZATION')
        with self.assertRaisesRegex(ValueError, 'must match the sealed artifact outer boundary'):
            self.bridge.recompute_sealed_buoyant_stability_at_pose(
                result_b.artifact, displacement_vertices=box((-1.,)*3, (1.,)*3)[0].tolist(),
                displacement_faces=box((-1.,)*3, (1.,)*3)[1].tolist(), rotation=pitch,
                translation=(0., 0., 0.), reference_origin=(0., 0., 0.),
                waterline_normal=(0., 0., 1.), liquid_density=b_target_mass / 51.2)
        with self.assertRaisesRegex(ValueError, 'rotation must be orthonormal'):
            self.bridge.recompute_sealed_buoyant_stability_at_pose(
                result_b.artifact, displacement_vertices=outer[0].tolist(),
                displacement_faces=outer[1].tolist(), rotation=2. * pitch,
                translation=(0., 0., 0.), reference_origin=(0., 0., 0.),
                waterline_normal=(0., 0., 1.), liquid_density=b_target_mass / 51.2)
        with self.assertRaisesRegex(ValueError, 'target displaced volume must be strictly interior'):
            self.bridge.recompute_sealed_buoyant_stability_at_pose(
                result_b.artifact, displacement_vertices=outer[0].tolist(),
                displacement_faces=outer[1].tolist(), rotation=pitch,
                translation=(0., 0., 0.), reference_origin=(0., 0., 0.),
                waterline_normal=(0., 0., 1.), liquid_density=0.1)

        reverse = np.diag((1., -1., -1.))
        reversed_stability = self.bridge.recompute_sealed_buoyant_stability_at_pose(
            result_b.artifact, displacement_vertices=outer[0].tolist(),
            displacement_faces=outer[1].tolist(), rotation=reverse,
            translation=(0., 0., 0.), reference_origin=(0., 0., 0.),
            waterline_normal=(0., 0., 1.), liquid_density=b_target_mass / 51.2)
        self.assertEqual(reversed_stability['artifact_scientific_digest'],
                         result_b.artifact.scientific_digest)
        self.assertEqual(reversed_stability['stable_at_pose'],
                         reversed_stability['minimum_metacentric_margin'] > 0.)

    def test_converged_buoyant_plane_a_is_sealed_from_actual_complement_geometry(self):
        def pair(mesh):
            return (mesh[0].tolist(), mesh[1].tolist())

        outer = box((-2.,) * 3, (2.,) * 3)
        design = box((-1.9,) * 3, (1.9,) * 3)
        shell = combine(outer, (design[0], design[1][:, ::-1]))
        fixed_regions = [
            ('shell_x_minus', box((-2., -2., -2.), (-1.9, 2., 2.)), 0.1),
            ('shell_x_plus', box((1.9, -2., -2.), (2., 2., 2.)), 0.1),
            ('shell_y_minus', box((-1.9, -2., -2.), (1.9, -1.9, 2.)), 0.1),
            ('shell_y_plus', box((-1.9, 1.9, -2.), (1.9, 2., 2.)), 0.1),
            ('shell_z_minus', box((-1.9, -1.9, -2.), (1.9, 1.9, -1.9)), 0.1),
            ('shell_z_plus', box((-1.9, -1.9, 1.9), (1.9, 1.9, 2.)), 0.1),
        ]
        problem_a = SpinArtifactProblem(
            original=outer,
            fixed_regions=fixed_regions,
            design_regions=[('dense_insert', design, 1.0)],
            length_unit='m', mass_unit='kg',
            target={'mode': 'BUOYANT_STABILITY', 'axis': [0., 0., 1.],
                    'origin': [0., 0., 0.], 'support_points': []},
            minimum_outer_wall=0.05,
        )
        request_a = spin.SpinItFasterRequest(
            *pair(outer), *pair(design), material_density=1.0,
            minimum_wall_thickness=0.05,
        )
        request_a = spin.with_measured_buoyancy_geometry(
            request_a, *pair(outer), *pair(shell),
            waterline_origin=(0., 0., 0.), waterline_normal=(0., 0., 1.),
            liquid_density=1.0, fixed_shell_density=0.1,
        )
        solved_a = spin.solve_spin_it_faster_buoyant_plane(
            request_a, initial_parameters=(0., 0., 0.),
            maximum_iterations=60, constraint_tolerance=1e-10,
        )

        sealed_a = self.bridge.seal_buoyant_plane_result_as_artifact(
            problem_a, request_a, solved_a,
            displacement_vertices=outer[0].tolist(),
            displacement_faces=outer[1].tolist(),
            waterline_origin=(0., 0., 0.),
            waterline_normal=(0., 0., 1.),
        )

        self.assertEqual(sealed_a.status, 'CONVERGED')
        self.assertIsNotNone(sealed_a.artifact)
        self.assertEqual(sealed_a.identity['operation'], 'SEAL_BUOYANT_PLANE_A')
        self.assertTrue(sealed_a.identity['prescribed_waterline_stable'])
        self.assertTrue(sealed_a.identity['plane_partition']['complete_partition'])
        self.assertTrue(sealed_a.identity['plane_interface_stitch']['exact_coordinate_matching'])
        self.assertFalse(sealed_a.identity['plane_interface_stitch']['vertex_coordinates_moved'])
        self.assertLessEqual(sealed_a.identity['plane_mesh_independent_moment_error'], 1e-10)
        self.assertLessEqual(sealed_a.identity['plane_partition_independent_moment_error'], 1e-10)
        self.assertLessEqual(sealed_a.identity['independent_moment_error'], 1e-10)
        self.assertEqual(sealed_a.artifact.problem.input_digest, problem_a.input_digest)
        self.assertEqual(len(sealed_a.artifact.design_results), 1)
        self.assertFalse(sealed_a.artifact.design_results[0][1].explicit_empty)
        self.assertFalse(sealed_a.artifact.design_results[0][2].explicit_empty)

        # The actual converged, prescribed-waterline solver output must pass
        # through the same sealed multi-material CAD path as a plane candidate.
        occt = Path(__file__).resolve().parents[2] / 'build' / 'occt-step-roundtrip-20260922a' / 'occt_step_roundtrip'
        self.assertTrue(occt.is_file(), 'pinned OCCT method-validation probe is missing')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cad = self.bridge.export_spin_it_faster_quadric_material_step(
                sealed_a, occt_executable=occt,
                brep_path=root / 'buoyant-a.brep', step_path=root / 'buoyant-a.step')
            self.assertEqual(cad['status'], 'FASTER_3D_MATERIAL_STEP_BREP_VALIDATED')
            self.assertEqual(cad['quadric_type'], 'PLANE')
            self.assertLessEqual(cad['maximum_normalized_material_moment_error'], 1e-6)
            sidecar = json.loads((root / 'buoyant-a.step.identity.json').read_text(
                encoding='ascii'))
            self.assertEqual(sidecar['artifact_scientific_digest'],
                             sealed_a.artifact.scientific_digest)
            self.assertEqual(sidecar['quadric_type'], 'PLANE')

        unstable = dict(solved_a, stable_at_prescribed_waterline=False)
        with self.assertRaisesRegex(ValueError, 'positive prescribed-waterline stability'):
            self.bridge.seal_buoyant_plane_result_as_artifact(
                problem_a, request_a, unstable,
                displacement_vertices=outer[0].tolist(),
                displacement_faces=outer[1].tolist(),
                waterline_origin=(0., 0., 0.),
                waterline_normal=(0., 0., 1.),
            )

    def test_floating_island_transfer_rejects_nonrigid_and_incompatible_geometry(self):
        source = self._sealed_source_material()
        target = self._floating_b_target()
        with self.assertRaisesRegex(ValueError, 'proper rigid rotation'):
            self.bridge.compose_sealed_faster_material_as_fixed_region(
                target, source, region_prefix='island_a',
                source_to_request_rotation=np.diag((1.,1.,-1.)),
                source_to_request_translation=(0.,0.,0.), source_length_to_request_unit=1.,
                reference_density=1., source_length_unit='m', target_length_unit='m',
                source_mass_unit='kg', target_mass_unit='kg')
        with self.assertRaisesRegex(ValueError, 'proper rigid rotation'):
            self.bridge.compose_sealed_faster_material_as_fixed_region(
                target, source, region_prefix='island_a',
                source_to_request_rotation=np.diag((2.,1.,1.)),
                source_to_request_translation=(0.,0.,0.), source_length_to_request_unit=1.,
                reference_density=1., source_length_unit='m', target_length_unit='m',
                source_mass_unit='kg', target_mass_unit='kg')
        common = dict(source_to_request_rotation=np.eye(3),
            source_length_to_request_unit=1., reference_density=1.,
            source_length_unit='m', target_length_unit='m',
            source_mass_unit='kg', target_mass_unit='kg')
        with self.assertRaisesRegex(ValueError, 'exact geometry check'):
            self.bridge.compose_sealed_faster_material_as_fixed_region(
                target, source, region_prefix='island_a',
                source_to_request_translation=(1.,0.,0.), **common)
        with self.assertRaisesRegex(ValueError, 'units must match'):
            self.bridge.compose_sealed_faster_material_as_fixed_region(
                target, source, region_prefix='island_a',
                source_to_request_translation=(0.,0.,0.), **{**common,'target_length_unit':'cm'})
        with self.assertRaisesRegex(ValueError, 'existing fixed geometry does not match'):
            self.bridge.compose_sealed_faster_material_as_fixed_region(
                replace(target, fixed_shell_moments=(1.,)+(0.,)*9), source,
                region_prefix='island_a', source_to_request_translation=(0.,0.,0.), **common)

    def test_floating_island_b_solver_consumes_a_and_rechecks_final_material(self):
        """A is frozen into B, then B's solved mesh closes the whole partition."""
        def shifted_box(lo, hi):
            vertices, faces = box(lo, hi)
            return vertices, faces

        def pair(mesh):
            return (mesh[0].tolist(), mesh[1].tolist())

        coords = (-2., -1., 0., 1., 2.)
        cells = []
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    cells.append(((i, j, k), shifted_box(
                        (coords[i], coords[j], coords[k]),
                        (coords[i + 1], coords[j + 1], coords[k + 1]))))
        source_mesh = next(mesh for index, mesh in cells if index == (2, 2, 1))
        design_mesh = next(mesh for index, mesh in cells if index == (2, 2, 2))
        fixed = [mesh for index, mesh in cells if index not in {(2, 2, 1), (2, 2, 2)}]

        source_problem = SpinArtifactProblem(
            original=source_mesh,
            fixed_regions=[],
            design_regions=[('editable', source_mesh, 1.5)],
            length_unit='m', mass_unit='kg',
            target={'mode': 'YO_YO', 'axis': [0., 0., 1.],
                    'origin': [0., 0., 0.], 'support_points': []},
            minimum_outer_wall=.05)
        source_v, source_f = source_problem.design_regions[0].mesh.arrays()
        frame = source_problem.normalization
        source_claim = spin.integrate_raw_mass_moments(
            ((source_v - np.asarray(frame['origin'])) / frame['L']).tolist(),
            source_f.tolist(), density=1.5 / frame['R'])
        source = SealedSpinArtifact.create(
            source_problem, outer=source_problem.original,
            design_results=[('editable', source_mesh, None)],
            final_material=source_mesh, newly_carved=None,
            solver_claim={'frame': FRAME, 'order': list(ORDER),
                          'moments': list(source_claim)})

        fixed_moments = np.zeros(10)
        for vertices, faces in fixed:
            props = compute_mass_properties(
                TriangleMesh(tuple(map(tuple, vertices.tolist())),
                             tuple(map(tuple, faces.tolist()))), 1.0)
            q = props.second_moment
            fixed_moments += np.asarray((props.mass, *props.first_moment,
                q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1]))

        outer = shifted_box((-2.,) * 3, (2.,) * 3)
        target_axes = np.asarray(((0., -1., 0.), (1., 0., 0.), (0., 0., 1.)))
        target_origin = np.asarray((0., 0., .25))
        request = spin.SpinItFasterRequest(
            *pair(outer), *pair(design_mesh),
            fixed_shell_moments=tuple(fixed_moments),
            target_axes=target_axes.tolist(), target_origin=target_origin.tolist(),
            mode='STATIC_BALANCE', static_objective='MAXIMUM_MASS',
            material_density=1.0, minimum_wall_thickness=.05,
            settings=spin.SpinSolverConfig(
                maximum_iterations=10, maximum_restarts=2,
                maximum_function_evaluations_per_restart=20,
                deterministic_seed=0, equation_profile='DUAL_MONOMIAL_RAW_V1'))
        existing = [(f'fixed_{i}', *pair(mesh), 1.0) for i, mesh in enumerate(fixed)]
        chain = self.bridge.solve_floating_island_b_from_sealed_a(
            request, source, region_prefix='island_a',
            source_to_request_rotation=np.eye(3),
            source_to_request_translation=(0., 0., 0.),
            source_length_to_request_unit=1.0, reference_density=1.0,
            source_length_unit='m', target_length_unit='m',
            source_mass_unit='kg', target_mass_unit='kg',
            existing_fixed_regions=existing,
            initial_coefficients=(1., -5., 0., 0.),
            maximum_iterations=10, maximum_line_search_steps=12,
            constraint_tolerance=1e-8)

        self.assertEqual(chain.status, 'CONVERGED')
        self.assertGreater(chain.solver_result['accepted_steps'], 0)
        self.assertTrue(chain.identity['final_material_domain_audit']['b_material_subset_design'])
        self.assertEqual(chain.identity['source_artifact_scientific_digest'],
                         source.scientific_digest)
        self.assertLessEqual(chain.identity['independent_normalized_constraints_norm'], 1e-8)
        b_v = np.asarray(chain.solver_result['variable_material_vertices'])
        b_f = np.asarray(chain.solver_result['variable_material_faces'])
        b_props = compute_mass_properties(
            TriangleMesh(tuple(map(tuple, b_v.tolist())),
                         tuple(map(tuple, b_f.tolist()))), 1.0)
        q = b_props.second_moment
        b_moments = np.asarray((b_props.mass, *b_props.first_moment,
            q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1]))
        model_to_target = target_axes.T
        translation = -model_to_target @ target_origin
        composed_fixed = np.asarray(chain.composition.request.fixed_shell_moments)
        m, lx, ly, lz, xx, yy, zz, yz, xz, xy = composed_fixed
        first = model_to_target @ np.asarray((lx, ly, lz))
        second = model_to_target @ np.asarray(((xx, xy, xz), (xy, yy, yz),
                                                (xz, yz, zz))) @ model_to_target.T
        second += np.outer(first, translation) + np.outer(translation, first) \
            + m * np.outer(translation, translation)
        first += m * translation
        fixed_target_moments = np.asarray((m, *first, second[0, 0], second[1, 1],
            second[2, 2], second[1, 2], second[0, 2], second[0, 1]))
        np.testing.assert_allclose(
            b_moments + fixed_target_moments,
            chain.solver_result['total_raw_moments'], rtol=1e-10, atol=1e-10)
        self.assertLessEqual(chain.identity['independent_total_moment_residual'], 1e-10)
        if os.environ.get('MIB_PRINT_FLOATING_CHAIN_IDENTITY') == '1':
            print('FLOATING_CHAIN_IDENTITY=' + json.dumps(
                chain.identity, sort_keys=True, separators=(',', ':')))

    def test_floating_island_chain_does_not_promote_unverified_solver_status(self):
        source = self._sealed_source_material()
        transfer = dict(
            region_prefix='island_a', source_to_request_rotation=np.eye(3),
            source_to_request_translation=(0., 0., 0.),
            source_length_to_request_unit=1., reference_density=1.,
            source_length_unit='m', target_length_unit='m',
            source_mass_unit='kg', target_mass_unit='kg')
        for solver_value, expected in (
                ({'status': 'CONVERGED',
                  'paper_optimality_conditions_verified': False},
                 'OPTIMALITY_NOT_VERIFIED'),
                ({'status': 'LINE_SEARCH_EXHAUSTED'}, 'LINE_SEARCH_EXHAUSTED')):
            with self.subTest(status=solver_value['status']), patch.object(
                    self.bridge, 'solve_spin_it_faster_static_plane_family',
                    return_value=solver_value):
                result = self.bridge.solve_floating_island_b_from_sealed_a(
                    self._floating_b_target(), source, **transfer,
                    initial_coefficients=(1., 0., 0., 0.))
                self.assertEqual(result.status, expected)
                self.assertFalse(result.identity['final_material_domain_audit'][
                    'b_material_subset_design'])

    def test_complete_material_includes_actual_fixed_shell_and_roundtrips(self):
        p = problem()
        result = self.evaluate(p)
        artifact = result.artifact
        self.assertEqual(result.status, 'NOT_OPTIMIZED_CANDIDATE')
        self.assertEqual(result.identity['solver_target_feasibility'], 'UNKNOWN')
        self.assertEqual(result.identity['project_status'], 'PARTIAL_BLOCKED')
        self.assertEqual(result.identity['physical_validation'], 'NOT_IN_SCOPE')
        self.assertEqual(result.identity['problem_digest'], p.input_digest)
        self.assertEqual(artifact.outer, p.original)
        self.assertEqual(artifact.facts['geometry_status'], 'MESH_MATERIAL_VALIDATED')
        np.testing.assert_allclose(result.identity['fixed_moments_input_normalized'],
            [7/8,0,0,0,31/384,31/384,31/384,0,0,0], atol=1e-14, rtol=0)
        regions = artifact.facts['regions']
        fixed = next(r for r in regions if r['kind']=='fixed')
        self.assertAlmostEqual(fixed['physical']['mass'], 95.2)
        self.assertGreater(artifact.facts['physical']['mass'], 95.2)
        self.assertGreater(artifact.final_material.face_count, p.original.face_count)
        self.assertLess(max(artifact.facts['moment_comparison']['absolute_errors']), 1e-6)
        repeat = self.evaluate(p)
        self.assertEqual(result.scientific_digest, repeat.scientific_digest)
        identity = result.identity
        identity['coefficients'][0] = 77.
        self.assertEqual(result.identity['coefficients'][0], COEFFICIENTS[0])
        with self.assertRaises((FrozenInstanceError, AttributeError)):
            result.artifact = None
        from make_it_stand.verification.spin_bundle import save_spin_artifact, load_spin_artifact
        with tempfile.TemporaryDirectory() as root:
            save_spin_artifact(artifact, root+'/candidate')
            restored = load_spin_artifact(root+'/candidate', expected_problem_digest=p.input_digest)
            self.assertEqual(restored.scientific_digest, artifact.scientific_digest)

    def test_density_once_and_translated_scaled_rotated_target(self):
        base = self.evaluate(problem())
        changed = self.evaluate(problem(fixed_density=.85))
        self.assertAlmostEqual(base.artifact.facts['physical']['mass'] - changed.artifact.facts['physical']['mass'], 47.6)
        moved = self.evaluate(problem(scale=2., offset=(2**20, -2**19, 2**18), axis=(1.,0.,0.)))
        self.assertEqual(moved.artifact.outer, moved.artifact.problem.original)
        self.assertAlmostEqual(moved.artifact.facts['physical']['mass'], 8*base.artifact.facts['physical']['mass'], places=6)
        self.assertNotEqual(base.scientific_digest, moved.scientific_digest)

    def test_multiple_actual_fixed_regions_keep_individual_densities(self):
        p = problem()
        slabs = [((-2.,-2.,-2.), (2.,2.,-1.), .85),
                 ((-2.,-2.,1.), (2.,2.,2.), 1.7),
                 ((-2.,-2.,-1.), (-1.,2.,1.), 1.7),
                 ((1.,-2.,-1.), (2.,2.,1.), 1.7),
                 ((-1.,-2.,-1.), (1.,-1.,1.), 1.7),
                 ((-1.,1.,-1.), (1.,2.,1.), 1.7)]
        many = SpinArtifactProblem(original=p.original,
            fixed_regions=[(str(i), box(lo,hi), rho) for i,(lo,hi,rho) in enumerate(slabs)],
            design_regions=[('design',p.design_regions[0].mesh,1.7)],
            length_unit='m', mass_unit='kg', target=p.target, minimum_outer_wall=.25)
        result = self.evaluate(many)
        fixed = [r for r in result.artifact.facts['regions'] if r['kind']=='fixed']
        self.assertEqual(len(fixed), 6)
        self.assertAlmostEqual(sum(r['physical']['mass'] for r in fixed), 81.6)
        self.assertAlmostEqual(result.identity['fixed_moments_input_normalized'][0], .75)
        self.assertAlmostEqual(result.identity['fixed_moments_input_normalized'][3], 3/64)

    def test_nonzero_target_origin_rotates_full_actual_native_moments(self):
        result = self.evaluate(problem(axis=(1.,0.,0.), target_origin=(.125,-.125,.25)))
        # Literal cyclic target axes: target (x,y,z) -> normalized (z,x,y).
        native = result.identity['native_evaluations'][0]['total_raw_moments']
        m, x, y, z, xx, yy, zz, yz, xz, xy = native
        tx, ty, tz = 1/32, -1/32, 1/16
        expected = [m, z+m*tx, x+m*ty, y+m*tz,
                    zz+2*tx*z+m*tx*tx, xx+2*ty*x+m*ty*ty, yy+2*tz*y+m*tz*tz,
                    xy+ty*y+tz*x+m*ty*tz, yz+tx*y+tz*z+m*tx*tz,
                    xz+tx*x+ty*z+m*tx*ty]
        independent = result.artifact.facts['normalized']['moments']
        np.testing.assert_allclose(independent, expected, atol=3e-14, rtol=0)
        self.assertGreater(max(abs(x) for x in independent[7:]), 1e-6)

    def test_non_dyadic_frame_rounding_is_rejected_without_repair(self):
        with self.assertRaisesRegex(ValueError, 'partition|boundary') as failure:
            self.evaluate(problem(axis=(1.,2.,3.)))
        print('Expected non-dyadic frame rejection:', failure.exception)

    def test_actual_unsuccessful_solve_keeps_events_and_has_no_seal(self):
        p = problem()
        result = self.bridge.solve_and_seal_faster(p, initial_quadrics=(COEFFICIENTS, COEFFICIENTS),
            evaluator_settings=self.config, solver_settings=self.solver)
        self.assertNotEqual(result.status, 'CONVERGED')
        self.assertIsNone(result.artifact)
        self.assertGreater(len(result.identity['native_run']['events']), 0)
        self.assertEqual(result.identity['native_run']['solve_function_evaluations'], 2)
        self.assertEqual(result.identity['native_run']['solve_restart_function_evaluations'], [1,1])

    def test_fullrank_solver_result_is_sealed_and_cad_moment_mismatch_fails_closed(self):
        solve_fullrank = getattr(self.bridge, 'solve_fullrank_and_seal_faster', None)
        self.assertTrue(callable(solve_fullrank),
                        'full-rank solver output is not connected to the locked material seal')
        def shifted(mesh, dz):
            vertices, faces = mesh
            return vertices + np.asarray((0., 0., dz)), faces

        outer = shifted(box((-1.,) * 3, (1.,) * 3), 2.)
        design = shifted(box((-.8,) * 3, (.8,) * 3), 2.)
        fixed_regions = [
            ('shell_x_minus', shifted(box((-1., -1., -1.), (-.8, 1., 1.)), 2.), 1.0),
            ('shell_x_plus', shifted(box((.8, -1., -1.), (1., 1., 1.)), 2.), 1.0),
            ('shell_y_minus', shifted(box((-.8, -1., -1.), (.8, -.8, 1.)), 2.), 1.0),
            ('shell_y_plus', shifted(box((-.8, .8, -1.), (.8, 1., 1.)), 2.), 1.0),
            ('shell_z_minus', shifted(box((-.8, -.8, -1.), (.8, .8, -.8)), 2.), 1.0),
            ('shell_z_plus', shifted(box((-.8, -.8, .8), (.8, .8, 1.)), 2.), 1.0),
        ]
        locked = SpinArtifactProblem(
            original=outer, fixed_regions=fixed_regions,
            design_regions=[('design', design, 1.0)],
            length_unit='m', mass_unit='kg',
            target={'mode': 'SPINNING_TOP', 'axis': [0., 0., 1.],
                    'origin': [0., 0., 2.], 'support_points': []},
            minimum_outer_wall=.15)
        seed = np.asarray((-.31, .07, -.09, .11, .83, 1.17, -.72,
                           .06, -.04, .03), dtype=np.float64)
        seed[0] += 1e-4
        seed[6] -= 1e-4
        config = spin.SpinEvaluatorConfig(
            derivative_method='DISCRETE_FULL_RANK_QUADRIC',
            icosphere_subdivisions=1)
        solver = spin.SpinSolverConfig(
            equation_profile='DUAL_MONOMIAL_RAW_V1', maximum_iterations=60,
            maximum_restarts=1, maximum_function_evaluations_per_restart=60)

        result = solve_fullrank(
            locked, initial_quadrics=(seed,), evaluator_settings=config,
            solver_settings=solver, center_of_mass_weight=1.0)

        self.assertEqual(result.status, 'CONVERGED')
        self.assertIsNotNone(result.artifact)
        self.assertEqual(result.identity['operation'], 'SOLVE_FULL_RANK')
        self.assertEqual(result.identity['request_overrides']['center_of_mass_weight'], 1.0)
        run = result.identity['native_run']
        self.assertEqual(run['solver_kind'], 'FULL_RANK_DISCRETE')
        self.assertGreater(len(run['trace']), 1)
        self.assertTrue(run['final_geometry_verified'])
        self.assertTrue(run['final_derivative_defined'])
        self.assertEqual(len(result.artifact.facts['normalized']['moments']), 10)
        occt = Path(__file__).resolve().parents[2] / 'build' / 'occt-step-roundtrip-20260922a' / 'occt_step_roundtrip'
        self.assertTrue(occt.is_file(), 'pinned OCCT method-validation probe is missing')
        with tempfile.TemporaryDirectory() as directory:
            brep_path = Path(directory) / 'optimized.brep'
            step_path = Path(directory) / 'optimized.step'
            with self.assertRaisesRegex(ValueError,
                                        'analytic OCCT material differs from the sealed triangulated material'):
                self.bridge.export_spin_it_faster_quadric_material_step(
                    result, occt_executable=occt,
                    brep_path=brep_path, step_path=step_path)
            self.assertFalse(brep_path.exists(), 'rejected B-rep must not escape validation')
            self.assertFalse(step_path.exists(), 'rejected STEP must not escape validation')
            self.assertFalse(step_path.with_suffix('.step.identity.json').exists(),
                             'rejected CAD output must not leave a trusted sidecar')

    def test_sealed_plane_candidate_roundtrips_as_exact_multi_material_cad(self):
        """The dedicated plane-family application must reach the sealed CAD path."""
        from make_it_stand.verification.faster_material import (
            export_spin_it_faster_quadric_material_step,
        )
        outer = box((-2.,) * 3, (2.,) * 3)
        design = box((-1.,) * 3, (1.,) * 3)
        fixed = [
            ('shell_x_minus', box((-2., -2., -2.), (-1., 2., 2.)), 0.8),
            ('shell_x_plus', box((1., -2., -2.), (2., 2., 2.)), 0.8),
            ('shell_y_minus', box((-1., -2., -2.), (1., -1., 2.)), 0.8),
            ('shell_y_plus', box((-1., 1., -2.), (1., 2., 2.)), 0.8),
            ('shell_z_minus', box((-1., -1., -2.), (1., 1., -1.)), 0.8),
            ('shell_z_plus', box((-1., -1., 1.), (1., 1., 2.)), 0.8),
        ]
        p = SpinArtifactProblem(
            original=outer, fixed_regions=fixed,
            design_regions=[('design', design, 1.4)],
            length_unit='m', mass_unit='kg',
            target={'mode': 'YO_YO', 'axis': [0., 0., 1.],
                    'origin': [0., 0., 0.], 'support_points': []},
            minimum_outer_wall=0.1)
        coefficients = np.asarray((0., 0., 0., 1., 0., 0., 0., 0., 0., 0.))
        candidate = self.bridge.evaluate_and_seal_faster_candidate(
            p, coefficients,
            evaluator_settings=spin.SpinEvaluatorConfig(
                derivative_method='DISCRETE_REDUCED_STRATUM',
                truncation_radius=2.0),
            solver_settings=self.solver)
        self.assertEqual(candidate.status, 'NOT_OPTIMIZED_CANDIDATE')
        self.assertIsNotNone(candidate.artifact)
        self.assertIn('PLANE', candidate.identity['native_evaluations'][0]['interface_type'])
        occt = Path(__file__).resolve().parents[2] / 'build' / 'occt-step-roundtrip-20260922a' / 'occt_step_roundtrip'
        self.assertTrue(occt.is_file(), 'pinned OCCT method-validation probe is missing')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cad = export_spin_it_faster_quadric_material_step(
                candidate, occt_executable=occt,
                brep_path=root / 'plane.brep', step_path=root / 'plane.step')
            self.assertEqual(cad['status'], 'FASTER_3D_MATERIAL_STEP_BREP_VALIDATED')
            self.assertEqual(cad['quadric_type'], 'PLANE')
            self.assertLessEqual(cad['maximum_normalized_material_moment_error'], 1e-6)
            self.assertTrue((root / 'plane.brep').is_file())
            self.assertTrue((root / 'plane.step').is_file())
            self.assertEqual(len(cad['material_region_ids']), 7)
            self.assertEqual(cad['entity_roles'], ['design'] + ['fixed'] * 6)
            sidecar = json.loads((root / 'plane.step.identity.json').read_text(
                encoding='ascii'))
            self.assertEqual(sidecar['cad_input_digest'], cad['cad_input_digest'])
            self.assertEqual(sidecar['quadric_type'], 'PLANE')
            self.assertEqual(len(sidecar['material_region_ids']), 7)
            self.assertRegex(sidecar['occt_executable_sha256'], r'^[0-9a-f]{64}$')
            self.assertRegex(sidecar['brep_sha256'], r'^[0-9a-f]{64}$')
            self.assertRegex(sidecar['step_sha256'], r'^[0-9a-f]{64}$')
            print(json.dumps({key: sidecar[key] for key in (
                'cad_input_digest', 'artifact_scientific_digest',
                'maximum_normalized_material_moment_error', 'brep_sha256',
                'step_sha256', 'entity_mass_moments', 'raw_mass_moments')},
                sort_keys=True))

    def test_sealed_paraboloid_candidate_rejects_cad_discretization_mismatch(self):
        """A coarse sealed paraboloid must not publish mismatched analytic CAD."""
        from make_it_stand.verification.faster_material import (
            evaluate_and_seal_faster_candidate,
            export_spin_it_faster_quadric_material_step,
        )
        outer = box((-2.,) * 3, (2.,) * 3)
        design = box((-1.,) * 3, (1.,) * 3)
        fixed = [
            ('shell_x_minus', box((-2., -2., -2.), (-1., 2., 2.)), 0.8),
            ('shell_x_plus', box((1., -2., -2.), (2., 2., 2.)), 0.8),
            ('shell_y_minus', box((-1., -2., -2.), (1., -1., 2.)), 0.8),
            ('shell_y_plus', box((-1., 1., -2.), (1., 2., 2.)), 0.8),
            ('shell_z_minus', box((-1., -1., -2.), (1., 1., -1.)), 0.8),
            ('shell_z_plus', box((-1., -1., 1.), (1., 1., 2.)), 0.8),
        ]
        p = SpinArtifactProblem(
            original=outer, fixed_regions=fixed,
            design_regions=[('design', design, 1.4)],
            length_unit='m', mass_unit='kg',
            target={'mode': 'YO_YO', 'axis': [0., 0., 1.],
                    'origin': [0., 0., 0.], 'support_points': []},
            minimum_outer_wall=0.1)
        coefficients = np.asarray((0., 0., 0., -1., 1., 1., 0., 0., 0., 0.))
        coefficients /= np.linalg.norm(coefficients)
        candidate = evaluate_and_seal_faster_candidate(
            p, coefficients,
            evaluator_settings=spin.SpinEvaluatorConfig(
                derivative_method='DISCRETE_REDUCED_STRATUM',
                paraboloid_grid_resolution=32),
            solver_settings=self.solver)
        self.assertEqual(candidate.status, 'NOT_OPTIMIZED_CANDIDATE')
        self.assertIsNotNone(candidate.artifact)
        self.assertEqual(candidate.identity['native_evaluations'][0]['interface_type'],
                         'PARABOLOID')
        occt = Path(__file__).resolve().parents[2] / 'build' / 'occt-step-roundtrip-20260922a' / 'occt_step_roundtrip'
        self.assertTrue(occt.is_file(), 'pinned OCCT method-validation probe is missing')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(ValueError, 'normalized moment error'):
                export_spin_it_faster_quadric_material_step(
                    candidate, occt_executable=occt,
                    brep_path=root / 'paraboloid.brep', step_path=root / 'paraboloid.step')
            self.assertFalse((root / 'paraboloid.brep').exists())
            self.assertFalse((root / 'paraboloid.step').exists())
            self.assertFalse((root / 'paraboloid.step.identity.json').exists())

    def test_sealed_hyperbolic_paraboloid_candidate_uses_saddle_cad_carrier(self):
        """A rank-two indefinite paraboloid must not enter the elliptic CAD path."""
        from make_it_stand.verification.faster_material import (
            evaluate_and_seal_faster_candidate,
            export_spin_it_faster_quadric_material_step,
        )
        outer = box((-2.,) * 3, (2.,) * 3)
        design = box((-1.,) * 3, (1.,) * 3)
        fixed = [
            ('shell_x_minus', box((-2., -2., -2.), (-1., 2., 2.)), 0.8),
            ('shell_x_plus', box((1., -2., -2.), (2., 2., 2.)), 0.8),
            ('shell_y_minus', box((-1., -2., -2.), (1., -1., 2.)), 0.8),
            ('shell_y_plus', box((-1., 1., -2.), (1., 2., 2.)), 0.8),
            ('shell_z_minus', box((-1., -1., -2.), (1., 1., -1.)), 0.8),
            ('shell_z_plus', box((-1., -1., 1.), (1., 1., 2.)), 0.8),
        ]
        problem = SpinArtifactProblem(
            original=outer, fixed_regions=fixed,
            design_regions=[('design', design, 1.4)],
            length_unit='m', mass_unit='kg',
            target={'mode': 'YO_YO', 'axis': [0., 0., 1.],
                    'origin': [0., 0., 0.], 'support_points': []},
            minimum_outer_wall=0.1)
        # Rotate the small saddle away from every coordinate axis. The weak
        # opposite-sign transverse eigenvalues exercise the general frame
        # transform while keeping one polarity inside the existing strict CAD
        # moment gate; the reversed polarity is retained as a fail-closed
        # diagnostic at two mesh resolutions.
        axis = np.asarray((1., 2., 3.), dtype=np.float64) / np.sqrt(14.)
        radial_x = np.cross(axis, np.asarray((0., 0., 1.)))
        radial_x /= np.linalg.norm(radial_x)
        radial_y = np.cross(axis, radial_x)
        quadratic = 1e-4 * (np.outer(radial_x, radial_x) -
                             np.outer(radial_y, radial_y))
        coefficients = np.asarray((0., *(-axis), quadratic[0, 0],
                                   quadratic[1, 1], quadratic[2, 2],
                                   quadratic[1, 2], quadratic[0, 2],
                                   quadratic[0, 1]), dtype=np.float64)
        coefficients /= np.linalg.norm(coefficients)
        candidate = evaluate_and_seal_faster_candidate(
            problem, coefficients,
            evaluator_settings=spin.SpinEvaluatorConfig(
                derivative_method='DISCRETE_REDUCED_STRATUM',
                paraboloid_grid_resolution=64),
            solver_settings=self.solver)
        self.assertEqual(candidate.status, 'NOT_OPTIMIZED_CANDIDATE')
        self.assertIsNotNone(candidate.artifact)
        self.assertEqual(candidate.identity['native_evaluations'][0]['interface_type'],
                         'PARABOLOID')
        occt = Path(__file__).resolve().parents[2] / 'build' / 'occt-step-roundtrip-20260922a' / 'occt_step_roundtrip'
        self.assertTrue(occt.is_file(), 'pinned OCCT method-validation probe is missing')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cad = export_spin_it_faster_quadric_material_step(
                candidate, occt_executable=occt,
                brep_path=root / 'saddle.brep', step_path=root / 'saddle.step')
            self.assertEqual(cad['quadric_type'], 'HYPERBOLIC_PARABOLOID')
            self.assertLessEqual(cad['maximum_normalized_material_moment_error'], 1e-6)
            self.assertTrue((root / 'saddle.brep').is_file())
            self.assertTrue((root / 'saddle.step').is_file())
            reversed_candidate = evaluate_and_seal_faster_candidate(
                problem, -coefficients,
                evaluator_settings=spin.SpinEvaluatorConfig(
                    derivative_method='DISCRETE_REDUCED_STRATUM',
                    paraboloid_grid_resolution=64),
                solver_settings=self.solver)
            self.assertEqual(reversed_candidate.status, 'NOT_OPTIMIZED_CANDIDATE')
            reversed_cad = export_spin_it_faster_quadric_material_step(
                reversed_candidate, occt_executable=occt,
                brep_path=root / 'saddle-reversed.brep',
                step_path=root / 'saddle-reversed.step')
            self.assertLessEqual(reversed_cad['maximum_normalized_material_moment_error'], 1e-6)
            self.assertTrue((root / 'saddle-reversed.brep').is_file())
            self.assertTrue((root / 'saddle-reversed.step').is_file())
            self.assertTrue((root / 'saddle-reversed.step.identity.json').is_file())
            whole_design = np.asarray((0.125, 0., 0., 0.,
                                       0.0026041666666666665,
                                       0.0026041666666666665,
                                       0.0026041666666666665, 0., 0., 0.))
            np.testing.assert_allclose(
                np.asarray(cad['entity_mass_moments'][0]) +
                np.asarray(reversed_cad['entity_mass_moments'][0]),
                whole_design, rtol=0., atol=1e-9)
            refined_candidate = evaluate_and_seal_faster_candidate(
                problem, -coefficients,
                evaluator_settings=spin.SpinEvaluatorConfig(
                    derivative_method='DISCRETE_REDUCED_STRATUM',
                    paraboloid_grid_resolution=128),
                solver_settings=self.solver)
            refined_cad = export_spin_it_faster_quadric_material_step(
                refined_candidate, occt_executable=occt,
                brep_path=root / 'saddle-reversed-r128.brep',
                step_path=root / 'saddle-reversed-r128.step')
            self.assertLessEqual(refined_cad['maximum_normalized_material_moment_error'], 1e-6)
            self.assertTrue((root / 'saddle-reversed-r128.brep').is_file())
            self.assertTrue((root / 'saddle-reversed-r128.step').is_file())
            self.assertTrue((root / 'saddle-reversed-r128.step.identity.json').is_file())

    def test_full_second_moment_transform_uses_offdiagonals_and_translation(self):
        # Independent Fraction-derived anchor supplied by controller, not native/bridge helpers.
        actual = self.bridge._transform_moments([2.,1.,-2.,3.,4.,5.,6.,.75,-.25,.5],
            [[.6,-.8,0.],[.8,.6,0.],[0.,0.,1.]], [.25,-.5,.75])
        np.testing.assert_allclose(actual,
            [2.,2.7,-1.4,4.5,1077/200,287/50,93/8,-23/10,81/40,-207/100], atol=2e-15, rtol=0)

    def test_actual_native_settings_cannot_share_effective_run_identity(self):
        p = problem()
        first = self.bridge.solve_and_seal_faster(p, initial_quadrics=(COEFFICIENTS,),
            evaluator_settings=self.config, solver_settings=self.solver)
        second = self.bridge.solve_and_seal_faster(p, initial_quadrics=(COEFFICIENTS,),
            evaluator_settings=replace(self.config, minimum_interface_gradient=2e-12), solver_settings=self.solver)
        self.assertNotEqual(first.identity['native_run']['input_hash'], second.identity['native_run']['input_hash'])
        self.assertNotEqual(first.scientific_digest, second.scientific_digest)

    def test_seal_rejects_tampered_moments_outer_cavity_and_fixed_partition(self):
        p = problem()
        artifact = self.evaluate(p).artifact
        from make_it_stand.verification.spin_artifact import SealedSpinArtifact
        kwargs = dict(outer=artifact.outer, design_results=artifact.design_results,
            final_material=artifact.final_material, newly_carved=artifact.newly_carved,
            solver_claim=artifact.solver_claim)
        claim = artifact.solver_claim
        claim['moments'][7] += .01
        with self.assertRaisesRegex(ValueError, 'moment claim mismatch'):
            SealedSpinArtifact.create(p, **{**kwargs, 'solver_claim': claim})
        with self.assertRaisesRegex(ValueError, 'outer canonical geometry bytes changed'):
            SealedSpinArtifact.create(p, **{**kwargs, 'outer': box((-3.,)*3,(3.,)*3)})
        with self.assertRaisesRegex(ValueError, 'partition'):
            SealedSpinArtifact.create(p, **{**kwargs, 'newly_carved': None})
        with self.assertRaisesRegex(ValueError, 'partition'):
            SpinArtifactProblem(original=p.original, fixed_regions=[('shell', p.original, 1.7)],
                design_regions=[('design', p.design_regions[0].mesh, 1.7)], length_unit='m', mass_unit='kg',
                target=p.target, minimum_outer_wall=.25)

    def test_scope_checks_reject_multiple_designs_buoyancy_and_unrepresented_support(self):
        p = problem()
        def target_problem(target):
            return SpinArtifactProblem(original=p.original,
                fixed_regions=[(r.region_id,r.mesh,r.density) for r in p.fixed_regions],
                design_regions=[(r.region_id,r.mesh,r.density) for r in p.design_regions],
                length_unit='m', mass_unit='kg', target=target, minimum_outer_wall=.25)
        for target in ({**p.target, 'mode':'BUOYANT_STABILITY'},
                       {**p.target, 'support_points':[[0.,0.,0.]]},
                       {**p.target, 'mode':'STATIC_BALANCE', 'support_points':[[1.,0.,0.]]}):
            with self.assertRaisesRegex(ValueError, 'NOT_SUPPORTED'):
                self.evaluate(target_problem(target))
        static = target_problem({**p.target, 'mode':'STATIC_BALANCE', 'support_points':[[0.,0.,0.]]})
        result = self.evaluate(static)
        self.assertEqual(result.status, 'NOT_OPTIMIZED_CANDIDATE')
        self.assertEqual(result.identity['solver_target_feasibility'], 'UNKNOWN')
        outer = box((-2.,)*3, (2.,)*3)
        two = SpinArtifactProblem(original=outer, fixed_regions=[],
            design_regions=[('a', box((-2.,)*3, (0.,2.,2.)), 1.),
                            ('b', box((0.,-2.,-2.), (2.,)*3), 1.)],
            length_unit='m', mass_unit='kg', target=p.target, minimum_outer_wall=.25)
        with self.assertRaisesRegex(ValueError, 'exactly one design'):
            self.evaluate(two)

    def test_unsupported_derivative_and_scope_cannot_be_sealed(self):
        with self.assertRaisesRegex(ValueError, 'derivative'):
            self.evaluate(problem(), (.01,1.,0.,0.,0.,0.,0.,0.,0.,0.))
        with self.assertRaisesRegex(ValueError, 'DISCRETE_FULL_RANK_QUADRIC'):
            self.bridge.evaluate_and_seal_faster_candidate(problem(), COEFFICIENTS,
                evaluator_settings=spin.SpinEvaluatorConfig(), solver_settings=self.solver)
        with self.assertRaisesRegex(ValueError, 'DUAL_MONOMIAL_RAW_V1'):
            self.bridge.evaluate_and_seal_faster_candidate(problem(), COEFFICIENTS,
                evaluator_settings=self.config, solver_settings=spin.SpinSolverConfig())


if __name__ == '__main__':
    unittest.main()
