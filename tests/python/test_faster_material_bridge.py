"""Actual native Faster geometry to independent complete material sealing."""
from dataclasses import FrozenInstanceError, replace
import importlib
import tempfile
import unittest

import numpy as np

from test_exact_material_region import box, combine
from make_it_stand import spin
from make_it_stand.verification.spin_problem import SpinArtifactProblem


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


class FasterSettingsTests(unittest.TestCase):
    def test_solver_exposes_immutable_actual_export_and_derivative_metadata(self):
        self.assertTrue(hasattr(spin.SpinItFasterResult, 'derivative_metadata'))
        v, f = box((-2.,)*3, (2.,)*3)
        d, g = box((-1.,)*3, (1.,)*3)
        request = spin.SpinItFasterRequest(v.tolist(), f.tolist(), d.tolist(), g.tolist(), mode='YO_YO',
            settings=spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1',
                maximum_restarts=1, maximum_function_evaluations_per_restart=1))
        run = spin.solve_spin_it_faster(request, initial_quadrics=(COEFFICIENTS,),
            evaluator_settings=spin.SpinEvaluatorConfig(derivative_method='DISCRETE_ELLIPSOID', icosphere_subdivisions=1))
        self.assertEqual(run.derivative_metadata['method'], 'discrete-ellipsoid-boolean-transport-v1')
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
        self.assertEqual(e.derivative_method, 'discrete-ellipsoid-boolean-transport-v1')
        self.assertEqual(e.interface_type, 'ELLIPSOID')
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
        with self.assertRaises(FrozenInstanceError):
            config.icosphere_subdivisions = 2
        for kw in ({'derivative_method': 'UNKNOWN'}, {'icosphere_subdivisions': True},
                   {'icosphere_subdivisions': 8}, {'minimum_interface_gradient': 0.},
                   {'retain_fixed_normal_provenance': 1}):
            with self.assertRaises(ValueError):
                replace(config, **kw)
        plane = spin.evaluate_spin_it_faster_quadric(request, (.01,1.,0.,0.,0.,0.,0.,0.,0.,0.),
            evaluator_settings=config)
        self.assertFalse(plane.derivative_defined)
        self.assertEqual(plane.derivative_method, 'discrete-ellipsoid-boolean-transport-v1')


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
        with self.assertRaisesRegex(ValueError, 'DISCRETE_ELLIPSOID'):
            self.bridge.evaluate_and_seal_faster_candidate(problem(), COEFFICIENTS,
                evaluator_settings=spin.SpinEvaluatorConfig(), solver_settings=self.solver)
        with self.assertRaisesRegex(ValueError, 'DUAL_MONOMIAL_RAW_V1'):
            self.bridge.evaluate_and_seal_faster_candidate(problem(), COEFFICIENTS,
                evaluator_settings=self.config, solver_settings=spin.SpinSolverConfig())


if __name__ == '__main__':
    unittest.main()
