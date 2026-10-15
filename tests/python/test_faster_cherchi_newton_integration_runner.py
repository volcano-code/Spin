import importlib.util
import hashlib
import json
from pathlib import Path
import signal
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RUNNER_PATH = ROOT / "experiments/spin_it_reproduction/strict/run_faster_cherchi_newton_integration.py"
sys.path.insert(0, str(RUNNER_PATH.parent))
SPEC = importlib.util.spec_from_file_location(
    "run_faster_cherchi_newton_integration", RUNNER_PATH
)
RUNNER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = RUNNER
SPEC.loader.exec_module(RUNNER)


class _Settings:
    def __init__(self):
        self.icosphere_subdivisions = None
        self.maximum_radius_doublings = None


class _Extension:
    SpinTriangulationSettings = _Settings


class CherchiNewtonRunnerConfigurationTests(unittest.TestCase):
    def test_paper_projective_cli_and_real_source_mapping(self):
        from make_it_stand import spin
        n = spin._native_module()
        parser = RUNNER._argument_parser()
        self.assertFalse(parser.parse_args(['--output', 'run']).paper_projective_hyperboloid)
        self.assertTrue(parser.parse_args([
            '--output', 'run', '--paper-projective-hyperboloid']).paper_projective_hyperboloid)
        settings, evaluator = RUNNER._make_triangulation_configuration(
            n, icosphere_subdivisions=2, maximum_radius_doublings=8,
            paper_projective_hyperboloid=True)
        self.assertTrue(evaluator.paper_projective_hyperboloid)
        self.assertEqual(evaluator.derivative_method, 'DISCRETE_FULL_RANK_QUADRIC')
        q = np.array((-1., .1, -.12, .08, 1., 1.3, -.9, .03, -.04, .02))
        domain = np.array(((-.2, -.2, -.2), (.2, .2, .2)))
        source = n.evaluate_spin_it_faster_source_motion(q, domain, settings)
        self.assertTrue(source['defined'], source['diagnostic'])
        other = n.evaluate_spin_it_faster_source_motion(q, domain,
            spin._native_evaluator_settings(n, evaluator).triangulation)
        for key in ('vertices', 'faces', 'vertex_velocity', 'source_projective_rays'):
            np.testing.assert_array_equal(source[key], other[key])
        branch, features = RUNNER._source_derivative_branch(
            source, paper_projective_hyperboloid=True)
        self.assertGreater(sum(row[0] >= 0 for row in features['projective_rays']), 0)
        self.assertEqual(branch['paper_source_features_sha256'], RUNNER.digest(features))
        self.assertEqual(features['schema'], 'paper-projective-source-features-v1')
        # Same array sizes but a different ideal root must invalidate the stratum.
        changed = dict(source)
        rays = np.asarray(source['source_projective_rays']).copy()
        i = np.flatnonzero(rays[:, 0] >= 0)[0]
        rays[i, 3] = 1 - rays[i, 3]
        changed['source_projective_rays'] = rays
        changed_branch, _ = RUNNER._source_derivative_branch(
            changed, paper_projective_hyperboloid=True)
        self.assertNotEqual(
            RUNNER._derivative_stratum_identity('intersection', 'b'*64, branch),
            RUNNER._derivative_stratum_identity('intersection', 'b'*64, changed_branch))
        with self.assertRaisesRegex(ValueError, 'paper.*source'):
            RUNNER._source_derivative_branch(
                dict(source, paper_projective_hyperboloid=False),
                paper_projective_hyperboloid=True)

    def test_paper_projective_rejects_mutable_normals_bad_flag_and_old_runtime(self):
        from make_it_stand import spin
        n = spin._native_module()
        with self.assertRaisesRegex(ValueError, 'fixed.*normal'):
            RUNNER._make_triangulation_configuration(n, icosphere_subdivisions=2,
                maximum_icosphere_subdivisions=3, paper_projective_hyperboloid=True)
        for bad in (1, 0, None, 'true'):
            with self.subTest(bad=bad), self.assertRaisesRegex(ValueError, 'boolean'):
                RUNNER._make_triangulation_configuration(n, icosphere_subdivisions=2,
                    paper_projective_hyperboloid=bad)
        with self.assertRaisesRegex(RuntimeError, 'unavailable'):
            RUNNER._make_triangulation_configuration(_Extension, icosphere_subdivisions=2,
                paper_projective_hyperboloid=True)
        legacy, config = RUNNER._make_triangulation_configuration(_Extension,
            icosphere_subdivisions=2)
        self.assertFalse(config.paper_projective_hyperboloid)
        self.assertFalse(hasattr(legacy, 'paper_projective_hyperboloid'))

    def test_full_bridge_identity_is_preserved_and_bound_without_relabelling(self):
        helper = getattr(RUNNER, '_bridge_identity_attachment', None)
        self.assertTrue(callable(helper), 'full bridge identity persistence is missing')
        identity = {'operation': 'SOLVE_CHERCHI_EXTERNAL_GEOMETRY',
                    'status': 'RESOURCE_LIMIT_REACHED', 'seal_error': None,
                    'problem_digest': 'locked', 'artifact_scientific_digest': None,
                    'native_run': {'status': 'RESOURCE_LIMIT_REACHED', 'events': []}}
        context = {'bridge_identity': identity,
                   'bridge_scientific_digest': RUNNER.digest(identity),
                   'bridge_status': identity['status'], 'bridge_seal_error': None,
                   'native_result': identity['native_run'],
                   'locked_problem_digest': 'locked', 'expected_artifact_digest': None}
        attachment = helper(context)
        self.assertEqual(attachment['identity'], identity)
        self.assertEqual(attachment['status'], 'BOUND_IDENTITY')
        self.assertEqual(attachment['digest'], context['bridge_scientific_digest'])
        for key, value in [('bridge_scientific_digest', 'wrong'),
                           ('bridge_status', 'CONVERGED'),
                           ('bridge_seal_error', 'different'),
                           ('native_result', {'status': 'CONVERGED'}),
                           ('locked_problem_digest', 'other'),
                           ('expected_artifact_digest', 'other')]:
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'bridge identity'):
                helper({**context, key: value})
        legacy = {k: v for k, v in context.items() if k != 'bridge_identity'}
        self.assertEqual(helper(legacy)['status'], 'LEGACY_IDENTITY_NOT_SAVED')
        self.assertIsNone(helper(legacy)['identity'])
        with self.assertRaisesRegex(ValueError, 'bridge identity'):
            helper({**context, 'bridge_identity': None})
        for key in identity:
            missing = {k:v for k,v in identity.items() if k != key}
            with self.subTest(missing=key), self.assertRaisesRegex(ValueError, 'bridge identity'):
                helper({**context, 'bridge_identity':missing,
                        'bridge_scientific_digest':RUNNER.digest(missing)})
        relabelled = {**identity, 'operation':'SOLVE_FULL_RANK'}
        with self.assertRaisesRegex(ValueError, 'bridge identity'):
            helper({**context, 'bridge_identity':relabelled,
                    'bridge_scientific_digest':RUNNER.digest(relabelled)})

    def test_directional_bridge_binding_uses_versioned_science_not_runtime_cost(self):
        from make_it_stand.verification.faster_material import _result
        result = _result('RESOURCE_LIMIT_REACHED', {
            'operation': 'SOLVE_CHERCHI_EXTERNAL_GEOMETRY', 'problem_digest': 'locked',
            'directional_mapping': {'schema': 'faster-material-supported-directions-v1'},
            'native_run': {'events': [{'error': 'rejected trial'}],
                           'evaluation_timings': [{'wall_seconds': .1}]}})
        identity = result.identity
        context = dict(bridge_identity=identity, bridge_scientific_digest=result.scientific_digest,
            bridge_status=result.status, bridge_seal_error=None, native_result=identity['native_run'],
            locked_problem_digest='locked', expected_artifact_digest=None)
        self.assertEqual(RUNNER._bridge_identity_attachment(context)['identity'], identity)
        identity['native_run']['events'][0]['error'] = 'concealed failure'
        with self.assertRaisesRegex(ValueError, 'bridge identity'):
            RUNNER._bridge_identity_attachment(context)

    def test_exact_motion_material_binds_raw_export_but_integrates_rational_geometry(self):
        from fractions import Fraction
        from test_exact_boolean_coordinate_reader import CoordinateReaderTests
        helper = getattr(RUNNER, '_prepare_exact_motion_material', None)
        self.assertTrue(callable(helper), 'exact moving-material preparation is missing')
        fixture = CoordinateReaderTests()
        rows, points, faces = fixture.fixture(Fraction(1, 2**60))
        # Keep a bounded, authenticated raw OBJ different from the canonical
        # rational display. Do not impose a new exact-float-export contract.
        points = points.copy()
        points[0, 0] = np.nextafter(0., 1.)
        with tempfile.TemporaryDirectory() as directory:
            sidecar = Path(directory) / 'material.exact.tsv'
            fixture.write(sidecar, rows, points)
            original = points.copy()
            snapshot, derivation, orientation = helper(points, faces, sidecar)
            np.testing.assert_array_equal(points, original)
            canonical, oriented = snapshot.arrays()
            self.assertNotEqual(canonical[0, 0], points[0, 0])
            self.assertTrue(snapshot.exact_region().facts['valid'])
            self.assertEqual(orientation['coordinate_authority'], 'RATIONAL_BOUNDARY')
            self.assertEqual(derivation['status'], 'BOUNDED_FLOAT_DERIVATION')
            self.assertEqual(derivation['rational_geometry_identity'], snapshot.identity)
            raw = RUNNER.integrate_rational_polyhedral_raw_moments(
                canonical, oriented, density=2.5, rational_vertices=snapshot.rational_vertices)
            independent = np.asarray(RUNNER._raw_moments(snapshot, 2.5))
            np.testing.assert_array_equal(raw, independent)
            self.assertGreater(raw[0], 0.)
            with self.assertRaises(RUNNER.CherchiDifferentialUnavailable):
                RUNNER.normalize_closed_component_orientation(points, faces)

    def test_exact_motion_material_rejects_intersecting_positive_components(self):
        from fractions import Fraction
        from test_exact_boolean_coordinate_reader import CoordinateReaderTests
        helper = getattr(RUNNER, '_prepare_exact_motion_material', None)
        self.assertTrue(callable(helper), 'exact moving-material preparation is missing')
        rows = tuple(tuple(Fraction(value) for value in point) for point in
                     ((0,0,0), (1,0,0), (0,1,0), (0,0,1)))
        faces = np.asarray(((0,2,1),(0,1,3),(0,3,2),(1,2,3)))
        rows += tuple(tuple(value + Fraction(1,8) for value in point) for point in rows)
        points = np.asarray(rows, dtype=float)
        faces = np.vstack((faces, faces + 4))
        with tempfile.TemporaryDirectory() as directory:
            sidecar = Path(directory) / 'material.exact.tsv'
            CoordinateReaderTests().write(sidecar, rows, points)
            with self.assertRaisesRegex(ValueError, 'intersecting_pairs'):
                helper(points, faces, sidecar)

    def test_source_side_selector_does_not_mistake_cross_product_cancellation_for_degeneracy(self):
        from fractions import Fraction as F
        q = np.asarray((-1.4960632120915978e-35, 3.99789656144784e-17,
            -1.5664854529958856e-17, -2.0735808749445873e-19,
            .39824375457455563, .3982438519330574, -.8263193972913151,
            3.172499238676896e-18, 2.966716193512342e-17, 3.9190872484896026e-19))
        design = np.asarray(((-.25,-.25,-.25),(-.25,-.25,.25),(.25,-.25,.25),
            (.25,-.25,-.25),(-.25,.25,-.25),(-.25,.25,.25),(.25,.25,.25),(.25,.25,-.25)))
        native = RUNNER.spin._native_module()
        settings, _ = RUNNER._make_triangulation_configuration(
            native, icosphere_subdivisions=5, maximum_radius_doublings=8)
        source = native.evaluate_spin_it_faster_source_motion(q, design, settings)
        self.assertTrue(source['defined'])
        vertices = [[F.from_float(float(x)) for x in v] for v in source['vertices']]
        rq = [F.from_float(float(x)) for x in q]
        A = ((rq[4],rq[9],rq[8]), (rq[9],rq[5],rq[7]), (rq[8],rq[7],rq[6]))
        positive = 0
        for indices, cap in zip(source['faces'], source['face_is_cap']):
            if cap: continue
            a,b,c = [vertices[i] for i in indices]
            u,v = [b[i]-a[i] for i in range(3)], [c[i]-a[i] for i in range(3)]
            area = [u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
            center = [(a[i]+b[i]+c[i])/3 for i in range(3)]
            gradient = [rq[i+1]+2*sum(A[i][j]*center[j] for j in range(3)) for i in range(3)]
            dot = sum(area[i]*gradient[i] for i in range(3))
            area2,gradient2 = sum(x*x for x in area),sum(x*x for x in gradient)
            self.assertGreater(area2, 0)
            self.assertGreater(gradient2, 0)
            self.assertGreater(dot, 0)
            self.assertGreater(dot*dot, F.from_float(1e-10)**2*area2*gradient2)
            positive += 1
        self.assertGreater(positive, 0)
        self.assertTrue(source['side_selection_defined'], source['side_selection_diagnostic'])
        self.assertFalse(source['interior_is_positive'])

    def test_top_objective_weight_reaches_actual_cherchi_evaluation_and_input_evidence(self):
        from make_it_stand.verification.faster_material import _transform_moments
        fixture = ROOT / 'experiments/spin_it_reproduction/faster/fixtures/cherchi-centered-axis-fixed-shell-v1'
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            shutil.copytree(fixture, root/'input')
            manifest_path = root/'input/problem.json'
            manifest = json.loads(manifest_path.read_text())
            source = (fixture/manifest['outer_mesh']).resolve()
            shutil.copy2(source, root/'input/outer.obj')
            manifest['outer_mesh'] = 'outer.obj'
            for kind in ('fixed_regions', 'design_regions'):
                for index, region in enumerate(manifest[kind]):
                    name = f'{kind}-{index}.obj'
                    shutil.copy2((fixture/region['mesh']).resolve(), root/'input'/name)
                    region['mesh'] = name
            manifest['target'] = {'mode': 'SPINNING_TOP', 'axis': [0., 0., 1.],
                                  'origin': [0., 0., -2.], 'support_points': []}
            # Frozen method fixture, not answer-independent paper initialization:
            # a radius-1/8 sphere centered at target-frame z=1/2.
            q = np.asarray((.25-.125**2, 0., 0., -1., 1., 1., 1., 0., 0., 0.))
            q /= np.linalg.norm(q)
            manifest['preflight_quadric_target_normalized'] = q.tolist()
            manifest_path.write_text(json.dumps(manifest))
            build = Path(RUNNER.spin._native_module().__file__).resolve().parent.parent
            probes = sorted((RUNNER.ROOT/'build').glob('cherchi-reference-singlethread-*/cherchi_probe'))
            self.assertEqual(len(probes), 1)
            RUNNER.run(root/'run', build=build, probe=probes[0],
                       problem_manifest=manifest_path, center_of_mass_weight=1.25,
                       maximum_function_evaluations=4, maximum_restarts=1,
                       icosphere_subdivisions=4, defer_finalization=True)
            inputs = json.loads((root/'run/inputs.json').read_text())
            self.assertEqual(inputs['center_of_mass_weight'], 1.25)
            payload = RUNNER._load_terminal_audit_checkpoint(root/'run/terminal-audit-checkpoint.json')
            bridge = RUNNER._bridge_identity_attachment(payload['run_finalization'])
            self.assertEqual(bridge['status'], 'BOUND_IDENTITY')
            self.assertEqual(bridge['identity']['operation'], 'SOLVE_CHERCHI_EXTERNAL_GEOMETRY')
            evaluated = [row for row in payload['run_finalization']['native_result']['events']
                         if row['type'] == 'evaluation' and row.get('objective') is not None]
            self.assertTrue(evaluated, 'the actual weighted objective must be evaluated')
            problem, _, _ = RUNNER._load_problem_manifest(manifest_path)
            _, request, _ = RUNNER._prepare(
                problem, RUNNER.spin.SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC'),
                RUNNER.spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1'))
            axes, origin = np.asarray(request.target_axes), np.asarray(request.target_origin)
            fixed = np.asarray(_transform_moments(request.fixed_shell_moments,
                                                  axes.T, -axes.T@origin))
            for row in evaluated:
                evaluation = next(e for e in payload['evaluations']
                                  if e['quadric_hash'] == row['quadric_hash'])
                s = np.asarray(evaluation['moment_map']['raw_moments']) + fixed
                mass, x, y, z, xx, yy, zz, yz, xz, xy = s
                i11 = yy+zz-(y*y+z*z)/mass
                i22 = xx+zz-(x*x+z*z)/mass
                i33 = xx+yy-(x*x+y*y)/mass
                i12 = -xy+x*y/mass
                expected = (i11*i11+i22*i22+2*i12*i12)/(i33*i33) + 1.25*z/mass
                self.assertAlmostEqual(row['objective'], expected, places=10)

    def test_cli_preserves_explicit_center_of_mass_objective_weight(self):
        args = RUNNER._argument_parser().parse_args((
            '--output', '/unused/run', '--center-of-mass-weight', '1.25',
        ))
        self.assertEqual(args.center_of_mass_weight, 1.25)

    def test_invalid_objective_weight_rejected_before_feasibility_or_output(self):
        manifest = ROOT / 'experiments/spin_it_reproduction/faster/fixtures/cherchi-rotated-box-fixed-shell-v2/problem.json'
        with tempfile.TemporaryDirectory() as temporary:
            for weight in (-1., float('nan'), float('inf'), True):
                with self.subTest(weight=weight):
                    output = Path(temporary) / 'run'
                    with self.assertRaisesRegex(ValueError, 'center.of.mass.weight'):
                        RUNNER.run(output, build=Path(temporary)/'missing-build',
                                   probe=Path(temporary)/'missing-probe',
                                   problem_manifest=manifest,
                                   center_of_mass_weight=weight)
                    self.assertFalse(output.exists())

    def test_terminal_complement_uses_exact_boundary_not_collapsed_display(self):
        from fractions import Fraction
        from make_it_stand.verification.rational_text import fraction_to_rational_text
        helper = getattr(RUNNER, '_terminal_exact_material', None)
        self.assertTrue(callable(helper), 'terminal material must use its authenticated exact boundary')
        rows = ((Fraction(0), Fraction(0), Fraction(1)),
                (Fraction(1, 2), Fraction(0), Fraction(1) + Fraction(1, 2**60)),
                (Fraction(1), Fraction(0), Fraction(1)),
                (Fraction(0), Fraction(1), Fraction(0)))
        vertices = np.asarray([[float(x) for x in row] for row in rows])
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)), dtype=np.int32)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'material.exact.tsv'
            header = ('vertex_index\tstatus\texact_x\texact_y\texact_z\t'
                      'memory_x\tmemory_y\tmemory_z\tobj_x\tobj_y\tobj_z\n')
            records = []
            for index, row in enumerate(rows):
                exact = [fraction_to_rational_text(x) for x in row]
                display = [format(float(x), '.17g') for x in row]
                records.append('\t'.join([str(index), 'MATCHED', *exact, *display, *display])+'\n')
            path.write_text(header+''.join(records), encoding='ascii')
            snapshot, derivation, orientation, moments = helper(vertices, faces, path, 2.5)
            self.assertTrue(snapshot.exact_region().facts['valid'])
            np.testing.assert_array_equal(snapshot.arrays()[1], faces)
            self.assertEqual(orientation['authority'], 'EXACT_RATIONAL_BOUNDARY')
            self.assertEqual(orientation['face_order_changed'], False)
            self.assertEqual(orientation['floating_orientation_check']['status'], 'REJECTED')
            self.assertIn('zero-area face', orientation['floating_orientation_check']['diagnostic'])
            self.assertEqual(derivation['status'], 'EXACT_FLOAT_MATCH')
            self.assertAlmostEqual(moments[0]/float(Fraction(1, 2**60)*Fraction(5, 2)/6), 1.0)
            with self.assertRaisesRegex(ValueError, 'material density'):
                helper(vertices, faces, path, 0.)
            # Correctly bound coordinates must not excuse a genuinely exact
            # zero-area surface: collapse the rational sliver as well.
            records[1] = records[1].replace(fraction_to_rational_text(rows[1][2]), '1/1', 1)
            path.write_text(header+''.join(records), encoding='ascii')
            with self.assertRaisesRegex(ValueError, 'invalid material boundary'):
                helper(vertices, faces, path, 2.5)

    def test_terminal_audit_uses_bound_exact_material_without_approving_collapsed_display(self):
        from fractions import Fraction as Q
        from make_it_stand.verification.rational_text import fraction_to_rational_text
        rows = ((Q(0), Q(0), Q(1)),
                (Q(1, 2), Q(0), Q(1) + Q(1, 2**60)),
                (Q(1), Q(0), Q(1)), (Q(0), Q(1), Q(0)))
        vertices = np.asarray([[float(x) for x in row] for row in rows])
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)))
        mass = Q(5, 2) * Q(1, 2**60) / 6
        sums = [sum(row[axis] for row in rows) for axis in range(3)]
        # Analytic affine-simplex averages, independently of either integrator.
        moments = [mass, *(mass*s/4 for s in sums)]
        moments += [mass*(sums[i]**2 + sum(row[i]**2 for row in rows))/20 for i in range(3)]
        moments += [mass*(sums[i]*sums[j] + sum(row[i]*row[j] for row in rows))/20
                    for i, j in ((1, 2), (0, 2), (0, 1))]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'terminal.obj'
            RUNNER._write_obj(path, vertices, faces)
            sidecar = Path(str(path)+'.exact.tsv')
            text = ('vertex_index\tstatus\texact_x\texact_y\texact_z\t'
                    'memory_x\tmemory_y\tmemory_z\tobj_x\tobj_y\tobj_z\n')
            for index, row in enumerate(rows):
                exact = [fraction_to_rational_text(x) for x in row]
                display = [format(float(x), '.17g') for x in row]
                text += '\t'.join([str(index), 'MATCHED', *exact, *display, *display])+'\n'
            sidecar.write_text(text, encoding='ascii')
            snapshot, derivation, _, _ = RUNNER._terminal_exact_material(vertices, faces, sidecar, 2.5)
            q = np.asarray((1., 0., 0., 0., 0., 0., 0., 0., 0., 0.))
            qhash = RUNNER.spin._native_module().spin_it_faster_quadric_hash(q)
            mesh_digest = hashlib.sha256(path.read_bytes()).hexdigest()
            evaluation = {'solver_evaluation_index': 1, 'quadric_coefficients': q.tolist(),
                'quadric_hash': qhash, 'geometry_hash': mesh_digest, 'topology_hash': 'a'*64,
                'oriented_output_path': str(path), 'output_path': str(path),
                'exact_float_derivation_audit': derivation,
                'exact_geometry_identity': RUNNER.digest({'rational_geometry': snapshot.identity,
                                                         'float_derivation_audit': derivation}),
                'moment_map': {'raw_moments': [float(x) for x in moments]},
                'independent_moment_audit': {'status': 'PASSED', 'density': 2.5, 'root_side': 4.}}
            result = {'solver_status': 'CONVERGED',
                'evaluation_map_binding_audit': {'status': 'CONSISTENT'},
                'evaluations': [evaluation], 'design_domain_faces': faces.tolist(),
                'design_domain_vertices': [(0.,0.,0.),(4.,0.,0.),(0.,4.,0.),(0.,0.,4.)],
                'trace': [{'function_evaluations': 1, 'quadric_hash': qhash,
                          'material_geometry_hash': mesh_digest, 'material_topology_hash': 'a'*64}]}
            try:
                audit = RUNNER._audit_terminal_geometry_binding(result,
                    quadric_hash=RUNNER.spin._native_module().spin_it_faster_quadric_hash)
            except ValueError as error:
                self.fail('authenticated rational material was rejected through its display: '+str(error))
            self.assertEqual(audit['status'], 'PASSED')
            self.assertEqual(audit['material_authority'], 'EXACT_RATIONAL_BOUNDARY')
            self.assertEqual(audit['float_validation_status'], 'REJECTED')
            self.assertFalse(audit['floating_boundary_audit']['boundary_valid'])
            self.assertEqual(audit['material_domain_status'], 'VERIFIED_SUBSET')
            self.assertTrue(audit['exact_boundary_audit']['region_facts']['valid'])
            self.assertEqual(audit['independent_moment_audit']['status'], 'PASSED')
            # Every displayed vertex lies in this unit cube, but the exact
            # sliver rises above z=1. Only the rational subset query sees it.
            outside = {**result,
                'design_domain_vertices': [(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
                                          (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.)],
                'design_domain_faces': [(0,2,1),(0,3,2),(4,5,6),(4,6,7),
                                        (0,1,5),(0,5,4),(3,7,6),(3,6,2),
                                        (0,4,7),(0,7,3),(1,2,6),(1,6,5)]}
            outside_audit = RUNNER._audit_terminal_geometry_binding(outside,
                quadric_hash=RUNNER.spin._native_module().spin_it_faster_quadric_hash)
            self.assertEqual(outside_audit['status'], 'FAILED')
            self.assertEqual(outside_audit['material_domain_status'], 'OUTSIDE_DESIGN_DOMAIN')
            dropped = {**result, 'evaluations': [dict(evaluation)]}
            dropped['evaluations'][0].pop('exact_geometry_identity')
            dropped['evaluations'][0].pop('exact_float_derivation_audit')
            with self.assertRaisesRegex(ValueError, 'sidecar.*identity'):
                RUNNER._audit_terminal_geometry_binding(dropped,
                    quadric_hash=RUNNER.spin._native_module().spin_it_faster_quadric_hash)
            # Preserve the already approved bounded-conversion contract: an
            # OBJ need not equal the canonical rounded rational coordinates.
            bounded_path = Path(temporary)/'bounded.obj'
            bounded_vertices = vertices.copy()
            bounded_vertices[0,0] = np.nextafter(0.,1.)
            RUNNER._write_obj(bounded_path, bounded_vertices, faces)
            bounded_sidecar = Path(str(bounded_path)+'.exact.tsv')
            bounded_text = text.replace('0\tMATCHED\t0/1\t0/1\t1/1\t0\t0\t1\t0\t0\t1\n',
                '0\tAPPROXIMATION_DISAGREES\t0/1\t0/1\t1/1\t'+
                '\t'.join([format(bounded_vertices[0,0],'.17g'),'0','1']*2)+'\n', 1)
            self.assertNotEqual(bounded_text, text)
            bounded_sidecar.write_text(bounded_text, encoding='ascii')
            bounded_snapshot, bounded_derivation, _, _ = RUNNER._terminal_exact_material(
                bounded_vertices, faces, bounded_sidecar, 2.5)
            self.assertEqual(bounded_derivation['status'],'BOUNDED_FLOAT_DERIVATION')
            bounded_hash = hashlib.sha256(bounded_path.read_bytes()).hexdigest()
            bounded_evaluation = {**evaluation, 'geometry_hash':bounded_hash,
                'output_path':str(bounded_path), 'oriented_output_path':str(bounded_path),
                'exact_float_derivation_audit':bounded_derivation,
                'exact_geometry_identity':RUNNER.digest({'rational_geometry':bounded_snapshot.identity,
                                                       'float_derivation_audit':bounded_derivation})}
            bounded_result = {**result, 'evaluations':[bounded_evaluation],
                'trace':[{**result['trace'][0], 'material_geometry_hash':bounded_hash}]}
            try:
                bounded_audit = RUNNER._audit_terminal_geometry_binding(bounded_result,
                    quadric_hash=RUNNER.spin._native_module().spin_it_faster_quadric_hash)
            except ValueError as error:
                self.fail('authenticated bounded conversion was rejected as exact inequality: '+str(error))
            self.assertEqual(bounded_audit['status'],'PASSED')
            self.assertEqual(bounded_audit['float_validation_status'],'REJECTED')
            # Rationals differing below binary64 resolution must still be
            # covered by the saved sidecar digest; no silent float fallback.
            sidecar.write_text(text.replace(fraction_to_rational_text(rows[1][2]),
                fraction_to_rational_text(Q(1)+Q(1,2**61)), 1), encoding='ascii')
            with self.assertRaisesRegex(ValueError, 'sidecar.*identity'):
                RUNNER._audit_terminal_geometry_binding(result,
                    quadric_hash=RUNNER.spin._native_module().spin_it_faster_quadric_hash)

    def test_full_finalization_rejects_audit_only_checkpoint_without_run_context(self):
        self.assertTrue(callable(getattr(RUNNER, "resume_run_finalization", None)),
                        "completed evidence needs a checked full-run finalizer")
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            RUNNER._save_terminal_audit_checkpoint(
                output, {"solver_status": "CONVERGED"}, dependencies=[])
            with self.assertRaisesRegex(ValueError, "full-run finalization context"):
                RUNNER.resume_run_finalization(output)
            self.assertFalse((output / "result.json").exists())
            self.assertFalse((output / "evidence.json").exists())

    def test_terminal_audit_checkpoint_rejects_changed_inputs_and_payload(self):
        self.assertTrue(callable(getattr(RUNNER, "_save_terminal_audit_checkpoint", None)),
                        "terminal audit must have a durable checkpoint writer")
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            dependency = output / "inputs.json"
            dependency.write_text('{"density":2}\n')
            payload = {"solver_status": "CONVERGED", "evaluations": [], "trace": []}
            checkpoint = RUNNER._save_terminal_audit_checkpoint(
                output, payload, dependencies=[dependency])
            self.assertEqual(RUNNER._load_terminal_audit_checkpoint(checkpoint), payload)
            dependency.write_text('{"density":3}\n')
            with self.assertRaisesRegex(ValueError, "dependency changed"):
                RUNNER._load_terminal_audit_checkpoint(checkpoint)
            dependency.write_text('{"density":2}\n')
            record = json.loads(checkpoint.read_text())
            record["payload"]["solver_status"] = "CONVERGED_FAKE"
            checkpoint.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, "checkpoint digest"):
                RUNNER._load_terminal_audit_checkpoint(checkpoint)

    def test_terminal_audit_checkpoint_rejects_different_loaded_checker(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            checkpoint = RUNNER._save_terminal_audit_checkpoint(
                output, {"solver_status": "CONVERGED"}, dependencies=[])
            record = json.loads(checkpoint.read_text())
            record["checker_identity"] = []
            record.pop("checkpoint_digest")
            record["checkpoint_digest"] = RUNNER.digest(record)
            checkpoint.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, "checker identity"):
                RUNNER._load_terminal_audit_checkpoint(checkpoint)

    def test_exact_newton_state_replay_preserves_the_frozen_solver_state(self):
        coefficients = np.asarray((0.31, -0.28, 0.17, 0.11, -0.43,
                                   -0.36, -0.52, 0.22, -0.13, 0.19))
        coefficients /= np.linalg.norm(coefficients)
        starts = RUNNER._solver_initial_quadrics(
            "exact_newton_state_replay", np.eye(10)[0],
            replay_coefficients=coefficients,
        )
        np.testing.assert_array_equal(np.asarray(starts[0]), coefficients)

    def test_certified_infeasible_manifest_stops_before_loading_build_or_starting_solver(self):
        manifest = ROOT / "experiments/spin_it_reproduction/faster/fixtures/cherchi-rotated-box-fixed-shell-v2/problem.json"
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "run"
            result = RUNNER.run(
                output, build=Path(temporary) / "missing-build",
                probe=Path(temporary) / "missing-probe", problem_manifest=manifest,
                maximum_function_evaluations=60, maximum_restarts=8,
                maximum_radius_doublings=8, deterministic_seed=1,
            )

            self.assertEqual(result["status"], "INFEASIBLE_CERTIFIED")
            self.assertEqual(result["solver_status"], "NOT_RUN_PREFLIGHT_CERTIFIED")
            self.assertEqual(result["function_evaluations"], 0)
            certificate_path = output / "preflight-feasibility-report.json"
            self.assertTrue(certificate_path.is_file())
            certificate = json.loads(certificate_path.read_text())
            self.assertEqual(certificate["status"], "INFEASIBLE_CERTIFIED")
            self.assertFalse((output / "evaluation-events.jsonl").exists())
            self.assertFalse((output / "sealed-material-artifact").exists())

    def test_run_diagnostic_stage_times_and_returns_the_wrapped_operation(self):
        with tempfile.TemporaryDirectory() as temporary:
            journal = Path(temporary) / "stage-timings.jsonl"
            RUNNER._initialize_stage_timing_journal(journal)

            result = RUNNER._run_diagnostic_stage(
                journal, -1, "newton_bridge", lambda: {"status": "RESOURCE_LIMIT_REACHED"},
            )

            self.assertEqual(result, {"status": "RESOURCE_LIMIT_REACHED"})
            records = RUNNER._read_stage_timing_journal(journal)
            self.assertEqual([record["event"] for record in records], [
                "stage_started", "stage_finished",
            ])
            self.assertEqual(records[0]["stage"], "newton_bridge")
            self.assertEqual(records[1]["outcome"], "succeeded")

    def test_stage_timing_journal_records_completed_and_failed_stage_durations(self):
        with tempfile.TemporaryDirectory() as temporary:
            journal = Path(temporary) / "stage-timings.jsonl"
            RUNNER._initialize_stage_timing_journal(journal)

            with RUNNER._diagnostic_stage(journal, 3, "source_motion"):
                pass
            with self.assertRaisesRegex(RuntimeError, "boolean fault"):
                with RUNNER._diagnostic_stage(journal, 3, "cherchi_boolean"):
                    raise RuntimeError("boolean fault")

            records = RUNNER._read_stage_timing_journal(journal)
            self.assertEqual(records[0]["event"], "stage_started")
            self.assertEqual(records[1]["event"], "stage_finished")
            self.assertEqual(records[1]["outcome"], "succeeded")
            self.assertGreaterEqual(records[1]["elapsed_ns"], 0)
            self.assertEqual(records[2]["stage"], "cherchi_boolean")
            self.assertEqual(records[3]["outcome"], "failed")
            self.assertIn("RuntimeError: boolean fault", records[3]["diagnostic"])

    def test_interruption_receipt_lists_unfinished_diagnostic_stage(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            journal = output / "stage-timings.jsonl"
            RUNNER._initialize_stage_timing_journal(journal)
            RUNNER._append_stage_timing_record(journal, {
                "event": "stage_started", "stage_id": "7:cherchi_boolean",
                "evaluation_index": 7, "stage": "cherchi_boolean",
            })
            solver_journal = output / "solver-events.jsonl"
            RUNNER._initialize_solver_event_journal(
                solver_journal, inputs_sha256="a" * 64,
            )
            RUNNER._append_solver_journal_event(
                solver_journal, "geometry_evaluation_started", {
                    "evaluation_index": 7, "quadric_coefficients": [0.0] * 10,
                },
            )

            receipt = RUNNER._write_interruption_receipt(output, signal.SIGTERM)

            self.assertEqual(receipt["open_stages"], ["7:cherchi_boolean"])
            self.assertEqual(receipt["open_solver_evaluation_indices"], [7])
            self.assertRegex(receipt["stage_timing_journal_sha256"], r"^[0-9a-f]{64}$")
            self.assertRegex(receipt["solver_event_journal_sha256"], r"^[0-9a-f]{64}$")

    def test_streamed_native_events_are_ordered_and_audited_against_terminal_events(self):
        with tempfile.TemporaryDirectory() as temporary:
            journal = Path(temporary) / "solver-events.jsonl"
            RUNNER._initialize_solver_event_journal(journal, inputs_sha256="b" * 64)
            events = [
                {"type": "attempt_start", "restart": 0, "outcome": "started"},
                {"type": "evaluation", "restart": 0, "outcome": "accepted"},
                {"type": "attempt_terminal", "restart": 0, "outcome": "terminated"},
            ]
            RUNNER._append_solver_journal_event(journal, "solver_event", events[0])
            RUNNER._append_solver_journal_event(journal, "geometry_evaluation_started", {
                "evaluation_index": 0, "quadric_coefficients": [0.0] * 10,
            })
            RUNNER._append_solver_journal_event(journal, "geometry_evaluation_finished", {
                "evaluation_index": 0, "outcome": "returned",
            })
            for event in events[1:]:
                RUNNER._append_solver_journal_event(journal, "solver_event", event)

            RUNNER._audit_streamed_solver_events(journal, events)
            with self.assertRaisesRegex(ValueError, "differ from finalized native events"):
                RUNNER._audit_streamed_solver_events(journal, events[:-1])

    def test_legacy_streamed_event_encoding_remains_readable_without_rewriting(self):
        with tempfile.TemporaryDirectory() as temporary:
            journal = Path(temporary) / "solver-events.jsonl"
            original = ('{"schema":"faster-streamed-solver-events-v1"}\n'
                        '{"record_type":"solver_event","record":{"residual":"+inf"}}\n')
            journal.write_text(original, encoding="utf-8")
            audit = RUNNER._audit_streamed_solver_events(journal, [{"residual": float("inf")}])
            self.assertEqual(audit["status"], "PASSED")
            self.assertEqual(journal.read_text(encoding="utf-8"), original)

    def test_typed_streamed_nonfinite_diagnostics_do_not_alias_strings(self):
        with tempfile.TemporaryDirectory() as temporary:
            journal = Path(temporary) / "solver-events.jsonl"
            RUNNER._initialize_solver_event_journal(journal, inputs_sha256="b" * 64)
            event = {"values": np.asarray([np.inf, -np.inf, np.nan]), "detail": "+inf"}
            RUNNER._append_solver_journal_event(journal, "solver_event", event)
            saved = RUNNER._read_streamed_solver_events(journal)
            self.assertEqual(saved, [{"values": [{"$nonfinite_float": "+inf"},
                                                  {"$nonfinite_float": "-inf"},
                                                  {"$nonfinite_float": "nan"}],
                                      "detail": "+inf"}])
            RUNNER._audit_streamed_solver_events(journal, [event])
            with self.assertRaisesRegex(ValueError, "differ from finalized native events"):
                RUNNER._audit_streamed_solver_events(journal, [{
                    "values": ["+inf", "-inf", "nan"], "detail": "+inf"}])

    def test_streamed_rank_failure_retains_actual_newton_system_for_replay(self):
        from make_it_stand.spin import (
            SpinItFasterRequest,
            SpinSolverConfig,
            solve_spin_it_faster_with_external_moment_evaluator,
        )

        vertices = np.asarray((
            (-3., -3., -3.), (3., -3., -3.), (3., 3., -3.), (-3., 3., -3.),
            (-3., -3., 3.), (3., -3., 3.), (3., 3., 3.), (-3., 3., 3.),
        ))
        faces = tuple(tuple(int(index) for index in face) for face in (
            (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 2, 6),
            (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
        ))
        design_vertices = vertices / 3.0
        request = SpinItFasterRequest(
            outer_vertices=vertices,
            outer_faces=faces,
            design_vertices=design_vertices,
            design_faces=faces,
            settings=SpinSolverConfig(
                equation_profile="DUAL_MONOMIAL_RAW_V1",
                maximum_restarts=1,
                maximum_function_evaluations_per_restart=1,
                capture_evaluation_payload=True,
            ),
        )
        initial = ((0.31, -0.28, 0.17, 0.11, -0.43,
                    -0.36, -0.52, 0.22, -0.13, 0.19),)

        def evaluator(_coefficients):
            return {
                "raw_moments": np.arange(1.0, 11.0),
                "moment_jacobian": np.diag((1.0,) + (0.0,) * 9),
                "derivative_defined": True,
                "geometry_hash": "persisted-rank-failure-geometry-v1",
                "topology_hash": "persisted-rank-failure-topology-v1",
                "evaluation_map_identity": "e" * 64,
                "boolean_backend_identity": "persisted-rank-failure-backend-v1",
                "derivative_identity": "persisted-rank-failure-derivative-v1",
            }

        with tempfile.TemporaryDirectory() as temporary:
            journal = Path(temporary) / "solver-events.jsonl"
            RUNNER._initialize_solver_event_journal(
                journal, inputs_sha256="f" * 64,
            )

            def persist(event):
                RUNNER._append_solver_journal_event(journal, "solver_event", event)

            result = solve_spin_it_faster_with_external_moment_evaluator(
                request, evaluator, initial_quadrics=initial, event_sink=persist,
            )
            streamed = RUNNER._read_streamed_solver_events(journal)
            RUNNER._audit_streamed_solver_events(journal, result["events"])
            # The real material bridge persists the native result as strict JSON.
            # Rank failure diagnostics contain infinities; that round trip must
            # preserve event identity instead of discarding the failed attempt.
            from make_it_stand.verification.faster_material import _json_value
            from make_it_stand.verification.spin_common import canonical, parse
            bridged_events = parse(canonical(_json_value(result["events"])))
            RUNNER._audit_streamed_solver_events(journal, bridged_events)
            changed = parse(canonical(bridged_events))
            failed = next(event for event in changed
                          if event.get("reason") == "bordered_system_rank_deficient")
            failed["newton_system_payload"]["linear_backward_error"] = 0.0
            with self.assertRaisesRegex(ValueError, "differ from finalized native events"):
                RUNNER._audit_streamed_solver_events(journal, changed)

        terminal = next(event for event in streamed
                        if event.get("reason") == "bordered_system_rank_deficient")
        payload = terminal.get("newton_system_payload")
        self.assertIsInstance(payload, dict)
        self.assertEqual(payload["schema"], "spin-newton-system-payload-v2")
        self.assertEqual(len(payload["matrix"]), len(payload["right_hand_side"]))
        self.assertLess(payload["equilibrated_rank"], len(payload["matrix"]), payload)
        self.assertIsNotNone(payload["multiplier_mu"])
        finalized = next(event for event in result["events"]
                         if event.get("reason") == "bordered_system_rank_deficient")
        finalized_payload = finalized["newton_system_payload"]
        np.testing.assert_array_equal(
            payload["matrix"], np.asarray(finalized_payload["matrix"]),
        )
        np.testing.assert_array_equal(
            payload["right_hand_side"],
            np.asarray(finalized_payload["right_hand_side"]),
        )
        self.assertEqual(payload["multiplier_mu"], finalized_payload["multiplier_mu"])

    def test_interruption_receipt_preserves_partial_ledger_without_forging_result(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            inputs = output / "inputs.json"
            inputs.write_text('{"schema":"run-inputs-v1"}\n', encoding="utf-8")
            journal = output / "evaluation-events.jsonl"
            journal.write_text(
                '{"schema":"event-ledger-v1"}\n'
                '{"record_type":"evaluation","record":{"evaluation_index":4}}\n'
                '{truncated',
                encoding="utf-8",
            )

            receipt = RUNNER._write_interruption_receipt(output, signal.SIGTERM)

            self.assertEqual(receipt["schema"], "faster-cherchi-interrupted-run-v1")
            self.assertEqual(receipt["status"], "INTERRUPTED")
            self.assertEqual(receipt["termination_signal"], "SIGTERM")
            self.assertEqual(receipt["partial_evaluation_count"], 1)
            self.assertEqual(receipt["last_evaluation_index"], 4)
            self.assertEqual(receipt["solver_terminal_status"], "NOT_AVAILABLE")
            self.assertEqual(receipt["freshness"], "NOT_GENERATED_INCOMPLETE_RUN")
            self.assertFalse(receipt["result_json_present"])
            self.assertEqual(json.loads((output / "interrupted-run.json").read_text()), receipt)

    def test_main_records_sigterm_and_exits_nonzero_after_interrupted_run(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            output.mkdir(exist_ok=True)
            (output / "inputs.json").write_text('{"schema":"run-inputs-v1"}\n', encoding="utf-8")
            argv = [
                str(RUNNER_PATH), "--output", str(output), "--build", str(output),
                "--cherchi-probe", str(output / "probe"),
            ]
            with patch.object(sys, "argv", argv), patch.object(
                RUNNER, "run",
                side_effect=lambda *_args, **_kwargs: signal.raise_signal(signal.SIGTERM),
            ):
                with self.assertRaises(SystemExit) as raised:
                    RUNNER.main()

            self.assertEqual(raised.exception.code, 128 + signal.SIGTERM)
            receipt = json.loads((output / "interrupted-run.json").read_text())
            self.assertEqual(receipt["termination_signal"], "SIGTERM")
            self.assertEqual(receipt["solver_terminal_status"], "NOT_AVAILABLE")

    def test_mesh_target_frame_transform_applies_origin_and_rotation(self):
        request = type("TargetFrame", (), {
            "target_origin": (1., 1., 1.),
            "target_axes": ((0., 1., 0.), (-1., 0., 0.), (0., 0., 1.)),
        })()

        result = RUNNER._mesh_vertices_in_target_frame(((2., 3., 4.),), request)

        np.testing.assert_array_equal(result, ((-2., 1., 3.),))

    def test_target_frame_round_trip_and_normalized_region_context(self):
        request = type("TargetFrame", (), {
            "target_origin": (0.25, -0.5, 0.75),
            "target_axes": ((0., 1., 0.), (-1., 0., 0.), (0., 0., 1.)),
        })()
        normalization = {"origin": (2., -1., 4.), "L": 2.0}
        model_vertices = np.asarray(((4., 1., 8.), (2., 3., 6.)))

        target_vertices = RUNNER._normalized_mesh_vertices_in_target_frame(
            model_vertices, normalization, request,
        )
        recovered_model = RUNNER._mesh_vertices_from_target_frame(
            target_vertices, request,
        )
        expected_normalized_model = (
            model_vertices - np.asarray(normalization["origin"])
        ) / normalization["L"]

        np.testing.assert_allclose(recovered_model, expected_normalized_model,
                                   rtol=0., atol=1e-15)

    def test_problem_manifest_loads_real_fixed_and_design_regions(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            outer_vertices = np.asarray((
                (-2., -2., -2.), (2., -2., -2.), (2., 2., -2.), (-2., 2., -2.),
                (-2., -2., 2.), (2., -2., 2.), (2., 2., 2.), (-2., 2., 2.),
            ))
            outer_faces = np.asarray((
                (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
                (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
                (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
            ), dtype=np.int32)
            design_vertices = outer_vertices / 2.
            fixed_vertices = np.vstack((outer_vertices, design_vertices))
            fixed_faces = np.vstack((outer_faces, outer_faces[:, ::-1] + 8))
            RUNNER._write_obj(root / "outer.obj", outer_vertices, outer_faces)
            RUNNER._write_obj(root / "design.obj", design_vertices, outer_faces)
            RUNNER._write_obj(root / "shell.obj", fixed_vertices, fixed_faces)
            manifest = {
                "schema": "faster-cherchi-problem-manifest-v1",
                "outer_mesh": "outer.obj",
                "fixed_regions": [{
                    "id": "shell", "mesh": "shell.obj", "density": 0.75,
                }],
                "design_regions": [{
                    "id": "core", "mesh": "design.obj", "density": 1.0,
                }],
                "length_unit": "m",
                "mass_unit": "kg",
                "target": {
                    "mode": "YO_YO", "axis": [0., 0., 1.],
                    "origin": [0., 0., 0.], "support_points": [],
                },
                "minimum_outer_wall": 0.1,
                "preflight_quadric_target_normalized": [
                    0., 0., 0., 0., 0.5773502691896258,
                    0.5773502691896258, 0.5773502691896258, 0., 0., 0.,
                ],
            }
            manifest_path = root / "problem.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            problem, assets, preflight = RUNNER._load_problem_manifest(manifest_path)

        self.assertEqual(len(problem.fixed_regions), 1)
        self.assertEqual(len(problem.design_regions), 1)
        self.assertEqual(problem.fixed_regions[0].density, 0.75)
        self.assertEqual(problem.design_regions[0].density, 1.0)
        self.assertEqual(problem.target["mode"], "YO_YO")
        self.assertEqual(problem.minimum_outer_wall, 0.1)
        np.testing.assert_allclose(preflight[4:7], [1. / np.sqrt(3.)] * 3,
                                   rtol=0., atol=1e-15)
        self.assertEqual({path.name for path in assets}, {
            "problem.json", "outer.obj", "design.obj", "shell.obj",
        })
        evaluator = RUNNER.spin.SpinEvaluatorConfig(
            derivative_method="DISCRETE_ELLIPSOID", icosphere_subdivisions=1,
        )
        solver = RUNNER.spin.SpinSolverConfig(
            maximum_iterations=2, maximum_restarts=1,
            maximum_function_evaluations_per_restart=2,
            equation_profile="DUAL_MONOMIAL_RAW_V1",
        )
        _, request, _ = RUNNER._prepare(problem, evaluator, solver)
        self.assertAlmostEqual(request.fixed_shell_moments[0], 0.65625, places=14)
        self.assertEqual(request.material_density, 1.0)
        self.assertEqual(request.original_units_per_normalized_unit, 4.0)

    def test_problem_manifest_rejects_unrecognized_top_level_fields(self):
        with tempfile.TemporaryDirectory() as temporary:
            manifest_path = Path(temporary) / "problem.json"
            manifest_path.write_text(json.dumps({
            "schema": "faster-cherchi-problem-manifest-v1",
            "paper_case": "invented",
            }), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "manifest fields"):
                RUNNER._load_problem_manifest(manifest_path)

    def test_run_uses_manifest_geometry_through_cherchi_moment_bridge(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            outer_vertices = np.asarray((
                (-2., -2., -2.), (2., -2., -2.), (2., 2., -2.), (-2., 2., -2.),
                (-2., -2., 2.), (2., -2., 2.), (2., 2., 2.), (-2., 2., 2.),
            ))
            outer_faces = np.asarray((
                (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
                (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
                (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
            ), dtype=np.int32)
            design_vertices = outer_vertices / 2.
            fixed_vertices = np.vstack((outer_vertices, design_vertices))
            fixed_faces = np.vstack((outer_faces, outer_faces[:, ::-1] + 8))
            RUNNER._write_obj(root / "outer.obj", outer_vertices, outer_faces)
            RUNNER._write_obj(root / "design.obj", design_vertices, outer_faces)
            RUNNER._write_obj(root / "shell.obj", fixed_vertices, fixed_faces)
            preflight = np.asarray((
                -0.125 ** 2, 0., 0., 0., 1., 1., 1., 0., 0., 0.,
            ))
            preflight /= np.linalg.norm(preflight)
            manifest = {
                "schema": "faster-cherchi-problem-manifest-v1",
                "outer_mesh": "outer.obj",
                "fixed_regions": [{"id": "shell", "mesh": "shell.obj", "density": 0.75}],
                "design_regions": [{"id": "core", "mesh": "design.obj", "density": 2.0}],
                "length_unit": "m", "mass_unit": "kg",
                "target": {
                    "mode": "YO_YO", "axis": [1., 0., 0.],
                    # Keep a nonidentity target frame while retaining a
                    # positive first-moment capacity in the editable box.
                    # A 0.25 x-offset makes this shell/core combination
                    # continuously infeasible and is covered by the separate
                    # certified-infeasibility test.
                    "origin": [0.01, -0.125, 0.125], "support_points": [],
                },
                "minimum_outer_wall": 0.1,
                "preflight_quadric_target_normalized": preflight.tolist(),
            }
            manifest_path = root / "problem.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            build = Path(RUNNER.spin._native_module().__file__).resolve().parent.parent
            probes = sorted((RUNNER.ROOT / "build").glob(
                "cherchi-reference-singlethread-*/cherchi_probe"
            ))
            self.assertEqual(len(probes), 1)

            result = RUNNER.run(
                root / "run", build=build, probe=probes[0],
                problem_manifest=manifest_path,
                maximum_function_evaluations=4, maximum_restarts=1,
                icosphere_subdivisions=4,
                initialization_mode="deterministic_mode_seed",
                deterministic_seed=0,
                defer_finalization=True,
            )

            self.assertEqual(result["status"], "PENDING_FINALIZATION")
            begin = json.loads((root / "run" / "begin.json").read_text())
            expected_configuration = RUNNER.native_build_configuration_files(
                build, RUNNER.spin._native_module().__file__)
            self.assertEqual(set(begin["request"]["build_configuration_files"]),
                             {str(RUNNER.ROOT / "CMakeLists.txt"),
                              *(str(path) for path in expected_configuration)})
            self.assertFalse((root / "run" / "result.json").exists())
            self.assertFalse((root / "run" / "evidence.json").exists())
            # Even a real, valid bundle from a different locked problem cannot
            # complete this run. No hand-written success flag is accepted.
            from test_sealed_spin_artifact import SpinSealTests
            from test_exact_material_region import box
            fixture = SpinSealTests()
            fixture.setUp()
            mesh = box()
            unrelated_problem = fixture.problem(mesh)
            unrelated_artifact = fixture.seal(
                unrelated_problem, mesh, mesh, None, mesh, None)
            wrong_output = root / "wrong-finalization"
            wrong_output.mkdir()
            RUNNER.save_spin_artifact(unrelated_artifact, wrong_output / "sealed-material-artifact")
            payload = RUNNER._load_terminal_audit_checkpoint(
                root / "run" / "terminal-audit-checkpoint.json")
            wrong_payload = json.loads(json.dumps(payload))
            wrong_payload["run_finalization"]["output_directory"] = str(wrong_output)
            wrong_payload["run_finalization"]["expected_artifact_digest"] = unrelated_artifact.scientific_digest
            RUNNER._save_terminal_audit_checkpoint(wrong_output, wrong_payload, dependencies=[])
            with self.assertRaisesRegex(ValueError, "wrong expected problem binding"):
                RUNNER.resume_run_finalization(wrong_output)
            self.assertFalse((wrong_output / "result.json").exists())
            self.assertFalse((wrong_output / "evidence.json").exists())
            process = subprocess.run(
                [sys.executable, str(RUNNER_PATH), "--output", str(root / "run"),
                 "--resume-finalization"], capture_output=True, text=True,
                timeout=120, check=True,
            )
            self.assertEqual(json.loads(process.stdout)["status"], "RESOURCE_LIMIT_REACHED")
            result = json.loads((root / "run" / "result.json").read_text())
            self.assertEqual(result["execution_protocol"], "CHECKPOINT_RESUMED_FINALIZATION_V1")
            self.assertEqual(result['bridge_identity_attachment'],
                             RUNNER._bridge_identity_attachment(payload['run_finalization']))
            self.assertEqual(result['bridge_identity_attachment']['status'], 'BOUND_IDENTITY')
            with self.assertRaisesRegex(ValueError, "existing completed outputs"):
                RUNNER.resume_run_finalization(root / "run")

            self.assertEqual(result["classification"], "RECONSTRUCTED_METHOD_VALIDATION")
            self.assertEqual(result["function_evaluations"], 4)
            self.assertEqual(result["event_stream_version"], 6)
            self.assertTrue(json.loads((root / "run" / "inputs.json").read_text())[
                "capture_evaluation_payload"
            ])
            self.assertTrue(any(
                item.get("independent_moment_audit", {}).get("status") == "PASSED"
                for item in result["evaluations"]
            ))
            self.assertEqual(result["status"], "RESOURCE_LIMIT_REACHED")
            self.assertEqual(result["restart_count"], 1)
            self.assertTrue(result["diagnostic"])
            self.assertEqual(json.loads((root / "run" / "freshness.json").read_text())["status"], "CURRENT")
            self.assertTrue(any(
                record.get("moment_map", {}).get("sha256")
                for record in result["evaluations"]
            ), "non-identity target frames must reach an actual Cherchi moment evaluation")

    def test_locked_synthetic_problem_contains_nonzero_integrated_fixed_shell(self):
        outer_vertices = np.asarray((
            (-2., -2., -2.), (2., -2., -2.), (2., 2., -2.), (-2., 2., -2.),
            (-2., -2., 2.), (2., -2., 2.), (2., 2., 2.), (-2., 2., 2.),
        ))
        outer_faces = np.asarray((
            (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
            (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
        ), dtype=np.int32)
        design_vertices = outer_vertices / 2.0
        problem = RUNNER._synthetic_locked_problem(
            outer_vertices, outer_faces, design_vertices, outer_faces,
        )

        self.assertEqual(len(problem.fixed_regions), 1)
        self.assertEqual(len(problem.design_regions), 1)
        self.assertEqual(problem.fixed_regions[0].density, 0.75)
        self.assertEqual(problem.design_regions[0].density, 1.0)
        self.assertEqual(problem.target["mode"], "YO_YO")
        self.assertEqual(problem.minimum_outer_wall, 0.1)
        self.assertNotEqual(problem.input_digest, "")
        evaluator = RUNNER.spin.SpinEvaluatorConfig(
            derivative_method="DISCRETE_ELLIPSOID", icosphere_subdivisions=1,
        )
        solver = RUNNER.spin.SpinSolverConfig(
            maximum_iterations=2, maximum_restarts=1,
            maximum_function_evaluations_per_restart=2,
            equation_profile="DUAL_MONOMIAL_RAW_V1",
        )
        problem, request, _ = RUNNER._prepare(problem, evaluator, solver)
        self.assertAlmostEqual(request.fixed_shell_moments[0], 0.65625, places=14)
        self.assertEqual(request.material_density, 1.0)
        self.assertEqual(request.original_units_per_normalized_unit, 4.0)
        self.assertEqual(np.max(np.ptp(request.outer_vertices, axis=0)), 1.0)

    def test_uniform_request_scaling_transforms_quadric_coefficients(self):
        source = np.arange(1.0, 11.0)
        scale = 0.25
        actual = RUNNER._quadric_coefficients_after_coordinate_scale(source, scale)
        expected = source * np.asarray((1.0, 4.0, 4.0, 4.0,
                                       16.0, 16.0, 16.0, 16.0, 16.0, 16.0))
        expected /= np.linalg.norm(expected)
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-15)

    def test_same_quadric_complement_uses_opposite_boolean_side(self):
        self.assertEqual(RUNNER._opposite_boolean_operation("intersection"), "subtraction")
        self.assertEqual(RUNNER._opposite_boolean_operation("subtraction"), "intersection")
        with self.assertRaisesRegex(ValueError, "unsupported Cherchi Boolean side"):
            RUNNER._opposite_boolean_operation("union")

    def test_canonical_event_digest_ignores_only_python_traceback_line_numbers(self):
        first = [{
            "type": "evaluation",
            "detail": (
                "_ClosedSourceSelectorError: native quadric-side witness unavailable"
                "\n\nAt:\n  /repo/runner.py(41): evaluate\n"
                "  /repo/spin.py(203): solve"
            ),
        }]
        shifted = [{
            "type": "evaluation",
            "detail": (
                "_ClosedSourceSelectorError: native quadric-side witness unavailable"
                "\n\nAt:\n  /repo/runner.py(48): evaluate\n"
                "  /repo/spin.py(211): solve"
            ),
        }]
        changed_failure = [{
            "type": "evaluation",
            "detail": (
                "_ClosedSourceSelectorError: source-side witness was ambiguous"
                "\n\nAt:\n  /repo/runner.py(48): evaluate\n"
                "  /repo/spin.py(211): solve"
            ),
        }]

        digest = RUNNER._canonical_event_digest(first)

        self.assertEqual(digest, RUNNER._canonical_event_digest(shifted))
        self.assertNotEqual(digest, RUNNER._canonical_event_digest(changed_failure))
        self.assertIn("runner.py(41)", first[0]["detail"],
                      "canonical hashing must not mutate the source evidence")

    def test_scientific_event_digest_excludes_timing_but_keeps_solver_values(self):
        first = [{
            "type": "evaluation", "objective": 1.25,
            "evaluation_payload": {
                "evaluator_wall_seconds": 0.2,
                "residual_evaluation_seconds": 0.03,
                "derivative_evaluation_seconds": 0.04,
                "solver_raw_moments": [2.0, 0.5],
                "evaluation_map_identity": "a" * 64,
            },
        }]
        shifted_timing = json.loads(json.dumps(first))
        shifted_timing[0]["evaluation_payload"].update({
            "evaluator_wall_seconds": 0.9,
            "residual_evaluation_seconds": 0.13,
            "derivative_evaluation_seconds": 0.24,
        })
        shifted_science = json.loads(json.dumps(first))
        shifted_science[0]["evaluation_payload"]["solver_raw_moments"][1] = 0.6

        digest = RUNNER._scientific_event_digest_v2(first)

        self.assertEqual(digest, RUNNER._scientific_event_digest_v2(shifted_timing))
        self.assertNotEqual(digest, RUNNER._scientific_event_digest_v2(shifted_science))
        self.assertNotEqual(digest, RUNNER._scientific_event_digest_v2(
            [{"type": "evaluation", "objective": 1.5,
              "evaluation_payload": first[0]["evaluation_payload"]}]))
        self.assertEqual(first[0]["evaluation_payload"]["evaluator_wall_seconds"], 0.2)

    def test_evidence_event_builder_preserves_success_and_failed_newton_trials(self):
        result = {
            "solver_status": "RESOURCE_LIMIT_REACHED",
            "evaluations": [
                {
                    "evaluation_index": 0,
                    "quadric_coefficients": [0.0] * 10,
                    "boolean_operation": "intersection",
                    "geometry_hash": "a" * 64,
                    "topology_hash": "b" * 64,
                    "moment_map": {"sha256": "c" * 64},
                },
                {
                    "evaluation_index": 1,
                    "quadric_coefficients": [1.0] * 10,
                    "stage": "cherchi_boolean",
                    "diagnostic": "arrangement failure",
                },
            ],
        }

        events = RUNNER._evidence_events(result)

        self.assertEqual([event["iteration"] for event in events], [0, 1])
        self.assertEqual([event["status"] for event in events], ["ok", "failed"])
        self.assertEqual(events[0]["data"]["moment_map_sha256"], "c" * 64)
        self.assertEqual(events[1]["data"]["diagnostic"], "arrangement failure")

    def test_durable_evaluation_journal_recovers_success_and_failure_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            journal = Path(directory) / "evaluation-events.jsonl"
            RUNNER._initialize_evaluation_journal(journal)
            success = {"evaluation_index": 0, "moment_map": {"sha256": "a" * 64}}
            failure = {"evaluation_index": 1, "stage": "cherchi_boolean", "diagnostic": "timeout"}
            RUNNER._append_evaluation_journal(journal, success)
            RUNNER._append_evaluation_journal(journal, failure)

            self.assertEqual(RUNNER._read_evaluation_journal(journal), [success, failure])

    def test_durable_evaluation_journal_fails_closed_on_truncated_tail(self):
        with tempfile.TemporaryDirectory() as directory:
            journal = Path(directory) / "evaluation-events.jsonl"
            journal.write_text(
                '{"schema":"faster-cherchi-evaluation-journal-v1"}\n'
                '{"record_type":"evaluation","record":',
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "journal row"):
                RUNNER._read_evaluation_journal(journal)

    def test_initialization_protocol_distinguishes_fixture_from_native_seeded_starts(self):
        fixture = np.arange(10, dtype=np.float64)

        fixed = RUNNER._solver_initial_quadrics("fixed_fixture", fixture)
        seeded = RUNNER._solver_initial_quadrics("deterministic_mode_seed", fixture)

        self.assertEqual(fixed, (tuple(float(x) for x in fixture),))
        self.assertEqual(seeded, ())
        with self.assertRaisesRegex(ValueError, "unsupported initialization mode"):
            RUNNER._solver_initial_quadrics("known_solution", fixture)

    def test_exact_line_search_replay_requires_and_preserves_explicit_unit_quadric(self):
        fixture = np.zeros(10, dtype=np.float64)
        replay = (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        self.assertEqual(
            RUNNER._solver_initial_quadrics(
                "exact_line_search_replay", fixture, replay_coefficients=replay,
            ),
            (replay,),
        )
        with self.assertRaisesRegex(ValueError, "requires explicit replay coefficients"):
            RUNNER._solver_initial_quadrics("exact_line_search_replay", fixture)
        with self.assertRaisesRegex(ValueError, "unit norm"):
            RUNNER._solver_initial_quadrics(
                "exact_line_search_replay", fixture,
                replay_coefficients=(2.0,) + (0.0,) * 9,
            )

    def test_exact_replay_keeps_cgal_preflight_fixture_separate_from_solver_start(self):
        fixture = (1.0,) + (0.0,) * 9
        replay = (0.0, 1.0) + (0.0,) * 8

        preflight, initial_quadrics = RUNNER._preflight_and_solver_initial_quadrics(
            "exact_line_search_replay", fixture, replay_coefficients=replay,
        )

        self.assertEqual(preflight, fixture)
        self.assertEqual(initial_quadrics, (replay,))

    def test_cli_requires_an_explicit_replay_quadric(self):
        args = RUNNER._argument_parser().parse_args((
            "--output", "run", "--build", "build", "--cherchi-probe", "probe",
            "--initialization-mode", "exact_line_search_replay",
            "--initial-quadric-coefficients", "1", "0", "0", "0", "0",
            "0", "0", "0", "0", "0",
        ))

        self.assertEqual(args.initialization_mode, "exact_line_search_replay")
        self.assertEqual(args.initial_quadric_coefficients, [1.0] + [0.0] * 9)

    def test_cli_accepts_a_locked_problem_manifest(self):
        args = RUNNER._argument_parser().parse_args((
            "--output", "run", "--build", "build", "--cherchi-probe", "probe",
            "--problem-manifest", "locked-problem.json",
        ))

        self.assertEqual(args.problem_manifest, Path("locked-problem.json"))

    def test_run_rejects_replay_coefficients_outside_diagnostic_mode_before_side_effects(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "should-not-be-created"
            with self.assertRaisesRegex(ValueError, "requires explicit replay coefficients"):
                RUNNER.run(
                    output, build=Path("missing-build"), probe=Path("missing-probe"),
                    initialization_mode="exact_line_search_replay",
                )
            self.assertFalse(output.exists())
            with self.assertRaisesRegex(ValueError, "only valid for exact line-search replay"):
                RUNNER.run(
                    output, build=Path("missing-build"), probe=Path("missing-probe"),
                    initialization_mode="fixed_fixture",
                    initial_quadric_coefficients=[1.0] + [0.0] * 9,
                )
            self.assertFalse(output.exists())

    def test_cli_records_native_deterministic_seed_without_changing_fixture_default(self):
        default = RUNNER._argument_parser().parse_args((
            "--output", "run", "--build", "build", "--cherchi-probe", "probe",
        ))
        seeded = RUNNER._argument_parser().parse_args((
            "--output", "run", "--build", "build", "--cherchi-probe", "probe",
            "--initialization-mode", "deterministic_mode_seed", "--deterministic-seed", "1",
        ))

        self.assertEqual(default.initialization_mode, "fixed_fixture")
        self.assertEqual(default.deterministic_seed, 0)
        self.assertEqual(seeded.initialization_mode, "deterministic_mode_seed")
        self.assertEqual(seeded.deterministic_seed, 1)

    def test_cli_exposes_a_bounded_restart_budget_without_changing_legacy_default(self):
        default = RUNNER._argument_parser().parse_args((
            "--output", "run", "--build", "build", "--cherchi-probe", "probe",
        ))
        seeded = RUNNER._argument_parser().parse_args((
            "--output", "run", "--build", "build", "--cherchi-probe", "probe",
            "--maximum-restarts", "8",
        ))

        self.assertEqual(default.maximum_restarts, 1)
        self.assertEqual(seeded.maximum_restarts, 8)

    def test_strict_default_icosphere_level_covers_captured_near_degenerate_trial(self):
        args = RUNNER._argument_parser().parse_args((
            "--output", "run",
            "--build", "build",
            "--cherchi-probe", "probe",
        ))

        self.assertEqual(args.icosphere_subdivisions, 5)

    def test_independent_polyhedral_raw_moments_match_analytic_density_cube(self):
        vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
        ))
        faces = np.asarray((
            (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
            (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)

        actual = RUNNER._independent_polyhedral_raw_moments(vertices, faces, 2.0)

        np.testing.assert_allclose(
            actual,
            (2.0, 1.0, 1.0, 1.0, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0,
             0.5, 0.5, 0.5),
            rtol=0.0, atol=1e-14,
        )

    def test_independent_moment_agreement_uses_fixed_input_scales_and_fails_closed(self):
        expected = np.asarray((2., 1., 1., 1., 2/3, 2/3, 2/3, .5, .5, .5))
        audit = RUNNER._audit_independent_moment_agreement(
            expected, expected, density=2.0, root_side=1.0,
        )
        self.assertEqual(audit["status"], "PASSED")
        self.assertEqual(audit["normalized_inf_residual"], 0.0)

        wrong = expected.copy()
        wrong[4] += 1e-3
        with self.assertRaisesRegex(ValueError, "independent raw-moment residual"):
            RUNNER._audit_independent_moment_agreement(
                expected, wrong, density=2.0, root_side=1.0,
            )

    def test_directional_map_identity_binds_both_real_basis_and_support(self):
        u=np.eye(10)[:,[0,1,2,3,7,8]]
        q=np.zeros(10);q[0]=1
        b=u[:,1:];j=np.arange(50,dtype=float).reshape(10,5)
        arguments=dict(quadric_coefficients=q, backend_identity='backend',
            derivative_identity='directional',source_identity='source',cap_identity='cap',
            native_tangent_basis=b,support_basis=u)
        value=RUNNER._moment_map_binding_sha256('a'*64,'b'*64,'intersection',np.arange(10),j,**arguments)
        self.assertNotEqual(value,RUNNER._moment_map_binding_sha256(
            'a'*64,'b'*64,'intersection',np.arange(10),-j,
            **{**arguments,'native_tangent_basis':-b}))
        # Same represented subspace but a different encoded support must not
        # alias its initialization/coordinate identity.
        rotated=u.copy();rotated[:,0]*=-1
        self.assertNotEqual(value,RUNNER._moment_map_binding_sha256(
            'a'*64,'b'*64,'intersection',np.arange(10),j,
            **{**arguments,'support_basis':rotated}))
        for changes in ({'native_tangent_basis':None}, {'native_tangent_basis':np.zeros((10,5))}):
            with self.subTest(changes=changes),self.assertRaises(ValueError):
                RUNNER._moment_map_binding_sha256('a'*64,'b'*64,'intersection',np.arange(10),j,
                    **{**arguments,**changes})

    def test_pure_height_source_rebuilds_genuine_motion_at_each_trial(self):
        from make_it_stand import spin
        native=spin._native_module();settings=native.SpinTriangulationSettings()
        settings.paraboloid_grid_resolution=8;settings.truncation_radius=4.
        settings.maximum_radius_doublings=4
        q=np.array([.125,.3,.7,-1,0,0,0,.2,.4,0.]);q/=np.linalg.norm(q)
        domain=np.array([[-.687,-.871,-.743],[.913,.729,.857]])
        first=RUNNER._pure_height_source_motion(native,q,domain,settings)
        q2=q.copy();q2[1]+=.003;q2/=np.linalg.norm(q2)
        second=RUNNER._pure_height_source_motion(native,q2,domain,settings)
        for coefficients,source in ((q,first),(q2,second)):
            b=source['native_tangent_basis'];motion=np.asarray(source['vertex_velocity'])
            self.assertEqual(b.shape,(10,5));self.assertEqual(motion.shape,(len(source['vertices'])*3,5))
            np.testing.assert_allclose(b.T@coefficients,0,atol=1e-12)
            np.testing.assert_allclose(b.T@b,np.eye(5),atol=1e-12)
            reference=native.evaluate_spin_it_faster_stratum_source_motion(coefficients,domain,settings)
            np.testing.assert_allclose(motion,np.asarray(reference['vertex_velocity'])@
                (np.asarray(reference['coefficient_basis']).T@b),atol=1e-12)
        self.assertFalse(np.array_equal(first['native_tangent_basis'],second['native_tangent_basis']))
        off=q.copy();off[4]=1e-3
        with self.assertRaisesRegex(ValueError,'support'):
            RUNNER._pure_height_source_motion(native,off,domain,settings)

    def test_pure_height_seed_protocol_has_no_fixture_or_answer_input(self):
        first=RUNNER._pure_height_initial_quadrics(0,8)
        second=RUNNER._pure_height_initial_quadrics(0,8)
        other=RUNNER._pure_height_initial_quadrics(1,8)
        np.testing.assert_array_equal(first,second)
        self.assertFalse(np.array_equal(first,other))
        self.assertEqual(np.asarray(first).shape,(8,10))
        np.testing.assert_allclose(np.linalg.norm(first,axis=1),1,atol=1e-15)
        np.testing.assert_array_equal(np.asarray(first)[:,[4,5,6,9]],np.zeros((8,4)))
        for seed,count in ((True,1),(-1,1),(0,0),(0,True)):
            with self.subTest(seed=seed,count=count),self.assertRaises(ValueError):
                RUNNER._pure_height_initial_quadrics(seed,count)

    def test_grid_stratum_uses_not_applicable_normals_without_admitting_mixed_levels(self):
        branch=dict(surface_type='PARABOLOID',truncation_radius_hex=(4.).hex(),
            initial_icosphere_subdivisions=-1,effective_icosphere_subdivisions=-1,
            directional_profile='PURE_HEIGHT_ACTIVE_STRATUM_V1',paraboloid_grid_resolution=8)
        self.assertEqual(len(RUNNER._derivative_stratum_identity('intersection','a'*64,branch)),64)
        for change in ({'initial_icosphere_subdivisions':0},{'effective_icosphere_subdivisions':0}):
            with self.assertRaises(ValueError):
                RUNNER._derivative_stratum_identity('intersection','a'*64,{**branch,**change})

    def test_pure_height_profile_rejects_wrong_mode_before_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'never-created'
            with self.assertRaisesRegex(ValueError,'pure-height'):
                RUNNER.run(output,build=Path('/unavailable'),probe=Path('/unavailable'),
                    directional_profile='PURE_HEIGHT_ACTIVE_STRATUM_V1')
            self.assertFalse(output.exists())

    def test_actual_directional_runner_records_seeded_cherchi_map_and_native_dimensions(self):
        manifest=ROOT/'experiments/spin_it_reproduction/faster/fixtures/cherchi-centered-axis-pure-height-method-v1/problem.json'
        build=Path(RUNNER.spin._native_module().__file__).resolve().parent.parent
        probe=ROOT/'build/cherchi-reference-singlethread-20260917/cherchi_probe'
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'run'
            result=RUNNER.run(output,build=build,probe=probe,problem_manifest=manifest,
                maximum_function_evaluations=1,maximum_restarts=1,maximum_radius_doublings=4,
                spinning_top_objective='LOWER_CENTER_OF_MASS',
                directional_profile='PURE_HEIGHT_ACTIVE_STRATUM_V1',
                initialization_mode='supported_mode_seed',deterministic_seed=0,
                defer_finalization=True)
            inputs=json.loads((output/'inputs.json').read_text())
            np.testing.assert_array_equal(inputs['initial_quadrics'],RUNNER._pure_height_initial_quadrics(0,1))
            self.assertNotEqual(inputs['initial_coefficients'],inputs['preflight_coefficients'])
            payload=RUNNER._load_terminal_audit_checkpoint(output/'terminal-audit-checkpoint.json')
            native=payload['run_finalization']['native_result']
            self.assertEqual(native['function_evaluations'],1)
            records=native['supported_space']['evaluations']
            self.assertEqual(len(records),1)
            self.assertIsNotNone(records[0]['evaluation'],records[0]['error'])
            maps=payload['evaluations'];self.assertEqual(len(maps),1)
            self.assertEqual(maps[0]['jacobian_shape'],[10,5])
            self.assertEqual(np.asarray(maps[0]['moment_map']['native_tangent_basis']).shape,(10,5))
            checked=RUNNER._audit_newton_evaluation_bindings(maps,native['events'],
                quadric_hash=RUNNER.spin._native_module().spin_it_faster_quadric_hash)
            self.assertEqual(checked['callback_map_count'],1)
            self.assertNotEqual(native['status'],'CONVERGED')

    def test_moment_map_binding_hash_covers_geometry_operation_moments_and_jacobian(self):
        moments = np.arange(10, dtype=np.float64)
        jacobian = np.arange(100, dtype=np.float64).reshape((10, 10))
        identity = {
            "quadric_coefficients": np.arange(10, dtype=np.float64),
            "backend_identity": "cherchi-reference:" + "c" * 64,
            "derivative_identity": "cherchi-provenance-reynolds-v1",
            "source_identity": "source-geometry-and-triangulation-v1",
            "cap_identity": "hyperboloid-truncation-and-tessellation-v1",
            "exact_geometry_identity": "rational-region-sha256-v1",
        }

        digest = RUNNER._moment_map_binding_sha256(
            "0" * 64, "1" * 64, "subtraction", moments, jacobian, **identity,
        )
        self.assertRegex(digest, r"^[0-9a-f]{64}$")
        self.assertNotEqual(digest, RUNNER._moment_map_binding_sha256(
            "2" * 64, "1" * 64, "subtraction", moments, jacobian, **identity,
        ))
        self.assertNotEqual(digest, RUNNER._moment_map_binding_sha256(
            "0" * 64, "1" * 64, "intersection", moments, jacobian, **identity,
        ))
        changed_moments = moments.copy()
        changed_moments[0] = np.nextafter(changed_moments[0], np.inf)
        self.assertNotEqual(digest, RUNNER._moment_map_binding_sha256(
            "0" * 64, "1" * 64, "subtraction", changed_moments, jacobian, **identity,
        ))
        changed_jacobian = jacobian.copy()
        changed_jacobian[9, 9] = np.nextafter(changed_jacobian[9, 9], np.inf)
        self.assertNotEqual(digest, RUNNER._moment_map_binding_sha256(
            "0" * 64, "1" * 64, "subtraction", moments, changed_jacobian, **identity,
        ))
        for key, value in identity.items():
            altered = dict(identity)
            altered[key] = (value + 1.0) if key == "quadric_coefficients" else value + "-changed"
            self.assertNotEqual(digest, RUNNER._moment_map_binding_sha256(
                "0" * 64, "1" * 64, "subtraction", moments, jacobian, **altered,
            ), key)

    def test_newton_event_geometry_must_match_its_returned_moment_map(self):
        moments = np.arange(10, dtype=np.float64)
        jacobian = np.eye(10, dtype=np.float64)
        geometry_hash, operation = "a" * 64, "intersection"
        construction_topology_hash = "b" * 64
        backend_identity = "cherchi-reference:" + "c" * 64
        derivative_identity = "cherchi-provenance-reynolds-v1"
        source_geometry_hash = "c" * 64
        source_triangulation_hash = "d" * 64
        cap_payload = {
            "surface_type": "HYPERBOLOID",
            "truncation_radius_hex": float(2.0).hex(),
            "initial_icosphere_subdivisions": 3,
            "effective_icosphere_subdivisions": 4,
            "maximum_radius_doublings": 4,
            "triangulation_hash": source_triangulation_hash,
        }
        derivative_branch = {
            key: cap_payload[key] for key in (
                "surface_type", "truncation_radius_hex",
                "initial_icosphere_subdivisions", "effective_icosphere_subdivisions",
            )
        }
        topology_hash = hashlib.sha256(json.dumps({
            "schema": "faster-cherchi-derivative-stratum-v2",
            "operation": operation,
            "construction_topology": construction_topology_hash,
            "derivative_branch": derivative_branch,
        }, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("ascii")).hexdigest()
        source_identity = hashlib.sha256(
            (source_geometry_hash + "\0" + source_triangulation_hash).encode("ascii")
        ).hexdigest()
        cap_identity = hashlib.sha256(json.dumps(
            cap_payload, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode("ascii")).hexdigest()
        identity = {
            "quadric_coefficients": np.arange(10, dtype=np.float64),
            "boolean_backend_identity": backend_identity,
            "derivative_identity": derivative_identity,
            "source_identity": source_identity,
            "cap_identity": cap_identity,
        }
        map_hash = RUNNER._moment_map_binding_sha256(
            geometry_hash, topology_hash, operation, moments, jacobian,
            quadric_coefficients=identity["quadric_coefficients"],
            backend_identity=identity["boolean_backend_identity"],
            derivative_identity=identity["derivative_identity"],
            source_identity=identity["source_identity"],
            cap_identity=identity["cap_identity"],
        )
        evaluation = {
            "evaluation_index": 4,
            "solver_evaluation_index": 5,
            "geometry_hash": geometry_hash,
            "topology_hash": topology_hash,
            "construction_topology_hash": construction_topology_hash,
            "derivative_branch": derivative_branch,
            "boolean_operation": operation,
            "source_geometry_hash": source_geometry_hash,
            "source_triangulation_hash": source_triangulation_hash,
            "cap_identity_payload": cap_payload,
            "quadric_hash": hashlib.sha256(identity["quadric_coefficients"].tobytes()).hexdigest(),
            **identity,
            "evaluation_map_identity": map_hash,
            "moment_map": {
                "raw_moments": moments.tolist(),
                "moment_jacobian": jacobian.tolist(),
                "sha256": map_hash,
            },
        }
        event = {
            "evaluation_index": 5,
            "material_geometry_hash": geometry_hash,
            "material_topology_hash": topology_hash,
            "quadric_coefficients": identity["quadric_coefficients"].tolist(),
            "quadric_hash": evaluation["quadric_hash"],
            "detail": (
                "evaluation_succeeded;evaluation_map=" + map_hash +
                ";backend=" + backend_identity + ";derivative=" + derivative_identity
            ),
        }

        quadric_hash = lambda values: hashlib.sha256(
            np.asarray(values, dtype=np.float64).tobytes()
        ).hexdigest()
        audit = RUNNER._audit_newton_evaluation_bindings(
            [evaluation], [event], quadric_hash=quadric_hash,
        )

        self.assertEqual(audit["status"], "CONSISTENT")
        self.assertEqual(audit["bound_evaluation_count"], 1)
        algebraic_branch=dict(derivative_branch,finite_point_map="ALGEBRAIC_V1")
        algebraic_topology=hashlib.sha256(json.dumps({
            "schema":"faster-cherchi-derivative-stratum-v2", "operation":operation,
            "construction_topology":construction_topology_hash,
            "derivative_branch":algebraic_branch,
        },sort_keys=True,separators=(',',':'),allow_nan=False).encode('ascii')).hexdigest()
        algebraic_cap=dict(cap_payload,finite_point_map="ALGEBRAIC_V1")
        algebraic_cap_hash=hashlib.sha256(json.dumps(algebraic_cap,sort_keys=True,
            separators=(',',':'),allow_nan=False).encode('ascii')).hexdigest()
        algebraic_map=RUNNER._moment_map_binding_sha256(
            geometry_hash,algebraic_topology,operation,moments,jacobian,
            quadric_coefficients=identity['quadric_coefficients'],backend_identity=backend_identity,
            derivative_identity=derivative_identity,source_identity=source_identity,
            cap_identity=algebraic_cap_hash)
        algebraic_evaluation=dict(evaluation,topology_hash=algebraic_topology,
            derivative_branch=algebraic_branch,cap_identity_payload=algebraic_cap,
            cap_identity=algebraic_cap_hash,evaluation_map_identity=algebraic_map,
            moment_map=dict(evaluation['moment_map'],sha256=algebraic_map))
        algebraic_event=dict(event,material_topology_hash=algebraic_topology,
            detail=event['detail'].replace(map_hash,algebraic_map))
        self.assertEqual(RUNNER._audit_newton_evaluation_bindings(
            [algebraic_evaluation],[algebraic_event],quadric_hash=quadric_hash)['status'],'CONSISTENT')
        with self.assertRaisesRegex(ValueError,'finite|map'):
            RUNNER._derivative_stratum_identity(operation,construction_topology_hash,
                dict(algebraic_branch,finite_point_map='UNSUPPORTED_MAP'))
        features = {
            'schema': 'paper-projective-source-features-v1',
            'surface_type': 'HYPERBOLOID', 'vertex_count': 3,
            'normal_vertices': [0, -1, -1],
            'normal_edges': [[-1, -1]] * 3,
            'projective_rays': [[-1]*4, [0, 1, 2, 0], [0, 2, 3, 0]],
            'faces': [[0, 1, 2]],
        }
        paper_branch = dict(derivative_branch, paper_source_features_sha256=RUNNER.digest(features))
        paper_branch['effective_icosphere_subdivisions'] = 3
        paper_topology = RUNNER._derivative_stratum_identity(
            operation, construction_topology_hash, paper_branch)
        paper_cap = dict(cap_payload, **paper_branch)
        paper_cap_hash = hashlib.sha256(json.dumps(paper_cap, sort_keys=True,
            separators=(',', ':'), allow_nan=False).encode('ascii')).hexdigest()
        paper_map = RUNNER._moment_map_binding_sha256(
            geometry_hash, paper_topology, operation, moments, jacobian,
            quadric_coefficients=identity['quadric_coefficients'],
            backend_identity=backend_identity, derivative_identity=derivative_identity,
            source_identity=source_identity, cap_identity=paper_cap_hash)
        paper_eval = dict(evaluation, topology_hash=paper_topology,
            derivative_branch=paper_branch, cap_identity_payload=paper_cap,
            cap_identity=paper_cap_hash, source_feature_payload=features,
            evaluation_map_identity=paper_map,
            moment_map=dict(evaluation['moment_map'], sha256=paper_map))
        paper_event = dict(event, material_topology_hash=paper_topology,
            detail=event['detail'].replace(map_hash, paper_map))
        self.assertEqual(RUNNER._audit_newton_evaluation_bindings(
            [paper_eval], [paper_event], quadric_hash=quadric_hash)['status'], 'CONSISTENT')
        for payload in (None, dict(features, projective_rays=[[-1]*4]*3)):
            with self.subTest(payload=payload), self.assertRaisesRegex(ValueError, 'source feature'):
                RUNNER._audit_newton_evaluation_bindings(
                    [dict(paper_eval, source_feature_payload=payload)], [paper_event],
                    quadric_hash=quadric_hash)
        returned_but_rejected = dict(evaluation, outcome="returned")
        rejected_event = {
            "type": "evaluation", "evaluation_index": 5,
            "outcome": "failed", "reason": "evaluation_exception",
            "detail": "invalid exact solver candidate boundary",
            "material_geometry_hash": None, "material_topology_hash": None,
        }
        rejected_audit = RUNNER._audit_newton_evaluation_bindings(
            [returned_but_rejected], [rejected_event], quadric_hash=quadric_hash,
        )
        self.assertEqual(rejected_audit["status"], "CONSISTENT")
        self.assertEqual(rejected_audit["bound_evaluation_count"], 0)
        self.assertEqual(rejected_audit["solver_rejected_geometry_evaluations"], [5])
        self.assertEqual(audit["map_identities"], [{
            "evaluation_index": 5,
            "evaluation_map_identity": map_hash,
            "boolean_backend_identity": backend_identity,
            "derivative_identity": derivative_identity,
            "source_identity": identity["source_identity"],
            "cap_identity": identity["cap_identity"],
        }])
        changed_event = dict(event, material_geometry_hash="c" * 64)
        with self.assertRaisesRegex(ValueError, "Newton event geometry identity mismatch"):
            RUNNER._audit_newton_evaluation_bindings(
                [evaluation], [changed_event], quadric_hash=quadric_hash,
            )
        changed_mapping = dict(evaluation)
        changed_mapping["moment_map"] = dict(evaluation["moment_map"])
        changed_mapping["moment_map"]["raw_moments"] = [2.0] + moments[1:].tolist()
        with self.assertRaisesRegex(ValueError, "moment-map binding hash mismatch"):
            RUNNER._audit_newton_evaluation_bindings(
                [changed_mapping], [event], quadric_hash=quadric_hash,
            )
        changed_backend = dict(event, detail=event["detail"].replace(backend_identity, "other-backend"))
        with self.assertRaisesRegex(ValueError, "Newton event backend identity mismatch"):
            RUNNER._audit_newton_evaluation_bindings(
                [evaluation], [changed_backend], quadric_hash=quadric_hash,
            )
        changed_quadric = dict(event, quadric_coefficients=[9.0] + [0.0] * 9)
        with self.assertRaisesRegex(ValueError, "Newton event quadric coefficients mismatch"):
            RUNNER._audit_newton_evaluation_bindings(
                [evaluation], [changed_quadric], quadric_hash=quadric_hash,
            )
        missing_quadric_hash = dict(event)
        missing_quadric_hash.pop("quadric_hash")
        with self.assertRaisesRegex(ValueError, "Newton event omitted quadric hash"):
            RUNNER._audit_newton_evaluation_bindings(
                [evaluation], [missing_quadric_hash], quadric_hash=quadric_hash,
            )
        wrong_quadric_hash = dict(event, quadric_hash="f" * 64)
        with self.assertRaisesRegex(ValueError, "Newton event quadric hash mismatch"):
            RUNNER._audit_newton_evaluation_bindings(
                [evaluation], [wrong_quadric_hash], quadric_hash=quadric_hash,
            )
        forged_quadric_hash = dict(event, quadric_hash="e" * 64)
        forged_evaluation = dict(evaluation, quadric_hash="e" * 64)
        with self.assertRaisesRegex(ValueError, "does not match native coefficient digest"):
            RUNNER._audit_newton_evaluation_bindings(
                [forged_evaluation], [forged_quadric_hash], quadric_hash=quadric_hash,
            )
        changed_source = dict(evaluation, source_geometry_hash="d" * 64)
        with self.assertRaisesRegex(ValueError, "source identity does not match source geometry"):
            RUNNER._audit_newton_evaluation_bindings(
                [changed_source], [event], quadric_hash=quadric_hash,
            )
        changed_cap = dict(evaluation)
        changed_cap["cap_identity_payload"] = dict(cap_payload, surface_type="ELLIPSOID")
        with self.assertRaisesRegex(ValueError, "cap identity does not match cap payload"):
            RUNNER._audit_newton_evaluation_bindings(
                [changed_cap], [event], quadric_hash=quadric_hash,
            )
        changed_cap_branch = dict(cap_payload, truncation_radius_hex=float(4.0).hex())
        changed_cap_identity = hashlib.sha256(json.dumps(
            changed_cap_branch, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode("ascii")).hexdigest()
        changed_branch_map_hash = RUNNER._moment_map_binding_sha256(
            geometry_hash, topology_hash, operation, moments, jacobian,
            quadric_coefficients=identity["quadric_coefficients"],
            backend_identity=backend_identity, derivative_identity=derivative_identity,
            source_identity=source_identity, cap_identity=changed_cap_identity,
        )
        stale_stratum = dict(evaluation)
        stale_stratum["cap_identity_payload"] = changed_cap_branch
        stale_stratum["cap_identity"] = changed_cap_identity
        stale_stratum["derivative_branch"] = dict(
            derivative_branch, truncation_radius_hex=float(4.0).hex(),
        )
        stale_stratum["evaluation_map_identity"] = changed_branch_map_hash
        stale_stratum["moment_map"] = dict(evaluation["moment_map"], sha256=changed_branch_map_hash)
        stale_event = dict(event, detail=event["detail"].replace(map_hash, changed_branch_map_hash))
        with self.assertRaisesRegex(ValueError, "derivative stratum identity mismatch"):
            RUNNER._audit_newton_evaluation_bindings(
                [stale_stratum], [stale_event], quadric_hash=quadric_hash,
            )
        incomplete_payload = dict(cap_payload)
        incomplete_payload.pop("maximum_radius_doublings")
        incomplete_cap_identity = hashlib.sha256(json.dumps(
            incomplete_payload, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode("ascii")).hexdigest()
        incomplete_map_hash = RUNNER._moment_map_binding_sha256(
            geometry_hash, topology_hash, operation, moments, jacobian,
            quadric_coefficients=identity["quadric_coefficients"],
            backend_identity=backend_identity, derivative_identity=derivative_identity,
            source_identity=source_identity, cap_identity=incomplete_cap_identity,
        )
        incomplete = dict(evaluation)
        incomplete["cap_identity_payload"] = incomplete_payload
        incomplete["cap_identity"] = incomplete_cap_identity
        incomplete["evaluation_map_identity"] = incomplete_map_hash
        incomplete["moment_map"] = dict(evaluation["moment_map"], sha256=incomplete_map_hash)
        incomplete_event = dict(event, detail=event["detail"].replace(map_hash, incomplete_map_hash))
        with self.assertRaisesRegex(ValueError, "cap identity payload schema is invalid"):
            RUNNER._audit_newton_evaluation_bindings(
                [incomplete], [incomplete_event], quadric_hash=quadric_hash,
            )

    def test_converged_terminal_geometry_audit_binds_the_accepted_native_state(self):
        vertices = np.asarray(((0., 0., 0.), (1., 0., 0.),
                               (0., 1., 0.), (0., 0., 1.)))
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)))
        with tempfile.TemporaryDirectory() as directory:
            mesh_path = Path(directory) / "terminal-oriented.obj"
            RUNNER._write_obj(mesh_path, vertices, faces)
            mesh_digest = hashlib.sha256(mesh_path.read_bytes()).hexdigest()
            coefficients = np.zeros(10, dtype=np.float64)
            coefficients[0] = 1.0
            native = RUNNER.spin._native_module()
            qhash = native.spin_it_faster_quadric_hash(coefficients)
            moments = np.asarray(RUNNER.spin.integrate_raw_mass_moments(
                vertices.tolist(), faces.tolist(), 1.0), dtype=np.float64)
            independent = RUNNER._independent_polyhedral_raw_moments(
                vertices, faces, 1.0)
            independent_audit = RUNNER._audit_independent_moment_agreement(
                moments, independent, density=1.0, root_side=1.0)
            evaluation = {
                # The solver's trace counter is local to a restart.  A prior
                # restart can therefore have the same per-restart value as a
                # different global evaluation index.
                "solver_evaluation_index": 9,
                "quadric_coefficients": coefficients.tolist(),
                "quadric_hash": qhash,
                "geometry_hash": mesh_digest,
                "topology_hash": "a" * 64,
                "oriented_output_path": str(mesh_path),
                "moment_map": {"raw_moments": moments.tolist()},
                "independent_moment_audit": independent_audit,
            }
            previous_restart_coefficients = np.zeros(10, dtype=np.float64)
            previous_restart_coefficients[1] = 1.0
            previous_restart_evaluation = {
                "solver_evaluation_index": 4,
                "quadric_coefficients": previous_restart_coefficients.tolist(),
                "quadric_hash": native.spin_it_faster_quadric_hash(
                    previous_restart_coefficients),
                "geometry_hash": "c" * 64,
                "topology_hash": "d" * 64,
                "oriented_output_path": str(mesh_path),
            }
            result = {
                "solver_status": "CONVERGED",
                "evaluation_map_binding_audit": {"status": "CONSISTENT"},
                "evaluations": [previous_restart_evaluation, evaluation],
                "design_domain_vertices": [
                    (0., 0., 0.), (2., 0., 0.), (0., 2., 0.), (0., 0., 2.),
                ],
                "design_domain_faces": faces.tolist(),
                "trace": [{
                    "function_evaluations": 4,
                    "restart": 2,
                    "quadric_hash": qhash,
                    "material_geometry_hash": mesh_digest,
                    "material_topology_hash": "a" * 64,
                }],
            }

            audit = RUNNER._audit_terminal_geometry_binding(
                result, quadric_hash=native.spin_it_faster_quadric_hash,
            )

            self.assertEqual(audit["status"], "PASSED")
            self.assertTrue(audit["geometry_verified"])
            self.assertEqual(audit["terminal_solver_evaluation_index"], 9)
            self.assertEqual(audit["terminal_restart_function_evaluations"], 4)
            self.assertEqual(audit["terminal_material_sha256"], mesh_digest)
            self.assertTrue(audit["exact_boundary_audit"]["self_intersection_free"])
            self.assertEqual(audit["material_domain_status"], "VERIFIED_SUBSET")
            self.assertTrue(audit["material_domain_relation"]["material_subset_of_design_domain"])
            self.assertFalse(result.get("final_geometry_verified", False),
                             "external audit must not rewrite the native solver claim")

            checkpoint = RUNNER._save_terminal_audit_checkpoint(
                Path(directory), result, dependencies=[mesh_path])
            resumed = RUNNER.resume_terminal_audit(Path(directory))
            self.assertEqual(resumed["audit"], audit)
            self.assertEqual(resumed["status"], "TERMINAL_AUDIT_PASSED")
            self.assertFalse(resumed["complete_run_evidence"])
            self.assertTrue(checkpoint.is_file())
            self.assertFalse((Path(directory) / "evidence.json").exists())
            for index in range(2):
                isolated = Path(directory) / f"independent-{index}"
                isolated.mkdir()
                (isolated / checkpoint.name).write_bytes(checkpoint.read_bytes())
                process = subprocess.run(
                    [sys.executable, str(RUNNER_PATH), "--output", str(isolated),
                     "--resume-terminal-audit"], capture_output=True, text=True,
                    timeout=60, check=True,
                )
                summary = json.loads(process.stdout)
                self.assertEqual(summary["status"], "TERMINAL_AUDIT_PASSED")
                saved = json.loads((isolated / "resumed-terminal-audit.json").read_text())
                self.assertEqual(saved["audit"], audit)
                self.assertFalse(saved["complete_run_evidence"])

            outside_domain = {**result, "design_domain_vertices": [
                (0., 0., 0.), (0.5, 0., 0.), (0., 0.5, 0.), (0., 0., 0.5),
            ]}
            failed_audit = RUNNER._audit_terminal_geometry_binding(
                outside_domain, quadric_hash=native.spin_it_faster_quadric_hash,
            )
            self.assertEqual(failed_audit["status"], "FAILED")
            self.assertFalse(failed_audit["geometry_verified"])
            self.assertEqual(failed_audit["material_domain_status"], "OUTSIDE_DESIGN_DOMAIN")
            self.assertTrue(failed_audit["material_domain_relation"]["completed"])
            self.assertFalse(
                failed_audit["material_domain_relation"]["material_subset_of_design_domain"]
            )

            mismatched = {**result, "trace": [dict(
                result["trace"][0], material_geometry_hash="b" * 64)]}
            with self.assertRaisesRegex(ValueError, "no unique material identity match"):
                RUNNER._audit_terminal_geometry_binding(
                    mismatched, quadric_hash=native.spin_it_faster_quadric_hash,
                )

    def test_exact_material_domain_relation_rejects_a_disjoint_terminal_solid(self):
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)))
        material = np.asarray(((0., 0., 0.), (1., 0., 0.),
                               (0., 1., 0.), (0., 0., 1.)))
        design_domain = np.asarray(((0., 0., 0.), (2., 0., 0.),
                                   (0., 2., 0.), (0., 0., 2.)))

        relation = RUNNER._audit_material_domain_subset(
            material + np.asarray((4.0, 0.0, 0.0)), faces,
            design_domain, faces,
        )

        self.assertTrue(relation["completed"])
        self.assertEqual(relation["status"], "OUTSIDE_DESIGN_DOMAIN")
        self.assertFalse(relation["material_subset_of_design_domain"])
        self.assertGreater(relation["aabb_diagnostic_only"][
            "vertices_strictly_outside_design_aabb"], 0)
        self.assertTrue(relation["aabb_diagnostic_only"]["not_a_containment_proof"])

    def test_newton_map_audit_rejects_empty_or_all_failed_evidence(self):
        with self.assertRaisesRegex(ValueError, "no successful evaluation maps to audit"):
            RUNNER._audit_newton_evaluation_bindings(
                [], [], quadric_hash=lambda values: "0" * 64,
            )
        with self.assertRaisesRegex(ValueError, "no successful evaluation maps to audit"):
            RUNNER._audit_newton_evaluation_bindings(
                [{"outcome": "failed", "evaluation_index": 0}], [],
                quadric_hash=lambda values: "0" * 64,
            )

    def test_raw_boolean_topology_audit_preserves_open_plane_boundary(self):
        vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
        ))
        # A cube missing its z=0 cap: the diagnostic must report an open
        # boundary, not silently close it or call the raw output valid.
        faces = np.asarray((
            (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
            (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6),
            (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)

        audit = RUNNER._audit_raw_boolean_topology(vertices, faces)

        self.assertEqual(audit["status"], "OPEN_BOUNDARY")
        self.assertEqual(audit["boundary_edge_count"], 4)
        self.assertEqual(audit["boundary_loop_sizes"], [4])
        self.assertEqual(audit["boundary_vertices_not_degree_two"], 0)
        self.assertTrue(any(
            plane["axis"] == "z" and plane["maximum_coordinate_spread"] == 0.0
            for plane in audit["near_axis_aligned_boundary_planes"]
        ))
        self.assertFalse(audit["closed_and_manifold"])

    def test_raw_boolean_topology_audit_detects_nonmanifold_edges(self):
        vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.),
            (0., -1., 0.),
        ))
        faces = np.asarray(((0, 1, 2), (1, 0, 3), (0, 1, 4)), dtype=np.int64)

        audit = RUNNER._audit_raw_boolean_topology(vertices, faces)

        self.assertEqual(audit["status"], "NONMANIFOLD_EDGE_INCIDENCE")
        self.assertEqual(audit["nonmanifold_edge_count"], 1)
        self.assertFalse(audit["closed_and_manifold"])

    def test_raw_boundary_audit_distinguishes_open_from_reorientable_faces(self):
        vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
        ))
        faces = np.asarray((
            (4, 5, 6), (4, 6, 7), (0, 5, 1), (0, 5, 4),
            (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6),
            (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)

        audit = RUNNER._audit_raw_boolean_topology(vertices, faces)

        self.assertEqual(audit["status"], "OPEN_BOUNDARY")
        self.assertFalse(audit["orientation_consistent_as_emitted"])
        self.assertEqual(audit["orientation_conflict_edge_count"], 0)
        self.assertEqual(audit["orientable_face_component_count"], 1)
        self.assertEqual(audit["suggested_face_flip_count"], 1)
        self.assertEqual(audit["directed_boundary_cycle_count_after_consistent_reorientation"], 1)
        self.assertFalse(audit["repair_performed"])

    def test_exact_boundary_provenance_joins_sidecar_by_coordinate(self):
        vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
        ))
        faces = np.asarray((
            (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
            (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6),
            (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)
        exact_header = (
            "vertex_index\tstatus\texact_x\texact_y\texact_z\t"
            "memory_x\tmemory_y\tmemory_z\tobj_x\tobj_y\tobj_z"
        )
        exact_rows = [exact_header]
        # Deliberately reverse the sidecar order; joining by row number would
        # misidentify the open boundary plane.
        for index in reversed(range(len(vertices))):
            point = vertices[index]
            exact = tuple(f"{int(value)}/1" for value in point)
            fields = [str(index), "MATCHED", *exact, *(repr(float(x)) for x in point),
                      *(repr(float(x)) for x in point)]
            exact_rows.append("\t".join(fields))
        label_rows = ["face_index\tinput_object_mask"] + [
            f"{index}\t0x2" for index in range(len(faces))
        ]
        with tempfile.TemporaryDirectory() as directory:
            exact_path = Path(directory) / "mesh.exact.tsv"
            labels_path = Path(directory) / "mesh.faces.tsv"
            exact_path.write_text("\n".join(exact_rows) + "\n", encoding="ascii")
            labels_path.write_text("\n".join(label_rows) + "\n", encoding="ascii")

            audit = RUNNER._audit_exact_boolean_boundary(
                vertices, faces, exact_path, labels_path,
            )

        self.assertEqual(audit["boundary_edge_count"], 4)
        self.assertEqual(audit["components"][0]["exact_axis_aligned_planes"], [
            {"axis": "z", "coordinate_exact": "0/1"},
        ])
        self.assertEqual(
            audit["components"][0]["adjacent_face_input_object_mask_counts"], {"0x2": 4}
        )
        self.assertTrue(audit["components"][0]["simple_planar_polygon_proven"])

    def test_exact_planar_cycle_audit_distinguishes_simple_and_self_crossing(self):
        from fractions import Fraction as Q

        square = [(Q(0), Q(0)), (Q(2), Q(0)), (Q(2), Q(2)), (Q(0), Q(2))]
        bow_tie = [(Q(0), Q(0)), (Q(2), Q(2)), (Q(0), Q(2)), (Q(2), Q(0))]

        simple = RUNNER._audit_exact_planar_cycle(square)
        crossing = RUNNER._audit_exact_planar_cycle(bow_tie)

        self.assertEqual(simple["status"], "CERTIFIED_SIMPLE_POLYGON")
        self.assertEqual(simple["intersection_count"], 0)
        self.assertEqual(crossing["status"], "SELF_INTERSECTION_OR_TOUCH")
        self.assertGreater(crossing["intersection_count"], 0)

    def test_exact_mesh_point_parity_requires_closed_mesh_and_rejects_edge_rays(self):
        from fractions import Fraction as Q

        vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
        ))
        faces = np.asarray((
            (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
            (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)

        inside = RUNNER._classify_exact_mesh_point_parity(vertices, faces, (Q(1, 4), Q(1, 3), Q(1, 5)))
        outside = RUNNER._classify_exact_mesh_point_parity(vertices, faces, (Q(2), Q(1, 3), Q(1, 5)))
        edge_ray = RUNNER._classify_exact_mesh_point_parity(vertices, faces, (Q(1, 2), Q(1, 3), Q(1, 2)))
        open_mesh = RUNNER._classify_exact_mesh_point_parity(vertices, faces[:-2], (Q(1, 4), Q(1, 3), Q(1, 5)))

        self.assertEqual(inside["status"], "INSIDE_ODD_PARITY")
        self.assertEqual(inside["positive_ray_intersection_count"], 1)
        self.assertEqual(outside["status"], "OUTSIDE_EVEN_PARITY")
        self.assertEqual(outside["positive_ray_intersection_count"], 0)
        self.assertEqual(edge_ray["status"], "UNKNOWN_RAY_EDGE_OR_VERTEX_HIT")
        self.assertEqual(open_mesh["status"], "UNKNOWN_SOURCE_NOT_CLOSED_EDGE_MANIFOLD")

    def test_exact_point_in_polygon_classifies_interior_exterior_and_boundary(self):
        from fractions import Fraction as Q

        square = [(Q(0), Q(0)), (Q(2), Q(0)), (Q(2), Q(2)), (Q(0), Q(2))]

        self.assertEqual(RUNNER._exact_point_in_polygon(square, (Q(1), Q(1))), "INSIDE")
        self.assertEqual(RUNNER._exact_point_in_polygon(square, (Q(3), Q(1))), "OUTSIDE")
        self.assertEqual(RUNNER._exact_point_in_polygon(square, (Q(0), Q(1))), "ON_BOUNDARY")

    def test_exact_mesh_plane_section_returns_closed_cycles_or_unknown(self):
        from fractions import Fraction as Q

        vertices = np.asarray((
            (0., -2., 0.), (-1., 0., -1.), (1., 0., -1.), (0., 0., 1.),
        ))
        faces = np.asarray(((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)), dtype=np.int64)

        section = RUNNER._exact_mesh_plane_section_cycles(vertices, faces, axis="y", plane=Q(-1))
        open_section = RUNNER._exact_mesh_plane_section_cycles(vertices, faces[:-1], axis="y", plane=Q(-1))

        self.assertEqual(section["status"], "CLOSED_SECTION_CYCLES")
        self.assertEqual(section["cycle_count"], 1)
        self.assertEqual(len(section["cycles"][0]), 3)
        self.assertEqual(open_section["status"], "UNKNOWN_NON_CYCLE_SECTION")

    def test_exact_cycle_subdivision_matches_only_collinear_boundary_refinement(self):
        from fractions import Fraction as Q

        source = [(Q(0), Q(0), Q(0)), (Q(2), Q(0), Q(0)),
                  (Q(2), Q(0), Q(2)), (Q(0), Q(0), Q(2))]
        refined = [(Q(0), Q(0), Q(0)), (Q(1), Q(0), Q(0)),
                   (Q(2), Q(0), Q(0)), (Q(2), Q(0), Q(2)),
                   (Q(0), Q(0), Q(2))]
        detoured = [(Q(0), Q(0), Q(0)), (Q(1), Q(0), Q(1, 3)),
                    (Q(2), Q(0), Q(0)), (Q(2), Q(0), Q(2)),
                    (Q(0), Q(0), Q(2))]

        matched = RUNNER._exact_cycle_matches_up_to_collinear_subdivision(
            refined, source, axis="y")
        rejected = RUNNER._exact_cycle_matches_up_to_collinear_subdivision(
            detoured, source, axis="y")

        self.assertTrue(matched["matches"])
        self.assertEqual(matched["reference_segment_count"], 4)
        self.assertEqual(matched["extra_collinear_vertices"], 1)
        self.assertFalse(rejected["matches"])

    def test_exact_planar_fan_requires_a_rounded_kernel_point(self):
        from fractions import Fraction as Q

        convex = [(Q(0), Q(0), Q(0)), (Q(2), Q(0), Q(0)),
                  (Q(2), Q(0), Q(2)), (Q(0), Q(0), Q(2))]
        concave_star = [(Q(0), Q(0), Q(0)), (Q(4), Q(0), Q(0)),
                        (Q(4), Q(0), Q(4)), (Q(2), Q(0), Q(39, 10)),
                        (Q(0), Q(0), Q(4))]

        square_fan = RUNNER._triangulate_planar_loop_from_kernel_point(
            convex, axis="y", plane=Q(0), desired_normal_sign=-1)
        star_fan = RUNNER._triangulate_planar_loop_from_kernel_point(
            concave_star, axis="y", plane=Q(0), desired_normal_sign=-1)

        self.assertEqual(square_fan["status"], "CERTIFIED_KERNEL_FAN")
        self.assertEqual(len(square_fan["triangles"]), 4)
        self.assertEqual(star_fan["status"], "CERTIFIED_KERNEL_FAN")
        self.assertEqual(len(star_fan["triangles"]), 5)

    def test_exact_planar_fan_compacts_large_rational_centroid_evidence(self):
        from fractions import Fraction as Q

        tiny_offset = Q(1, 1 << 15000)
        loop = [(tiny_offset, tiny_offset, Q(0)),
                (tiny_offset + 2, tiny_offset, Q(0)),
                (tiny_offset + 2, tiny_offset + 2, Q(0)),
                (tiny_offset, tiny_offset + 2, Q(0))]

        result = RUNNER._triangulate_planar_loop_from_kernel_point(
            loop, axis="z", plane=Q(0), desired_normal_sign=1,
        )

        self.assertEqual(result["status"], "CERTIFIED_KERNEL_FAN")
        self.assertEqual(result["center_binary64"], [1.0, 1.0, 0.0])
        self.assertEqual(len(result["center_exact_sha256"]), 64)

    def test_diagnostic_cap_candidate_closes_only_witness_supported_planar_boundary(self):
        vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
        ))
        # An outward cube without its z=0 cap.  The patch audit is deliberately
        # explicit: this helper builds a diagnostic candidate, not a Boolean.
        faces = np.asarray((
            (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
            (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6),
            (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)
        patch_audit = {
            "status": "WITNESS_SUPPORTED_INTERIOR_DISK",
            "selected_patch_candidate": "INTERIOR_DISK",
            "axis": "z",
            "plane_coordinate_exact": "0/1",
            "output_boundary_matches_unique_source_section": True,
        }

        candidate = RUNNER._build_diagnostic_cap_candidate(
            vertices, faces, patch_audit=patch_audit,
            source_self_intersection_status="INTERSECTION_FREE",
            density=1.0, root_side=1.0,
        )

        self.assertEqual(candidate["status"], "CAP_CANDIDATE_GEOMETRY_AND_MOMENTS_PASSED")
        self.assertTrue(candidate["exact_boundary_audit"]["self_intersection_free"])
        self.assertEqual(candidate["topology_audit"]["status"], "CLOSED_EDGE_MANIFOLD")
        self.assertEqual(candidate["topology_audit"]["boundary_edge_count"], 0)
        self.assertEqual(candidate["topology_audit"]["orientation_conflict_edge_count"], 0)
        self.assertFalse(candidate["boolean_equivalence_proven"])
        self.assertEqual(candidate["source_self_intersection_status"], "INTERSECTION_FREE")
        self.assertEqual(candidate["independent_moment_audit"]["status"], "PASSED")
        expected = np.asarray((1., .5, .5, .5, 1/3, 1/3, 1/3, .25, .25, .25))
        np.testing.assert_allclose(candidate["evaluator_raw_moments"], expected, atol=1e-14)
        np.testing.assert_allclose(candidate["independent_raw_moments"], expected, atol=1e-14)

    def test_diagnostic_cap_candidate_restores_exact_sidecar_boundary_coordinates(self):
        from fractions import Fraction as Q

        vertices = np.asarray((
            (0., 0., 1e-16), (1., 0., 1e-16), (1., 1., 1e-16), (0., 1., 1e-16),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
        ))
        faces = np.asarray((
            (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
            (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6),
            (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)
        patch_audit = {
            "status": "WITNESS_SUPPORTED_INTERIOR_DISK", "selected_patch_candidate": "INTERIOR_DISK",
            "axis": "z", "plane_coordinate_exact": "0/1",
            "output_boundary_matches_unique_source_section": True,
        }
        exact_boundary = {
            index: (Q.from_float(float(vertices[index, 0])),
                    Q.from_float(float(vertices[index, 1])), Q(0))
            for index in range(4)
        }

        candidate = RUNNER._build_diagnostic_cap_candidate(
            vertices, faces, patch_audit=patch_audit,
            source_self_intersection_status="INTERSECTION_FREE",
            exact_boundary_coordinates=exact_boundary,
            density=1.0, root_side=1.0,
        )

        self.assertEqual(candidate["status"], "CAP_CANDIDATE_GEOMETRY_AND_MOMENTS_PASSED")
        self.assertTrue(candidate["exact_boundary_audit"]["self_intersection_free"])
        self.assertEqual(candidate["exact_boundary_coordinate_replacements"], 4)
        self.assertGreater(candidate["maximum_boundary_coordinate_displacement"], 0.0)
        self.assertTrue(np.all(candidate["vertices"][:4, 2] == 0.0))

    def test_diagnostic_cap_candidate_requires_exact_self_intersection_audit(self):
        vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
        ))
        faces = np.asarray((
            (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
            (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6),
            (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)
        patch_audit = {
            "status": "WITNESS_SUPPORTED_INTERIOR_DISK",
            "selected_patch_candidate": "INTERIOR_DISK",
            "axis": "z", "plane_coordinate_exact": "0/1",
            "output_boundary_matches_unique_source_section": True,
        }

        candidate = RUNNER._build_diagnostic_cap_candidate(
            vertices, faces, patch_audit=patch_audit,
            source_self_intersection_status="INTERSECTION_FREE",
            density=1.0, root_side=1.0,
        )

        self.assertEqual(candidate["status"], "CAP_CANDIDATE_GEOMETRY_AND_MOMENTS_PASSED")
        self.assertEqual(candidate["candidate_self_intersection_status"], "INTERSECTION_FREE")
        self.assertEqual(candidate["candidate_exact_audit_face_count"], len(candidate["faces"]))
        self.assertEqual(candidate["boolean_equivalence_proven"], False)

    def test_exact_candidate_audit_rejects_intersecting_closed_components(self):
        native = RUNNER.spin._native_module()
        cube_vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
        ))
        cube_faces = np.asarray((
            (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
            (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)
        vertices = np.vstack((cube_vertices, cube_vertices + 0.5))
        faces = np.vstack((cube_faces, cube_faces + len(cube_vertices)))

        audit = native.audit_exact_triangle_boundary(vertices, faces, "AnyNonzero")

        self.assertFalse(audit["self_intersection_free"])
        self.assertGreater(audit["intersecting_pair_count"], 0)
        self.assertEqual(len(set(audit["face_components"])), 2)

    def test_diagnostic_cap_candidate_fails_closed_without_exact_source_clearance(self):
        vertices = np.asarray(((0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
                               (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.)))
        faces = np.asarray(((4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
                            (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6),
                            (3, 0, 4), (3, 4, 7)), dtype=np.int64)
        patch_audit = {
            "status": "WITNESS_SUPPORTED_INTERIOR_DISK", "selected_patch_candidate": "INTERIOR_DISK",
            "axis": "z", "plane_coordinate_exact": "0/1",
            "output_boundary_matches_unique_source_section": True,
        }

        candidate = RUNNER._build_diagnostic_cap_candidate(
            vertices, faces, patch_audit=patch_audit,
            source_self_intersection_status="UNKNOWN",
            density=1.0, root_side=1.0,
        )

        self.assertEqual(candidate["status"], "UNKNOWN_SOURCE_SELF_INTERSECTION")
        self.assertIsNone(candidate["vertices"])
        self.assertIsNone(candidate["faces"])

    def test_diagnostic_cap_candidate_does_not_mistake_face_complement_for_disk(self):
        vertices = np.asarray(((0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
                               (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.)))
        faces = np.asarray(((4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
                            (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6),
                            (3, 0, 4), (3, 4, 7)), dtype=np.int64)
        patch_audit = {
            "status": "WITNESS_SUPPORTED_FACE_COMPLEMENT", "selected_patch_candidate": "FACE_COMPLEMENT",
            "axis": "z", "plane_coordinate_exact": "0/1",
            "output_boundary_matches_unique_source_section": True,
        }

        candidate = RUNNER._build_diagnostic_cap_candidate(
            vertices, faces, patch_audit=patch_audit,
            source_self_intersection_status="INTERSECTION_FREE",
            density=1.0, root_side=1.0,
        )

        self.assertEqual(candidate["status"], "UNKNOWN_COMPLEMENT_PATCH_REQUIRES_HOLE_TRIANGULATION")
        self.assertIsNone(candidate["vertices"])
        self.assertIsNone(candidate["faces"])

    def test_exact_radial_annulus_triangulation_covers_square_minus_diamond(self):
        from fractions import Fraction as Q

        outer = [(Q(-2), Q(-2)), (Q(2), Q(-2)), (Q(2), Q(2)), (Q(-2), Q(2))]
        inner = [(Q(1), Q(0)), (Q(0), Q(1)), (Q(-1), Q(0)), (Q(0), Q(-1))]
        result = RUNNER._triangulate_exact_star_shaped_polygon_difference(
            outer, inner, kernel_point=(Q(0), Q(0)),
        )

        self.assertEqual(result["status"], "CERTIFIED_RADIAL_ANNULUS")
        self.assertEqual(result["triangle_count"], 16)
        self.assertEqual(result["exact_area_encoding"], "14/1")
        self.assertEqual(result["boundary_edge_count"], 16)

    def test_exact_radial_annulus_compacts_large_rational_area_evidence(self):
        from fractions import Fraction as Q

        tiny = Q(1, 1 << 15000)
        outer = [(Q(0), Q(0)), (Q(2) + tiny, Q(0)),
                 (Q(2) + tiny, Q(2)), (Q(0), Q(2))]
        inner = [(Q(3, 2), Q(1)), (Q(1), Q(3, 2)),
                 (Q(1, 2), Q(1)), (Q(1), Q(1, 2))]
        result = RUNNER._triangulate_exact_star_shaped_polygon_difference(
            outer, inner, kernel_point=(Q(1), Q(1)),
        )

        self.assertEqual(result["status"], "CERTIFIED_RADIAL_ANNULUS")
        self.assertTrue(result["exact_area_encoding"].startswith("sha256:"))
        self.assertAlmostEqual(result["exact_area_float"], 3.5, places=12)

    def test_exact_radial_annulus_rejects_non_star_shaped_inner_loop(self):
        from fractions import Fraction as Q

        outer = [(Q(-1), Q(-1)), (Q(4), Q(-1)), (Q(4), Q(4)), (Q(-1), Q(4))]
        inner = [(Q(0), Q(0)), (Q(3), Q(0)), (Q(3), Q(1)), (Q(1), Q(1)),
                 (Q(1), Q(2)), (Q(3), Q(2)), (Q(3), Q(3)), (Q(0), Q(3))]
        result = RUNNER._triangulate_exact_star_shaped_polygon_difference(
            outer, inner, kernel_point=(Q(3, 2), Q(3, 2)),
        )

        self.assertNotEqual(result["status"], "CERTIFIED_RADIAL_ANNULUS")
        self.assertEqual(result["triangles"], [])

    def test_boolean_face_patch_selection_requires_separating_source_witnesses(self):
        self.assertEqual(
            RUNNER._classify_boolean_face_patch("intersection", "INSIDE_ODD_PARITY", "OUTSIDE_EVEN_PARITY"),
            "INTERIOR_DISK",
        )
        self.assertEqual(
            RUNNER._classify_boolean_face_patch("subtraction", "OUTSIDE_EVEN_PARITY", "INSIDE_ODD_PARITY"),
            "INTERIOR_DISK",
        )
        self.assertEqual(
            RUNNER._classify_boolean_face_patch("intersection", "OUTSIDE_EVEN_PARITY", "INSIDE_ODD_PARITY"),
            "FACE_COMPLEMENT",
        )
        self.assertEqual(
            RUNNER._classify_boolean_face_patch("intersection", "INSIDE_ODD_PARITY", "INSIDE_ODD_PARITY"),
            "UNKNOWN_FACE_REGION",
        )

    def test_planar_patch_audit_finds_intersection_disk_from_mesh_witnesses(self):
        from fractions import Fraction as Q

        design_vertices = np.asarray((
            (-2., -1., -2.), (2., -1., -2.), (2., -1., 2.), (-2., -1., 2.),
            (-2., 1., -2.), (2., 1., -2.), (2., 1., 2.), (-2., 1., 2.),
        ))
        design_faces = np.asarray((
            (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
            (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)
        source_vertices = np.asarray((
            (0., -2., 0.), (-1., 0., -1.), (1., 0., -1.), (0., 0., 1.),
        ))
        source_faces = np.asarray(((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)), dtype=np.int64)
        section = [(Q(-1, 2), Q(-1, 1), Q(-1, 2)),
                   (Q(1, 2), Q(-1, 1), Q(-1, 2)),
                   (Q(0), Q(-1, 1), Q(1, 2))]

        audit = RUNNER._audit_planar_boolean_patch_region(
            section, axis="y", design_vertices=design_vertices, design_faces=design_faces,
            source_vertices=source_vertices, source_faces=source_faces, operation="intersection",
        )

        self.assertEqual(audit["status"], "WITNESS_SUPPORTED_INTERIOR_DISK")
        self.assertEqual(audit["selected_patch_candidate"], "INTERIOR_DISK")
        self.assertEqual(audit["loop_interior_source_classification"]["status"], "INSIDE_ODD_PARITY")
        self.assertEqual(audit["face_complement_source_classification"]["status"], "OUTSIDE_EVEN_PARITY")
        self.assertEqual(audit["source_plane_section_audit"]["status"], "CLOSED_SECTION_CYCLES")
        self.assertTrue(audit["output_boundary_matches_unique_source_section"])
        self.assertFalse(audit["cap_generated"])

    def test_level5_is_explicitly_passed_to_source_and_native_evaluators(self):
        source_settings, evaluator_settings = RUNNER._make_triangulation_configuration(
            _Extension, icosphere_subdivisions=5
        )

        self.assertEqual(source_settings.icosphere_subdivisions, 5)
        self.assertEqual(source_settings.maximum_icosphere_subdivisions, -1)
        self.assertEqual(source_settings.maximum_radius_doublings, 6)
        self.assertEqual(evaluator_settings.icosphere_subdivisions, 5)
        self.assertEqual(evaluator_settings.maximum_icosphere_subdivisions, -1)
        self.assertEqual(evaluator_settings.maximum_radius_doublings, 6)
        self.assertEqual(
            evaluator_settings.derivative_method, "DISCRETE_ELLIPSOID"
        )

    def test_cap_radius_doubling_budget_is_explicit_and_bound_to_both_evaluators(self):
        source_settings, evaluator_settings = RUNNER._make_triangulation_configuration(
            _Extension, icosphere_subdivisions=4, maximum_radius_doublings=8
        )

        self.assertEqual(source_settings.maximum_radius_doublings, 8)
        self.assertEqual(evaluator_settings.maximum_radius_doublings, 8)

    def test_cap_radius_doubling_budget_rejects_unbounded_or_negative_values(self):
        for invalid in (-1, 17, True):
            with self.subTest(invalid=invalid):
                with self.assertRaisesRegex(ValueError, "radius doublings"):
                    RUNNER._make_triangulation_configuration(
                        _Extension, icosphere_subdivisions=4,
                        maximum_radius_doublings=invalid,
                    )

    def test_strict_profile_defaults_to_certified_cap_budget_and_accepts_override(self):
        default = RUNNER._argument_parser().parse_args((
            "--output", "run", "--build", "build", "--cherchi-probe", "probe",
        ))
        refined = RUNNER._argument_parser().parse_args((
            "--output", "run", "--build", "build", "--cherchi-probe", "probe",
            "--maximum-radius-doublings", "8",
        ))

        self.assertEqual(default.maximum_radius_doublings, 6)
        self.assertEqual(refined.maximum_radius_doublings, 8)

    def test_bounded_cap_refinement_is_explicit_and_cannot_lower_resolution(self):
        source_settings, evaluator_settings = RUNNER._make_triangulation_configuration(
            _Extension, icosphere_subdivisions=3,
            maximum_icosphere_subdivisions=5,
        )
        self.assertEqual(source_settings.icosphere_subdivisions, 3)
        self.assertEqual(source_settings.maximum_icosphere_subdivisions, 5)
        self.assertEqual(evaluator_settings.maximum_icosphere_subdivisions, 5)
        for maximum in (-2, 2, 8, True):
            with self.subTest(maximum=maximum):
                with self.assertRaises(ValueError):
                    RUNNER._make_triangulation_configuration(
                        _Extension, icosphere_subdivisions=3,
                        maximum_icosphere_subdivisions=maximum,
                    )

    def test_cli_records_requested_icosphere_level(self):
        args = RUNNER._argument_parser().parse_args((
            "--output", "run",
            "--build", "build",
            "--cherchi-probe", "probe",
            "--icosphere-subdivisions", "5",
            "--maximum-function-evaluations", "1",
        ))

        self.assertEqual(args.icosphere_subdivisions, 5)
        self.assertEqual(args.maximum_icosphere_subdivisions, -1)
        self.assertEqual(args.maximum_function_evaluations, 1)

        adaptive = RUNNER._argument_parser().parse_args((
            "--output", "run", "--build", "build", "--cherchi-probe", "probe",
            "--icosphere-subdivisions", "3",
            "--maximum-icosphere-subdivisions", "5",
        ))
        self.assertEqual(adaptive.icosphere_subdivisions, 3)
        self.assertEqual(adaptive.maximum_icosphere_subdivisions, 5)

    def test_cherchi_boolean_timeout_is_an_explicit_resource_setting(self):
        args = RUNNER._argument_parser().parse_args((
            "--output", "run",
            "--build", "build",
            "--cherchi-probe", "probe",
        ))
        self.assertEqual(args.cherchi_timeout_seconds, 180.0)

        args = RUNNER._argument_parser().parse_args((
            "--output", "run",
            "--build", "build",
            "--cherchi-probe", "probe",
            "--cherchi-timeout-seconds", "240",
        ))
        self.assertEqual(args.cherchi_timeout_seconds, 240.0)

    def test_configuration_rejects_unsupported_levels(self):
        for level in (-1, 8):
            with self.subTest(level=level):
                with self.assertRaises(ValueError):
                    RUNNER._make_triangulation_configuration(
                        _Extension,
                        icosphere_subdivisions=level,
                    )

    def test_captured_one_sheet_trial_stays_on_negative_subtraction_branch(self):
        coefficients = (
            -0.1998242186450357, 0.02382394696275607,
            0.005243276046278874, 0.0015024879022398209,
            0.5734432875092278, 0.7938771401867838,
            -0.0000336431262224418, 0.0007933293534953297,
            -0.0019952147391396123, 0.01975405514590195,
        )

        self.assertEqual(
            RUNNER._supported_quadric_operation(coefficients), "subtraction"
        )

    def test_captured_two_sheet_trial_uses_center_positive_intersection_branch(self):
        coefficients = (
            -0.00010189795769451712, 0.0011279444835807586,
            0.0010546441324340969, -0.02356662364217407,
            0.021837996108963594, 0.021768397637837767,
            -0.9988986803954175, 0.018233537792266833,
            0.018955367064424793, -0.0010566870211034748,
        )

        self.assertEqual(
            RUNNER._supported_quadric_operation(coefficients), "intersection"
        )

    def test_captured_near_cone_trial_uses_exact_nonzero_center_sign(self):
        coefficients = (
            9.383970775630317e-13, -1.2425364931169101e-8,
            1.1871684596774663e-8, -9.086419627621456e-9,
            0.4031536826392824, 0.40315352542109545,
            -0.8215438777781616, -3.8394933782711187e-7,
            1.8571867393924378e-8, -3.025244045465824e-10,
        )

        self.assertEqual(
            RUNNER._exact_quadric_center_level_sign(coefficients), 1
        )
        self.assertEqual(
            RUNNER._supported_quadric_operation(coefficients), "intersection"
        )

    def test_exact_cone_center_level_remains_rejected(self):
        cone = (0., 0., 0., 0., 1., 1., -1., 0., 0., 0.)

        self.assertEqual(RUNNER._exact_quadric_center_level_sign(cone), 0)
        with self.assertRaisesRegex(ValueError, "exactly zero center level"):
            RUNNER._supported_quadric_operation(cone)

    @staticmethod
    def _cube_with_interface_faces(interface_faces):
        vertices = np.asarray((
            (-1., -1., -1.), (1., -1., -1.), (1., 1., -1.), (-1., 1., -1.),
            (-1., -1., 1.), (1., -1., 1.), (1., 1., 1.), (-1., 1., 1.),
        ))
        faces = np.asarray((
            (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
            (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
        ), dtype=np.int64)
        caps = np.ones(len(faces), dtype=np.int8)
        caps[list(interface_faces)] = 0
        return vertices, faces, caps

    def test_boolean_side_uses_oriented_interface_not_centroid(self):
        coefficients = (1., 1., 0., 0., 0., 0., 0., 0., 0., 0.)
        vertices, faces, caps = self._cube_with_interface_faces((8, 9))
        details = RUNNER._closed_source_selector_details(
            coefficients, vertices, faces, caps,
        )
        self.assertEqual(details["operation"], "intersection")

    def test_strict_runner_requires_native_side_witness(self):
        selected = RUNNER._source_side_selection_details({
            "defined": True,
            "side_selection_defined": True,
            "interior_is_positive": False,
            "side_selection_method": "interior-centroid-winding-or-oriented-interface-v1",
        })
        self.assertEqual(selected["operation"], "subtraction")
        self.assertEqual(selected["method"], "interior-centroid-winding-or-oriented-interface-v1")
        with self.assertRaisesRegex(ValueError, "native quadric-side witness"):
            RUNNER._source_side_selection_details({
                "defined": True,
                "side_selection_defined": False,
                "side_selection_diagnostic": "mixed oriented interface signs",
            })

    def test_closed_source_selector_reports_value_scale_and_threshold_margin(self):
        vertices, faces, caps = self._cube_with_interface_faces((8, 9))
        details = RUNNER._closed_source_selector_details(
            (1., 1., 0., 0., 0., 0., 0., 0., 0., 0.),
            vertices, faces, caps,
        )
        self.assertEqual(details["operation"], "intersection")
        self.assertEqual(details["negative_exterior_face_count"], 2)
        self.assertGreater(details["minimum_normalized_alignment"], 0.99)
        self.assertAlmostEqual(details["threshold"], 1e-10, places=24)

    def test_closed_source_selector_reports_negative_branch_and_rejects_small_margin(self):
        vertices, faces, caps = self._cube_with_interface_faces((10, 11))
        details = RUNNER._closed_source_selector_details(
            (-1., 1., 0., 0., 0., 0., 0., 0., 0., 0.),
            vertices, faces, caps,
        )
        self.assertEqual(details["operation"], "subtraction")
        self.assertFalse(details["interior_is_positive"])

    def test_closed_source_side_rejects_ambiguous_or_inconsistent_interface(self):
        vertices, faces, caps = self._cube_with_interface_faces((8, 9))
        with self.assertRaisesRegex(ValueError, "ambiguous"):
            RUNNER._closed_source_selector_details(
                (0., 0., 1., 0., 0., 0., 0., 0., 0., 0.),
                vertices, faces, caps,
            )

        caps[10] = 0
        caps[11] = 0
        with self.assertRaisesRegex(ValueError, "disagree"):
            RUNNER._closed_source_selector_details(
                (1., 1., 0., 0., 0., 0., 0., 0., 0., 0.),
                vertices, faces, caps,
            )

        caps[8] = 0
        caps[9] = 0
        with self.assertRaisesRegex(ValueError, "disagree"):
            RUNNER._closed_source_selector_details(
                (1., 1., 0., 0., 0., 0., 0., 0., 0., 0.),
                vertices, faces, caps,
            )

    def test_failed_source_motion_and_center_stratum_trials_have_distinct_stages(self):
        coefficients = np.zeros(10)

        source_motion = RUNNER._failed_evaluation_record(
            1, coefficients,
            ValueError("quadric source motion unavailable: fixed-normal bordered differential solve did not certify a finite 1e-10 relative residual"),
        )
        center_stratum = RUNNER._failed_evaluation_record(
            2, coefficients,
            ValueError("quadric has exactly zero center level (singular cone)"),
        )

        self.assertEqual(source_motion["failure_stage"], "source_motion_differential")
        self.assertEqual(center_stratum["failure_stage"], "quadric_stratum_classification")

    def test_material_topology_identity_binds_boolean_and_cap_branches(self):
        construction = {"records": [
            {"vertex_index": i, "point_type": "SOURCE_VERTEX",
             "operands": [{"status": "UNIQUE", "source_vertices": [
                 {"input_object": 0, "input_vertex": i}]}],
             "surface_planes": []}
            for i in range(3)
        ]}
        faces = np.asarray(((0, 1, 2),), dtype=np.int64)
        cap_branch = {
            "surface_type": "ONE_SHEET_HYPERBOLOID",
            "truncation_radius_hex": float(2.0).hex(),
            "initial_icosphere_subdivisions": 4,
            "effective_icosphere_subdivisions": 4,
        }
        subtraction = RUNNER._material_topology_identity(
            "subtraction", construction, faces, derivative_branch=cap_branch,
        )
        intersection = RUNNER._material_topology_identity(
            "intersection", construction, faces, derivative_branch=cap_branch,
        )
        self.assertEqual(len(subtraction), 64)
        self.assertEqual(subtraction, RUNNER._material_topology_identity(
            "subtraction", construction, faces, derivative_branch=cap_branch,
        ))
        self.assertNotEqual(subtraction, intersection)
        larger_cap = {**cap_branch, "truncation_radius_hex": float(4.0).hex()}
        self.assertNotEqual(subtraction, RUNNER._material_topology_identity(
            "subtraction", construction, faces, derivative_branch=larger_cap,
        ))
        refined = {**cap_branch, "effective_icosphere_subdivisions": 5}
        self.assertNotEqual(subtraction, RUNNER._material_topology_identity(
            "subtraction", construction, faces, derivative_branch=refined,
        ))

    def test_rank_deficient_trial_is_rejected_before_geometry_callback(self):
        plane = (0., 1., 0., 0., 0., 0., 0., 0., 0., 0.)

        with self.assertRaisesRegex(ValueError, "unsupported quadratic signature"):
            RUNNER._supported_quadric_operation(plane)

    def test_failed_evaluation_record_binds_exception_to_index_and_coefficients(self):
        coefficients = (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0)
        with tempfile.TemporaryDirectory() as temporary:
            artifact_dir = Path(temporary)
            source_bytes = b"v 0 0 0\n"
            material_bytes = b"v 1 0 0\n"
            construction_bytes = b"face,point,type\n"
            (artifact_dir / "source-017.obj").write_bytes(source_bytes)
            (artifact_dir / "material-017.obj").write_bytes(material_bytes)
            (artifact_dir / "material-017.obj.construction.tsv").write_bytes(construction_bytes)
            record = RUNNER._failed_evaluation_record(
                17, coefficients,
                ValueError("orientation input contains a degenerate face"),
                artifact_directory=artifact_dir,
            )

        self.assertEqual(record["evaluation_index"], 17)
        self.assertEqual(record["quadric_coefficients"], list(coefficients))
        self.assertEqual(record["outcome"], "failed")
        self.assertEqual(record["failure_type"], "ValueError")
        self.assertEqual(
            record["failure_message"], "orientation input contains a degenerate face"
        )
        self.assertEqual(record["failure_stage"], "orientation_normalization")
        artifacts = record["available_artifacts"]
        self.assertEqual(artifacts["source_geometry"]["path"], "source-017.obj")
        self.assertEqual(
            artifacts["source_geometry"]["sha256"], hashlib.sha256(source_bytes).hexdigest()
        )
        self.assertEqual(artifacts["raw_boolean_geometry"]["path"], "material-017.obj")
        self.assertEqual(
            artifacts["raw_boolean_geometry"]["sha256"],
            hashlib.sha256(material_bytes).hexdigest(),
        )
        self.assertEqual(
            artifacts["construction_provenance"]["sha256"],
            hashlib.sha256(construction_bytes).hexdigest(),
        )

    def test_failed_evaluation_record_preserves_known_boolean_operation(self):
        record = RUNNER._failed_evaluation_record(
            13, np.zeros(10), RuntimeError("orientation normalization failed"),
            operation="intersection",
        )

        self.assertEqual(record["boolean_operation"], "intersection")

    def test_failed_evaluation_record_without_artifact_directory_is_safe(self):
        record = RUNNER._failed_evaluation_record(
            3, np.zeros(10), RuntimeError("no geometry artifacts were written"),
        )

        self.assertEqual(record["failure_stage"], "source_quadric_triangulation")
        self.assertEqual(record["available_artifacts"], {})
        self.assertNotIn("boolean_operation", record)
        self.assertIsNone(record["raw_boolean_topology_diagnostic"])
        self.assertIsNone(record["exact_boolean_boundary_provenance"])


if __name__ == "__main__":
    unittest.main()
