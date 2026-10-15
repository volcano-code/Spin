"""Independent analytic and file-backed tests for the metric-only checker."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / 'experiments/teapot_spin2014/effect_audit.py'
# A missing implementation is an assertion failure, not a collection error.
spec = importlib.util.spec_from_file_location('teapot_effect_audit', MODULE) if MODULE.exists() else None
audit = importlib.util.module_from_spec(spec) if spec else None
if spec:
    spec.loader.exec_module(audit)

# Tetrahedron with axis intercepts 1, 2, 3. Volume = 1; all cross terms differ.
TETRA = [1, .25, .5, .75, .1, .4, .9, .3, .15, .1]
OBJ = 'v 0 0 0\nv 1 0 0\nv 0 2 0\nv 0 0 3\nf 1 3 2\nf 1 2 4\nf 1 4 3\nf 2 3 4\n'


def dump(path, value):
    path.write_text(json.dumps(value, allow_nan=False))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_fixture(base):
    """Tiny files with the actual stream_probe manifest/header schema."""
    inputs = base / 'inputs'
    inputs.mkdir()
    (inputs / 'outer.obj').write_text(OBJ)
    (inputs / 'design.obj').write_text(OBJ)
    geometry = {str(inputs / name): sha(inputs / name) for name in ('outer.obj', 'design.obj')}
    dump(inputs / 'manifest.json', {'geometry_hashes': geometry})
    reference = base / 'reference.json'
    cfg = {'input': str(inputs), 'meters_per_model_unit': .1, 'density_kg_m3': 1000,
           'geometry': {'root_side': 3, 'root_origin': [0, 0, 0]},
           'solver': {'initial_fill': .25}, 'boundary': {}, 'minimum_wall': .01}
    dump(reference, cfg)
    settings = {'source_settings': cfg['geometry'],
                'solver_settings': {**cfg['solver'], 'termination_policy': 'PAPER_POST_TRANSITION_V1'},
                'boundary_settings': {}, 'density': 1., 'minimum_outer_wall': .01}
    payload = {'settings': settings,
               'source': {'outer_surface_moments': TETRA,
                          'carveable_moments': [v * .5 for v in TETRA]},
               'solve': {'status': 'RELAXATION_UNKNOWN', 'trace_hash': 'trace', 'output_hash': 'output',
                         'final_partition': [{'level': 0, 'coordinates': [0, 0, 0],
                                              'moments': [v * .5 for v in TETRA]}],
                         'relaxed_fill': [.5], 'binary_fill': [], 'binary_moments': None,
                         'binary_check': None}}
    run = base / 'run'
    run.mkdir()
    dump(run / 'header.json', {'kind': 'header', 'trace_counts': [1], 'payload': payload})
    dump(run / 'config.json', {**settings, 'reference_config': str(reference),
                              'reference_sha256': sha(reference), 'reference_settings': cfg})
    tracked = [reference, inputs / 'manifest.json', inputs / 'outer.obj', inputs / 'design.obj']
    dump(run / 'source.json', {'files': {str(p): sha(p) for p in tracked}, 'commit': 'fixture'})
    dump(run / 'summary.json', {k: payload['solve'][k] for k in ('status', 'trace_hash', 'output_hash', 'binary_check')})
    # The fake digest is deliberately invalid: this checker must not claim to check the trace.
    refresh_outputs(run)
    return run


def refresh_outputs(run):
    files = {str(run / name): sha(run / name)
             for name in ('header.json', 'config.json', 'source.json', 'summary.json')}
    files[str(run / 'native-run.jsonl')] = 'separate-parent-owned-trace-gate'
    dump(run / 'outputs.json', files)


class EffectAuditTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(audit, 'missing independent effect checker implementation')

    def test_tetrahedron_cross_terms_parallel_axis_and_si(self):
        result = audit.moment_metrics(TETRA, .1, 1000)
        native, si = result['model_units'], result['si']
        self.assertAlmostEqual(native['mass'], 1.)
        np.testing.assert_allclose(native['com'], [.25, .5, .75])
        np.testing.assert_allclose(native['pivot_inertia'], [[1.3, -.1, -.15], [-.1, 1., -.3], [-.15, -.3, .5]])
        np.testing.assert_allclose(native['central_inertia'], [[.4875, .025, .0375], [.025, .375, .075], [.0375, .075, .1875]])
        self.assertAlmostEqual(si['mass'], 1.)
        np.testing.assert_allclose(si['com'], [.025, .05, .075])
        self.assertAlmostEqual(si['central_inertia'][0][0], .004875)
        self.assertAlmostEqual(native['com_to_axis_distance'], np.sqrt(.3125))

    def test_positive_pivot_gap_uses_pivot_not_central_denominator(self):
        # Mass=1, COM=(0,0,1), central I=(4,5,7), pivot I=(5,6,7).
        result = audit.moment_metrics([1, 0, 0, 1, 4, 3, 2, 0, 0, 0], .1, 1000)
        self.assertEqual(result['fixed_point_bound']['status'], 'CONDITIONS_SATISFIED_FOR_FORMAL_BOUND')
        self.assertAlmostEqual(result['fixed_point_bound']['omega_rad_s'], np.sqrt(98.0665))
        self.assertAlmostEqual(result['formal_com_paper_score']['omega_rad_s'], np.sqrt(49.03325))
        self.assertAlmostEqual(result['model_units']['central_lateral_squared_ratio_sum'], 41 / 49)

    def test_nonpositive_pivot_gap_is_not_applicable_not_infeasible(self):
        # Translate previous shape to z=2; central gap stays positive, pivot gap becomes -2.
        result = audit.moment_metrics([1, 0, 0, 2, 4, 3, 5, 0, 0, 0], .1, 1000)
        self.assertEqual(result['fixed_point_bound']['status'], 'NOT_APPLICABLE')
        self.assertIn('NONPOSITIVE_PIVOT_GAP', result['fixed_point_bound']['reasons'])
        self.assertIsNone(result['fixed_point_bound']['omega_rad_s'])
        self.assertGreater(result['formal_com_paper_score']['omega_rad_s'], 0)

    def test_off_axis_or_coupled_body_never_admits_fixed_point_bound(self):
        for raw, reason in [([1, .1, 0, 1, 4.01, 3, 2, 0, .1, 0], 'COM_OFF_AXIS'),
                            ([1, 0, 0, 1, 4, 3, 2, .1, 0, 0], 'AXIS_NOT_PRINCIPAL')]:
            result = audit.moment_metrics(raw, 1, 1)
            self.assertEqual(result['fixed_point_bound']['status'], 'NOT_APPLICABLE')
            self.assertIn(reason, result['fixed_point_bound']['reasons'])
            self.assertIsNone(result['fixed_point_bound']['omega_rad_s'])

    def test_repeated_eigenvalues_report_subspace_without_arbitrary_axis(self):
        # Central eigenvalues 2, 3, 3, with target z inside the dominant plane.
        r = audit.moment_metrics([1, 0, 0, 0, 2, 1, 1, 0, 0, 0], 1, 1)
        self.assertEqual(r['model_units']['dominant_axis_status'], 'NONUNIQUE')
        self.assertIsNone(r['model_units']['dominant_axis'])
        self.assertIsNone(r['model_units']['dominant_axis_angle_degrees'])
        self.assertAlmostEqual(r['model_units']['dominant_subspace_angle_degrees'], 0)
        # Repeated lower eigenvalues also preclude the strict fixed-point condition.
        r = audit.moment_metrics([1, 0, 0, 1, 3, 3, 2, 0, 0, 0], 1, 1)
        self.assertIn('DEGENERATE_PIVOT_EIGENVALUES', r['fixed_point_bound']['reasons'])

    def test_rejects_invalid_moments_and_density_conventions(self):
        for index, value in [(0, 0), (0, -1), (4, float('nan')), (7, float('inf')), (3, None)]:
            raw = TETRA.copy()
            raw[index] = value
            with self.assertRaises(ValueError):
                audit.moment_metrics(raw, .1, 1000)
        for raw in [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0], TETRA[:9], [1, 0, 0, 0, -1, 1, 1, 0, 0, 0]]:
            with self.assertRaises(ValueError):
                audit.moment_metrics(raw, 1, 1)
        for scale, density, native in [(0, 1, 1), (1, float('inf'), 1), (1, 1000, 1000)]:
            with self.assertRaises(ValueError):
                audit.moment_metrics(TETRA, scale, density, native_density=native)

    def test_extracts_configured_initial_and_relaxed_states_without_binary_invention(self):
        with tempfile.TemporaryDirectory() as d:
            run = run_fixture(Path(d))
            result = audit.audit_run(run)
            self.assertEqual(result['native_status'], 'RELAXATION_UNKNOWN')
            self.assertEqual(result['scope'], 'METRIC_ONLY_NO_ADMISSION')
            self.assertEqual(result['full_trace_verification'], 'NOT_PERFORMED_SEPARATE_PARENT_GATE')
            self.assertIsNone(result['states']['native_binary'])
            self.assertAlmostEqual(result['states']['initial_mixture']['model_units']['mass'], .875)
            self.assertAlmostEqual(result['states']['terminal_relaxed']['model_units']['mass'], .75)
            self.assertLessEqual(result['checks']['original_normalized_error'], 1e-10)
            self.assertIn(str(MODULE), result['checker_files'])
            self.assertTrue(any(p.endswith('/mass_moments.py') for p in result['checker_files']))
            self.assertEqual(result['geometry_integration']['vertices_dtype'], 'float64')
            self.assertEqual(result['geometry_integration']['faces_dtype'], 'int64')

    def test_binary_requires_native_presence_and_matching_binary_beta(self):
        with tempfile.TemporaryDirectory() as d:
            run = run_fixture(Path(d))
            h = json.loads((run / 'header.json').read_text())
            s = h['payload']['solve']
            s.update(binary_fill=[1], binary_moments=[v * .5 for v in TETRA])
            dump(run / 'header.json', h)
            refresh_outputs(run)
            r = audit.audit_run(run)
            self.assertAlmostEqual(r['states']['native_binary']['model_units']['mass'], .5)
            self.assertEqual(r['native_status'], 'RELAXATION_UNKNOWN')
            for fill, moments in [([.5], s['binary_moments']), ([0], s['binary_moments']),
                                  ([1], None), ([], s['binary_moments'])]:
                bad = copy.deepcopy(h)
                bad['payload']['solve'].update(binary_fill=fill, binary_moments=moments)
                dump(run / 'header.json', bad)
                refresh_outputs(run)
                with self.assertRaises(ValueError):
                    audit.audit_run(run)

    def test_rejects_bad_beta_partition_original_and_config_density(self):
        for mutation in ('dimension', 'range', 'partition', 'original', 'density', 'null_beta'):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as d:
                run = run_fixture(Path(d))
                h = json.loads((run / 'header.json').read_text())
                if mutation == 'dimension': h['payload']['solve']['relaxed_fill'] = []
                if mutation == 'range': h['payload']['solve']['relaxed_fill'] = [1.01]
                if mutation == 'null_beta': h['payload']['solve']['relaxed_fill'] = [None]
                if mutation == 'partition': h['payload']['source']['carveable_moments'][9] += .01
                if mutation == 'original': h['payload']['source']['outer_surface_moments'][8] += .01
                if mutation == 'density': h['payload']['settings']['density'] = 1000
                dump(run / 'header.json', h)
                refresh_outputs(run)
                with self.assertRaises(ValueError): audit.audit_run(run)

    def test_stale_header_source_config_summary_and_upstream_rejected(self):
        for name in ('header.json', 'source.json', 'config.json', 'summary.json', '../inputs/outer.obj'):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as d:
                run = run_fixture(Path(d))
                with (run / name).open('a') as stream: stream.write('\n')
                with self.assertRaisesRegex(ValueError, 'STALE'):
                    audit.audit_run(run)

    def test_cli_creates_new_output_only_and_never_mutates_run(self):
        with tempfile.TemporaryDirectory() as d:
            run = run_fixture(Path(d))
            before = {p.name: sha(p) for p in run.iterdir()}
            output = Path(d) / 'audit'
            command = [sys.executable, str(MODULE), '--run', str(run), '--output', str(output)]
            p = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertTrue((output / 'effect.json').is_file())
            self.assertEqual(before, {p.name: sha(p) for p in run.iterdir()})
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)


if __name__ == '__main__':
    unittest.main()
