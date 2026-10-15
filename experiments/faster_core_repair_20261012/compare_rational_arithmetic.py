"""Compare kernels on saved actual Cherchi materials, with no new Boolean."""
from fractions import Fraction
from math import fsum
import json
from pathlib import Path
import sys
import time

import numpy as np
import _cherchi_rational_moments as arithmetic
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import exact_fraction
from make_it_stand.verification.cherchi_motion_differential import (
    recover_construction_vertex_motion, integrate_rational_polyhedral_raw_moments,
    differentiate_polyhedral_raw_moments)

HERE = Path(__file__).resolve().parent


def timed(operation):
    started = time.monotonic()
    value = operation()
    return value, time.monotonic() - started


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    fixed = HERE / 'projective-formal-pipeline-job-01/output'
    frozen = json.loads((fixed / 'registered-inputs.json').read_text())
    samples = [(name, fixed / 'fixed' / (name + '-material.json'), frozen['request'], frozen['settings'],
                dict(L=1., origin=[0., 0., 0.]), frozen['request']['design_vertices'])
               for name in ('ellipsoid', 'one_sheet', 'two_sheet')]
    run = HERE / 'core-benchmark-faces-12288-seed-0-job-01/output/run'
    config = json.loads((run / 'configuration.json').read_text())
    result = json.loads((run / 'result.json').read_text())
    physical = json.loads((HERE / 'core-benchmark-inputs-job-01/output/faces-12288.json').read_text())['request']['design_vertices']
    terminal = json.loads((run / 'terminal-material-reference.json').read_text())['evaluation_file']
    for name, filename in (('faces-12288-initial', 'evaluation-000001.json'), ('faces-12288-terminal', terminal)):
        samples.append((name, run / filename, config['request'], config['evaluator'], result['identity']['normalization'], physical))
    rows = []
    native = spin._native_module()
    for name, path, request, settings, norm, design in samples:
        material = json.loads(path.read_text())
        tri = native.SpinTriangulationSettings()
        for key in ('icosphere_subdivisions', 'maximum_icosphere_subdivisions', 'maximum_radius_doublings',
                    'retain_fixed_normal_provenance', 'truncation_radius', 'paper_projective_hyperboloid', 'algebraic_finite_map'):
            setattr(tri, key, settings[key])
        axes = np.asarray(request['target_axes'])
        origin = np.asarray(request['target_origin'])
        source = native.evaluate_spin_it_faster_source_motion(np.asarray(material['quadric_coefficients']),
            (np.asarray(request['design_vertices']) - origin) @ axes, tri)
        if not source['defined'] or not source['side_selection_defined']:
            raise ValueError('SOURCE_REPLAY_UNAVAILABLE')
        source_vertices = (np.asarray(source['vertices']) @ axes.T + origin) * norm['L'] + np.asarray(norm['origin'])
        if not np.array_equal(source_vertices, material['quadric_source_vertices']):
            raise ValueError('SOURCE_REPLAY_GEOMETRY_CHANGED')
        source_velocity = np.einsum('vkd,jk->vjd', np.asarray(source['vertex_velocity']).reshape(-1, 3, 10), axes) * norm['L']
        physical_vertices = np.array([[float(exact_fraction(x)) for x in p]
                                     for p in material['material_proof']['exact_vertices']])
        velocity = recover_construction_vertex_motion(material['material_proof']['resolved_construction_records'],
            design, source_vertices, source_velocity, physical_vertices)
        velocity = np.einsum('vkd,kj->vjd', velocity / norm['L'], axes)
        vertices = np.asarray(material['variable_vertices'])
        faces = np.asarray(material['variable_faces'], dtype=np.int64)
        rational = [[Fraction(x) for x in p] for p in material['variable_rational_vertices']]
        tokens = [[str(x) for x in p] for p in rational]
        density = request['material_density']
        old_moments, old_integral_seconds = timed(lambda: integrate_rational_polyhedral_raw_moments(
            vertices, faces, density, rational_vertices=rational))
        new_moments, new_integral_seconds = timed(lambda: arithmetic.integrate_moments(tokens, faces, density))
        old_jacobian, old_derivative_seconds = timed(lambda: differentiate_polyhedral_raw_moments(
            vertices, faces, velocity, density, rational_vertices=rational))

        def new_derivative():
            coefficients = arithmetic.transport_coefficients(tokens, faces, density)
            terms = np.einsum('fjmc,fjcd->fmd', coefficients, velocity[faces])
            return np.array([[fsum(terms[:, m, d]) for d in range(10)] for m in range(10)])

        new_jacobian, new_derivative_seconds = timed(new_derivative)
        row = dict(name=name, faces=len(faces), source_replayed=True, boolean_calls=0,
                   original_moments_match_saved=np.array_equal(old_moments, material['variable_raw_moments']),
                   original_jacobian_match_saved=np.array_equal(old_jacobian, material['moment_jacobian']),
                   moments_equal=np.array_equal(old_moments, new_moments),
                   jacobian_equal=np.array_equal(old_jacobian, new_jacobian),
                   maximum_jacobian_error=float(np.max(np.abs(old_jacobian - new_jacobian))),
                   old_integral_seconds=old_integral_seconds, new_integral_seconds=new_integral_seconds,
                   old_derivative_seconds=old_derivative_seconds, new_derivative_seconds=new_derivative_seconds)
        rows.append(row)
        (output / 'rows.json').write_text(json.dumps(rows))
    (output / 'report.json').write_text(json.dumps(dict(rows=rows,
        passed=all(r['original_moments_match_saved'] and r['original_jacobian_match_saved']
                   and r['moments_equal'] and r['jacobian_equal'] for r in rows),
        source_calls=len(rows), boolean_calls=0, core_admission=False,
        scope='Arithmetic kernels on saved actual material; full adapter integration pending')))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
