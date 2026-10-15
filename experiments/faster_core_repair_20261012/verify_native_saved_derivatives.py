"""Kernel equivalence on the registered 243 saved production samples."""
from fractions import Fraction
from math import fsum
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import _cherchi_rational_moments as arithmetic
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import exact_fraction
from make_it_stand.verification.cherchi_motion_differential import recover_construction_vertex_motion

HERE = Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    registered = HERE / 'projective-formal-pipeline-job-01/output/registered-inputs.json'
    frozen = json.loads(registered.read_text())
    request, settings = frozen['request'], frozen['settings']
    axes = np.asarray(request['target_axes']); origin = np.asarray(request['target_origin'])
    design = np.asarray(request['design_vertices'])
    native = spin._native_module()
    tri = native.SpinTriangulationSettings()
    for key in ('icosphere_subdivisions', 'maximum_icosphere_subdivisions', 'maximum_radius_doublings',
                'retain_fixed_normal_provenance', 'truncation_radius', 'paper_projective_hyperboloid', 'algebraic_finite_map'):
        setattr(tri, key, settings[key])
    reports = []
    for state in ('ellipsoid', 'one_sheet', 'two_sheet'):
        folder = HERE / 'projective-formal-pipeline-job-01/output/derivatives' / state
        manifest = json.loads((folder / 'manifest.json').read_text())
        certificate = json.loads((HERE / 'projective-fixed-branch-job-01/output' / (state + '.json')).read_text())
        if not certificate['finite_sample_certification']:
            raise ValueError('REGISTERED_SAMPLE_BRANCH_NOT_CERTIFIED')
        rows = []
        with (output / (state + '-events.jsonl')).open('x') as stream:
            for index, q in enumerate(manifest['samples']):
                path = folder / f'point-{index:03d}-material.json'
                if hashlib.sha256(path.read_bytes()).hexdigest() != certificate['rows'][index]['material_sha256']:
                    raise ValueError('SAVED_MATERIAL_CHANGED')
                material = json.loads(path.read_text())['same_call_material']
                if material['quadric_coefficients'] != q:
                    raise ValueError('SAVED_PARAMETER_CHANGED')
                source = native.evaluate_spin_it_faster_source_motion(np.asarray(q), (design - origin) @ axes, tri)
                if not source['defined'] or not source['side_selection_defined']:
                    raise ValueError('SOURCE_REPLAY_UNAVAILABLE')
                source_vertices = np.asarray(source['vertices']) @ axes.T + origin
                if not np.array_equal(source_vertices, material['quadric_source_vertices']) or not np.array_equal(source['faces'], material['quadric_source_faces']):
                    raise ValueError('SOURCE_REPLAY_CHANGED')
                source_velocity = np.einsum('vkd,jk->vjd', np.asarray(source['vertex_velocity']).reshape(-1, 3, 10), axes)
                physical_vertices = np.array([[float(exact_fraction(x)) for x in p]
                    for p in material['material_proof']['exact_vertices']])
                moving = recover_construction_vertex_motion(material['material_proof']['resolved_construction_records'],
                    design, source_vertices, source_velocity, physical_vertices)
                moving = np.einsum('vkd,kj->vjd', moving, axes)
                rational = [[Fraction(x) for x in p] for p in material['variable_rational_vertices']]
                tokens = [[str(x) for x in p] for p in rational]
                faces = np.asarray(material['variable_faces'], dtype=np.int64)
                moments = arithmetic.integrate_moments(tokens, faces, request['material_density'])
                coefficients = arithmetic.transport_coefficients(tokens, faces, request['material_density'])
                terms = np.einsum('fjmc,fjcd->fmd', coefficients, moving[faces])
                jacobian = np.array([[fsum(terms[:, m, d]) for d in range(10)] for m in range(10)])
                row = dict(index=index, material_sha256=certificate['rows'][index]['material_sha256'],
                    source_replay_equal=True, moments_equal=np.array_equal(moments, material['variable_raw_moments']),
                    jacobian_equal=np.array_equal(jacobian, material['moment_jacobian']),
                    maximum_moment_error=float(np.max(np.abs(moments - material['variable_raw_moments']))),
                    maximum_jacobian_error=float(np.max(np.abs(jacobian - material['moment_jacobian']))))
                rows.append(row)
                stream.write(json.dumps(row) + '\n'); stream.flush()
                if not row['moments_equal'] or not row['jacobian_equal']:
                    raise ValueError('NATIVE_KERNEL_CHANGED_SAVED_PRODUCTION_MAP')
        report = dict(state=state, sample_count=len(rows), passed=len(rows) == 81,
            rows=rows, original_branch_certificate_sha256=hashlib.sha256(
                (HERE / 'projective-fixed-branch-job-01/output' / (state + '.json')).read_bytes()).hexdigest())
        (output / (state + '.json')).write_text(json.dumps(report)); reports.append(report)
    (output / 'report.json').write_text(json.dumps(dict(states=reports,
        passed=all(r['passed'] for r in reports), source_calls=243, boolean_calls=0,
        arithmetic_module_sha256=hashlib.sha256(Path(arithmetic.__file__).read_bytes()).hexdigest(),
        core_admission=False, scope='Finite saved production map equivalence; no neighborhood proof')))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
