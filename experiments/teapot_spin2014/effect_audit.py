"""Read-only reconstructed Spin-It moment diagnostics, never an admission gate.

The native header contains density-one binary64 moments ordered
M,x,y,z,x2,y2,z2,yz,xz,xy. Beta is removal fraction. Full JSONL verification,
topology certification, physical validation and paper reproduction are outside
this check. No native solver is imported or called.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import trimesh

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification import mass_moments as independent

TOL = 1e-10
GRAVITY = 9.80665
ORDER = ['M', 'x', 'y', 'z', 'x2', 'y2', 'z2', 'yz', 'xz', 'xy']


def _finite(values, shape, name):
    try:
        array = np.asarray(values, dtype=np.float64)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(f'{name}: invalid numeric values') from error
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError(f'{name}: wrong dimension or nonfinite values')
    return array


def _positive(value, name):
    result = float(_finite(value, (), name))
    if result <= 0:
        raise ValueError(f'{name}: must be positive')
    return result


def _tensor_metrics(raw):
    mass = _positive(raw[0], 'mass')
    com = raw[1:4] / mass
    second = np.array([[raw[4], raw[9], raw[8]],
                       [raw[9], raw[5], raw[7]],
                       [raw[8], raw[7], raw[6]]])
    pivot = np.trace(second) * np.eye(3) - second
    centered_second = second - mass * np.outer(com, com)
    central = pivot - mass * (float(com @ com) * np.eye(3) - np.outer(com, com))
    if not np.isfinite(central).all() or not np.isfinite(pivot).all():
        raise ValueError('nonfinite derived inertia')
    eigenvalues, eigenvectors = np.linalg.eigh(central)
    if eigenvalues[0] <= 0 or np.linalg.eigvalsh(centered_second)[0] <= 0:
        raise ValueError('nonpositive central moments or inertia')
    pivot_eigenvalues = np.linalg.eigvalsh(pivot)
    lateral = np.linalg.eigvalsh(central[:2, :2])
    pivot_lateral = np.linalg.eigvalsh(pivot[:2, :2])
    dominant = eigenvalues[-1] - eigenvalues <= TOL * eigenvalues[-1]
    basis = eigenvectors[:, dominant]
    projection = min(1., float(np.linalg.norm(basis[2, :])))
    subspace_angle = math.degrees(math.acos(projection))
    unique = int(dominant.sum()) == 1
    ratios = lateral / central[2, 2]
    return dict(mass=mass, com=com.tolist(), mass_times_com_height=float(mass * com[2]),
                com_to_axis_distance=float(np.linalg.norm(com[:2])),
                central_inertia=central.tolist(), pivot_inertia=pivot.tolist(),
                central_eigenvalues=eigenvalues.tolist(), pivot_eigenvalues=pivot_eigenvalues.tolist(),
                dominant_axis_status='UNIQUE' if unique else 'NONUNIQUE',
                dominant_axis=eigenvectors[:, -1].tolist() if unique else None,
                dominant_axis_angle_degrees=subspace_angle if unique else None,
                dominant_subspace_angle_degrees=subspace_angle,
                dominant_subspace_dimension=int(dominant.sum()),
                central_lateral_eigenvalues=lateral.tolist(), central_Izz=float(central[2, 2]),
                central_lateral_ratios=ratios.tolist(),
                central_lateral_squared_ratio_sum=float(ratios @ ratios),
                central_transverse_coupling=float(np.linalg.norm(central[:2, 2])),
                normalized_central_transverse_coupling=float(np.linalg.norm(central[:2, 2]) / eigenvalues[-1]),
                pivot_transverse_coupling=float(np.linalg.norm(pivot[:2, 2])),
                normalized_pivot_transverse_coupling=float(np.linalg.norm(pivot[:2, 2]) / pivot_eigenvalues[-1]),
                pivot_lateral_eigenvalues=pivot_lateral.tolist(), pivot_Izz=float(pivot[2, 2]),
                central_target_minus_lateral_gap=float(central[2, 2] - lateral[-1]),
                pivot_target_minus_lateral_gap=float(pivot[2, 2] - pivot_lateral[-1]))


def _speed(mass, height, gap):
    omega = math.sqrt(mass * GRAVITY * height / gap)
    if not math.isfinite(omega):
        raise ValueError('nonfinite formal speed')
    return dict(omega_rad_s=omega, revolutions_per_second=omega / (2 * math.pi),
                rpm=omega * 60 / (2 * math.pi))


def moment_metrics(raw, meters_per_model_unit, density_kg_m3, *, native_density=1., root_scale=1.):
    """Compute tensors in the fixed origin/target-z frame and convert once to SI.

    Numerical alignment uses TOL relative to the declared input root scale;
    eigenspace and coupling checks use TOL relative to the largest eigenvalue.
    The bound describes an ideal fixed pivot only, never simulated stability.
    """
    raw = _finite(raw, (10,), 'raw moments')
    if _positive(native_density, 'native density') != 1.:
        raise ValueError('incompatible native density: density-one moments required')
    scale = _positive(meters_per_model_unit, 'meters_per_model_unit')
    density = _positive(density_kg_m3, 'density_kg_m3')
    root = _positive(root_scale, 'root scale')
    model = _tensor_metrics(raw)
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        factors = density * np.power(scale, np.array([3, 4, 4, 4, 5, 5, 5, 5, 5, 5]))
        converted = _finite(raw * factors, (10,), 'SI moments')
    si = _tensor_metrics(converted)
    no_speed = dict(omega_rad_s=None, revolutions_per_second=None, rpm=None)
    formal = dict(label='FORMAL_COM_BASED_PAPER_SCORE_SURROGATE_NOT_PHYSICAL_CRITICAL_SPEED',
                  status='NOT_APPLICABLE', **no_speed)
    if si['central_target_minus_lateral_gap'] > 0 and si['com'][2] > 0:
        formal.update(status='FORMAL_SURROGATE_ONLY',
                      **_speed(si['mass'], si['com'][2], si['central_target_minus_lateral_gap']))
    reasons = []
    if model['pivot_target_minus_lateral_gap'] <= 0:
        reasons.append('NONPOSITIVE_PIVOT_GAP')
    if model['com_to_axis_distance'] > TOL * root:
        reasons.append('COM_OFF_AXIS')
    if model['normalized_pivot_transverse_coupling'] > TOL:
        reasons.append('AXIS_NOT_PRINCIPAL')
    if np.min(np.diff(model['pivot_eigenvalues'])) <= TOL * model['pivot_eigenvalues'][-1]:
        reasons.append('DEGENERATE_PIVOT_EIGENVALUES')
    if model['com'][2] <= 0:
        reasons.append('NONPOSITIVE_COM_HEIGHT')
    fixed = dict(label='LEWIS1992_IDEAL_FIXED_POINT_SUFFICIENT_BOUND_PIVOT_INERTIA',
                 status='NOT_APPLICABLE' if reasons else 'CONDITIONS_SATISFIED_FOR_FORMAL_BOUND',
                 reasons=reasons, inequality='omega_squared_strictly_greater_than_bound_squared',
                 **no_speed)
    if not reasons:
        fixed.update(**_speed(si['mass'], si['com'][2], si['pivot_target_minus_lateral_gap']))
    return dict(raw_moments_model_units=raw.tolist(), model_units=model, si=si,
                formal_com_paper_score=formal, fixed_point_bound=fixed)


def normalized_error(left, right, root_scale):
    """Absolute component error normalized by root_side^(3+polynomial degree)."""
    root = _positive(root_scale, 'root scale')
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        factors = np.power(root, np.array([3, 4, 4, 4, 5, 5, 5, 5, 5, 5]))
    if not np.isfinite(factors).all() or np.any(factors <= 0):
        raise ValueError('unrepresentable root normalization')
    return float(np.max(np.abs(_finite(left, (10,), 'left moments') -
                               _finite(right, (10,), 'right moments')) / factors))


def extract_states(payload, initial_fill, root_scale):
    """Reconstruct represented material, never independently round a relaxation."""
    source, solve = payload['source'], payload['solve']
    outer = _finite(source['outer_surface_moments'], (10,), 'outer moments')
    carveable = _finite(source['carveable_moments'], (10,), 'carveable moments')
    count = len(solve['final_partition'])
    cells = _finite([leaf['moments'] for leaf in solve['final_partition']], (count, 10), 'partition')
    total = np.array([math.fsum(cells[:, j]) for j in range(10)])
    partition_error = normalized_error(total, carveable, root_scale)
    if partition_error > TOL:
        raise ValueError('partition does not match carveable moments')

    def material(beta, binary=False):
        values = _finite(beta, (count,), 'beta')
        if np.any((values < 0) | (values > 1)) or (binary and np.any((values != 0) & (values != 1))):
            raise ValueError('beta outside removal range or nonbinary native output')
        return outer - np.array([math.fsum(values * cells[:, j]) for j in range(10)])

    initial = float(_finite(initial_fill, (), 'initial fill'))
    states = dict(original=outer, initial_mixture=material(np.full(count, initial)),
                  terminal_relaxed=material(solve['relaxed_fill']), native_binary=None)
    binary_error = None
    fill, moments = solve.get('binary_fill', []), solve.get('binary_moments')
    if moments is not None or len(fill):
        if moments is None or not len(fill):
            raise ValueError('incomplete native binary output')
        binary = material(fill, binary=True)
        binary_error = normalized_error(binary, moments, root_scale)
        if binary_error > TOL:
            raise ValueError('native binary moments do not match binary_fill')
        states['native_binary'] = binary
    return states, dict(partition_normalized_error=partition_error, binary_normalized_error=binary_error)


def integrate_original(path):
    """Use the source runner's process=False OBJ ingestion and independent integral.

    Conversion with tolist preserves binary64 coordinates and integer indices;
    NumPy int64 faces must become Python ints for TriangleMesh's strict checker.
    This does not certify source self-intersection freedom or region topology.
    """
    mesh = trimesh.load(Path(path), force='mesh', process=False)
    vertices, faces = np.asarray(mesh.vertices), np.asarray(mesh.faces)
    if vertices.dtype != np.dtype('float64') or faces.dtype.kind not in 'iu':
        raise ValueError('source OBJ ingestion must preserve binary64 vertices and integer faces')
    result = independent.integrate_geometric_moments(independent.TriangleMesh(vertices.tolist(), faces.tolist()))
    q = result.second_moment
    raw = np.array([result.volume, *result.first_moment, q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1]])
    return raw, dict(method='PYTHON_ORIENTED_TETRAHEDRA_MASS_MOMENTS',
                     vertices_dtype=str(vertices.dtype), faces_dtype=str(faces.dtype),
                     vertices=len(vertices), triangles=len(faces), process=False,
                     full_source_topology_certification='NOT_PROVIDED')


def _sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _read(path):
    def invalid(token):
        raise ValueError('nonfinite JSON token: ' + token)
    return json.loads(Path(path).read_text(), parse_constant=invalid)


def _verify(files):
    for name, expected in files.items():
        if not Path(name).is_file() or _sha(name) != expected:
            raise ValueError('STALE_RESULT: ' + name)


def audit_run(run):
    """Verify small-file evidence and current upstream files; never read JSONL."""
    run = Path(run).resolve()
    outputs_path = run / 'outputs.json'
    outputs_hash = _sha(outputs_path)
    outputs = _read(outputs_path)
    required = [run / name for name in ('header.json', 'config.json', 'source.json', 'summary.json')]
    if any(str(path) not in outputs for path in required):
        raise ValueError('output manifest missing required file')
    trace = run / 'native-run.jsonl'
    # The archive digest is recorded but intentionally not verified here.
    checked_outputs = {p: h for p, h in outputs.items() if Path(p).resolve() != trace}
    if any(Path(p).suffix == '.jsonl' for p in checked_outputs):
        raise ValueError('unexpected trace path in output manifest')
    _verify(checked_outputs)
    header, config, source, summary = map(_read, required)
    if header.get('kind') != 'header':
        raise ValueError('not a stream_probe header')
    upstream = source['files']
    if any(Path(p).suffix == '.jsonl' for p in upstream):
        raise ValueError('trace not allowed in upstream manifest')
    _verify(upstream)
    reference = Path(config['reference_config']).resolve()
    cfg = config['reference_settings']
    inputs = Path(cfg['input']).resolve()
    manifest = inputs / 'manifest.json'
    for path in (reference, manifest, inputs / 'outer.obj', inputs / 'design.obj'):
        if str(path) not in upstream:
            raise ValueError('upstream manifest missing required input: ' + str(path))
    if _sha(reference) != config['reference_sha256'] or _read(reference) != cfg:
        raise ValueError('STALE_REFERENCE_CONFIG')
    geometry = _read(manifest)['geometry_hashes']
    if any(upstream.get(p) != h for p, h in geometry.items()):
        raise ValueError('geometry and upstream manifests disagree')
    payload = header['payload']
    for key in ('source_settings', 'solver_settings', 'boundary_settings', 'density', 'minimum_outer_wall'):
        if payload['settings'][key] != config[key]:
            raise ValueError('header/config settings mismatch: ' + key)
    if (config['source_settings'] != cfg['geometry'] or
            config['boundary_settings'] != cfg['boundary'] or
            config['minimum_outer_wall'] != cfg['minimum_wall'] or
            config['solver_settings'] != {**cfg['solver'], 'termination_policy': 'PAPER_POST_TRANSITION_V1'}):
        raise ValueError('reference/runner settings mismatch')
    for key in ('status', 'trace_hash', 'output_hash', 'binary_check'):
        if summary[key] != payload['solve'][key]:
            raise ValueError('header/summary mismatch: ' + key)
    root_scale = config['source_settings']['root_side']
    states, checks = extract_states(payload, config['solver_settings']['initial_fill'], root_scale)
    original, integration = integrate_original(inputs / 'outer.obj')
    checks['original_normalized_error'] = normalized_error(original, states['original'], root_scale)
    if checks['original_normalized_error'] > TOL:
        raise ValueError('original geometry moments mismatch')
    checker_files = {str(Path(p).resolve()): _sha(p) for p in (__file__, independent.__file__)}
    report = dict(schema='teapot-spin2014-moment-effect-v1', scope='METRIC_ONLY_NO_ADMISSION',
                  run=str(run), native_status=payload['solve']['status'], moment_order=ORDER,
                  full_trace_verification='NOT_PERFORMED_SEPARATE_PARENT_GATE',
                  recorded_trace_sha256=outputs.get(str(trace)),
                  original_paper_case=False, physical_validation='NOT_PERFORMED',
                  material_admission='NOT_GRANTED', geometry_integration=integration, checks=checks,
                  normalization='max_abs_component_error / root_side^(3+polynomial_degree)',
                  tolerance=TOL, gravity_m_s2=GRAVITY,
                  units=dict(meters_per_model_unit=cfg['meters_per_model_unit'],
                             density_kg_m3=cfg['density_kg_m3'], native_density=config['density'],
                             model_mass='unit_density_model_volume', model_inertia='model_length^5',
                             si_mass='kg', si_length='m', si_inertia='kg*m^2'),
                  checker_files=checker_files, library_versions=dict(numpy=np.__version__, trimesh=trimesh.__version__),
                  verified_output_files=checked_outputs, verified_upstream_files=upstream,
                  output_manifest_sha256=outputs_hash,
                  states={name: None if raw is None else moment_metrics(raw, cfg['meters_per_model_unit'],
                          cfg['density_kg_m3'], native_density=config['density'], root_scale=root_scale)
                          for name, raw in states.items()})
    # Catch modifications made during the computation without opening the trace.
    _verify({**checked_outputs, **upstream, **checker_files, str(outputs_path): outputs_hash})
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('--output must be a new directory')
    report = audit_run(args.run)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'effect.json').write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n')
    print('METRIC_ONLY_NO_ADMISSION', args.output / 'effect.json')


if __name__ == '__main__':
    main()
