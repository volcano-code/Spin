"""Bounded exact construction and admission of the preregistered Heart v2.

All Boolean operands and results retain rational coordinates. The connected
union uses De Morgan inside an enclosing box because the pinned rational
construction interface exposes intersection and difference, not union.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import gmpy2
import numpy as np
import trimesh

from rounded_input import SPEC, capacity_bounds, exact_moments, publish, write_mesh


def exact_rows(vertices):
    result = []
    for row in vertices:
        converted = []
        for value in row:
            q = gmpy2.mpq(float(value)) if isinstance(value, (float, np.floating)) else gmpy2.mpq(value)
            converted.append(f'{q.numerator}/{q.denominator}')
        result.append(converted)
    return result


def construct(native, a, b, operation):
    result = native.libigl_exact_boolean_construction_from_rationals(
        a[0], np.asarray(a[1], dtype=np.int32), b[0], np.asarray(b[1], dtype=np.int32), operation)
    if not result['succeeded']:
        raise ValueError('EXACT_CONSTRUCTION_FAILED: '+operation)
    return exact_rows(result['vertices']), np.asarray(result['faces'], dtype=np.int32).tolist()


def union(native, a, b, enclosure):
    rest = construct(native, enclosure, a, 'difference')
    rest = construct(native, rest, b, 'difference')
    return construct(native, enclosure, rest, 'difference')


def save_region(path, region):
    publish(path.with_suffix('.rational.json'), dict(schema='authoritative-rational-input-boundary-v1',
        vertices=region[0], faces=region[1], obj_role='FLOAT_DISPLAY_NOT_AUTHORITY'))
    vertices = np.array([[float(gmpy2.mpq(x)) for x in row] for row in region[0]])
    write_mesh(path, trimesh.Trimesh(vertices, np.array(region[1]), process=False))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    with (output/'events.jsonl').open('x') as events:
        def event(stage):
            value = dict(stage=stage, elapsed_seconds=time.monotonic()-started)
            events.write(json.dumps(value)+'\n'); events.flush(); os.fsync(events.fileno())
            print(json.dumps(value), flush=True)
        try:
            event('SOURCE_IDENTITY')
            prepared = json.loads((args.input/'preparation.json').read_text())
            for name, expected in prepared['sources'].items():
                if hashlib.sha256((args.input/name).read_bytes()).hexdigest() != expected:
                    raise ValueError('STALE_INPUT: '+name)
            if json.loads((args.input/'spec.json').read_text()) != SPEC:
                raise ValueError('SPEC_CHANGED_AFTER_PREPARATION')
            import _make_it_stand as native
            from make_it_stand.verification.spin_problem import MeshSnapshot, SpinArtifactProblem
            event('LOAD_BINARY64_SOURCE_AS_EXACT_RATIONALS')
            def load(name):
                mesh = trimesh.load(args.input/name, force='mesh', process=False)
                MeshSnapshot((mesh.vertices, mesh.faces))
                return exact_rows(mesh.vertices), mesh.faces.tolist()
            body, inner, pin = [load(name) for name in ('heart-body.obj', 'inner-source.obj', 'pin.obj')]
            all_vertices = np.array([[float(gmpy2.mpq(x)) for x in row] for row in body[0]+pin[0]])
            low, high = all_vertices.min(axis=0)-1, all_vertices.max(axis=0)+1
            box = trimesh.creation.box(extents=high-low)
            box.apply_translation((low+high)/2)
            enclosure = exact_rows(box.vertices), box.faces.tolist()
            event('EXACT_CONNECTED_UNION')
            outer = union(native, body, pin, enclosure)
            save_region(output/'outer.obj', outer)
            event('EXACT_DOMAIN_EXCLUDING_PIN')
            design = construct(native, inner, pin, 'difference')
            save_region(output/'design.obj', design)
            event('EXACT_FIXED_MATERIAL')
            fixed = construct(native, outer, design, 'difference')
            save_region(output/'fixed.obj', fixed)
            event('INDEPENDENT_SIGNED_MOMENTS_AND_CAPACITY')
            moments = {name: exact_moments(*region) for name, region in
                       (('outer', outer), ('fixed', fixed), ('design', design))}
            partition_delta = [a-b-c for a, b, c in zip(moments['outer'], moments['fixed'], moments['design'])]
            if any(x != 0 for x in partition_delta):
                raise ValueError('INDEPENDENT_PARTITION_MOMENTS_DIFFER')
            capacity = capacity_bounds(*fixed, *design)
            event('EXACT_REGION_AND_PARTITION_RELOAD')
            snapshots = {name: MeshSnapshot.from_rationals(*region) for name, region in
                         (('outer', outer), ('fixed', fixed), ('design', design))}
            problem = SpinArtifactProblem(original=snapshots['outer'],
                fixed_regions=[('fixed-shell-and-rounded-pin', snapshots['fixed'], 1.)],
                design_regions=[('heart-inner-excluding-pin', snapshots['design'], 1.)],
                length_unit='25mm', mass_unit='18.75g',
                target=dict(mode='SPINNING_TOP', axis=[0., 0., 1.], origin=[0., 0., 0.], support_points=[]),
                minimum_outer_wall=SPEC['minimum_outer_clearance_model'])
            capacity['partition_digest'] = problem.input_digest
            capacity['geometry_assumptions_rechecked'] = True
            publish(output/'capacity.json', capacity)
            event('EXACT_OUTER_DOMAIN_CLEARANCE')
            clearance = native.check_exact_boundary_clearance(snapshots['outer'].exact_region(),
                snapshots['design'].exact_region(), SPEC['minimum_outer_clearance_model'])
            publish(output/'clearance.json', clearance)
            s = np.array([float(x) for x in moments['outer']])
            q = np.array([[s[4], s[9], s[8]], [s[9], s[5], s[7]], [s[8], s[7], s[6]]])
            center = s[1:4]/s[0]
            inertia = np.trace(q)*np.eye(3)-q-s[0]*(center@center*np.eye(3)-np.outer(center, center))
            si_mass_scale = SPEC['density_kg_m3']*SPEC['model_unit_m']**3
            input_ok = (clearance['completed'] and clearance['at_least_minimum'] is True)
            status = ('INPUT_GEOMETRY_BLOCKED' if not input_ok else
                      'INPUT_CAPACITY_EXCLUDED' if capacity['status']=='EXCLUDED_BY_FIRST_MOMENT_BOUND' else
                      'INPUT_ADMITTED_FEASIBILITY_UNKNOWN')
            record = dict(schema='rounded-heart-input-admission-v1', status=status,
                problem_identity=problem.identity, problem_digest=problem.input_digest,
                material_facts={name: snapshot.exact_region().facts for name, snapshot in snapshots.items()},
                capacity_status=capacity['status'], mass_kg=s[0]*si_mass_scale,
                center_m=(center*SPEC['model_unit_m']).tolist(),
                inertia_center_kg_m2=(inertia*si_mass_scale*SPEC['model_unit_m']**2).tolist(),
                cross_axis_inertia_relative=float(np.linalg.norm(inertia[:2, 2])/np.linalg.norm(inertia)),
                independent_moments_exact={name: [str(x) for x in s] for name, s in moments.items()},
                partition_moment_delta_exact=[str(x) for x in partition_delta],
                native_path=native.__file__, native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
                files={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob('*')
                       if p.suffix in ('.json', '.obj')},
                optimization='NOT_RUN', simulation='NOT_RUN',
                full_manufacturing_wall_thickness='UNKNOWN')
            publish(output/'admission.json', record)
            event(status)
            publish(output/'terminal.json', dict(status=status, completed=True,
                optimization_invocations=0, simulation_invocations=0,
                elapsed_seconds=time.monotonic()-started))
        except Exception as exc:
            publish(output/'failure.json', dict(type=type(exc).__name__, message=str(exc),
                completed=False, elapsed_seconds=time.monotonic()-started))
            raise


if __name__ == '__main__':
    main()
