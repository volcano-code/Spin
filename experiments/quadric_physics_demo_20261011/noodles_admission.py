"""One bounded exact admission of frozen Noodles source parts; no optimizer."""
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT/'python'), str(ROOT/'experiments/faster_heart_contact_20261011')]
from noodles_input import SPEC, spec_identity, write_obj
from admit_rounded_input import exact_rows, construct, union, save_region
from rounded_input import exact_moments, capacity_bounds


def publish(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.flush(); os.fsync(stream.fileno())
    fd = os.open(Path(path).parent, os.O_RDONLY|os.O_DIRECTORY)
    try: os.fsync(fd)
    finally: os.close(fd)


def boundary_components(faces):
    """Face connectivity; combined with native valid embedded region evidence."""
    parent = list(range(len(faces))); first = {}
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for i, face in enumerate(faces):
        for a, b in zip(face, face[1:]+face[:1]):
            edge = tuple(sorted((int(a), int(b))))
            if edge in first: parent[root(i)] = root(first[edge])
            else: first[edge] = i
    return len({root(i) for i in range(len(faces))})


def translated(region, origin):
    q = [gmpy2.mpq(float(x)) for x in origin]
    return [[str(gmpy2.mpq(x)-q[k]) for k, x in enumerate(row)] for row in region[0]], region[1]


def measured_interface_clearance(clearance):
    """Derived positive interface bound; never a printing-thickness criterion."""
    if not clearance['completed'] or not clearance['applicable']:
        raise ValueError('OUTER_DOMAIN_DISTANCE_UNKNOWN')
    squared = gmpy2.mpq(clearance['squared_distance_exact'])
    if squared <= 0: raise ValueError('DESIGN_DOMAIN_CONTACTS_OUTER_BOUNDARY')
    candidate = float(gmpy2.sqrt(squared))/2
    if candidate <= 0 or not np.isfinite(candidate):
        raise ValueError('MEASURED_DISTANCE_NOT_REPRESENTABLE')
    while gmpy2.mpq(candidate)**2 > squared/4:
        candidate = np.nextafter(candidate, 0.)
    return candidate


def validate_manifest(manifest,input_version='v1'):
    if input_version=='v1':
        expected_spec,expected_identity=SPEC,spec_identity()
        tube_names=[f'tube-{i}-fixed' for i in range(len(SPEC['tube_angles_degrees']))]
        generator=HERE/'noodles_input.py'
    elif input_version=='dense-v2':
        from noodles_dense_input import SPEC as dense_spec,identity,tube_layout
        expected_spec,expected_identity=dense_spec,identity()
        tube_names=[f'tube-{i:03d}-fixed' for i in range(len(tube_layout()))]
        generator=HERE/'noodles_dense_input.py'
        if manifest.get('tube_mesh_helper_sha256')!=hashlib.sha256((HERE/'noodles_input.py').read_bytes()).hexdigest():
            raise ValueError('STALE_SOURCE_TUBE_HELPER')
    elif input_version=='dense-v3':
        from noodles_dense_v3_input import SPEC as dense_spec,identity,tube_layout
        expected_spec,expected_identity=dense_spec,identity()
        tube_names=[f'tube-{i:03d}-fixed' for i in range(len(tube_layout()))]
        generator=HERE/'noodles_dense_v3_input.py'
        if manifest.get('tube_mesh_helper_sha256')!=hashlib.sha256((HERE/'noodles_input.py').read_bytes()).hexdigest():
            raise ValueError('STALE_SOURCE_TUBE_HELPER')
    elif input_version in ('dense-v4','dense-v5'):
        if input_version=='dense-v4':
            from noodles_dense_v4_input import SPEC as dense_spec,identity,tube_layout
        else:
            from noodles_dense_v5_input import SPEC as dense_spec,identity,tube_layout
        expected_spec,expected_identity=dense_spec,identity()
        tube_names=[f'tube-{i:03d}-fixed' for i in range(len(tube_layout()))]
        generator=HERE/('noodles_dense_v4_input.py' if input_version=='dense-v4' else 'noodles_dense_v5_input.py')
        if manifest.get('tube_mesh_helper_sha256')!=hashlib.sha256((HERE/'noodles_input.py').read_bytes()).hexdigest():
            raise ValueError('STALE_SOURCE_TUBE_HELPER')
        if manifest.get('inherited_geometry_generator_sha256')!=hashlib.sha256((HERE/'noodles_dense_v3_input.py').read_bytes()).hexdigest():
            raise ValueError('STALE_INHERITED_GEOMETRY_GENERATOR')
    else:raise ValueError('UNSUPPORTED_SOURCE_INPUT_VERSION')
    if manifest['spec']!=expected_spec or manifest['spec_sha256']!=expected_identity:
        raise ValueError('STALE_FROZEN_SPEC')
    if manifest.get('creator_sha256')!=hashlib.sha256(generator.read_bytes()).hexdigest():
        raise ValueError('STALE_SOURCE_GENERATOR')
    if manifest['construction']['fixed_parts']!=['hub-fixed','pin-fixed']+tube_names:
        raise ValueError('FROZEN_FIXED_PARTS_IDENTITY_MISMATCH')
    if set(manifest['parts'])!=set(['outer-ellipsoid','inner-ellipsoid-candidate','hub-fixed','pin-fixed']+tube_names):
        raise ValueError('FROZEN_SOURCE_PART_SET_MISMATCH')
    return expected_spec


def admit(case, output,*,input_version='v1'):
    case, output = Path(case), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    with (output/'events.jsonl').open('x') as events:
        def event(stage, **details):
            record = dict(stage=stage, elapsed_seconds=time.monotonic()-started, **details)
            events.write(json.dumps(record, allow_nan=False)+'\n'); events.flush(); os.fsync(events.fileno())
            print(json.dumps(record), flush=True)
        try:
            event('FROZEN_SOURCE_IDENTITY')
            manifest = json.loads((case/'manifest.json').read_text())
            spec=validate_manifest(manifest,input_version)
            import _make_it_stand as native
            from make_it_stand.verification.spin_problem import MeshSnapshot, SpinArtifactProblem
            from make_it_stand.verification.spin_common import native as shared_native
            if Path(shared_native().__file__).resolve() != Path(native.__file__).resolve():
                raise ValueError('NATIVE_LOADING_IDENTITY_MISMATCH')
            event('SOURCE_PART_EXACT_SELF_INTERSECTION_AND_SI_VALIDITY')
            regions, part_facts = {}, {}
            for name, entry in manifest['parts'].items():
                for suffix in ('model', 'si'):
                    event('SOURCE_PART_CHECK_START',part=name,coordinate_unit=suffix)
                    path = case/entry[f'{suffix}_obj']
                    if hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256'][suffix]:
                        raise ValueError('STALE_SOURCE_OBJ:'+name+':'+suffix)
                    mesh = trimesh.load(path, force='mesh', process=False)
                    snapshot = MeshSnapshot((mesh.vertices, mesh.faces))
                    part_facts.setdefault(name, {})[suffix] = snapshot.exact_region().facts
                    if suffix == 'model': regions[name] = exact_rows(mesh.vertices), mesh.faces.tolist()
                event('SOURCE_PART_CHECKED', part=name)
            publish(output/'source-part-facts.json', part_facts)
            vertices = np.array([[float(gmpy2.mpq(x)) for x in row] for name, region in regions.items()
                                 if name != 'inner-ellipsoid-candidate' for row in region[0]])
            low, high = vertices.min(axis=0)-1, vertices.max(axis=0)+1
            box = trimesh.creation.box(extents=high-low); box.apply_translation((low+high)/2)
            enclosure = exact_rows(box.vertices), box.faces.tolist()
            fixed_names = manifest['construction']['fixed_parts']
            outer = regions['outer-ellipsoid']
            for name in fixed_names:
                event('EXACT_OUTER_UNION_START', part=name)
                outer = union(native, outer, regions[name], enclosure)
                event('EXACT_OUTER_UNION_END', part=name, triangles=len(outer[1]))
            event('EXACT_CONNECTED_OUTER_CHECK')
            outer_snapshot = MeshSnapshot.from_rationals(*outer)
            count = boundary_components(outer[1])
            if count != 1: raise ValueError('OUTER_BOUNDARY_NOT_SINGLE_CONNECTED:'+str(count))
            save_region(output/'outer.obj', outer)
            design = regions['inner-ellipsoid-candidate']
            for name in fixed_names:
                event('EXACT_DOMAIN_EXCLUSION_START', part=name)
                design = construct(native, design, regions[name], 'difference')
                event('EXACT_DOMAIN_EXCLUSION_END', part=name, triangles=len(design[1]))
            event('EXACT_FIXED_MATERIAL')
            fixed = construct(native, outer, design, 'difference')
            save_region(output/'design.obj', design); save_region(output/'fixed.obj', fixed)
            snapshots = dict(outer=outer_snapshot, fixed=MeshSnapshot.from_rationals(*fixed),
                             design=MeshSnapshot.from_rationals(*design))
            event('EXACT_POSITIVE_OUTER_DOMAIN_SEPARATION')
            clearance = native.check_exact_boundary_clearance(snapshots['outer'].exact_region(),
                snapshots['design'].exact_region(), 0.)
            publish(output/'clearance.json', clearance)
            measured = measured_interface_clearance(clearance)
            event('EXACT_PARTITION_AND_PROBLEM_BINDING')
            problem = SpinArtifactProblem(original=snapshots['outer'],
                fixed_regions=[('noodles-fixed-tubes-hub-pin-shell', snapshots['fixed'], 1.)],
                design_regions=[('noodles-inner-excluding-fixed-parts', snapshots['design'], 1.)],
                length_unit='25mm', mass_unit='18.75g',
                target=dict(mode='SPINNING_TOP', axis=spec['target_axis'],
                            origin=spec['target_origin_model'], support_points=[]), minimum_outer_wall=measured)
            event('INDEPENDENT_EXACT_MOMENTS_AND_TRANSLATED_CAPACITY')
            values = {name: exact_moments(*region) for name, region in
                      [('outer', outer), ('fixed', fixed), ('design', design)]}
            delta = [a-b-c for a, b, c in zip(values['outer'], values['fixed'], values['design'])]
            if any(x != 0 for x in delta): raise ValueError('PARTITION_INDEPENDENT_MOMENTS_MISMATCH')
            capacity = capacity_bounds(*translated(fixed, spec['target_origin_model']),
                                       *translated(design, spec['target_origin_model']))
            capacity['target_origin_model'] = spec['target_origin_model']
            capacity['coordinate_frame'] = 'TARGET_ORIGIN_TRANSLATED_MODEL'
            capacity['problem_digest'] = problem.input_digest
            publish(output/'capacity.json', capacity)
            event('SI_DERIVED_BOUNDARIES_AND_MOMENTS')
            scale = gmpy2.mpq(float(spec['model_unit_m']))
            si_audits = {}
            for name, region in [('outer', outer), ('fixed', fixed), ('design', design)]:
                rows = [[float(gmpy2.mpq(x)*scale) for x in row] for row in region[0]]
                mesh = trimesh.Trimesh(rows, np.asarray(region[1]), process=False)
                path = output/f'{name}-si.obj'; write_obj(path, mesh)
                loaded = trimesh.load(path, force='mesh', process=False)
                audit = native.make_exact_material_region(loaded.vertices, np.asarray(loaded.faces, dtype=np.int32)).facts
                floating = exact_moments(loaded.vertices.tolist(), loaded.faces.tolist())
                powers = [3]+[4]*3+[5]*6
                expected = [x*scale**power for x, power in zip(values[name], powers)]
                length = max(float(np.ptp(loaded.vertices, axis=0).max()), 1e-30)
                errors = [float(abs(a-b)/gmpy2.mpq(length)**p) for a, b, p in zip(floating, expected, powers)]
                si_audits[name] = dict(native_facts=audit, maximum_normalized_moment_error=max(errors),
                    moments_pass=max(errors)<=1e-6, file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    scope='SOURCE_EXPORT_VALIDITY_AND_MOMENTS_NOT_COMPLETE_DERIVED_PARTITION_ADMISSION')
            publish(output/'si-export-audit.json', si_audits)
            s = np.array([float(x) for x in values['outer']]); center = s[1:4]/s[0]
            q = np.array([[s[4],s[9],s[8]],[s[9],s[5],s[7]],[s[8],s[7],s[6]]])
            inertia = np.trace(q)*np.eye(3)-q-s[0]*(center.dot(center)*np.eye(3)-np.outer(center,center))
            mass_factor = spec['density_kg_m3']*spec['model_unit_m']**3
            status = ('INPUT_CAPACITY_EXCLUDED' if capacity['status']=='EXCLUDED_BY_FIRST_MOMENT_BOUND'
                      else 'INPUT_ADMITTED_FEASIBILITY_UNKNOWN')
            record = dict(schema='noodles-exact-input-admission-v1', status=status,
                input_version=input_version,classification=spec['classification'], problem_identity=problem.identity, problem_digest=problem.input_digest,
                exact_region_facts={name:snapshot.exact_region().facts for name,snapshot in snapshots.items()},
                outer_connected_boundary_components=count, independent_moments_exact={name:list(map(str,v)) for name,v in values.items()},
                partition_moment_delta_exact=list(map(str,delta)), capacity_status=capacity['status'],
                mass_kg=float(s[0]*mass_factor), center_m=(center*spec['model_unit_m']).tolist(),
                inertia_center_kg_m2=(inertia*mass_factor*spec['model_unit_m']**2).tolist(),
                measured_input_isolation_interface_bound_model=measured, printing_wall_thickness='UNKNOWN',
                si_exports_valid=all(x['native_facts']['completed'] and x['native_facts']['valid'] and x['moments_pass'] for x in si_audits.values()),
                collision_admission='NOT_GRANTED', native_path=native.__file__,
                native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
                optimization_invocations=0, simulation_invocations=0)
            publish(output/'admission.json', record)
            event(status)
            publish(output/'terminal.json', dict(status=status, completed=True, elapsed_seconds=time.monotonic()-started,
                                                optimization_invocations=0, simulation_invocations=0))
        except Exception as exc:
            publish(output/'failure.json', dict(type=type(exc).__name__, message=str(exc), completed=False,
                                               elapsed_seconds=time.monotonic()-started))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True); args=parser.parse_args()
    admit(args.input,args.output)
