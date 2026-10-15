"""Read-only pin-minus-fixed witness. It never changes or admits an input."""
import argparse
from fractions import Fraction as Q
import json
import math
from pathlib import Path
import time
import gmpy2

from admit_bounded_request import save,text,cross
from rational_request_view import sha


def load_pin_boundary(path):
    vertices=[];faces=[]
    for line in Path(path).read_text().splitlines():
        fields=line.split()
        if not fields or fields[0].startswith('#'):continue
        if fields[0]=='v':
            if len(fields)!=4:raise ValueError('INVALID_PIN_OBJ_VERTEX')
            values=[float(x) for x in fields[1:]]
            if not all(math.isfinite(x) for x in values):raise ValueError('NONFINITE_PIN_OBJ_VERTEX')
            vertices.append([text(Q.from_float(x)) for x in values])
        elif fields[0]=='f':
            if len(fields)!=4:raise ValueError('NONTRIANGULAR_PIN_OBJ')
            indices=[int(x.split('/')[0])-1 for x in fields[1:]]
            if any(i<0 or i>=len(vertices) for i in indices):raise ValueError('INVALID_PIN_OBJ_INDEX')
            faces.append(indices)
        elif fields[0] not in ('o','g','s','vn','vt','usemtl','mtllib'):
            raise ValueError('UNSUPPORTED_PIN_OBJ_RECORD')
    if not vertices or not faces:raise ValueError('EMPTY_PIN_SOURCE')
    return dict(vertices=vertices,faces=faces)


def coordinate_correspondence(pin,fixed):
    fixed_points={tuple(Q(x) for x in row) for row in fixed['vertices']}
    found=sum(tuple(Q(x) for x in row) in fixed_points for row in pin['vertices'])
    return dict(source_pin_vertex_count=len(pin['vertices']),source_pin_vertices_found_exactly=found,
                containment_proved=False,scope='VERTEX_IDENTITY_ONLY_NOT_SOLID_CONTAINMENT')


def corresponding_pin_planes(pin,authority,derived):
    """Affine plane values bound a whole corresponding triangle, not samples.

    This local witness does not decide which fixed solid bounds are active;
    only the independent exact set difference establishes containment failure.
    """
    if authority['faces']!=derived['faces']:raise ValueError('FIXED_FACE_CORRESPONDENCE_CHANGED')
    def exact_points(packet):return [tuple(gmpy2.mpq(x) for x in row) for row in packet['vertices']]
    def plane(points,face):
        a,b,c=[points[k] for k in face]
        n=cross(tuple(b[k]-a[k] for k in range(3)),tuple(c[k]-a[k] for k in range(3)))
        if not any(n):raise ValueError('ZERO_AREA_PLANE_WITNESS')
        p=(*n,-sum(n[k]*a[k] for k in range(3)))
        scale=next(x for x in p[:3] if x)
        return tuple(x/scale for x in p),p
    pv,av,dv=map(exact_points,(pin,authority,derived))
    source_planes={}
    for index,face in enumerate(pin['faces']):
        key,p=plane(pv,face)
        source_planes.setdefault(key,(index,p))
    matched=0;departed=[]
    for index,face in enumerate(authority['faces']):
        key,_=plane(av,face)
        if key not in source_planes:continue
        matched+=1;source_index,p=source_planes[key]
        values=[sum(p[k]*dv[v][k] for k in range(3))+p[3] for v in face]
        if any(values):
            departed.append(dict(authority_fixed_face_index=index,source_pin_face_index=source_index,
                signed_plane_vertex_values_exact=[text(x) for x in values],
                minimum_signed_plane_value_exact=text(min(values)),maximum_signed_plane_value_exact=text(max(values))))
    return dict(scope='EXACT_AFFINE_TRIANGLE_PLANE_VALUES_NOT_GLOBAL_CONTAINMENT',
                matched_authority_fixed_faces=matched,departed_faces=departed,containment_proved=False)


def run(failed_job,output,expected_native_sha256):
    import _make_it_stand as native
    import numpy as np
    from make_it_stand.verification.spin_problem import MeshSnapshot
    from rounded_input import exact_moments
    failed_job,output=Path(failed_job),Path(output)
    output.mkdir(parents=True,exist_ok=False)
    config=json.loads((failed_job/'config.json').read_text())
    identities=config['identity']
    keys=('authority-fixed','candidate-fixed.rational.json','pin','source_binding','original_admission','candidate','native')
    inputs={key:identities[key] for key in keys}
    began=time.monotonic();stage='SOURCE_IDENTITY'
    with (output/'events.jsonl').open('x') as stream:
        def event(name):
            nonlocal stage
            stage=name;stream.write(json.dumps(dict(stage=stage,elapsed_seconds=time.monotonic()-began))+'\n')
            stream.flush()
            import os
            os.fsync(stream.fileno());print(stage,flush=True)
        try:
            if sha(native.__file__)!=expected_native_sha256 or identities['native']['sha256']!=expected_native_sha256:
                raise ValueError('NATIVE_IDENTITY_MISMATCH')
            for info in inputs.values():
                if sha(info['path'])!=info['sha256']:raise ValueError('STALE_DIAGNOSTIC_SOURCE: '+info['path'])
            save(output/'source-identity.json',dict(schema='pin-role-witness-identity-v1',
                failed_config_sha256=sha(failed_job/'config.json'),inputs=inputs,native_sha256=sha(native.__file__),
                worker_sha256=sha(__file__),geometry_admission_granted=False,optimization_invocations=0))
            pin=load_pin_boundary(inputs['pin']['path'])
            # The diagnostic must match exactly the loader used by the failed gate.
            import trimesh
            mesh=trimesh.load(inputs['pin']['path'],force='mesh',process=False)
            loaded=dict(vertices=[[text(Q.from_float(float(x))) for x in row] for row in mesh.vertices],
                        faces=mesh.faces.tolist())
            if loaded!=pin:raise ValueError('FAILED_GATE_PIN_LOADER_COORDINATE_OR_INDEX_MISMATCH')
            fixed={name:json.loads(Path(inputs[key]['path']).read_text()) for name,key in
                   [('authority','authority-fixed'),('derived','candidate-fixed.rational.json')]}
            save(output/'pin-source.rational.json',pin)
            save(output/'coordinate-correspondence.json',{name:coordinate_correspondence(pin,p) for name,p in fixed.items()})
            event('EXACT_CORRESPONDING_PIN_PLANES')
            save(output/'corresponding-pin-plane-departures.json',corresponding_pin_planes(pin,fixed['authority'],fixed['derived']))
            event('PIN_AND_FIXED_EXACT_REGIONS')
            ps=MeshSnapshot.from_rationals(pin['vertices'],pin['faces'])
            snapshots={name:MeshSnapshot.from_rationals(p['vertices'],p['faces']) for name,p in fixed.items()}
            relations={name:native.relate_exact_material_regions(ps.exact_region(),snapshot.exact_region())
                       for name,snapshot in snapshots.items()}
            save(output/'pin-fixed-relations.json',relations)
            event('ACTUAL_PIN_MINUS_DERIVED_FIXED')
            result=native.libigl_exact_boolean_construction_from_rationals(pin['vertices'],np.asarray(pin['faces'],dtype=np.int32),
                fixed['derived']['vertices'],np.asarray(fixed['derived']['faces'],dtype=np.int32),'difference')
            if not result['succeeded']:raise ValueError('PIN_DIFFERENCE_CONSTRUCTION_FAILED')
            packet=dict(schema='pin-minus-derived-fixed-exact-witness-v1',
                vertices=[[text(gmpy2.mpq(str(x))) for x in row] for row in result['vertices']],
                faces=np.asarray(result['faces'],dtype=np.int32).tolist(),
                role='DIAGNOSTIC_DIFFERENCE_NOT_REPLACEMENT_PIN')
            save(output/'pin-minus-fixed.rational.json',packet)
            moments=exact_moments(packet['vertices'],packet['faces']) if packet['faces'] else [Q(0)]*10
            bounds=[[text(min(Q(v[k]) for v in packet['vertices'])),text(max(Q(v[k]) for v in packet['vertices']))]
                    for k in range(3)] if packet['vertices'] else []
            report=dict(schema='bounded-pin-role-diagnostic-v1',source_pin_loader_identical=True,
                authority_contains_source_pin=relations['authority']['a_subset_b'],
                derived_contains_source_pin=relations['derived']['a_subset_b'],
                exact_difference_moments=[text(x) for x in moments],difference_bounds_exact=bounds,
                difference_vertex_count=len(packet['vertices']),difference_face_count=len(packet['faces']),
                geometry_admission_granted=False,optimization_invocations=0,simulation_invocations=0,
                conclusion='NO_AUTOMATIC_REPAIR_OR_ENTITY_SUBSTITUTION',elapsed_seconds=time.monotonic()-began)
            for info in inputs.values():
                if sha(info['path'])!=info['sha256']:raise ValueError('SOURCE_CHANGED_DURING_WITNESS')
            save(output/'terminal.json',report)
            return report
        except Exception as exc:
            save(output/'failure.json',dict(stage=stage,type=type(exc).__name__,message=str(exc),inputs=inputs,
                failed_config_sha256=sha(failed_job/'config.json'),geometry_admission_granted=False,
                optimization_invocations=0,simulation_invocations=0,elapsed_seconds=time.monotonic()-began))
            raise


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--failed-job',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--expected-native-sha256',required=True)
    args=parser.parse_args();run(args.failed_job,args.output,args.expected_native_sha256)
