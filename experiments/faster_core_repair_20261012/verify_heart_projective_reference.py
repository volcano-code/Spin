"""Independent saved-state reference for the first failed Heart direction.

Frozen recipes only: no native calls, Boolean, new finite-difference samples,
or change to the failed five-point protocol. This is diagnostic evidence and
does not certify a parameter neighborhood or admit the complete case.
"""
import hashlib
import json
from pathlib import Path
import sys
import time
import mpmath as mp
from independent_projective_reference import projective_points
from independent_directional_reference import Dual,line_plane,tetra_moments,classify

HERE=Path(__file__).resolve().parent
load=lambda p:json.loads(p.read_text())


def main(output):
    capture=HERE/'heart-projective-reference-source-job-01/output'
    source=load(capture/'source.json'); replay=load(capture/'report.json')
    saved_path=HERE/'heart-cherchi-heart-first-small-face-rejection-derivative-job-01/output/point-000-material.json'
    registration=HERE/'heart-cherchi-coordinate-map-job-01/output/registration.json'
    if not replay['passed'] or replay['material_sha256']!=hashlib.sha256(saved_path.read_bytes()).hexdigest():
        raise ValueError('REPLAY_MATERIAL_IDENTITY_MISMATCH')
    if replay['registration_sha256']!=hashlib.sha256(registration.read_bytes()).hexdigest():
        raise ValueError('REPLAY_REGISTRATION_IDENTITY_MISMATCH')
    material=load(saved_path);config=load(registration);request=config['normalized']
    if not source['paper_projective_hyperboloid'] or source['effective_icosphere_subdivisions']!=2:
        raise ValueError('REFERENCE_SCOPE_PROJECTIVE_LEVEL2_ONLY')
    if request['target_axes']!=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]:
        raise ValueError('REFERENCE_FIXTURE_FRAME_NOT_IDENTITY')
    normal_path=HERE/'core-registered-failure-capture-job-01/output/source_side.json'
    normals=load(normal_path)['snapshots'][0]['source_fixed_normals'][:162]
    records=sorted(material['material_proof']['resolved_construction_records'],key=lambda r:r['vertex_index'])
    if [r['vertex_index'] for r in records]!=list(range(len(records))):
        raise ValueError('NONCONTIGUOUS_CONSTRUCTION_RECORDS')
    output.mkdir(parents=True,exist_ok=False)
    protocol=dict(precisions=[70,110],direction=0,old_five_point_gate_unchanged=True,
        geometry_tolerance=1e-11,derivative_tolerance=1e-6,precision_agreement=1e-20,
        source_sha256=hashlib.sha256((capture/'source.json').read_bytes()).hexdigest(),
        normals_sha256=hashlib.sha256(normal_path.read_bytes()).hexdigest(),
        material_sha256=hashlib.sha256(saved_path.read_bytes()).hexdigest(),case_admission=False)
    (output/'protocol.json').write_text(json.dumps(protocol))
    answers=[]
    for precision in protocol['precisions']:
        start=time.monotonic()
        with mp.workdps(precision):
            q=list(map(mp.mpf,material['quadric_coefficients']))
            def geometry(t):
                coefficients=[q[0]+t,*q[1:]]
                return mp.matrix(projective_points(coefficients,normals,source['source_normal_vertices'],
                    source['source_projective_rays'],source['interface_vertex_count'],
                    source['faces'],source['truncation_radius']))
            points=geometry(0);motion=mp.diff(geometry,mp.mpf(0))
            source_error=max(abs(points[i,a]-source['vertices'][i][a]) for i in range(points.rows) for a in range(3))
            if source_error>mp.mpf('1e-11'):
                (output/'failure.json').write_text(json.dumps(dict(stage='SOURCE_RECONSTRUCTION',error=str(source_error),dps=precision)))
                raise ValueError('INDEPENDENT_SOURCE_RECONSTRUCTION_MISMATCH')
            length=mp.mpf(config['length']);origin=list(map(mp.mpf,config['origin']))
            target=list(map(mp.mpf,request['target_origin']))
            physical=[[Dual((points[i,a]+target[a])*length+origin[a],motion[i,a]*length)
                       for a in range(3)] for i in range(points.rows)]
            objects=[[[Dual(x) for x in p] for p in config['original']['design_vertices']],physical]
            vertices=[]
            for record in records:
                operands=[]
                for operand in record['operands']:
                    if operand['status']!='UNIQUE' or len(operand['source_vertices'])!=1:
                        raise ValueError('AMBIGUOUS_CONSTRUCTION_SOURCE')
                    item=operand['source_vertices'][0]
                    operands.append(objects[item['input_object']][item['input_vertex']])
                if record['point_type']=='EXPLICIT3D' and len(operands)==1:p=operands[0]
                elif record['point_type']=='LPI' and len(operands)==5:p=line_plane(operands)
                else:raise ValueError('UNSUPPORTED_CONSTRUCTION_TYPE')
                vertices.append([(p[a]-origin[a])/length-target[a] for a in range(3)])
            material_error=max(abs(p[a].v-mp.mpf(material['variable_rational_vertices'][i][a]))
                               for i,p in enumerate(vertices) for a in range(3))
            if material_error>mp.mpf('1e-11'):raise ValueError('INDEPENDENT_MATERIAL_RECONSTRUCTION_MISMATCH')
            moments=tetra_moments(vertices,material['variable_faces'],mp.mpf(request['material_density']))
            answer=dict(dps=precision,source_error=str(source_error),material_error=str(material_error),
                derivative=[mp.nstr(x.d,precision) for x in moments],seconds=time.monotonic()-start)
        (output/f'reference-{precision}.json').write_text(json.dumps(answer));answers.append(answer)
        print('REFERENCE_FINISHED',precision,answer['seconds'],flush=True)
    with mp.workdps(110):
        production=[[row[0]] for row in material['moment_jacobian']]
        low=[[x] for x in answers[0]['derivative']];high=[[x] for x in answers[1]['derivative']]
        status=classify(production,low,high,True)
        wrong=[list(row) for row in production];wrong[0][0]+=1e-3
        negative=classify(wrong,low,high,True)
        error=max(abs(a-b) for a,b in zip(mp.matrix(production),mp.matrix(high)))
    report=dict(status=status,wrong_jacobian_status=negative,maximum_derivative_error=str(error),
        direction=0,old_five_point_gate_passed=False,case_admission=False,
        neighborhood_certification=False,native_calls=0,boolean_calls=0,
        scope='POINTWISE_FROZEN_PROJECTIVE_RECIPE_REFERENCE_NOT_BRANCH_CERTIFICATION')
    (output/'report.json').write_text(json.dumps(report));print(json.dumps(report),flush=True)
    if status!='SUPPORTED' or negative!='INCONSISTENT':raise SystemExit(1)


if __name__=='__main__':main(Path(sys.argv[1]))
