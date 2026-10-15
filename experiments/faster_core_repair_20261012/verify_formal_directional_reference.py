"""Saved-state multiprecision reference: no native/Boolean/Jacobian invocation."""
import hashlib
import json
from pathlib import Path
import sys
import time
import mpmath as mp
from independent_directional_reference import (Dual,finite_point,root_point_and_jacobian,
                                              line_plane,tetra_moments,cross,classify)

HERE=Path(__file__).resolve().parent
INPUT=HERE/'formal-cherchi-derivative-job-01/output'
load=lambda p:json.loads(p.read_text())
base=load(INPUT/'point-000-material.json');material=base['same_call_material'];branch=base['source_branch']
config=load(HERE/'formal-cherchi-small-closure-job-04/output/run/configuration.json')
result=load(HERE/'formal-cherchi-small-closure-job-04/output/run/result.json')
norm=result['identity']['normalization'];request=config['request']
physical_design=load(HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json')['request']['design_vertices']
# The first 162 immutable level2 normals are shared input geometry, not derivatives.
normal_tape=load(HERE/'core-registered-failure-capture-job-01/output/source_side.json')['snapshots'][0]['source_fixed_normals'][:162]
certificate=load(HERE/'formal-cherchi-branch-job-01/output/report.json')
if not certificate['finite_sample_certification']:raise ValueError('FINITE_SAMPLE_BRANCH_CERTIFICATE_REQUIRED')
if certificate['manifest_sha256']!=hashlib.sha256((INPUT/'manifest.json').read_bytes()).hexdigest():raise ValueError('MANIFEST_MISMATCH')
if branch['projective_ray_ids']:raise ValueError('REFERENCE_SCOPE_LEGACY_TRUNCATION_ONLY')
records=sorted(material['material_proof']['resolved_construction_records'],key=lambda r:r['vertex_index'])
if any(r['point_type'] not in ('EXPLICIT3D','LPI') for r in records):raise ValueError('REFERENCE_CONSTRUCTION_UNSUPPORTED')
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
with (out/'protocol.json').open('x') as stream:
    json.dump(dict(version='INDEPENDENT_IMPLICIT_TETRA_REFERENCE_V1',precisions=[70,110],
        directions=list(range(10)),jacobian_tolerance=1e-6,precision_agreement=1e-20,
        source_reconstruction_tolerance=1e-11,material_reconstruction_tolerance=1e-11,
        maximum_constraint_condition=1e12,old_five_point_gate_unchanged=True,
        scope='Saved normalized small legacy hyperboloid state only; not projective-path admission'),stream)

def compute(dps):
    with mp.workdps(dps):
        q=list(map(mp.mpf,base['actual_coefficients']))
        length=mp.mpf(norm['L']);origin=list(map(mp.mpf,norm['origin']))
        target=list(map(mp.mpf,request['target_origin']));axes=mp.matrix(request['target_axes'])
        source_actual=[list(mp.matrix([(mp.mpf(v[i])-origin[i])/length for i in range(3)]).T*axes)
                       for v in material['quadric_source_vertices']]
        source_actual=[[v[i]-target[i] for i in range(3)] for v in source_actual]
        # Above target subtraction is valid for this registered identity frame only.
        if axes!=mp.eye(3):raise ValueError('REFERENCE_FIXTURE_FRAME_NOT_IDENTITY')
        points={};jacobians={};max_condition=mp.mpf(0)
        cap=[]
        for i,(nid,edge) in enumerate(zip(branch['normal_vertices'],branch['normal_edges'])):
            if nid>=0:
                n=list(map(mp.mpf,normal_tape[nid]));points[i]=finite_point(q,n)
                jacobians[i]=mp.matrix([[mp.diff(lambda t:finite_point([q[k]+(t if k==j else 0) for k in range(10)],n)[a],0)
                                        for j in range(10)] for a in range(3)])
            elif edge[0]>=0:
                plane=cross(list(map(mp.mpf,normal_tape[edge[0]])),list(map(mp.mpf,normal_tape[edge[1]])))
                p,J,condition=root_point_and_jacobian(q,plane,mp.mpf(branch['radius']),source_actual[i])
                max_condition=max(max_condition,condition);points[i]=p;jacobians[i]=J
            else:cap.append(i)
        for i in cap:
            neighbors=sorted({v for f in branch['faces'] if i in f for v in f if v!=i})
            if len(neighbors)<3 or any(v not in points for v in neighbors):raise ValueError('UNSUPPORTED_CAP_RECIPE')
            points[i]=[sum(points[v][a] for v in neighbors)/len(neighbors) for a in range(3)]
            jacobians[i]=sum((jacobians[v] for v in neighbors),mp.zeros(3,10))/len(neighbors)
        source_error=max(abs(points[i][a]-source_actual[i][a]) for i in points for a in range(3))
        columns=[];material_error=mp.mpf(0)
        for direction in range(10):
            source_physical=[[Dual((points[i][a]+target[a])*length+origin[a],jacobians[i][a,direction]*length)
                              for a in range(3)] for i in range(len(points))]
            objects=[[[Dual(x) for x in p] for p in physical_design],source_physical]
            vertices=[]
            for record in records:
                operands=[]
                for item in record['operands']:
                    if item['status']!='UNIQUE' or len(item['source_vertices'])!=1:raise ValueError('AMBIGUOUS_SOURCE')
                    item=item['source_vertices'][0];operands.append(objects[item['input_object']][item['input_vertex']])
                p=operands[0] if record['point_type']=='EXPLICIT3D' else line_plane(operands)
                vertices.append([(p[a]-origin[a])/length-target[a] for a in range(3)])
            for i,p in enumerate(vertices):
                for a in range(3):
                    num,den=material['variable_rational_vertices'][i][a].split('/')
                    material_error=max(material_error,abs(p[a].v-mp.mpf(num)/mp.mpf(den)))
            moments=tetra_moments(vertices,material['variable_faces'],mp.mpf(request['material_density']))
            columns.append([mp.nstr(x.d,dps) for x in moments])
        return dict(dps=dps,jacobian=[[columns[j][i] for j in range(10)] for i in range(10)],
            source_error=mp.nstr(source_error,dps),material_error=mp.nstr(material_error,dps),
            maximum_constraint_condition=mp.nstr(max_condition,dps))

answers=[]
for precision in (70,110):
    start=time.monotonic();answer=compute(precision);answer['seconds']=time.monotonic()-start
    with (out/f'reference-{precision}.json').open('x') as stream:json.dump(answer,stream)
    answers.append(answer);print('REFERENCE_PRECISION_FINISHED',precision,answer['seconds'],flush=True)
with mp.workdps(110):
    valid=all(mp.mpf(a['source_error'])<=mp.mpf('1e-11') and mp.mpf(a['material_error'])<=mp.mpf('1e-11')
              and mp.mpf(a['maximum_constraint_condition'])<=mp.mpf('1e12') for a in answers)
    J=material['moment_jacobian'];low=answers[0]['jacobian'];high=answers[1]['jacobian']
    status=classify(J,low,high,valid)
    corrupted=[list(row) for row in J];corrupted[0][0]+=1e-3
    negative=classify(corrupted,low,high,valid)
    error=max(abs(a-b) for a,b in zip(mp.matrix(J),mp.matrix(high)))
    report=dict(status=status,geometry_and_branch_valid=valid,maximum_derivative_error=str(error),
        wrong_jacobian_status=negative,old_five_point_gate_passed=False,
        core_admission=False,scope='Independent implicit source + scalar dual tetrahedra, saved state only',
        native_calls=0,boolean_calls=0,production_derivative_calls=0,
        base_sha256=hashlib.sha256((INPUT/'point-000-material.json').read_bytes()).hexdigest())
with (out/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(report),flush=True)
if status!='SUPPORTED' or negative!='INCONSISTENT':raise SystemExit(1)
