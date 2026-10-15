"""Pre-registered five-point source-coordinate differential, not a moment gate."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
from make_it_stand import spin
from probe_inprocess_derivatives import weights
from run_small_inprocess_solver import safe

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def main(output):
    frozen=json.loads((ROOT/'experiments/faster_heart_contact_20261011'
        /'rounded-clearance-v3-faster-frozen-v2/request.json').read_text())
    r=frozen['requests'][0]
    design=(np.asarray(r['design_vertices'])-r['target_origin'])@np.asarray(r['target_axes'])
    settings=spin._native_module().SpinTriangulationSettings()
    for key in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
                'retain_fixed_normal_provenance','truncation_radius','algebraic_finite_map'):
        setattr(settings,key,frozen['evaluator_settings'][key])
    settings.paper_projective_hyperboloid=True
    states=json.loads((HERE/'frozen-heart-states.json').read_text())['records']
    output.mkdir(parents=True,exist_ok=False)
    manifest=dict(scope='SOURCE_VERTEX_DIFFERENTIAL_NOT_PRODUCTION_MOMENT_JACOBIAN',
        base_step=3e-6,scales=[3e-6,1.5e-6,.75e-6],layout='five-point',
        states={name:states[name]['coefficients'] for name in ('heart-first-cap-failure','heart-first-small-face-rejection')},
        native_sha256=hashlib.sha256(Path(spin._native_module().__file__).read_bytes()).hexdigest())
    with (output/'manifest.json').open('x') as stream:json.dump(manifest,stream)
    reports=[]
    with (output/'events.jsonl').open('x') as stream:
        def event(row):stream.write(json.dumps(safe(row),allow_nan=False)+'\n');stream.flush();os.fsync(stream.fileno())
        for name,raw in manifest['states'].items():
            base=np.asarray(raw);samples=[base.copy()];index={tuple(base):0};layouts=[]
            for direction in range(10):
                for H in manifest['scales']:
                    ids=[];positions=[]
                    for shift in (-H,-H/2,0,H/2,H):
                        q=base.copy();q[direction]+=shift;key=tuple(q)
                        if key not in index:index[key]=len(samples);samples.append(q)
                        ids.append(index[key]);positions.append(float(q[direction]-base[direction]))
                    layouts.append(dict(direction=direction,scale=H,indices=ids,positions=positions))
            if len(samples)!=81:raise ValueError('SOURCE_POINT_COUNT_MISMATCH')
            records=[];started=time.monotonic()
            for i,q in enumerate(samples):
                try:
                    result=spin._native_module().evaluate_spin_it_faster_source_motion(q,design,settings)
                    if not result['defined'] or not result['side_selection_defined']:raise ValueError(result['diagnostic'])
                    with (output/f'{name}-point-{i:03d}.json').open('x') as saved:
                        json.dump(safe(dict(coefficients=q,source=result)),saved,allow_nan=False)
                    records.append(result);event(dict(state=name,index=i,status='SOURCE_RETURNED',seconds=time.monotonic()-started))
                except Exception as exc:
                    event(dict(state=name,index=i,status='SOURCE_FAILED',exception=str(exc)));raise
            centre=records[0];jac=np.asarray(centre['vertex_velocity']).reshape(len(centre['vertices']),3,10)
            topology=True;cap_stable=True
            for record in records:
                if record['truncation_radius']!=centre['truncation_radius'] or \
                    record['source_projective_rays']!=centre['source_projective_rays'] or \
                    record['source_normal_vertices']!=centre['source_normal_vertices']:
                    topology=False
                interface=lambda row:np.asarray(row['faces'])[np.asarray(row['face_is_cap'])==0].tolist()
                topology&=interface(record)==interface(centre)
                cap_stable&=np.asarray(record['faces']).tolist()==np.asarray(centre['faces']).tolist()
            checks=[]
            for layout in layouts:
                w=weights(layout['positions'])
                numerical=np.tensordot(w,np.asarray([records[i]['vertices'] for i in layout['indices']]),axes=1)
                error=np.linalg.norm(numerical-jac[:,:,layout['direction']],axis=1)
                scale=np.maximum(1.,np.linalg.norm(jac[:,:,layout['direction']],axis=1))
                checks.append(dict(**layout,maximum_relative_error=float(np.max(error/scale)),
                                   passed=bool(np.all(error/scale<=1e-6))))
            reports.append(dict(state=name,point_count=81,checks=checks,interface_and_point_provenance_stable=topology,
                cap_triangulation_stable=cap_stable,passed=topology and all(row['passed'] for row in checks)))
    with (output/'report.json').open('x') as stream:
        json.dump(dict(records=reports,passed=all(r['passed'] for r in reports),
            production_moment_derivative_admission=False,scope=manifest['scope']),stream,indent=2)
    print(json.dumps([dict(state=r['state'],passed=r['passed'],
        maximum_relative_error=max(x['maximum_relative_error'] for x in r['checks'])) for r in reports]),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
