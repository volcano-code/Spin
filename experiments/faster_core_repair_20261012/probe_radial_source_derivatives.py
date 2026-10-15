"""Original 81-point source protocol for the registered seed1 failure."""
import json
from pathlib import Path
import sys
import numpy as np
import _make_it_stand as native
def weights(nodes):
    nodes=np.asarray(nodes,dtype=float);scale=max(abs(nodes))
    matrix=np.vstack([(nodes/scale)**degree for degree in range(5)])
    return np.linalg.solve(matrix,np.array([0.,1.,0.,0.,0.]))/scale


def main(output):
    root=Path(__file__).resolve().parent
    config=json.loads((root/'heart-certified-formal-seed-1-job-01/output/run/configuration.json').read_text())
    request=config['request'];settings=native.SpinTriangulationSettings()
    for key in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
                'retain_fixed_normal_provenance','truncation_radius','paper_projective_hyperboloid','algebraic_finite_map'):
        setattr(settings,key,config['evaluator'][key])
    design=(np.asarray(request['design_vertices'])-request['target_origin'])@np.asarray(request['target_axes']).T
    base=np.array(json.loads((root/'heart-seed1-source-failure-registration-01.json').read_text())['event']['quadric_coefficients'])
    samples=[base.copy()];lookup={tuple(base):0};layouts=[]
    for direction in range(10):
        for h in (3e-6,1.5e-6,.75e-6):
            indices=[];positions=[]
            for shift in (-h,-h/2,0,h/2,h):
                q=base.copy();q[direction]+=shift
                if tuple(q) not in lookup:lookup[tuple(q)]=len(samples);samples.append(q)
                indices.append(lookup[tuple(q)]);positions.append(q[direction]-base[direction])
            layouts.append((direction,h,indices,positions))
    assert len(samples)==81
    output=Path(output);output.mkdir(exist_ok=False,parents=True)
    records=[]
    for i,q in enumerate(samples):
        result=native.evaluate_spin_it_faster_source_motion(q,design,settings)
        if not result['defined'] or not result['side_selection_defined']:raise ValueError(result['diagnostic'])
        records.append(result)
        np.savez(output/f'point-{i:03d}.npz',coefficients=q,vertices=result['vertices'],faces=result['faces'],velocity=result['vertex_velocity'])
        print('SOURCE_POINT',i,flush=True)
    center=records[0];jac=np.asarray(center['vertex_velocity']).reshape(-1,3,10)
    stable=all(r['source_projective_rays']==center['source_projective_rays'] and
        r['truncation_radius']==center['truncation_radius'] and
        np.array_equal(r['faces'],center['faces']) for r in records)
    checks=[]
    for direction,h,indices,positions in layouts:
        fd=np.tensordot(weights(positions),np.array([records[i]['vertices'] for i in indices]),axes=1)
        error=np.linalg.norm(fd-jac[:,:,direction],axis=1)/np.maximum(1.,np.linalg.norm(jac[:,:,direction],axis=1))
        checks.append(dict(direction=direction,scale=h,maximum_relative_error=float(max(error))))
    report=dict(point_count=81,branch_stable=stable,checks=checks,
        passed=stable and all(r['maximum_relative_error']<=1e-6 for r in checks),
        production_material_derivative_admission=False)
    (output/'report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main(sys.argv[1])
