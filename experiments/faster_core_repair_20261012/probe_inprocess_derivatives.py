"""Pre-registered 81-point actual-map five-point derivative check."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from make_it_stand import spin
from cherchi_core_map import CherchiCoreMap
from independent_tetra_moments import integrate


def weights(nodes):
    nodes=np.asarray(nodes,dtype=float)
    scale=float(np.max(np.abs(nodes)))
    if scale<=0 or len(np.unique(nodes))!=5:raise ValueError('INVALID_FIVE_POINT_NODES')
    normalized=nodes/scale
    matrix=np.vstack([normalized**degree for degree in range(5)])
    unit=np.linalg.solve(matrix,np.array((0.,1.,0.,0.,0.)))
    if np.max(np.abs(matrix@unit-np.array((0.,1.,0.,0.,0.))))>1e-12:
        raise ValueError('DIFFERENCE_WEIGHT_POLYNOMIAL_CHECK_FAILED')
    return unit/scale


def main(run,output):
    config=json.loads((run/'configuration.json').read_text())
    payload=dict(config['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload);settings=spin.SpinEvaluatorConfig(**config['evaluator'])
    terminal=json.loads((run/'numerical-result.json').read_text())['result']
    if terminal['status']!='CONVERGED':raise ValueError('NO_CONVERGED_BASE_STATE')
    base=np.asarray(terminal['quadric_coefficients'])
    matrix=np.array(((base[4],base[9],base[8]),(base[9],base[5],base[7]),(base[8],base[7],base[6])))
    discriminant=float(base[1:4]@np.linalg.solve(matrix,base[1:4])-4*base[0])
    if abs(discriminant)<=1e-10:raise ValueError('NEAR_DEGENERATE_BASE_NOT_REGULAR_FULLRANK_PROBE')
    output.mkdir(parents=True,exist_ok=False)
    nodes={};samples=[]
    def add(q):
        key=tuple(q)
        if key not in nodes:nodes[key]=len(samples);samples.append(np.asarray(q).copy())
        return nodes[key]
    base_index=add(base);layouts=[]
    for direction in range(10):
        for level in range(3):
            H=1e-4/(2**level)
            indices=[];positions=[]
            for shift in (-H,-H/2,0.,H/2,H):
                q=base.copy();q[direction]+=shift
                indices.append(add(q));positions.append(float(q[direction]-base[direction]))
            layouts.append(dict(direction=direction,scale=H,indices=indices,positions=positions))
    if len(samples)!=81:raise ValueError('UNEXPECTED_UNIQUE_POINT_COUNT')
    with (output/'manifest.json').open('x') as stream:
        json.dump(dict(role='NEW_REGISTERED_INPROCESS_FULLRANK_PRODUCTION_MAP_CHECK',base=base.tolist(),
            base_center_level_discriminant=discriminant,base_step=1e-4,layouts=layouts,
            samples=[q.tolist() for q in samples],point_count=len(samples)),stream)
    mapping=CherchiCoreMap(request,settings);moments=[];jacobian=None;branches=[];integration_errors=[]
    fixture=json.loads((Path(__file__).parent/'same-input-job-03/output/frozen-inputs.json').read_text())['fixed_material']
    rotation=np.asarray(request.target_axes).T
    fixed_vertices=(np.asarray(fixture['vertices'])-np.asarray(request.target_origin))@rotation.T
    independent_fixed=integrate(fixed_vertices,fixture['faces'],fixture['density'])
    with (output/'events.jsonl').open('x') as events:
        for index,q in enumerate(samples):
            start=time.perf_counter()
            try:
                result=mapping(q)
                last=mapping.last
                independent=integrate(last['exact_vertices'],last['faces'],request.material_density)+independent_fixed
                integration_errors.append(float(np.max(np.abs(independent-np.asarray(result['raw_moments'])))))
                moments.append(independent)
                source=last['source']
                branch=dict(surface_type=source['surface_type'],radius=source['truncation_radius'],
                    cap_count=source['cap_face_count'],faces=np.asarray(source['faces']).tolist(),
                    normal_vertices=np.asarray(source['source_normal_vertices']).tolist(),
                    normal_edges=np.asarray(source['source_normal_edges']).tolist(),
                    projective_ray_ids=np.asarray(source['source_projective_rays']).tolist())
                branches.append(branch)
                with (output/f'point-{index:03d}-material.json').open('x') as geometry_stream:
                    json.dump(dict(actual_coefficients=q.tolist(),vertices=[
                        [str(x.numerator)+'/'+str(x.denominator) for x in point] for point in last['exact_vertices']],
                        faces=last['faces'].tolist(),source_branch=branch,geometry_admission=False),geometry_stream)
                if index==base_index:jacobian=np.asarray(result['moment_jacobian'])
                row=dict(index=index,status='RETURNED_ACTUAL_MAP',seconds=time.perf_counter()-start,
                    actual_coefficients=q.tolist(),raw_moments=moments[-1].tolist(),
                    geometry_hash=result['geometry_hash'],topology_hash=result['topology_hash'])
            except Exception as exc:
                row=dict(index=index,status='FAILED',exception_type=type(exc).__name__,exception=str(exc))
                events.write(json.dumps(row)+'\n');events.flush();raise
            events.write(json.dumps(row)+'\n');events.flush()
    checks=[]
    for layout in layouts:
        w=weights(layout['positions'])
        numerical=w@np.asarray([moments[index] for index in layout['indices']])
        error=np.abs(numerical-jacobian[:,layout['direction']])
        checks.append(dict(**layout,numerical=numerical.tolist(),errors=error.tolist(),
            maximum_error=float(np.max(error)),passed=bool(np.all(error<=1e-6))))
    with (output/'report.json').open('x') as stream:
        json.dump(dict(checks=checks,point_count=len(samples),passed=all(x['passed'] for x in checks),
            production_jacobian=jacobian.tolist(),independent_integration_max_error=max(integration_errors),
            independent_material_integration_admission=max(integration_errors)<=1e-10,
            source_branch_stable=all(x==branches[0] for x in branches),branch_certification=False,
            scope='independent exact tetrahedral moments; source branch checked; material connection changes still require classification'),stream)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();main(args.run,args.output)
