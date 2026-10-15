"""Registered fixed states, original three scales; independent saved-material integration."""
import json
import sys
from pathlib import Path
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator, _json_value
from independent_tetra_moments import integrate
from probe_formal_cherchi_derivatives import weights

HERE=Path(__file__).resolve().parent


def main(output, frozen=None, normalize=True):
    output.mkdir(parents=True,exist_ok=False)
    if frozen is None:
        frozen=json.loads((HERE/'same-input-job-03/output/frozen-inputs.json').read_text())
    payload=dict(frozen['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload)
    settings=spin.SpinEvaluatorConfig(**frozen['settings'])
    reports=[]
    for name,coefficients in frozen['states'].items():
        folder=output/name;folder.mkdir()
        base=np.asarray(coefficients,dtype=float)
        if normalize:base/=np.linalg.norm(base)
        samples=[base.copy()];lookup={tuple(base):0};layouts=[]
        for direction in range(10):
            for scale in (1e-4,5e-5,2.5e-5):
                indices=[];positions=[]
                for offset in (-scale,-scale/2,0.,scale/2,scale):
                    q=base.copy();q[direction]+=offset
                    key=tuple(q)
                    if key not in lookup:lookup[key]=len(samples);samples.append(q)
                    indices.append(lookup[key]);positions.append(float(q[direction]-base[direction]))
                layouts.append(dict(direction=direction,scale=scale,indices=indices,positions=positions))
        if len(samples)!=81:raise ValueError('POINT_COUNT')
        (folder/'manifest.json').write_text(json.dumps(dict(samples=[q.tolist() for q in samples],layouts=layouts)))
        adapter=InProcessCherchiEvaluator(request,settings)
        values=[];errors=[];branches=[];jac=None
        with (folder/'events.jsonl').open('x') as events:
            for index,q in enumerate(samples):
                try:
                    mapped=adapter(q);last=adapter.last
                    independent=integrate(last['exact_vertices'],last['faces'],request.material_density)
                    values.append(independent);errors.append(float(np.max(np.abs(independent-mapped['variable_raw_moments']))))
                    if index==0:jac=np.asarray(mapped['moment_jacobian'])
                    source=last['source']
                    branch={key:source[key] for key in ('surface_type','truncation_radius','cap_face_count','faces','source_normal_vertices','source_normal_edges','source_projective_rays')}
                    branch=json.loads(json.dumps(branch,default=_json_value))
                    branches.append(branch)
                    (folder/f'point-{index:03d}-material.json').write_text(json.dumps(dict(
                        same_call_material=adapter.material(q),source_branch=branch),default=_json_value,allow_nan=False))
                    row=dict(index=index,status='RETURNED',integration_error=errors[-1])
                except Exception as exc:
                    events.write(json.dumps(dict(index=index,status='FAILED',error=str(exc)))+'\n');events.flush();raise
                events.write(json.dumps(row)+'\n');events.flush()
        checks=[]
        for layout in layouts:
            numerical=weights(layout['positions'])@np.asarray([values[i] for i in layout['indices']])
            error=float(np.max(np.abs(numerical-jac[:,layout['direction']])))
            checks.append(dict(**layout,maximum_error=error,passed=error<=1e-6))
        report=dict(state=name,checks=checks,passed=all(x['passed'] for x in checks),
            independent_integration_max_error=max(errors),production_jacobian=jac.tolist(),
            source_branch_stable=all(b==branches[0] for b in branches),
            material_branch_certification=False,core_admission=False)
        (folder/'report.json').write_text(json.dumps(report));reports.append(report)
    (output/'report.json').write_text(json.dumps(dict(states=reports,core_admission=False)))


if __name__=='__main__':main(Path(sys.argv[1]))
