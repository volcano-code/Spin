"""Two frozen difficult states, same Cherchi map in both coordinate frames.

Diagnostic only: no seed initialization, solve, branch or neighborhood admission.
"""
from dataclasses import asdict,replace
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator,_json_value
from make_it_stand.verification.faster_material import _transform_moments
from case_coordinate_contract import normalize_target_moments,normalize_target_quadric,normalized_faster_weight

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main(output):
    output.mkdir(parents=True,exist_ok=False)
    heart=ROOT/'experiments/faster_heart_contact_20261011'
    frozen_path=heart/'rounded-clearance-v3-faster-frozen-v2/request.json'
    frozen=json.loads(frozen_path.read_text())
    for filename,expected in frozen['input_sources'].items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest()!=expected:
            raise ValueError('FROZEN_INPUT_CHANGED')
    payload=dict(frozen['requests'][0]);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    original=spin.SpinItFasterRequest(**payload)
    v=np.asarray(original.outer_vertices); low=v.min(axis=0); spans=v.max(axis=0)-low
    center=low+spans*.5; length=float(max(spans))
    fixed_shifted=_transform_moments(original.fixed_shell_moments,np.eye(3),-center)
    normalized=replace(original,
        outer_vertices=((v-center)/length).tolist(),
        design_vertices=((np.asarray(original.design_vertices)-center)/length).tolist(),
        target_origin=((np.asarray(original.target_origin)-center)/length).tolist(),
        fixed_shell_moments=normalize_target_moments(fixed_shifted,length,1.).tolist(),
        original_units_per_normalized_unit=length,
        center_of_mass_weight=normalized_faster_weight(original.center_of_mass_weight,length))
    sys.path.insert(0,str(heart))
    from core_faster_heart import paper_core_settings
    settings=spin.SpinEvaluatorConfig(**paper_core_settings(frozen['evaluator_settings']))
    (output/'registration.json').write_text(json.dumps(dict(original=asdict(original),
        normalized=asdict(normalized),settings=asdict(settings),length=length,origin=center.tolist(),
        request_sha256=hashlib.sha256(frozen_path.read_bytes()).hexdigest(),case_admission=False)))
    evaluators=[InProcessCherchiEvaluator(original,settings,rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1'),
        InProcessCherchiEvaluator(normalized,settings,normalization=dict(origin=center.tolist(),L=length),
            design_vertices_physical=original.design_vertices,rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1')]
    states=json.loads((HERE/'frozen-heart-states.json').read_text())['records']
    rows=[]
    for name in ('heart-first-cap-failure','heart-first-small-face-rejection'):
        q=np.asarray(states[name]['coefficients']);p,chain=normalize_target_quadric(q,length)
        maps=[]
        for label,evaluator,coeff in zip(('original','normalized'),evaluators,(q,p)):
            mapped=evaluator(coeff)
            (output/f'{name}-{label}.json').write_text(json.dumps(evaluator.material(coeff),default=_json_value))
            maps.append(mapped)
        scales=np.array([length**power for power in (3,4,4,4,5,5,5,5,5,5)])
        moment_error=float(np.max(np.abs(maps[0]['raw_moments']/scales-maps[1]['raw_moments'])))
        jac_error=float(np.max(np.abs(maps[1]['moment_jacobian']@chain-maps[0]['moment_jacobian']/scales[:,None])))
        row=dict(state=name,moment_error=moment_error,jacobian_chain_error=jac_error,
            passed=moment_error<=1e-10 and jac_error<=1e-6)
        rows.append(row)
        (output/'rows.json').write_text(json.dumps(rows))
        if not row['passed']:raise ValueError('HEART_CHERCHI_COORDINATE_MAP_MISMATCH')
    (output/'report.json').write_text(json.dumps(dict(passed=True,rows=rows,case_admission=False,
        scope='Two current Cherchi fixed-state coordinate comparisons, not production FD gate')))

if __name__=='__main__':main(Path(sys.argv[1]))
