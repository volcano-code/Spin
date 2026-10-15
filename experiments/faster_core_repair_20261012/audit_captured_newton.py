"""Three-evaluation registered prefix, native matrix vs independent symbolic Eq14."""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import numpy as np
from make_it_stand import spin
from cherchi_core_map import CherchiCoreMap
from independent_eq14 import formulas, dual
from run_small_inprocess_solver import safe


def reconstruct(event,functions):
    f,df,hf,g,dg,hg=functions
    payload=event['evaluation_payload'];r=np.asarray(payload['solver_raw_moments'])
    p=dual(event['quadric_coefficients']);J=np.asarray(payload['solver_moment_jacobian'])
    grad=np.asarray(df(*r)).ravel();G=np.asarray(dg(*r));constraints=np.asarray(g(*r)).ravel()
    coeff=np.linalg.lstsq(np.column_stack((-G.T,p)),-grad,rcond=None)[0]
    lam,mu=coeff[:-1],coeff[-1]
    H=np.asarray(hf(*r))-np.einsum('i,ijk->jk',lam,np.asarray(hg(*r)))
    A=np.block([[H@J+mu*np.eye(10),-G.T,p[:,None]],
                [G@J,np.zeros((4,4)),np.zeros((4,1))],
                [2*p[None,:],np.zeros((1,4)),np.zeros((1,1))]])
    rhs=np.r_[-(grad+mu*p-G.T@lam),-constraints,0.]
    return A,rhs,lam,mu


def main(output,saved=None):
    here=Path(__file__).parent
    config=json.loads((here/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json').read_text())
    payload=dict(config['request']);solver=dict(payload['settings'])
    solver.update(maximum_iterations=3,maximum_restarts=1,maximum_function_evaluations_per_restart=3)
    payload['settings']=spin.SpinSolverConfig(**solver)
    request=spin.SpinItFasterRequest(**payload);settings=spin.SpinEvaluatorConfig(**config['evaluator'])
    output.mkdir(parents=True,exist_ok=False)
    with (output/'configuration.json').open('x') as stream:
        json.dump(safe(dict(role='THREE_EVALUATION_NEWTON_MATRIX_DIAGNOSTIC_NOT_FULL_OPTIMIZATION',
                           request=asdict(request),evaluator=config['evaluator'])),stream)
    events=[]
    if saved is None:
        with (output/'events.jsonl').open('x') as stream:
            def sink(row):
                events.append(row);stream.write(json.dumps(safe(row),allow_nan=False)+'\n');stream.flush();os.fsync(stream.fileno())
            result=spin.solve_spin_it_faster_with_external_moment_evaluator(request,CherchiCoreMap(request,settings),event_sink=sink)
        evaluations=result['function_evaluations']
    else:
        events=[json.loads(line) for line in (saved/'events.jsonl').read_text().splitlines()]
        evaluations=json.loads((saved/'report.json').read_text())['evaluations']
    functions=formulas();checks=[];accepted=None
    for event in events:
        if event['type']=='evaluation' and event['outcome']=='accepted' and event.get('evaluation_payload'):
            accepted=event
        captured=event.get('newton_system_payload')
        if not captured:continue
        if accepted is None or accepted['evaluation_index']!=event['evaluation_index']:
            raise ValueError('NEWTON_PAYLOAD_STATE_MISMATCH')
        A,rhs,lam,mu=reconstruct(accepted,functions)
        errors=dict(matrix_error=float(np.max(np.abs(A-captured['matrix']))),
            rhs_error=float(np.max(np.abs(rhs-captured['right_hand_side']))),
            multiplier_error=float(np.max(np.abs(lam-captured['constraint_multipliers']))),
            mu_error=abs(float(mu)-captured['multiplier_mu']))
        checks.append(dict(evaluation_index=event['evaluation_index'],errors=errors,
                           passed=all(value<=1e-10 for value in errors.values())))
    baseline_path=here/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/events.jsonl'
    baseline=[json.loads(line) for line in baseline_path.read_text().splitlines()]
    old=[r for r in baseline if r['type']=='evaluation' and r['restart']==0][:3]
    new=[r for r in events if r['type']=='evaluation']
    numerical_same=len(new)==len(old)==3
    for a,b in zip(old,new):
        for k in ('quadric_coefficients','objective','raw_equation_residual','multiplier_mu'):
            numerical_same=numerical_same and bool(np.array_equal(np.asarray(a[k]),np.asarray(b[k])))
        for k in ('solver_coordinates','solver_raw_moments','variable_raw_moments','solver_moment_jacobian'):
            numerical_same=numerical_same and bool(np.array_equal(np.asarray(a['evaluation_payload'][k]),np.asarray(b['evaluation_payload'][k])))
    # The loop assembles the last accepted state's system before its next
    # backtracking loop notices the evaluation budget; three systems are valid.
    report=dict(checks=checks,evaluations=evaluations,prefix_scientific_values_equal=numerical_same,
                native_matrix_capture_count=len(checks),passed=len(checks)==3 and numerical_same and all(c['passed'] for c in checks),
                saved_event_recheck=str(saved) if saved else None,additional_mapping_calls=0 if saved else evaluations,
                complete_optimization_admission=False,scope='Native regular 15x15 systems from three-evaluation top prefix only')
    with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(report),flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--saved',type=Path)
    args=p.parse_args();main(args.output,args.saved)
