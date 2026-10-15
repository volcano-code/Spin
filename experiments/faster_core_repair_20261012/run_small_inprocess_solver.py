"""Seed-only numerical pilot using the actual in-process Cherchi map."""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import time
import numpy as np
from make_it_stand import spin
from cherchi_core_map import CherchiCoreMap


def safe(value):
    if isinstance(value,np.ndarray):return safe(value.tolist())
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {str(k):safe(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [safe(v) for v in value]
    if isinstance(value,float) and not np.isfinite(value):return dict(nonfinite=str(value))
    return value


def main(output,seed):
    here=Path(__file__).parent
    frozen=json.loads((here/'same-input-job-03/output/frozen-inputs.json').read_text())
    payload=dict(frozen['request'])
    # New method-pilot configuration fixed before its first run. This is not
    # the historical successful request or an author case.
    payload.update(mode='YO_YO',center_of_mass_weight=0.)
    settings=dict(payload['settings']);settings.update(equation_profile='DUAL_MONOMIAL_RAW_V1',
        deterministic_seed=seed,maximum_iterations=60,maximum_restarts=8,
        maximum_function_evaluations_per_restart=60,capture_evaluation_payload=True)
    payload['settings']=spin.SpinSolverConfig(**settings)
    request=spin.SpinItFasterRequest(**payload)
    evaluator=spin.SpinEvaluatorConfig(**frozen['settings'])
    output.mkdir(parents=True,exist_ok=False)
    with (output/'configuration.json').open('x') as stream:
        json.dump(safe(dict(request=asdict(request),evaluator=asdict(evaluator),
            role='NEW_SMALL_NONZERO_FIXED_MATERIAL_NUMERICAL_METHOD_PILOT',initialization='REQUEST_AND_SEED_ONLY')),
            stream,allow_nan=False)
    mapping=CherchiCoreMap(request,evaluator)
    with (output/'events.jsonl').open('x') as events:
        def sink(row):
            events.write(json.dumps(safe(row),allow_nan=False)+'\n');events.flush();os.fsync(events.fileno())
        started=time.perf_counter()
        result=spin.solve_spin_it_faster_with_external_moment_evaluator(request,mapping,event_sink=sink)
        elapsed=time.perf_counter()-started
    with (output/'numerical-result.json').open('x') as stream:
        json.dump(safe(dict(result=result,seconds=elapsed,independent_material_admission=False,
            derivative_admission=False,physics_admission=False)),stream,allow_nan=False)
    if result['status']=='CONVERGED' and mapping.last is not None:
        last=mapping.last
        with (output/'last-material.json').open('x') as stream:
            json.dump(safe(dict(vertices=last['vertices'],faces=last['faces'],
                exact_vertices=[[str(v.numerator)+'/'+str(v.denominator) for v in row] for row in last['exact_vertices']],
                geometry_hash=last['geometry_hash'],topology_hash=last['topology_hash'],
                material_admission=False)),stream,allow_nan=False)
    print(json.dumps(dict(seed=seed,status=result['status'],seconds=elapsed,
        function_evaluations=result['function_evaluations'],raw_residual=safe(result.get('raw_equation_residual')))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--seed',type=int,choices=(0,1),required=True)
    args=p.parse_args();main(args.output,args.seed)
