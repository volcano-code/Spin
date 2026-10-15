"""Registered box top pilot: support origin below the body, linear height."""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import time
from make_it_stand import spin
from cherchi_core_map import CherchiCoreMap
from run_small_inprocess_solver import safe


def main(output,seed):
    here=Path(__file__).parent
    frozen=json.loads((here/'same-input-job-03/output/frozen-inputs.json').read_text())
    payload=dict(frozen['request'])
    payload.update(mode='SPINNING_TOP',target_origin=(0.,0.,-1.1),center_of_mass_weight=1.)
    config=dict(payload['settings']);config.update(equation_profile='DUAL_MONOMIAL_RAW_V1',deterministic_seed=seed,
        maximum_iterations=60,maximum_restarts=8,maximum_function_evaluations_per_restart=60,capture_evaluation_payload=True)
    payload['settings']=spin.SpinSolverConfig(**config)
    request=spin.SpinItFasterRequest(**payload);settings=spin.SpinEvaluatorConfig(**frozen['settings'])
    output.mkdir(parents=True,exist_ok=False)
    with (output/'configuration.json').open('x') as stream:
        json.dump(safe(dict(request=asdict(request),evaluator=asdict(settings),
            role='NEW_SMALL_FIXED_MATERIAL_SPINNING_TOP_METHOD_PILOT',initialization='REQUEST_AND_SEED_ONLY')),stream,allow_nan=False)
    mapping=CherchiCoreMap(request,settings)
    with (output/'events.jsonl').open('x') as stream:
        def sink(row):
            stream.write(json.dumps(safe(row),allow_nan=False)+'\n');stream.flush();os.fsync(stream.fileno())
        start=time.perf_counter()
        result=spin.solve_spin_it_faster_with_external_moment_evaluator(request,mapping,event_sink=sink)
        elapsed=time.perf_counter()-start
    with (output/'numerical-result.json').open('x') as stream:
        json.dump(safe(dict(result=result,seconds=elapsed,independent_material_admission=False,
            derivative_admission=False,physics_admission=False)),stream,allow_nan=False)
    if result['status']=='CONVERGED':
        last=mapping.last
        with (output/'last-material.json').open('x') as stream:
            json.dump(safe(dict(vertices=last['vertices'],faces=last['faces'],
                exact_vertices=[[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in last['exact_vertices']],
                geometry_hash=last['geometry_hash'],topology_hash=last['topology_hash'],material_admission=False)),stream,allow_nan=False)
    print(json.dumps(dict(seed=seed,status=result['status'],seconds=elapsed,
        evaluations=result['function_evaluations'],raw_residual=safe(result.get('raw_equation_residual')))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--seed',type=int,choices=(0,1),required=True)
    args=p.parse_args();main(args.output,args.seed)
