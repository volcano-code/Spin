"""One preregistered hypothesis: heuristic material intersection false positive."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time


def main(run,output):
    from make_it_stand.spin import (SpinItFasterRequest,SpinSolverConfig,SpinEvaluatorConfig,
        evaluate_spin_it_faster_fullrank_quadric)
    from rounded_solver_journal import save
    run,output=Path(run),Path(output);output.mkdir(parents=True,exist_ok=False)
    config=json.loads((run/'mapping-config.json').read_text())
    rows=[json.loads(row) for row in (run/'mapping-events.jsonl').read_text().splitlines()]
    failed=next(row for row in rows if row['status']=='FAILED')
    payload=dict(config['request']);payload['settings']=SpinSolverConfig(**payload['settings'])
    request=SpinItFasterRequest(**payload)
    settings=dict(config['evaluator']);settings['exact_binary64_self_intersection']=True
    started=time.monotonic()
    save(output/'hypothesis.json',dict(hypothesis='HEURISTIC_INTERSECTION_FALSE_POSITIVE',
        coefficients=failed['coefficients'],old_error=failed['exception_message'],evaluator_settings=settings,
        optimizer_invocations=0))
    try:
        result=evaluate_spin_it_faster_fullrank_quadric(request,failed['coefficients'],
            evaluator_settings=SpinEvaluatorConfig(**settings))
        save(output/'evaluation.json',asdict(result))
        save(output/'result.json',dict(actual_mapping_returned=True,derivative_defined=result.derivative_defined,
            derivative_diagnostic=result.derivative_diagnostic,seconds=time.monotonic()-started,
            geometry_admission=False,optimization_admission=False))
    except Exception as error:
        save(output/'result.json',dict(actual_mapping_returned=False,exception_type=type(error).__name__,
            error=str(error),seconds=time.monotonic()-started,geometry_admission=False,optimization_admission=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    main(args.run,args.output)
