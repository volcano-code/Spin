"""One saved actual Heart state for native stage accounting, no restart."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from phase_timing import PhaseTiming

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


def main(output,repetitions=1,perturb_constant=False,replay_record=None):
    from make_it_stand import spin
    heart=ROOT/'experiments/faster_heart_contact_20261011'
    sys.path.insert(0,str(heart))
    from core_faster_heart import core_settings
    from rounded_solver_journal import json_value
    frozen=json.loads((heart/'rounded-clearance-v3-faster-frozen-v2/request.json').read_text())
    for path,digest in frozen['input_sources'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise ValueError('FROZEN_INPUT_CHANGED')
    payload=dict(frozen['requests'][0]);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload)
    settings=spin.SpinEvaluatorConfig(**core_settings(frozen['evaluator_settings']))
    # Diagnose the stalled first attempt, not the subsequent random restart.
    events_path=heart/'rounded-clearance-v3-core-faster-v3/seed-0-job-01/output/run/newton-events.jsonl'
    events=[json.loads(x) for x in events_path.read_text().splitlines()]
    rows=[x for x in events if x.get('restart')==0 and x.get('outcome')=='accepted' and (x.get('step') or 0)>0]
    state=rows[-1];q=np.asarray(state['quadric_coefficients'])
    if replay_record is not None:
        prior=json.loads(Path(replay_record).read_text())
        if prior.get('role')!='SAVED_STATE_TIMING_NOT_OPTIMIZATION':raise ValueError('NOT_A_DIAGNOSTIC_REPLAY_RECORD')
        q=np.asarray(prior['actual_coefficients'],dtype=float)
        if q.shape!=(10,) or not np.isfinite(q).all():raise ValueError('INVALID_REPLAY_COEFFICIENTS')
    output.mkdir(parents=True,exist_ok=False)
    identity=dict(role='SAVED_STATE_TIMING_NOT_OPTIMIZATION',state=state,
        native_path=spin._native_module().__file__,
        native_sha256=hashlib.sha256(Path(spin._native_module().__file__).read_bytes()).hexdigest(),
        event_source_sha256=hashlib.sha256(events_path.read_bytes()).hexdigest())
    for repetition in range(repetitions):
        timer=PhaseTiming();started=time.perf_counter();record=dict(identity,repetition=repetition)
        actual_q=q.copy()
        if perturb_constant and repetition>0:
            actual_q[0]+=0.001*repetition
            actual_q/=np.linalg.norm(actual_q)
        record['actual_coefficients']=actual_q.tolist()
        record['diagnostic_perturbation']='PRE_REGISTERED_CONSTANT_PLUS_0.001' if perturb_constant else 'NONE_SAME_STATE'
        try:
            evaluated=spin.evaluate_spin_it_faster_fullrank_quadric(request,actual_q,
                evaluator_settings=settings,diagnostic_stage_observer=timer)
            record.update(evaluation=asdict(evaluated),status='ACTUAL_MAPPING_RETURNED')
        except Exception as exc:
            record.update(status='ACTUAL_MAPPING_FAILED',exception_type=type(exc).__name__,exception=str(exc))
        record.update(seconds=time.perf_counter()-started,phases=timer.finish(),optimization_admission=False)
        name='heart-state-timing.json' if repetition==0 else f'heart-state-timing-{repetition}.json'
        with (output/name).open('x') as stream:json.dump(json_value(record),stream,allow_nan=False)
        print(json.dumps(dict(status=record['status'],seconds=record['seconds'],phases=record['phases']['stages'])),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--repetitions',type=int,choices=(1,2,3,4),default=1)
    p.add_argument('--perturb-constant',action='store_true')
    p.add_argument('--replay-record',type=Path)
    args=p.parse_args();main(args.output,args.repetitions,args.perturb_constant,args.replay_record)
