"""Fixed failed production point under existing exact material policy."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import time
from make_it_stand import spin
from phase_timing import PhaseTiming
from run_small_inprocess_solver import safe

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def main(output):
    heart=ROOT/'experiments/faster_heart_contact_20261011';sys.path.insert(0,str(heart))
    from core_faster_heart import core_settings
    request_file=heart/'rounded-clearance-v3-faster-frozen-v2/request.json'
    frozen=json.loads(request_file.read_text())
    for path,digest in frozen['input_sources'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise ValueError('FROZEN_INPUT_CHANGED')
    payload=dict(frozen['requests'][0]);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload)
    manifest_file=HERE/'repaired-heart-production-moments-job-02/output/manifest.json'
    manifest=json.loads(manifest_file.read_text());q=manifest['samples'][4]
    config=core_settings(frozen['evaluator_settings'])
    config.update(paper_projective_hyperboloid=True,exact_binary64_material_boundary=True)
    output.mkdir(parents=True,exist_ok=False)
    with (output/'manifest.json').open('x') as stream:
        json.dump(dict(role='FIXED_FAILED_POINT_EXACT_MATERIAL_POLICY_DIAGNOSIS',point_index=4,
            coefficients=q,settings=config,source_manifest_sha256=hashlib.sha256(manifest_file.read_bytes()).hexdigest()),stream)
    timer=PhaseTiming();start=time.monotonic();report=dict(optimization_admission=False)
    try:
        result=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,
            evaluator_settings=spin.SpinEvaluatorConfig(**config),diagnostic_stage_observer=timer)
        with (output/'evaluation.json').open('x') as stream:json.dump(safe(asdict(result)),stream,allow_nan=False)
        report.update(actual_mapping_returned=True,derivative_defined=result.derivative_defined,
                      diagnostic=result.derivative_diagnostic)
    except Exception as exc:report.update(actual_mapping_returned=False,exception_type=type(exc).__name__,exception=str(exc))
    report.update(seconds=time.monotonic()-start,phases=timer.finish())
    with (output/'report.json').open('x') as stream:json.dump(safe(report),stream,allow_nan=False,indent=2)
    print(json.dumps(safe(report)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
