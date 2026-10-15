"""Fixed failure regressions: actual source, actual CGAL material, never optimize."""
import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import numpy as np
from make_it_stand import spin
from phase_timing import PhaseTiming
from run_small_inprocess_solver import safe

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def main(output):
    heart=ROOT/'experiments/faster_heart_contact_20261011'
    sys.path.insert(0,str(heart))
    from core_faster_heart import core_settings
    frozen_path=heart/'rounded-clearance-v3-faster-frozen-v2/request.json'
    frozen=json.loads(frozen_path.read_text())
    for path,digest in frozen['input_sources'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise ValueError('FROZEN_INPUT_CHANGED')
    payload=dict(frozen['requests'][0]);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload)
    settings=core_settings(frozen['evaluator_settings']);settings['paper_projective_hyperboloid']=True
    evaluator=spin.SpinEvaluatorConfig(**settings)
    native=spin._native_module()
    source_settings=native.SpinTriangulationSettings()
    for attr in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
                 'retain_fixed_normal_provenance','truncation_radius','algebraic_finite_map','paper_projective_hyperboloid'):
        setattr(source_settings,attr,settings[attr])
    design=(np.asarray(request.design_vertices)-request.target_origin)@np.asarray(request.target_axes)
    states=json.loads((HERE/'frozen-heart-states.json').read_text())['records']
    output.mkdir(parents=True,exist_ok=False)
    with (output/'configuration.json').open('x') as stream:
        json.dump(dict(role='FIXED_FAILED_STATES_NOT_OPTIMIZATION',settings=settings,
            request_sha256=hashlib.sha256(frozen_path.read_bytes()).hexdigest(),
            native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
            states={name:states[name]['coefficients'] for name in ('heart-first-cap-failure','heart-first-small-face-rejection')}),stream)
    reports=[]
    with (output/'events.jsonl').open('x') as events:
        def event(row):events.write(json.dumps(safe(row),allow_nan=False)+'\n');events.flush();os.fsync(events.fileno())
        for name in ('heart-first-cap-failure','heart-first-small-face-rejection'):
            q=np.asarray(states[name]['coefficients']);start=time.monotonic();report=dict(state=name)
            event(dict(status='FIXED_STATE_STARTED',state=name))
            try:
                source=native.evaluate_spin_it_faster_source_motion(q,design,source_settings)
                audit=native.audit_exact_triangle_boundary(source['vertices'],source['faces'],'OutwardPositive')
                with (output/f'{name}-source.json').open('x') as stream:json.dump(safe(source),stream,allow_nan=False)
                with (output/f'{name}-source-audit.json').open('x') as stream:json.dump(safe(audit),stream,allow_nan=False)
                if not source['defined'] or not source['side_selection_defined'] or not audit['boundary_valid']:
                    raise ValueError('SOURCE_OR_SELECTOR_NOT_CERTIFIED')
                report['source_boundary_valid']=True
                event(dict(status='ACTUAL_CGAL_MAPPING_STARTED',state=name))
                timer=PhaseTiming()
                result=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=evaluator,
                    diagnostic_stage_observer=timer)
                with (output/f'{name}-evaluation.json').open('x') as stream:
                    json.dump(safe(asdict(result)),stream,allow_nan=False)
                report.update(actual_material_returned=True,derivative_defined=result.derivative_defined,
                    material_vertices=len(result.material_vertices),material_faces=len(result.material_faces),
                    phases=timer.finish(),source_derivative_diagnostic=source['diagnostic'])
            except Exception as exc:
                report.update(actual_material_returned=False,exception_type=type(exc).__name__,exception=str(exc))
            report.update(seconds=time.monotonic()-start,optimization_admission=False)
            reports.append(report);event(dict(status='FIXED_STATE_TERMINAL',**report))
    with (output/'report.json').open('x') as stream:json.dump(safe(dict(records=reports)),stream,indent=2,allow_nan=False)
    print(json.dumps(safe(reports)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
