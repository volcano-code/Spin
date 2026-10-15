"""Diagnostic timings of one real saved evaluation; no normal-solver timing claim."""
import cProfile
from dataclasses import replace
import hashlib
import io
import json
from pathlib import Path
import pstats
import sys
import time
import numpy as np
from make_it_stand import spin
from make_it_stand.verification import cherchi_inprocess as adapter

HERE=Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    core=json.loads((HERE/'core-benchmark-inputs-job-01/output/faces-12288.json').read_text())
    run=HERE/'layered-formal-benchmark-faces-12288-seed-0-job-04/output/run'
    config=json.loads((run/'configuration.json').read_text())
    payload=dict(config['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload)
    outer=np.asarray(core['request']['outer_vertices']);low=outer.min(0);span=outer.max(0)-low
    normalization=dict(origin=(low+span*.5).tolist(),L=float(max(span)))
    start=time.monotonic()
    producer=adapter.InProcessCherchiEvaluator(request,spin.SpinEvaluatorConfig(**config['evaluator']),
        normalization=normalization,design_vertices_physical=core['request']['design_vertices'],
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1',fixed_boundary_reuse=True)
    preparation=time.monotonic()-start
    times={};calls={}
    def timed(name, operation):
        def wrapper(*args,**kwargs):
            start=time.monotonic()
            try:return operation(*args,**kwargs)
            finally:
                times[name]=times.get(name,0.)+time.monotonic()-start
                calls[name]=calls.get(name,0)+1
        return wrapper
    adapter.normalize_closed_component_orientation=timed('orientation_normalization',adapter.normalize_closed_component_orientation)
    adapter.recover_construction_vertex_motion=timed('construction_vertex_motion',adapter.recover_construction_vertex_motion)
    producer._target_view=timed('exact_target_coordinate_view',producer._target_view)
    material=json.loads((run/'evaluation-000073.json').read_text())
    profiler=cProfile.Profile();start=time.monotonic()
    profiler.enable();actual=producer(material['quadric_coefficients']);profiler.disable()
    elapsed=time.monotonic()-start
    stream=io.StringIO();pstats.Stats(profiler,stream=stream).strip_dirs().sort_stats('cumulative').print_stats(35)
    (output/'profile.txt').write_text(stream.getvalue())
    checked={field:np.array_equal(np.asarray(actual[field]),np.asarray(material[field]))
             for field in ('raw_moments','variable_raw_moments','moment_jacobian')}
    checked.update(geometry_hash=actual['geometry_hash']==material['geometry_hash'],
                   topology_hash=actual['topology_hash']==material['topology_hash'])
    report=dict(passed=all(checked.values()),checks=checked,constructor_seconds=preparation,
        profiled_evaluation_seconds=elapsed,phase_seconds=times,phase_calls=calls,
        producer_stage_seconds=producer.stage_seconds,producer_calls=producer.call_counts,
        native_sha256=hashlib.sha256(Path(producer.native.__file__).read_bytes()).hexdigest(),
        core_admission=False,scope='One saved actual terminal state; cProfile adds overhead, not normal timing')
    (output/'report.json').write_text(json.dumps(report,sort_keys=True,allow_nan=False))
    print(json.dumps(report),flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main(Path(sys.argv[1]))
