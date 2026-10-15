"""Replay the registered source only; no Boolean, initialization or Newton."""
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import _make_it_stand as native

HERE = Path(__file__).resolve().parent


def main():
    registration = json.loads((HERE/'heart-seed1-source-failure-registration-01.json').read_text())
    run = HERE/'heart-certified-formal-seed-1-job-01/output/run'
    for name, filename in [('events', 'events.jsonl'), ('configuration', 'configuration.json')]:
        if hashlib.sha256((run/filename).read_bytes()).hexdigest() != registration[name+'_sha256']:
            raise ValueError('REGISTERED_SOURCE_RUN_CHANGED')
    config = json.loads((run/'configuration.json').read_text())
    request = config['request']
    design = (np.asarray(request['design_vertices'])-np.asarray(request['target_origin']))@np.asarray(request['target_axes']).T
    settings = native.SpinTriangulationSettings()
    for name in ('icosphere_subdivisions', 'maximum_icosphere_subdivisions',
                 'maximum_radius_doublings', 'retain_fixed_normal_provenance',
                 'truncation_radius', 'paper_projective_hyperboloid', 'algebraic_finite_map'):
        setattr(settings, name, config['evaluator'][name])
    start = time.monotonic()
    try:
        source = native.evaluate_spin_it_faster_source_motion(
            np.asarray(registration['event']['quadric_coefficients']), design, settings)
    except Exception as error:
        print(json.dumps(dict(status='SOURCE_EXCEPTION', error=str(error),
            source_seconds=time.monotonic()-start, optimization_invocations=0,
            native_module=str(native.__file__))), flush=True)
        return
    print(json.dumps(dict(status='SOURCE_RETURNED', defined=source['defined'],
        diagnostic=source['diagnostic'], source_seconds=time.monotonic()-start,
        optimization_invocations=0, native_module=str(native.__file__))), flush=True)


if __name__ == '__main__':
    main()
