"""Pre-registered legacy/paper source-only replay; no Boolean or optimization."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from make_it_stand import spin
from run_small_inprocess_solver import safe


def main(output):
    here=Path(__file__).parent;root=here.parents[1]
    frozen_path=root/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
    frozen=json.loads(frozen_path.read_text())
    for path,expected in frozen['input_sources'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected:raise ValueError('FROZEN_INPUT_CHANGED')
    request=frozen['requests'][0];rotation=np.asarray(request['target_axes']).T
    design=(np.asarray(request['design_vertices'])-np.asarray(request['target_origin']))@rotation.T
    states=json.loads((here/'frozen-heart-states.json').read_text())['records']
    native=spin._native_module();output.mkdir(parents=True,exist_ok=False)
    manifest=dict(role='FIXED_FAILED_SOURCE_ONLY_PROFILE_COMPARISON_NOT_OPTIMIZATION',
        states={k:states[k]['coefficients'] for k in ('heart-first-cap-failure','heart-first-small-face-rejection')},
        profiles=[dict(paper_projective_hyperboloid=False),dict(paper_projective_hyperboloid=True)],
        base_settings=frozen['evaluator_settings'],native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        request_sha256=hashlib.sha256(frozen_path.read_bytes()).hexdigest())
    with (output/'manifest.json').open('x') as stream:json.dump(manifest,stream)
    records=[]
    for name,q in manifest['states'].items():
        for paper in (False,True):
            settings=native.SpinTriangulationSettings()
            for attr in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
                         'retain_fixed_normal_provenance','truncation_radius','algebraic_finite_map'):
                setattr(settings,attr,frozen['evaluator_settings'][attr])
            settings.paper_projective_hyperboloid=paper
            row=dict(state=name,paper_projective_hyperboloid=paper,optimization_admission=False)
            start=time.perf_counter()
            try:
                source=native.evaluate_spin_it_faster_source_motion(np.asarray(q),design,settings)
                row.update(defined=bool(source['defined']),diagnostic=source['diagnostic'],
                           side_selection_defined=bool(source['side_selection_defined']))
                with (output/f'{name}-paper-{int(paper)}-source.json').open('x') as stream:
                    json.dump(safe(source),stream,allow_nan=False)
            except Exception as exc:
                row.update(defined=False,exception_type=type(exc).__name__,exception=str(exc))
            row['seconds']=time.perf_counter()-start;records.append(row)
    with (output/'report.json').open('x') as stream:json.dump(dict(records=records,source_only=True),stream,indent=2)
    print(json.dumps(records),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
