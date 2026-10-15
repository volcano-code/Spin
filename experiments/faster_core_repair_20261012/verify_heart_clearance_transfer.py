"""Reuse bound exact distance and bidirectional perturbation certificates."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import sys
from clearance_bound_transfer import transfer

HERE = Path(__file__).resolve().parent
HEART = HERE.parent / 'faster_heart_contact_20261011'


def main(output):
    load = lambda p: json.loads(p.read_text())
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    frozen_path = HEART / 'rounded-clearance-v3-faster-frozen-v2/request.json'
    frozen = load(frozen_path)
    for filename, expected in frozen['input_sources'].items():
        if sha(Path(filename)) != expected:
            raise ValueError('FROZEN_INPUT_CHANGED')
    derived_path = HEART / 'rounded-clearance-v3-bounded-admission-job-03/output/admission.json'
    derived = load(derived_path)
    for filename, expected in derived['source_files'].items():
        if sha(Path(filename)) != expected:
            raise ValueError('DERIVED_EVIDENCE_SOURCE_CHANGED')
    authority_path = HEART / 'rounded-clearance-v3-admission-job-01/output/admission.json'
    authority = load(authority_path)
    if sha(authority_path) != derived['authority_admission_sha256']:
        raise ValueError('AUTHORITY_NOT_BOUND')
    distance_path = authority_path.parent / 'clearance.json'
    if sha(distance_path) != authority['files']['clearance.json']:
        raise ValueError('DISTANCE_NOT_BOUND')
    distance = load(distance_path)
    if not distance['completed'] or not distance['applicable']:
        raise ValueError('DISTANCE_UNAVAILABLE')
    bounds = derived['bounds']
    if not derived['geometry_admission_granted'] or not all(
        bounds[role]['entity_bound_applied'] and bounds[role]['numerical_bounds_within_tolerance']
        for role in ('outer', 'design')):
        raise ValueError('BIDIRECTIONAL_BOUNDS_NOT_ADMITTED')
    result = transfer(distance['squared_distance_exact'],
        bounds['outer']['bidirectional_boundary_hausdorff_linf_upper_bound_exact'],
        bounds['design']['bidirectional_boundary_hausdorff_linf_upper_bound_exact'],
        Fraction(375, 1000000) / Fraction(25, 1000))
    report = dict(result, request_sha256=sha(frozen_path), derived_admission_sha256=sha(derived_path),
        authority_admission_sha256=sha(authority_path), distance_sha256=sha(distance_path),
        wall_m_exact='3/8000', model_unit_m_exact='1/40',
        full_manufacturing_wall_thickness='UNKNOWN', optimization_invocations=0,
        case_admission=False, scope='FROZEN_DERIVED_BOUNDARY_CLEARANCE_ONLY')
    with output.open('x') as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps(report), flush=True)
    if not result['at_least_required']:
        raise ValueError('BOUND_INSUFFICIENT')


if __name__ == '__main__':
    main(Path(sys.argv[1]))
