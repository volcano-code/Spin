"""Capture missing source recipes, binding replay to the saved actual material.

No Boolean, finite-difference rerun, or optimization. This is not a same-call
historical tape and cannot repair the failed five-point acceptance by itself.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import _json_value

HERE=Path(__file__).resolve().parent


def main(output):
    registration=HERE/'heart-cherchi-coordinate-map-job-01/output/registration.json'
    saved_path=HERE/'heart-cherchi-heart-first-small-face-rejection-derivative-job-01/output/point-000-material.json'
    config=json.loads(registration.read_text());saved=json.loads(saved_path.read_text())
    request=config['normalized'];settings=config['settings']
    target=(np.asarray(request['design_vertices'])-np.asarray(request['target_origin'])) @ np.asarray(request['target_axes']).T
    native=spin._native_module();tri=native.SpinTriangulationSettings()
    for name in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
                 'retain_fixed_normal_provenance','truncation_radius','paper_projective_hyperboloid','algebraic_finite_map'):
        setattr(tri,name,settings[name])
    source=native.evaluate_spin_it_faster_source_motion(np.asarray(saved['quadric_coefficients']),target,tri)
    if not source['defined'] or not source['side_selection_defined']:
        raise ValueError('SOURCE_REPLAY_UNDEFINED')
    physical=(np.asarray(source['vertices']) @ np.asarray(request['target_axes']).T+
              np.asarray(request['target_origin']))*config['length']+np.asarray(config['origin'])
    if not np.array_equal(physical,saved['quadric_source_vertices']) or not np.array_equal(source['faces'],saved['quadric_source_faces']):
        raise ValueError('SOURCE_REPLAY_NOT_IDENTICAL_TO_SAVED_MATERIAL')
    output.mkdir(parents=True,exist_ok=False)
    (output/'source.json').write_text(json.dumps(source,default=_json_value,allow_nan=False))
    (output/'report.json').write_text(json.dumps(dict(
        passed=True,source_calls=1,boolean_calls=0,optimization_calls=0,
        material_sha256=hashlib.sha256(saved_path.read_bytes()).hexdigest(),
        registration_sha256=hashlib.sha256(registration.read_bytes()).hexdigest(),
        native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        scope='SOURCE_REPLAY_IDENTICAL_TO_SAVED_GEOMETRY_NOT_HISTORICAL_SAME_CALL_TAPE',
        derivative_admission=False,case_admission=False)))


if __name__=='__main__':main(Path(sys.argv[1]))
