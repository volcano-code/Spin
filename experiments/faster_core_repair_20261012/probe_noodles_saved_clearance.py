"""Measure the frozen input boundary gap; never choose/change solver wall settings."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import trimesh
from make_it_stand import spin

HERE=Path(__file__).resolve().parent


def main(output):
    request_path=HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
    request=json.loads(request_path.read_text())
    for name,digest in request['input_sources'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=digest:
            raise ValueError('FROZEN_INPUT_CHANGED')
    paths={Path(p).name:Path(p) for p in request['input_sources']}
    output.mkdir(parents=True,exist_ok=False)
    native=spin._native_module();regions={};timings={}
    for role in ('outer','design'):
        started=time.monotonic()
        path=paths[role+'-si-candidate.obj']
        mesh=trimesh.load(path,force='mesh',process=False)
        region=native.make_exact_material_region(mesh.vertices,mesh.faces)
        (output/(role+'-facts.json')).write_text(json.dumps(region.facts))
        if not region.facts['valid'] or not region.facts['completed']:
            raise ValueError('INVALID_EXACT_INPUT_REGION:'+role)
        regions[role]=region;timings[role+'_region_seconds']=time.monotonic()-started
        print('REGION_READY',role,timings[role+'_region_seconds'],flush=True)
    started=time.monotonic()
    result=native.check_exact_boundary_clearance(regions['outer'],regions['design'],0.)
    timings['distance_seconds']=time.monotonic()-started
    (output/'distance.json').write_text(json.dumps(result))
    if not result['completed'] or not result['applicable']:
        raise ValueError('EXACT_INPUT_DISTANCE_UNAVAILABLE')
    squared=Fraction(result['squared_distance_exact'])
    report=dict(measurement_completed=True,strictly_positive_gap=squared>0,
        squared_distance_exact=str(squared),witness_faces=result['witness_faces'],
        request_wall_values=[r['minimum_wall_thickness'] for r in request['requests']],
        request_sha256=hashlib.sha256(request_path.read_bytes()).hexdigest(),
        native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        timings=timings,case_admission=False,printing_qualification=False,
        scope='FROZEN_INPUT_BOUNDARY_DISTANCE_ONLY_REQUEST_UNCHANGED')
    (output/'report.json').write_text(json.dumps(report));print(json.dumps(report),flush=True)


if __name__=='__main__':main(Path(sys.argv[1]))
