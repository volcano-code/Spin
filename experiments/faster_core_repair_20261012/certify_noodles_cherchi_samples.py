"""Offline exact classification of the registered Cherchi noodles samples.

No new source/Boolean calls. Real construction changes remain real changes;
finite classification alone does not certify a smooth parameter neighborhood.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent


def validate_gate(report,manifest,request_sha256):
    if report.get('numerical_passed') is not True or report.get('point_count')!=81:
        raise ValueError('OWN_PRODUCTION_DIFFERENCE_GATE_NOT_PASSED')
    if manifest.get('backend')!='CHERCHI' or manifest.get('point_count')!=81:
        raise ValueError('WRONG_PRODUCTION_BACKEND_OR_SAMPLE_COUNT')
    if manifest.get('request_sha256')!=request_sha256:
        raise ValueError('FROZEN_REQUEST_IDENTITY_MISMATCH')


def main(output):
    import numpy as np
    import trimesh
    from saved_connection_certificate import describe_material_connection
    inputs=HERE/'noodles-formal-cherchi-derivatives-job-01/output'
    request_path=HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
    request_raw=request_path.read_bytes(); request=json.loads(request_raw)
    manifest=json.loads((inputs/'manifest.json').read_text())
    report=json.loads((inputs/'report.json').read_text())
    validate_gate(report,manifest,hashlib.sha256(request_raw).hexdigest())
    samples=manifest['samples']
    if len(samples)!=81 or len({tuple(q) for q in samples})!=81:
        raise ValueError('INVALID_UNIQUE_SAMPLE_SET')
    for name,digest in request['input_sources'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=digest:
            raise ValueError('FROZEN_INPUT_CHANGED')
    paths=[Path(p) for p in request['input_sources'] if Path(p).name=='design-si-candidate.obj']
    if len(paths)!=1:raise ValueError('PHYSICAL_DESIGN_IDENTITY_AMBIGUOUS')
    mesh=trimesh.load(paths[0],force='mesh',process=False)
    if not np.array_equal(mesh.faces,request['requests'][0]['design_faces']):
        raise ValueError('PHYSICAL_DESIGN_FACE_IDENTITY_MISMATCH')
    base=json.loads((inputs/'point-000-material.json').read_text())
    output.mkdir(parents=True,exist_ok=False)
    rows=[]
    with (output/'events.jsonl').open('x') as events:
        for index,q in enumerate(samples):
            path=inputs/f'point-{index:03d}-material.json'
            raw=path.read_bytes(); material=json.loads(raw)
            digest=hashlib.sha256(raw).hexdigest()
            try:
                if material['quadric_coefficients']!=q:raise ValueError('PARAMETER_IDENTITY_MISMATCH')
                if material['evaluation_map_identity']!=base['evaluation_map_identity']:
                    raise ValueError('MIXED_PRODUCTION_MAPS')
                details=describe_material_connection(base,material,mesh.vertices)
                row=dict(index=index,status='FINITE_SAMPLE_CLASSIFIED',**details,material_sha256=digest)
                rows.append(row)
            except Exception as exc:
                events.write(json.dumps(dict(index=index,status='UNCLASSIFIED',
                    error=f'{type(exc).__name__}: {exc}',material_sha256=digest))+'\n')
                events.flush();os.fsync(events.fileno());raise
            events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno())
    (output/'report.json').write_text(json.dumps(dict(point_count=len(rows),rows=rows,
        finite_classification_complete=len(rows)==81,finite_material_connection_certified=False,
        neighborhood_certification=False,case_admission=False,source_calls=0,boolean_calls=0,
        request_sha256=hashlib.sha256(request_raw).hexdigest(),
        input_manifest_sha256=hashlib.sha256((inputs/'manifest.json').read_bytes()).hexdigest())))


if __name__=='__main__':main(Path(sys.argv[1]))
