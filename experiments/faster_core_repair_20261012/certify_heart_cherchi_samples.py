"""Exact offline replay of finite material samples; no new geometry evaluations."""
import hashlib
import json
import os
from pathlib import Path
import sys
from saved_connection_certificate import classify_material_connection,describe_material_connection
HERE=Path(__file__).resolve().parent

def main(output,state,mode='strict'):
    if mode not in ('strict','classify'):raise ValueError('UNREGISTERED_CLASSIFICATION_MODE')
    if state not in ('heart-first-cap-failure','heart-first-small-face-rejection'):
        raise ValueError('UNREGISTERED_HEART_STATE')
    inputs=HERE/f'heart-cherchi-{state}-derivative-job-01/output'
    report=json.loads((inputs/'report.json').read_text())
    if not report['numerical_passed'] or report['point_count']!=81:
        raise ValueError('PRODUCTION_DIFFERENCE_GATE_NOT_PASSED')
    manifest=json.loads((inputs/'manifest.json').read_text())
    base=json.loads((inputs/'point-000-material.json').read_text())
    registration=json.loads((HERE/'heart-cherchi-coordinate-map-job-01/output/registration.json').read_text())
    design=registration['original']['design_vertices']
    output.mkdir(parents=True,exist_ok=False)
    rows=[]
    with (output/'events.jsonl').open('x') as events:
        for index,q in enumerate(manifest['samples']):
            path=inputs/f'point-{index:03d}-material.json'
            raw=path.read_bytes(); material=json.loads(raw)
            try:
                if material['quadric_coefficients']!=q:raise ValueError('PARAMETER_IDENTITY_MISMATCH')
                if material['evaluation_map_identity']!=base['evaluation_map_identity']:
                    raise ValueError('MIXED_PRODUCTION_MAPS')
                if mode=='strict':
                    details=dict(classification=classify_material_connection(base,material,design))
                else:details=describe_material_connection(base,material,design)
                row=dict(index=index,status='FINITE_SAMPLE_CERTIFIED' if mode=='strict' else 'FINITE_SAMPLE_CLASSIFIED',
                    **details,material_sha256=hashlib.sha256(raw).hexdigest())
                rows.append(row)
            except Exception as exc:
                events.write(json.dumps(dict(index=index,status='UNCLASSIFIED',error=f'{type(exc).__name__}: {exc}',
                    material_sha256=hashlib.sha256(raw).hexdigest()))+'\n')
                events.flush();os.fsync(events.fileno());raise
            events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno())
    (output/'report.json').write_text(json.dumps(dict(state=state,point_count=len(rows),
        finite_material_connection_certified=mode=='strict' and len(rows)==81,
        finite_classification_complete=len(rows)==81,mode=mode,neighborhood_certification=False,
        source_calls=0,boolean_calls=0,case_admission=False,rows=rows)))

if __name__=='__main__':main(Path(sys.argv[1]),sys.argv[2],sys.argv[3] if len(sys.argv)>3 else 'strict')
