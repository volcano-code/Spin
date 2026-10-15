"""Classify saved Heart geometry without relabeling its failed FD protocol."""
import hashlib
import json
import os
from pathlib import Path
import sys
from saved_connection_certificate import describe_material_connection

HERE=Path(__file__).resolve().parent


def main(output):
    inputs=HERE/'heart-cherchi-heart-first-small-face-rejection-derivative-job-01/output'
    numeric=json.loads((inputs/'report.json').read_text())
    manifest=json.loads((inputs/'manifest.json').read_text())
    reference_dir=HERE/'heart-projective-independent-reference-job-01/output'
    reference=json.loads((reference_dir/'report.json').read_text())
    protocol=json.loads((reference_dir/'protocol.json').read_text())
    base_path=inputs/'point-000-material.json';base_raw=base_path.read_bytes();base=json.loads(base_raw)
    if numeric['point_count']!=81 or len(manifest['samples'])!=81:
        raise ValueError('INCOMPLETE_REGISTERED_SAMPLES')
    if (reference['status']!='SUPPORTED' or reference['wrong_jacobian_status']!='INCONSISTENT' or
        reference['direction']!=0 or protocol['material_sha256']!=hashlib.sha256(base_raw).hexdigest()):
        raise ValueError('POINTWISE_REFERENCE_NOT_BOUND')
    registration=json.loads((HERE/'heart-cherchi-coordinate-map-job-01/output/registration.json').read_text())
    design=registration['original']['design_vertices']
    output.mkdir(parents=True,exist_ok=False);rows=[]
    with (output/'events.jsonl').open('x') as events:
        for index,q in enumerate(manifest['samples']):
            raw=(inputs/f'point-{index:03d}-material.json').read_bytes();material=json.loads(raw)
            try:
                if material['quadric_coefficients']!=q:raise ValueError('SAMPLE_COEFFICIENT_CHANGED')
                if material['evaluation_map_identity']!=base['evaluation_map_identity']:
                    raise ValueError('MIXED_PRODUCTION_MAPS')
                row=dict(index=index,status='FINITE_SAMPLE_CLASSIFIED',
                    **describe_material_connection(base,material,design),
                    material_sha256=hashlib.sha256(raw).hexdigest())
                rows.append(row)
            except Exception as exc:
                events.write(json.dumps(dict(index=index,status='UNCLASSIFIED',error=f'{type(exc).__name__}: {exc}'))+'\n')
                events.flush();os.fsync(events.fileno());raise
            events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno())
    report=dict(point_count=len(rows),rows=rows,finite_classification_complete=len(rows)==81,
        original_five_point_numerical_passed=numeric['numerical_passed'],
        pointwise_reference_status=reference['status'],reference_direction=0,
        neighborhood_certification=False,case_admission=False,source_calls=0,boolean_calls=0,
        scope='GEOMETRY_CLASSIFICATION_WITH_FAILED_FD_RETAINED',
        manifest_sha256=hashlib.sha256((inputs/'manifest.json').read_bytes()).hexdigest(),
        reference_sha256=hashlib.sha256((reference_dir/'report.json').read_bytes()).hexdigest())
    (output/'report.json').write_text(json.dumps(report))


if __name__=='__main__':main(Path(sys.argv[1]))
