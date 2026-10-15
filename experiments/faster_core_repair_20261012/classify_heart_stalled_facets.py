"""Pilot exact supporting-facet evidence for two saved production samples."""
import hashlib
import json
from pathlib import Path
import sys
import time

from make_it_stand.verification.cherchi_inprocess import exact_fraction
from make_it_stand.verification.faster_facet_connection import facet_connection_witness

HERE = Path(__file__).resolve().parent


def main(output, *, all_samples=False):
    probe = HERE/'heart-stalled-production-81-job-01/output'
    report = json.loads((probe/'report.json').read_text())
    if report['point_count'] != 81:
        raise ValueError('INCOMPLETE_PRODUCTION_SAMPLE')
    frozen = json.loads((HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json').read_text())
    path = next(Path(name) for name in frozen['input_sources'] if Path(name).name == 'design.binary64.json')
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != frozen['input_sources'][str(path)]:
        raise ValueError('FROZEN_DESIGN_CHANGED')
    design = json.loads(raw)
    output.mkdir(parents=True, exist_ok=False)
    witnesses = []
    for index in (range(81) if all_samples else (0,16)):
        material_path = probe/f'point-{index:03d}-material.json'
        payload = material_path.read_bytes(); material = json.loads(payload)
        started = time.monotonic()
        try:
            witness = facet_connection_witness(
                [[exact_fraction(token) for token in row] for row in material['material_proof']['exact_vertices']],
                material['variable_faces'], [(design['vertices'],design['faces']),
                    (material['quadric_source_vertices'],material['quadric_source_faces'])])
        except ValueError as error:
            with (output/'report.json').open('x') as stream:
                json.dump(dict(schema='heart-stalled-facet-pilot-v1', completed=False,
                    failed_point=index, error=str(error), elapsed_seconds=time.monotonic()-started,
                    derivative_admission=False, case_complete=False), stream)
            raise
        witness.update(point=index, material_sha256=hashlib.sha256(payload).hexdigest(),
                       elapsed_seconds=time.monotonic()-started)
        with (output/f'witness-{index:03d}.json').open('x') as stream:
            json.dump(witness,stream,indent=2,allow_nan=False)
        witnesses.append(witness)
        print('FACET_WITNESS',index,witness['elapsed_seconds'],flush=True)
    classifications = [dict(point=w['point'], material_sha256=w['material_sha256'],
        same_oriented_supporting_facet_connections=w['connections']==witnesses[0]['connections'],
        elapsed_seconds=w['elapsed_seconds']) for w in witnesses]
    result = dict(schema='heart-stalled-facet-classification-v1' if all_samples else 'heart-stalled-facet-pilot-v1',
        completed=True, sample_count=len(witnesses), classifications=classifications,
        same_oriented_supporting_facet_connections=all(
            row['same_oriented_supporting_facet_connections'] for row in classifications),
        independently_audited_production_probe_numerical_passed=report['numerical_passed'],
        derivative_admission=False, neighborhood_certification=False, case_complete=False,
        finite_connection_classification_complete=all_samples and len(witnesses)==81,
        scope='Registered finite samples only; no transfer to unclassified samples or neighborhoods')
    with (output/'report.json').open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main(Path(sys.argv[1]),all_samples='--all-samples' in sys.argv[2:])
