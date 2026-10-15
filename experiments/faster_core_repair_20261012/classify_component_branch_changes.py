"""Diagnose saved failed correspondences; never grants a smooth derivative gate."""
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

from cherchi_core_map import replay_lpi_constraint
from exact_planar_patches import oriented_planar_chain_signature
from make_it_stand.verification.cherchi_inprocess import construction_constraint_key, exact_fraction

HERE = Path(__file__).resolve().parent


def replay(records, operands):
    records = sorted(records, key=lambda r: r['vertex_index'])
    if [r['vertex_index'] for r in records] != list(range(len(records))):
        raise ValueError('NONCONTIGUOUS_RECORDS')
    keys = [construction_constraint_key(r) for r in records]
    if len(set(keys)) != len(keys):
        raise ValueError('AMBIGUOUS_FEATURES')
    result = []
    for record in records:
        if record['point_type'] == 'EXPLICIT3D':
            vertex = record['operands'][0]['source_vertices'][0]
            point = operands[vertex['input_object']][vertex['input_vertex']]
        else:
            point = replay_lpi_constraint(record, operands)
        result.append(point)
    return result


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    registration = json.loads((HERE / 'noodles-first-core-failures-registration-01/registration.json').read_text())
    raw = Path(registration['mapping_config']).read_bytes()
    if hashlib.sha256(raw).hexdigest() != registration['mapping_config_sha256']:
        raise ValueError('CHANGED_CONFIG')
    design = json.loads(raw)['request']['design_vertices']
    folder = HERE / 'component-selector-derivatives-job-01/output/source_side'
    prior_path = HERE / 'component-selector-branch-job-01/output/source_side.json'
    prior = json.loads(prior_path.read_text())
    manifest = json.loads((folder / 'manifest.json').read_text())
    base_saved = json.loads((folder / 'point-000-material.json').read_text())
    base = base_saved['same_call_material']
    rows = []
    for old in prior['rows']:
        index = old['index']
        path = folder / f'point-{index:03d}-material.json'
        if hashlib.sha256(path.read_bytes()).hexdigest() != old['material_sha256']:
            raise ValueError('SAVED_MATERIAL_CHANGED')
        row = dict(old)
        if old['status'] == 'FINITE_SAMPLE_CERTIFIED':
            row['classification_resolved'] = True
            rows.append(row)
            continue
        try:
            saved = json.loads(path.read_text())
            current = saved['same_call_material']
            if current['quadric_coefficients'] != manifest['samples'][index]:
                raise ValueError('PARAMETER_MISMATCH')
            if saved['source_branch'] != base_saved['source_branch']:
                raise ValueError('SOURCE_BRANCH_CHANGED')
            operands = [[tuple(Fraction(float(x)) for x in p) for p in vertices]
                        for vertices in (design, current['quadric_source_vertices'])]
            reference = replay(base['material_proof']['resolved_construction_records'], operands)
            actual = [tuple(exact_fraction(x) for x in p)
                      for p in current['material_proof']['exact_vertices']]
            if replay(current['material_proof']['resolved_construction_records'], operands) != actual:
                raise ValueError('CURRENT_PROVENANCE_REPLAY_DIFFERS')
            a = oriented_planar_chain_signature(reference, base['variable_faces'])
            b = oriented_planar_chain_signature(actual, current['variable_faces'])
            if a['exact_null_triangle_indices'] or b['exact_null_triangle_indices']:
                raise ValueError('NULL_REFERENCE_OR_CURRENT_TRIANGLE')
            equal = a['signature'] == b['signature']
            row.update(classification_resolved=True, exact_surface_equal=equal,
                       classification='EXACT_COPLANAR_SUBDIVISION' if equal else 'ACTUAL_CONSTRUCTION_CHANGE',
                       current_provenance_replayed=True,
                       reference_vertices=len(reference), current_vertices=len(actual),
                       reference_plane_count=len(a['signature']), current_plane_count=len(b['signature']))
            if not equal:
                planes_a = {plane: lines for plane, lines in a['signature']}
                planes_b = {plane: lines for plane, lines in b['signature']}
                differing = [plane for plane in sorted(set(planes_a) | set(planes_b))
                             if planes_a.get(plane) != planes_b.get(plane)]
                row['differing_exact_planes'] = [[str(x) for x in plane] for plane in differing]
            row.pop('error', None)
        except Exception as exc:
            row.update(classification_resolved=False, diagnostic=f'{type(exc).__name__}: {exc}')
        rows.append(row)
        with (output / 'events.jsonl').open('a') as stream:
            stream.write(json.dumps(dict(index=index, resolved=row['classification_resolved'],
                                         classification=row.get('classification'))) + '\n')
    (output / 'report.json').write_text(json.dumps(dict(
        rows=rows, all_finite_samples_classified=len(rows) == 81 and all(r['classification_resolved'] for r in rows),
        smooth_derivative_gate=False, neighborhood_certification=False, core_admission=False,
        boolean_calls=0, source_calls=0, prior_certificate_sha256=hashlib.sha256(prior_path.read_bytes()).hexdigest())))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
