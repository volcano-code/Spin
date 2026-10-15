"""Describe saved exact facet differences without granting branch admission."""
import hashlib
import json
from pathlib import Path
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def audit(folder):
    from collections import defaultdict
    from make_it_stand.verification.faster_facet_connection import _plane, _points
    def source_groups(material):
        points=_points(material['quadric_source_vertices']);groups=defaultdict(list)
        for index,face in enumerate(material['quadric_source_faces']):groups[_plane(points,face)].append(index)
        return sorted(tuple(v) for v in groups.values())
    report_path = folder/'report.json'
    report = json.loads(report_path.read_text())
    base_path = folder/'witness-000.json'
    base = json.loads(base_path.read_text())
    base_facets = {key(row['facet']): row for row in base['connections']}
    materials = folder.parents[1]/'heart-stalled-production-81-job-01/output'
    base_material_path = materials/'point-000-material.json'
    if digest(base_material_path) != base['material_sha256']:
        raise ValueError('BASE_MATERIAL_CHANGED')
    base_material = json.loads(base_material_path.read_text())
    base_groups=source_groups(base_material)
    rows = []
    for sample in report['classifications']:
        path = folder/f"witness-{sample['point']:03d}.json"
        witness = json.loads(path.read_text())
        if witness['material_sha256'] != sample['material_sha256']:
            raise ValueError('WITNESS_MATERIAL_BINDING_CHANGED')
        facets = {key(row['facet']): row for row in witness['connections']}
        removed = sorted(base_facets.keys()-facets.keys())
        added = sorted(facets.keys()-base_facets.keys())
        changed = sorted(k for k in base_facets.keys() & facets.keys()
                         if base_facets[k] != facets[k])
        same = not (removed or added or changed)
        if same != sample['same_oriented_supporting_facet_connections']:
            raise ValueError('ORDER_OR_DUPLICATE_REQUIRES_SEPARATE_DIAGNOSIS')
        correspondence = None
        if not same:
            material_path = materials/f"point-{sample['point']:03d}-material.json"
            if digest(material_path) != sample['material_sha256']:
                raise ValueError('MATERIAL_CHANGED')
            material = json.loads(material_path.read_text())
            correspondence = dict(
                source_face_indices_identical=material['quadric_source_faces']==base_material['quadric_source_faces'],
                source_vertex_count_identical=len(material['quadric_source_vertices'])==len(base_material['quadric_source_vertices']),
                source_coordinates_identical=material['quadric_source_vertices']==base_material['quadric_source_vertices'])
            correspondence['source_coplanar_face_groups_identical']=source_groups(material)==base_groups
        rows.append(dict(point=sample['point'], witness_sha256=digest(path),
                         material_sha256=sample['material_sha256'],
                         source_index_correspondence=correspondence,
                         classification=('SAME_ORIENTED_SUPPORT_CONNECTIONS' if same
                                         else 'REAL_CONSTRUCTION_SUPPORT_CHANGE' if
                                         correspondence['source_face_indices_identical'] and
                                         correspondence['source_coplanar_face_groups_identical'] else
                                         'SUPPORT_CONNECTION_CHANGE_REQUIRES_SOURCE_CORRESPONDENCE'),
                         removed_facets=[json.loads(k) for k in removed],
                         added_facets=[json.loads(k) for k in added],
                         changed_cycles=[dict(facet=json.loads(k),
                             base=base_facets[k]['oriented_cycles'],
                             sample=facets[k]['oriented_cycles']) for k in changed]))
    return dict(schema='heart-exact-facet-deltas-v1', report_sha256=digest(report_path),
                base_witness_sha256=digest(base_path), samples=rows,
                unchanged=sum(x['classification']=='SAME_ORIENTED_SUPPORT_CONNECTIONS' for x in rows),
                changed=sum(x['classification']!='SAME_ORIENTED_SUPPORT_CONNECTIONS' for x in rows),
                derivative_admission=False, convergence_certificate=False,
                finite_support_classification_complete=all(x['classification'] in
                    ('SAME_ORIENTED_SUPPORT_CONNECTIONS','REAL_CONSTRUCTION_SUPPORT_CHANGE') for x in rows),
                scope='Exact saved support-cycle differences; moving source face correspondence still required')


if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    result = audit(here/'heart-stalled-facet-pilot-job-03/output')
    with Path(sys.argv[1]).open('x') as stream:
        json.dump(result, stream, separators=(',', ':'), allow_nan=False)
    print(json.dumps({k: v for k, v in result.items() if k != 'samples'}))
