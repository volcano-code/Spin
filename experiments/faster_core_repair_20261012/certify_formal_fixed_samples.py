"""Classify only saved finite samples using exact construction and oriented chains."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from cherchi_core_map import replay_lpi_constraint
from exact_planar_patches import oriented_planar_chain_signature
from make_it_stand.verification.cherchi_inprocess import construction_constraint_key,exact_fraction

HERE=Path(__file__).resolve().parent


def main(output, inputs=None, design=None, states=('ellipsoid','one_sheet','two_sheet')):
    output.mkdir(parents=True,exist_ok=False)
    if inputs is None:inputs=HERE/'formal-fixed-derivatives-job-01/output'
    if design is None:design=json.loads((HERE/'same-input-job-03/output/frozen-inputs.json').read_text())['request']['design_vertices']
    reports=[]
    for name in states:
        folder=inputs/name
        manifest=json.loads((folder/'manifest.json').read_text())
        base=json.loads((folder/'point-000-material.json').read_text())
        material=base['same_call_material']
        records=sorted(material['material_proof']['resolved_construction_records'],key=lambda r:r['vertex_index'])
        if [r['vertex_index'] for r in records]!=list(range(len(records))):raise ValueError('NONCONTIGUOUS_RECORDS')
        if len(set(construction_constraint_key(r) for r in records))!=len(records):raise ValueError('AMBIGUOUS_FEATURES')
        rows=[]
        for index,q in enumerate(manifest['samples']):
            path=folder/f'point-{index:03d}-material.json'
            saved=json.loads(path.read_text());m=saved['same_call_material']
            row=dict(index=index,material_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            try:
                if m['quadric_coefficients']!=q:raise ValueError('PARAMETER_MISMATCH')
                if saved['source_branch']!=base['source_branch']:raise ValueError('SOURCE_BRANCH_CHANGED')
                operands=[[tuple(Fraction(float(x)) for x in p) for p in vertices]
                    for vertices in (design,m['quadric_source_vertices'])]
                points=[]
                for record in records:
                    if record['point_type']=='EXPLICIT3D':
                        v=record['operands'][0]['source_vertices'][0]
                        point=operands[v['input_object']][v['input_vertex']]
                    else:point=replay_lpi_constraint(record,operands)
                    points.append(point)
                actual=[tuple(exact_fraction(x) for x in p) for p in m['material_proof']['exact_vertices']]
                if len(set(points))!=len(points) or len(points)!=len(actual) or set(points)!=set(actual):
                    raise ValueError('RECIPE_NOT_BIJECTIVE')
                correspondence={p:i for i,p in enumerate(points)}
                current=[[correspondence[actual[v]] for v in face] for face in m['variable_faces']]
                reference=material['variable_faces']
                a=oriented_planar_chain_signature(points,reference)
                b=oriented_planar_chain_signature(points,current)
                if a['signature']!=b['signature']:raise ValueError('ORIENTED_SURFACE_CHAIN_DIFFERS')
                row.update(status='FINITE_SAMPLE_CERTIFIED',classification='IDENTICAL_OR_RENUMBERED'
                    if current==reference else 'EXACT_COPLANAR_RETRIANGULATION')
            except Exception as exc:
                row.update(status='NOT_CERTIFIED',error=f'{type(exc).__name__}: {exc}')
            rows.append(row)
        report=dict(state=name,rows=rows,finite_sample_certification=len(rows)==81 and all(
            r['status']=='FINITE_SAMPLE_CERTIFIED' for r in rows),neighborhood_certification=False,
            boolean_calls=0,source_calls=0,manifest_sha256=hashlib.sha256((folder/'manifest.json').read_bytes()).hexdigest())
        (output/(name+'.json')).write_text(json.dumps(report));reports.append(report)
    (output/'report.json').write_text(json.dumps(dict(states=reports,core_admission=False)))


if __name__=='__main__':main(Path(sys.argv[1]))
