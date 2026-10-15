"""Finite-sample exact oriented-chain classification, not mesh admission."""
from fractions import Fraction
from cherchi_core_map import replay_lpi_constraint
from make_it_stand.verification.cherchi_inprocess import construction_constraint_key,exact_fraction
from exact_planar_patches import oriented_planar_chain_signature


def _replayed_boundary(reference,current,design):
    if (current['boolean_operation']!=reference['boolean_operation'] or
        current['quadric_source_faces']!=reference['quadric_source_faces']):
        raise ValueError('SOURCE_CONNECTIVITY_OR_OPERATION_CHANGED')
    records=sorted(reference['material_proof']['resolved_construction_records'],key=lambda r:r['vertex_index'])
    if [r['vertex_index'] for r in records]!=list(range(len(records))):
        raise ValueError('NONCONTIGUOUS_RECORDS')
    if len(set(construction_constraint_key(r) for r in records))!=len(records):
        raise ValueError('AMBIGUOUS_FEATURES')
    operands=[[tuple(Fraction(float(x)) for x in row) for row in vertices]
              for vertices in (design,current['quadric_source_vertices'])]
    points=[]
    for record in records:
        if record['point_type']=='EXPLICIT3D':
            vertex=record['operands'][0]['source_vertices'][0]
            point=operands[vertex['input_object']][vertex['input_vertex']]
        else:point=replay_lpi_constraint(record,operands)
        points.append(point)
    actual=[tuple(exact_fraction(x) if isinstance(x,str) else Fraction(x) for x in row)
            for row in current['material_proof']['exact_vertices']]
    return points,actual


def classify_material_connection(reference,current,design):
    points,actual=_replayed_boundary(reference,current,design)
    if len(set(points))!=len(points) or len(set(actual))!=len(actual) or set(points)!=set(actual):
        raise ValueError('RECIPE_NOT_BIJECTIVE')
    correspondence={p:i for i,p in enumerate(points)}
    current_faces=[[correspondence[actual[v]] for v in face] for face in current['variable_faces']]
    original_faces=reference['variable_faces']
    if oriented_planar_chain_signature(points,original_faces)['signature']!=oriented_planar_chain_signature(points,current_faces)['signature']:
        raise ValueError('ORIENTED_SURFACE_CHAIN_DIFFERS')
    return 'IDENTICAL_OR_RENUMBERED' if current_faces==original_faces else 'EXACT_COPLANAR_RETRIANGULATION'


def describe_material_connection(reference,current,design):
    """Classify finite geometry without granting smooth-branch admission."""
    own,actual=_replayed_boundary(current,current,design)
    if own!=actual or len(set(actual))!=len(actual):
        raise ValueError('CURRENT_CONSTRUCTION_REPLAY_MISMATCH')
    predicted,actual=_replayed_boundary(reference,current,design)
    same=oriented_planar_chain_signature(predicted,reference['variable_faces'])['signature']==oriented_planar_chain_signature(actual,current['variable_faces'])['signature']
    if same:
        try:classification=classify_material_connection(reference,current,design)
        except ValueError as exc:
            if str(exc)!='RECIPE_NOT_BIJECTIVE':raise
            classification='EXACT_COPLANAR_REFINEMENT'
    else:
        keys=lambda m:{construction_constraint_key(r) for r in m['material_proof']['resolved_construction_records']}
        if keys(reference)==keys(current):
            raise ValueError('SURFACE_CHANGED_WITHOUT_CONSTRUCTION_SWITCH')
        classification='REAL_CONSTRUCTION_BRANCH_CHANGE'
    return dict(classification=classification,own_construction_replayed=True,
        same_reference_surface=same,smooth_branch_certified=False,
        neighborhood_certification=False)
