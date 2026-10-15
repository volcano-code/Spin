from copy import deepcopy
import pytest
from saved_connection_certificate import classify_material_connection


def fixture():
    vertices=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
    faces=[[1,2,3],[0,3,2],[0,1,3],[0,2,1]]
    records=[dict(vertex_index=i,point_type='EXPLICIT3D',operands=[dict(status='UNIQUE',
        source_vertices=[dict(input_object=0,input_vertex=i)])]) for i in range(4)]
    material=dict(quadric_source_faces=[],quadric_source_vertices=[],boolean_operation='intersection',
        variable_faces=faces,material_proof=dict(exact_vertices=vertices,resolved_construction_records=records))
    return vertices,material


def test_exact_reindexing_is_classified_without_coordinate_nearest_matching():
    design,base=fixture();current=deepcopy(base)
    current['material_proof']['exact_vertices']=list(reversed(design))
    current['variable_faces']=[[3-i for i in face] for face in base['variable_faces']]
    assert classify_material_connection(base,current,design)=='IDENTICAL_OR_RENUMBERED'


def test_orientation_or_unmatched_vertex_is_not_certified():
    design,base=fixture();current=deepcopy(base)
    current['variable_faces'][0].reverse()
    with pytest.raises(ValueError):classify_material_connection(base,current,design)


def test_exact_coplanar_diagonal_change_preserves_oriented_surface():
    design,base=fixture()
    design[:]=[[0,0,0],[1,0,0],[1,1,0],[0,1,0],[.5,.5,1]]
    base['material_proof']['exact_vertices']=design
    base['material_proof']['resolved_construction_records'].append(dict(vertex_index=4,
        point_type='EXPLICIT3D',operands=[dict(status='UNIQUE',source_vertices=[dict(input_object=0,input_vertex=4)])]))
    base['variable_faces']=[[0,2,1],[0,3,2],[0,1,4],[1,2,4],[2,3,4],[3,0,4]]
    current=deepcopy(base)
    current['variable_faces'][:2]=[[0,3,1],[1,3,2]]
    assert classify_material_connection(base,current,design)=='EXACT_COPLANAR_RETRIANGULATION'


def test_real_construction_switch_is_recorded_not_promoted_to_smooth_branch():
    from saved_connection_certificate import describe_material_connection
    design,base=fixture();design=design+[[0,0,2]]
    current=deepcopy(base)
    current['material_proof']['exact_vertices'][3]=[0,0,2]
    current['material_proof']['resolved_construction_records'][3]['operands'][0]['source_vertices'][0]['input_vertex']=4
    result=describe_material_connection(base,current,design)
    assert result['classification']=='REAL_CONSTRUCTION_BRANCH_CHANGE'
    assert result['own_construction_replayed']
    assert not result['same_reference_surface']
    assert not result['smooth_branch_certified']
    current['material_proof']['exact_vertices'][3]=[0,0,3]
    with pytest.raises(ValueError):describe_material_connection(base,current,design)
    current=deepcopy(base);current['material_proof']['exact_vertices'][0]=[1,1,1]
    with pytest.raises(ValueError):classify_material_connection(base,current,design)
