"""Exact supporting-facet witnesses distinguish storage from real changes."""
from fractions import Fraction
import importlib

import pytest


VERTICES = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1)]
FACES = [(0,2,1),(0,3,2),(4,5,6),(4,6,7),
         (0,1,5),(0,5,4),(1,2,6),(1,6,5),
         (2,3,7),(2,7,6),(3,0,4),(3,4,7)]


def witness(vertices, faces, inputs=None):
    api = importlib.import_module('make_it_stand.verification.faster_facet_connection')
    return api.facet_connection_witness(vertices, faces, inputs or [(VERTICES, FACES)])


def test_exact_coplanar_diagonal_change_has_same_oriented_facet_connections():
    altered = list(FACES)
    altered[:2] = [(0,3,1),(1,3,2)]
    a, b = witness(VERTICES,FACES), witness(VERTICES,altered)
    assert a['connections'] == b['connections']
    assert a['material_boundary_sha256'] != b['material_boundary_sha256']
    assert not a['neighborhood_certification']


def test_exact_collinear_edge_subdivision_preserves_facet_connections():
    vertices = VERTICES+[(Fraction(1,2),0,0)]
    faces = list(FACES)
    faces[0:1] = [(0,2,8),(8,2,1)]
    index = faces.index((0,1,5))
    faces[index:index+1] = [(0,8,5),(8,1,5)]
    assert witness(vertices,faces)['connections'] == witness(VERTICES,FACES)['connections']


def test_material_side_reversal_changes_oriented_connections():
    reversed_faces = [(a,c,b) for a,b,c in FACES]
    assert witness(VERTICES,reversed_faces)['connections'] != witness(VERTICES,FACES)['connections']


def test_support_identity_survives_continuous_input_geometry_change():
    vertices = [(2*x,y,z) for x,y,z in VERTICES]
    a, b = witness(VERTICES,FACES), witness(vertices,FACES,[(vertices,FACES)])
    assert a['connections'] == b['connections']
    assert a['input_boundary_sha256'] != b['input_boundary_sha256']


def test_an_unexplained_material_plane_is_rejected():
    vertices = [(Fraction(x,2),y,z) for x,y,z in VERTICES]
    with pytest.raises(ValueError,match='UNEXPLAINED_MATERIAL_PLANE'):
        witness(vertices,FACES)


def test_open_material_boundary_is_rejected():
    with pytest.raises(ValueError,match='OPEN_OR_INCONSISTENT_BOUNDARY'):
        witness(VERTICES,FACES[:-1])


def test_coplanar_patches_at_a_saddle_vertex_follow_actual_triangle_fans():
    ring = [(1,0,0),(1,1,0),(0,1,1),(-1,0,0),(-1,-1,0),(0,-1,-1)]
    vertices = [(0,0,0)]+ring+[(x,y,-2) for x,y,z in ring]
    faces = [(0,i+1,(i+1)%6+1) for i in range(6)]
    faces += [(7,9,8),(7,10,9),(7,11,10),(7,12,11)]
    for i in range(6):
        a,b,c,d = i+1,(i+1)%6+1,i+7,(i+1)%6+7
        faces += [(a,c,d),(a,d,b)]
    result = witness(vertices,faces,[(vertices,faces)])
    # Top triangles 0 and 3 are distinct z=0 patches meeting at the center.
    facet = next(row for row in result['connections'] if row['facet']==((0,(0,3)),))
    assert len(facet['oriented_cycles'])==2
