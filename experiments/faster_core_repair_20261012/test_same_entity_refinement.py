"""The performance fixture must preserve exact planes and shared edges."""
import importlib.util
import json
from collections import defaultdict
import math
from pathlib import Path

import pytest
from exact_planar_patches import oriented_planar_chain_signature

HERE = Path(__file__).resolve().parent


def api():
    spec = importlib.util.find_spec('same_entity_refinement')
    assert spec is not None, 'same-entity refinement helper missing'
    import same_entity_refinement
    return same_entity_refinement


def cube():
    request = json.loads((HERE / 'same-input-job-03/output/frozen-inputs.json').read_text())['request']
    return request['design_vertices'], request['design_faces']


@pytest.mark.parametrize('longest,count', [(False, 48), (True, 96)])
def test_refinement_preserves_exact_surface_and_closed_shared_edges(longest, count):
    vertices, faces = cube()
    refined, triangles = api().refine_mesh(vertices, faces, rounds=1, split_longest=longest)
    assert len(triangles) == count
    assert len(set(map(tuple, refined))) == len(refined)
    edges = defaultdict(list)
    for face in triangles:
        for u, v in zip(face, face[1:] + face[:1]):
            edges[tuple(sorted((u, v)))].append((u, v))
    assert all(len(uses) == 2 and uses[0] == uses[1][::-1] for uses in edges.values())
    assert oriented_planar_chain_signature(vertices, faces)['signature'] == oriented_planar_chain_signature(refined, triangles)['signature']


def test_nonrepresentable_midpoint_is_rejected():
    with pytest.raises(ValueError, match='MIDPOINT_NOT_EXACT_BINARY64'):
        api().refine_mesh([[1., 0., 0.], [math.nextafter(1., math.inf), 1., 0.], [0., 0., 0.]], [[0, 1, 2]], rounds=1)


def test_no_refinement_does_not_alias_input():
    vertices, faces = cube()
    refined, triangles = api().refine_mesh(vertices, faces, rounds=0)
    refined[0][0] = 99.
    triangles[0][0] = 7
    assert vertices[0][0] == -1.
    assert faces[0][0] == 1
