"""Three registered mesh scales of one exact planar design entity."""
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from exact_planar_patches import oriented_planar_chain_signature
from independent_tetra_moments import integrate
from same_entity_refinement import refine_mesh

HERE = Path(__file__).resolve().parent


def digest(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def closed_edges(faces):
    edges = defaultdict(list)
    for face in faces:
        for u, v in zip(face, face[1:] + face[:1]):
            edges[tuple(sorted((u, v)))].append((u, v))
    if not all(len(uses) == 2 and uses[0] == uses[1][::-1] for uses in edges.values()):
        raise ValueError('NONMANIFOLD_OR_OPEN_OR_MISORIENTED_INPUT')
    return len(edges)


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    config = json.loads((HERE / 'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json').read_text())
    frozen = json.loads((HERE / 'same-input-job-03/output/frozen-inputs.json').read_text())
    vertices = config['request']['design_vertices']
    faces = config['request']['design_faces']
    original = oriented_planar_chain_signature(vertices, faces)
    original_moments = integrate([[Fraction(float(x)) for x in p] for p in vertices], faces, config['request']['material_density'])
    rows = []
    for name, rounds, longest in (('small', 0, False), ('faces-12288', 5, False), ('faces-98304', 6, True)):
        started = time.monotonic()
        v, f = refine_mesh(vertices, faces, rounds=rounds, split_longest=longest)
        edge_count = closed_edges(f)
        exact = oriented_planar_chain_signature(v, f)
        if exact['exact_null_triangle_indices'] or exact['signature'] != original['signature']:
            raise ValueError('REFINEMENT_CHANGED_EXACT_ORIENTED_SURFACE')
        moments = integrate([[Fraction(float(x)) for x in p] for p in v], f, config['request']['material_density'])
        error = float(np.max(np.abs(moments - original_moments)))
        if error > 1e-10:
            raise ValueError('REFINEMENT_CHANGED_INDEPENDENT_MOMENTS')
        case = json.loads(json.dumps(config))
        case['request']['design_vertices'] = v
        case['request']['design_faces'] = f
        case['request']['settings']['deterministic_seed'] = 0
        case['evaluator']['paper_projective_hyperboloid'] = True
        case['fixed_material'] = frozen['fixed_material']
        case['role'] = 'SAME_ENTITY_CORE_BENCHMARK'
        case['initialization'] = 'REQUEST_AND_SEED_ONLY'
        path = output / (name + '.json')
        payload = json.dumps(case, sort_keys=True, allow_nan=False).encode()
        path.write_bytes(payload)
        row = dict(name=name, configuration_file=path.name, sha256=hashlib.sha256(payload).hexdigest(),
                   design_vertices=len(v), design_faces=len(f), closed_edge_count=edge_count,
                   exact_oriented_surface_equal=True, independent_moment_error=error,
                   original_surface_digest=hashlib.sha256(repr(original['signature']).encode()).hexdigest(),
                   fixed_material_digest=digest(case['fixed_material']),
                   outer_digest=digest([case['request']['outer_vertices'], case['request']['outer_faces']]),
                   seconds=time.monotonic() - started, core_admission=False)
        rows.append(row)
        with (output / 'events.jsonl').open('a') as stream:
            stream.write(json.dumps(row) + '\n')
    (output / 'manifest.json').write_text(json.dumps(dict(
        rows=rows, all_inputs_same_entity=True, refinement_scope='DESIGN_BOUNDARY_ONLY',
        seeds=list(range(10)), maximum_attempts=8, maximum_evaluations_per_attempt=60,
        raw_residual_threshold=1e-12, alpha=0.5, t_max=0.5, core_admission=False)))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
