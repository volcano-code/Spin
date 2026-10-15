"""Offline two-state exact recipe/chain pilot. No Boolean or new source calls."""
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from birth_face_recipes import recover_vertex_recipe, replay_vertex_recipe
from exact_planar_patches import oriented_planar_chain_signature

HERE = Path(__file__).resolve().parent
out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=False)
start = time.monotonic()
frozen = json.loads((HERE / 'noodles-balanced-request-freeze-job-01/output/request.json').read_text())
request = frozen['requests'][0]
design = (np.asarray(request['design_vertices']) - request['target_origin']) @ np.asarray(request['target_axes'])
design = [tuple(Fraction(float(x)) for x in row) for row in design]
captures = HERE / 'noodles-same-call-proof-capture-job-01/output'
sources = HERE / 'noodles-interval-algebraic-source-derivative-job-01/output'
saved = HERE / 'noodles-algebraic-production-moments-job-01/output'
states = []
report = dict(scope='TWO_CAPTURED_STATES_EXACT_RECIPE_CHAIN_PILOT',
              material_boolean_calls=0, source_replay_calls=0,
              complete_frozen_sample_certification=False, optimization_admission=False)
try:
    for index in (0, 1):
        capture_path = captures / f'point-{index:03d}.json'
        capture = json.loads(capture_path.read_text())
        source = json.loads((sources / f'point-{index:03d}.json').read_text())
        if hashlib.sha256((saved / f'point-{index:03d}-material.json').read_bytes()).hexdigest() != capture['saved_sha256']:
            raise ValueError('SAVED_CAPTURE_IDENTITY_MISMATCH')
        if source['coefficients'] != capture['coefficients'] or source['source']['scientific_hash'] != capture['evaluation']['source_scientific_hash']:
            raise ValueError('SAVED_SOURCE_IDENTITY_MISMATCH')
        proof = capture['proof']; material = capture['evaluation']
        if proof['source_faces'] != source['source']['faces'] or proof['geometry_hash'] != material['geometry_hash']:
            raise ValueError('SAME_CALL_TAPE_IDENTITY_MISMATCH')
        operands = [design, [tuple(Fraction(float(x)) for x in row) for row in source['source']['vertices']]]
        points = [tuple(Fraction(x) for x in row) for row in proof['exact_material_vertices']]
        incident = [set() for _ in points]
        for face, birth in zip(material['material_faces'], proof['material_birth_faces'], strict=True):
            for v in face: incident[v].add(birth)
        recipes = []
        for vertex, point in enumerate(points):
            try:
                recipe = recover_vertex_recipe(point, incident[vertex], operands,
                    [request['design_faces'], proof['source_faces']])
                if replay_vertex_recipe(recipe, operands) != point:
                    raise ValueError('RECIPE_REPLAY_MISMATCH')
                recipes.append(recipe)
            except ValueError as exc:
                report['first_failure'] = dict(index=index, vertex=vertex, reason=str(exc),
                    point=[str(x) for x in point], birth_faces=sorted(incident[vertex]),
                    resolved_prefix=len(recipes))
                raise
        states.append(dict(points=points, operands=operands, faces=material['material_faces'], recipes=recipes))
        with (out / f'point-{index:03d}-recipes.json').open('x') as stream:
            json.dump(dict(index=index, recipes=recipes, capture_sha256=hashlib.sha256(capture_path.read_bytes()).hexdigest()), stream)
        print(json.dumps(dict(index=index, resolved=len(recipes), kinds=dict(Counter(r[0] for r in recipes)), seconds=time.monotonic()-start)), flush=True)
    base = states[0]
    report['states'] = []
    for index, state in enumerate(states):
        replay = [replay_vertex_recipe(r, state['operands']) for r in base['recipes']]
        if len(set(replay)) != len(replay) or set(replay) != set(state['points']):
            raise ValueError('REPLAY_NOT_EXACT_BIJECTION')
        lookup = {p: i for i, p in enumerate(replay)}
        current = [[lookup[state['points'][v]] for v in face] for face in state['faces']]
        a = oriented_planar_chain_signature(replay, base['faces'])
        b = oriented_planar_chain_signature(replay, current)
        if a['signature'] != b['signature']:
            raise ValueError('EXACT_ORIENTED_SURFACE_CHAIN_DIFFERS')
        report['states'].append(dict(index=index, exact_bijection=True, oriented_chain_equal=True,
            proof_only_null_triangles=[len(a['exact_null_triangle_indices']), len(b['exact_null_triangle_indices'])]))
    report['pair_certified'] = True
except Exception as exc:
    report.update(pair_certified=False, exception_type=type(exc).__name__, exception=str(exc))
finally:
    report['seconds'] = time.monotonic() - start
    with (out / 'report.json').open('x') as stream: json.dump(report, stream, indent=2)
    print(json.dumps(report), flush=True)
if not report['pair_certified']: sys.exit(1)
