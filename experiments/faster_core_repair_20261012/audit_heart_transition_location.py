"""Auxiliary distance diagnostics for the isolated cone limit, not a certificate."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import trimesh


def main(output):
    here = Path(__file__).resolve().parent
    run = here/'heart-guarded-formal-seed-0-job-01/output/run'
    proof_path = here/'heart-exact-quadric-transition-01.json'
    proof = json.loads(proof_path.read_text())
    config_path = run/'configuration.json'
    request = json.loads(config_path.read_text())['request']
    q = []
    for state in proof['states']:
        path = run/f"evaluation-{state['evaluation']:06d}.json"
        if hashlib.sha256(path.read_bytes()).hexdigest() != state['material_sha256']:
            raise ValueError('SAVED_MATERIAL_CHANGED')
        q.append(np.asarray(json.loads(path.read_text())['quadric_coefficients']))
    from fractions import Fraction
    interval = proof['degeneracy_parameter_intervals'][0]
    t = float((Fraction(interval[0])+Fraction(interval[1]))/2)
    c = (1-t)*q[0]+t*q[1]
    A = np.diag(c[4:7]); A[0,1]=A[1,0]=c[9]/2
    A[0,2]=A[2,0]=c[8]/2; A[1,2]=A[2,1]=c[7]/2
    center = -np.linalg.solve(A, c[1:4])/2
    point = center@np.asarray(request['target_axes']).T+np.asarray(request['target_origin'])
    queries = {}
    for name in ('design', 'outer'):
        mesh = trimesh.Trimesh(request[name+'_vertices'], request[name+'_faces'], process=False)
        _, distances, faces = trimesh.proximity.closest_point(mesh, [point])
        queries[name] = dict(approximate_contains=bool(mesh.contains([point])[0]),
            approximate_distance=float(distances[0]), closest_face=int(faces[0]))
    report = dict(schema='heart-cone-location-auxiliary-v1',
        proof_sha256=hashlib.sha256(proof_path.read_bytes()).hexdigest(),
        configuration_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),
        parameter_midpoint=t, center_target=center.tolist(), center_request=point.tolist(),
        queries=queries, domain_intersection_certified=False,
        inference='Auxiliary queries suggest cone apex lies in fixed material rather than design domain; exact interval location remains required',
        stagnation_root_cause_certified=False, case_complete=False)
    with output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(report))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
