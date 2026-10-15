"""One preregistered Heart internal-domain construction, no optimizer invoked."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
import trimesh
from clearance_domain import extract_distance_domain
from rounded_input import SPEC, publish, write_mesh, preview

CONSTRUCTION = dict(schema='rounded-heart-distance-domain-input-v1',
    minimum_clearance_model=.015, candidate_grid_spacing_model=.02,
    interpolation_guard='TWO_ACTUAL_GRID_CELL_DIAGONALS',
    distance_query='TRIMESH_SIGNED_DISTANCE_BINARY64_CANDIDATE_ONLY',
    fixed_pin_exclusion=True, optimization_invocations=0, simulation_invocations=0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    old = json.loads((args.source/'preparation.json').read_text())
    if json.loads((args.source/'spec.json').read_text()) != SPEC:
        raise ValueError('SOURCE_SPEC_REVISION_MISMATCH')
    for name, sha in old['sources'].items():
        if hashlib.sha256((args.source/name).read_bytes()).hexdigest() != sha:
            raise ValueError('STALE_FROZEN_HEART_SOURCE: '+name)
    body, pin = [trimesh.load(args.source/name, force='mesh', process=False)
                 for name in ('heart-body.obj', 'pin.obj')]
    source_hashes = {name: hashlib.sha256((args.source/name).read_bytes()).hexdigest()
                     for name in ('heart-body.obj', 'pin.obj', 'spec.json')}
    prereg = dict(CONSTRUCTION, frozen_outer_source=source_hashes,
                  inherited_physics_and_solver=SPEC,
                  status='PREREGISTERED_SINGLE_CONSTRUCTION_NOT_ADMITTED')
    publish(args.output/'construction.json', prereg)
    began = time.monotonic()
    with (args.output/'events.jsonl').open('x') as events:
        last = -1
        def event(row):
            nonlocal last
            bucket = int(100*row.get('completed', 0)/max(1, row.get('total', 1)))
            if row.get('stage') == 'DISTANCE_GRID' and bucket == last:
                return
            last = bucket
            value = dict(row, elapsed_seconds=time.monotonic()-began)
            events.write(json.dumps(value)+'\n'); events.flush(); os.fsync(events.fileno())
            print(json.dumps(value), flush=True)
        try:
            def distance(points):
                return np.minimum(trimesh.proximity.signed_distance(body, points),
                                  -trimesh.proximity.signed_distance(pin, points))
            domain, report = extract_distance_domain(distance,
                body.bounds[0]-.04, body.bounds[1]+.04,
                spacing=CONSTRUCTION['candidate_grid_spacing_model'],
                clearance=CONSTRUCTION['minimum_clearance_model'], event=event)
            for name, mesh in (('heart-body.obj', body), ('pin.obj', pin), ('inner-source.obj', domain)):
                write_mesh(args.output/name, mesh)
            publish(args.output/'candidate-report.json', report)
            preview(args.output, body, domain, pin, filename='distance-domain-preview.png')
            publish(args.output/'terminal.json', dict(status='CANDIDATE_PREPARED_NOT_ADMITTED',
                completed=True, optimization_invocations=0, simulation_invocations=0,
                elapsed_seconds=time.monotonic()-began,
                files={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in args.output.glob('*.obj')}))
        except Exception as exc:
            publish(args.output/'failure.json', dict(type=type(exc).__name__, message=str(exc),
                completed=False, elapsed_seconds=time.monotonic()-began))
            raise


if __name__ == '__main__':
    main()
