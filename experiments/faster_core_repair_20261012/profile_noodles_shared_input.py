"""Actual frozen noodles shared-operand comparison; no optimization claim."""
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import trimesh
from make_it_stand.verification.spin_problem import MeshSnapshot, partition
from make_it_stand.verification.cherchi_run_journal import RunJournal

HERE = Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    frozen = json.loads((HERE/'noodles-balanced-request-freeze-job-01/output/request.json').read_text())
    for name, sha in frozen['input_sources'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != sha:
            raise ValueError('FROZEN_INPUT_CHANGED')
    paths = {Path(p).name: Path(p) for p in frozen['input_sources']}
    with RunJournal(output, durability_profile='PER_EVENT_FSYNC_V1') as journal:
        times = {}
        def stage(name, fn):
            start = time.monotonic()
            journal.emit(dict(kind='STAGE_STARTED', stage=name))
            value = fn()
            times[name] = time.monotonic()-start
            journal.emit(dict(kind='STAGE_RETURNED', stage=name, seconds=times[name]))
            return value
        meshes = {r: trimesh.load(paths[r+'-si-candidate.obj'], force='mesh', process=False)
                  for r in ('outer', 'design', 'fixed')}
        outer = stage('outer', lambda: MeshSnapshot((meshes['outer'].vertices, meshes['outer'].faces)))
        design = stage('design', lambda: MeshSnapshot((meshes['design'].vertices, meshes['design'].faces)))
        fixed = stage('owned_fixed', lambda: MeshSnapshot.interior_difference(outer, design))
        # Frozen fixed boundary must be literally the shared boundary, before
        # accepting a newly identified owned representation of this input.
        vertices = np.vstack((meshes['outer'].vertices, meshes['design'].vertices))
        faces = np.vstack((meshes['outer'].faces, meshes['design'].faces[:, ::-1]+len(meshes['outer'].vertices)))
        if not (np.array_equal(vertices, meshes['fixed'].vertices) and
                np.array_equal(faces, meshes['fixed'].faces)):
            raise ValueError('FROZEN_FIXED_NOT_SHARED_BOUNDARY')
        facts = stage('independent_partition', lambda: partition(outer.exact_region(),
                    [fixed.exact_region(), design.exact_region()]))
        journal.publish('report.json', dict(schema='noodles-owned-input-probe-v1',
            fixed_identity=fixed.identity, partition=dict(facts), stage_seconds=times,
            frozen_shared_boundary_equal=True, optimization_invocations=0,
            case_complete=False, performance_gate_passed=False))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
