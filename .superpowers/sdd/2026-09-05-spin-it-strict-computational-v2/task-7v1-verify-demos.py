"""Controller verification of actual runs, using locked analytical expectations."""
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
from make_it_stand.verification.spin_bundle import load_spin_artifact
from make_it_stand.verification.spin_common import native
from make_it_stand.verification.spin_view import load_spin_scene

native()
EXPECTED_PROBLEM = '2f45a44a9c74fd3d02bc68a7c36e11b8d7fdee0a0fe22a9857236338d6ed3929'
roots = [Path(f'build/task-7v1-source-run-{i}') for i in (1, 2)]
world = [503/64, 0, 0, 2003/256, 32687/12288, 32687/12288, 128363/12288, 0, 0, 0]
normalized = [503/512, 0, 0, -9/4096, 32687/393216, 32687/393216, 32651/393216, 0, 0, 0]
inertia = np.diag([16432021/3090432]*2+[32687/6144])
artifacts = []
for root in roots:
    scientific = json.loads((root/'scientific.json').read_text())
    identity = scientific['identity']
    fresh = check_freshness(PaperRunEvidenceV2.from_json((root/'evidence.json').read_text()))
    assert fresh.status == 'CURRENT', asdict(fresh)
    assert identity['problem_digest'] == EXPECTED_PROBLEM
    artifact = load_spin_artifact(root/'sealed', expected_problem_digest=EXPECTED_PROBLEM)
    assert artifact.scientific_digest == identity['artifact_scientific_digest']
    facts = artifact.facts
    np.testing.assert_allclose(facts['physical']['moments'], world, rtol=0, atol=1e-12)
    np.testing.assert_allclose(facts['normalized']['moments'], normalized, rtol=0, atol=1e-14)
    np.testing.assert_allclose(facts['physical']['central_inertia'], inertia, rtol=0, atol=1e-12)
    assert facts['clearance']['squared_distance_exact'] == '25/64'
    assert facts['bounded_void_components'] == facts['newly_carved_components'] == 1
    for kind, a in zip(('vertices', 'faces'), artifact.outer.arrays()):
        original = np.load(root/f'outer_{kind}.npy', allow_pickle=False)
        assert original.dtype == a.dtype and original.shape == a.shape and original.tobytes() == a.tobytes()
    actual = identity['native_run']
    solve = actual['solve']
    assert sum(len(p['trace']) for p in solve['passes']) == 39
    assert actual['geometry']['union']['input_cuboids'] == 72
    assert solve['input_hash'] == '83359bc0abc68ff056f49e18ec52ad6bb34d1121e0743f28d074cbab22d7f7a9'
    assert solve['trace_hash'] == '72002479f41c4a8a4ef71731e1ff8bc41c42537cd2990e8dbd3417130838d246'
    assert solve['output_hash'] == '39b18d241e42e8978295e6bb97c7e2e0c743c59893cee9ca609dbc355b6d9283'
    assert actual['geometry']['material_mesh_hash'] == '920dd72259050025e6ed85417082de89d43008c31d99674169c6b42d58238a9b'
    scene = load_spin_scene(root/'sealed', expected_problem_digest=EXPECTED_PROBLEM)
    assert scene['artifact_digest'] == artifact.scientific_digest
    np.testing.assert_allclose(scene['physical']['central_inertia'], inertia, rtol=0, atol=1e-12)
    assert sorted(scene['principal_axes']['multiplicities']) == [1, 2]
    assert len(scene['principal_axes']['unique_axes']) == 1
    np.testing.assert_array_equal(scene['principal_axes']['unique_axes'][0]['direction'], [0, 0, 1])
    assert scene['support_display'] == 'TARGET_POINTS_ONLY_NO_CERTIFIED_HULL'
    assert sorted(m['kind'] for m in scene['meshes']) == ['cavity','design','fixed','material','original','outer','retained']
    for m in scene['meshes']:
        expected_counts = (8, 12) if m['kind'] in ('original','outer','design','cavity') else (16, 24)
        assert (len(m['vertices']), len(m['faces'])) == expected_counts
    artifacts.append(artifact)
    print(json.dumps({'directory': str(root), 'freshness': asdict(fresh),
        'reintegrated_artifact': artifact.scientific_digest, 'scientific_digest': scientific['scientific_digest'],
        'max_normalized_error': max(facts['moment_comparison']['absolute_errors']),
        'ordered_records': 39, 'input_cuboids': 72, 'actual_cavities': facts['bounded_void_components'],
        'principal_axes': scene['principal_axes'], 'scene_layer_count': len(scene['meshes'])}, sort_keys=True))
for name in ('inputs.json', 'scientific.json', 'native-run.json'):
    a,b = [(r/name).read_bytes() for r in roots]
    assert a == b, name
    print(name, hashlib.sha256(a).hexdigest(), 'BYTE_IDENTICAL')
for attr in ('outer', 'final_material', 'newly_carved'):
    assert getattr(artifacts[0], attr) == getattr(artifacts[1], attr), attr
assert artifacts[0].scientific_digest == artifacts[1].scientific_digest
print('PASS: actual source solve pair, reintegration, analytic moments, outer bytes, trace determinism and independently admitted viewer layers')
