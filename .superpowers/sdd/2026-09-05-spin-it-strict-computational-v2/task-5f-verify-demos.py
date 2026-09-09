"""Controller-only independent saved-result checks; no optimization or edits."""
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
from make_it_stand.verification.spin_bundle import load_spin_artifact
from make_it_stand.verification.spin_common import native

native()  # Freshness observes the actually loaded extension, never imports it implicitly.
roots = [Path(f'build/task-5f-demo-run-{i}') for i in (1, 2)]
world = [503/64, 0, 0, 2003/256, 32687/12288, 32687/12288, 128363/12288, 0, 0, 0]
normalized = [503/512, 0, 0, -9/4096, 32687/393216, 32687/393216, 32651/393216, 0, 0, 0]
artifacts = []
for root in roots:
    scientific = json.loads((root/'scientific.json').read_text())
    identity = scientific['identity']
    evidence = PaperRunEvidenceV2.from_json((root/'evidence.json').read_text())
    fresh = check_freshness(evidence)
    assert fresh.status == 'CURRENT', asdict(fresh)
    artifact = load_spin_artifact(root/'sealed', expected_problem_digest=identity['problem_digest'])
    assert artifact.scientific_digest == identity['artifact_scientific_digest']
    facts = artifact.facts
    np.testing.assert_allclose(facts['physical']['moments'], world, rtol=0, atol=1e-12)
    np.testing.assert_allclose(facts['normalized']['moments'], normalized, rtol=0, atol=1e-14)
    np.testing.assert_allclose(facts['physical']['central_inertia'], np.diag([16432021/3090432]*2+[32687/6144]), rtol=0, atol=1e-12)
    assert facts['clearance']['squared_distance_exact'] == '25/64'
    assert facts['bounded_void_components'] == facts['newly_carved_components'] == 1
    for kind, a in zip(('vertices', 'faces'), artifact.outer.arrays()):
        original = np.load(root/f'outer_{kind}.npy', allow_pickle=False)
        assert original.dtype == a.dtype and original.shape == a.shape and original.tobytes() == a.tobytes()
    native = identity['native_run']
    solve = native['solve']
    assert sum(len(p['trace']) for p in solve['passes']) == 39
    assert native['geometry']['union']['input_cuboids'] == 72
    assert solve['input_hash'] == '83359bc0abc68ff056f49e18ec52ad6bb34d1121e0743f28d074cbab22d7f7a9'
    assert solve['trace_hash'] == '72002479f41c4a8a4ef71731e1ff8bc41c42537cd2990e8dbd3417130838d246'
    assert solve['output_hash'] == '39b18d241e42e8978295e6bb97c7e2e0c743c59893cee9ca609dbc355b6d9283'
    assert native['geometry']['material_mesh_hash'] == '920dd72259050025e6ed85417082de89d43008c31d99674169c6b42d58238a9b'
    artifacts.append(artifact)
    print(json.dumps({'directory': str(root), 'freshness': asdict(fresh), 'reintegrated_artifact': artifact.scientific_digest,
        'scientific_digest': scientific['scientific_digest'], 'max_normalized_error': max(facts['moment_comparison']['absolute_errors']),
        'ordered_records': 39, 'input_cuboids': 72, 'actual_cavities': facts['bounded_void_components']}, sort_keys=True))
for name in ('inputs.json', 'scientific.json', 'native-run.json'):
    a,b = [(r/name).read_bytes() for r in roots]
    assert a == b, name
    print(name, hashlib.sha256(a).hexdigest(), 'BYTE_IDENTICAL')
for attr in ('outer', 'final_material', 'newly_carved'):
    assert getattr(artifacts[0], attr) == getattr(artifacts[1], attr), attr
assert artifacts[0].scientific_digest == artifacts[1].scientific_digest
print('PASS: two actual runs, reloaded independent seals, analytic moments, full inertia, outer bytes, geometry and trace determinism')
