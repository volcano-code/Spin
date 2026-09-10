"""Independently integrate actual public-cage output using Python tetrahedra."""
import hashlib
import json
import math
from pathlib import Path
import sys

from make_it_stand.verification.mass_moments import TriangleMesh, compute_mass_properties

paths = [Path(p) for p in sys.argv[1:]]
assert len(paths) == 2, 'two actual probe outputs are required'
expected_q = [[513/4096 if i == j else 27/256 for j in range(3)] for i in range(3)]
for path in paths:
    p = json.loads(path.read_text())
    assert p['scope'] == 'SYNTHETIC_CAGE_GEOMETRY_ONLY_NOT_MATERIAL'
    assert p['selected_count'] == 216 and p['handle_count'] == 218
    assert len(p['vertices']) == 218 and len(p['faces']) == 432
    assert len({tuple(v) for v in p['vertices']}) == 218
    for v in p['vertices']:
        assert all(.125 <= x <= .875 for x in v)
        assert any(x == .125 or x == .875 for x in v)
    # This checker does not call the native source-tree/cage integrator.
    actual = compute_mass_properties(TriangleMesh(p['vertices'], p['faces']))
    pairs = [(actual.mass,27/64)]
    pairs += list(zip(actual.first_moment,[27/128]*3))
    pairs += list(zip(actual.center_of_mass,[.5]*3))
    for i in range(3):
        for j in range(3):
            pairs += [(actual.second_moment[i][j], expected_q[i][j]),
                      (actual.central_inertia[i][j],81/2048 if i == j else 0)]
    assert all(math.isfinite(a) and abs(a-b) <= 1e-12 for a,b in pairs)
    print(json.dumps({'file':str(path),'file_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'output_sha256':p['output_sha256'],'cells':216,'handles':218,'triangles':432,
        'independent_envelope_volume':actual.mass,'maximum_absolute_error':max(abs(a-b) for a,b in pairs),
        'scope':'CAGE_GEOMETRY_ONLY_NOT_MATERIAL_NOT_PAPER_EXPERIMENT'},sort_keys=True))
assert paths[0].read_bytes() == paths[1].read_bytes(), 'actual repeated cage output differs'
print('PASS: actual C++ cage pair byte-identical, independent Python tetrahedron integration matches literal envelope moments')
