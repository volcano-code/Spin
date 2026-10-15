"""Existing saved-input regression, not a Cherchi derivative admission.

Catches lost translation, repeated scaling/density and a mismatched mesh when
preparing the physical Boolean view. No Boolean or optimization is run here.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import trimesh

HERE = Path(__file__).resolve().parent


def test_saved_noodles_request_matches_actual_si_input():
    frozen = json.loads((HERE / 'noodles-balanced-request-freeze-job-01/output/request.json').read_text())
    for name, digest in frozen['input_sources'].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == digest
    source = HERE / 'noodles-balanced-shared-export-job-02/output'
    manifest = json.loads((source / 'manifest.json').read_text())
    model = json.loads((HERE.parents[1] / 'experiments/quadric_physics_demo_20261011/noodles-dense-v5/manifest.json').read_text())
    spec = model['spec']
    normalization = frozen['normalization']
    length = normalization['L']
    origin = np.asarray(normalization['origin'])
    assert normalization['mass_unit_kg'] == normalization['R'] * length ** 3
    assert normalization['R'] == spec['density_kg_m3']
    meshes = {name: trimesh.load(source / manifest['files'][name]['name'], force='mesh', process=False)
              for name in ('outer', 'design', 'fixed')}
    assert length == float(np.max(np.ptp(meshes['outer'].vertices, axis=0)))
    for request in frozen['requests']:
        for name in ('outer', 'design'):
            np.testing.assert_array_equal(request[name + '_faces'], meshes[name].faces)
            np.testing.assert_allclose(request[name + '_vertices'],
                (meshes[name].vertices - origin) / length, atol=1e-14, rtol=0)
        np.testing.assert_allclose(request['target_origin'],
            (np.asarray(spec['target_origin_model']) * spec['model_unit_m'] - origin) / length,
            atol=1e-14, rtol=0)
        assert request['material_density'] == 1.
        assert request['original_units_per_normalized_unit'] == length
        assert request['center_of_mass_weight'] == spec['center_of_mass_weight']
        assert spec['weight_coordinate_contract'] == 'FASTER_NORMALIZED_REQUEST'
    fixed = json.loads((HERE / 'noodles-balanced-request-freeze-job-01/output/fixed.binary64.json').read_text())
    np.testing.assert_array_equal(fixed['faces'], meshes['fixed'].faces)
    np.testing.assert_allclose(fixed['vertices'], (meshes['fixed'].vertices-origin)/length, atol=1e-14, rtol=0)
