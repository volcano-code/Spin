"""Persist/replay actual owned set operations without changing scientific facts."""
import json
from unittest.mock import patch

import numpy as np
import pytest

from make_it_stand.verification.spin_problem import MeshSnapshot, SpinArtifactProblem
from make_it_stand.verification.spin_artifact import SealedSpinArtifact
from make_it_stand.verification.spin_common import native, FRAME, ORDER
from make_it_stand.verification.spin_bundle import save_spin_artifact, load_spin_artifact
from test_exact_material_region import box, combine


def artifact(owned=True):
    outer, cavity = box(), box((.25,)*3, (.75,)*3)
    outer[0][0, 0] = -0.
    retained = combine(outer, (cavity[0], cavity[1][:, ::-1]))
    problem = SpinArtifactProblem(original=outer, fixed_regions=[], design_regions=[('d', outer, 1.)],
        length_unit='m', mass_unit='kg', minimum_outer_wall=.25,
        target=dict(mode='SPINNING_TOP', axis=[0., 0., 1.], origin=[0., 0., 0.], support_points=[]))
    material = retained
    if owned:
        region = native().make_exact_interior_difference_region(problem.original.exact_region(),
                                                                 MeshSnapshot(cavity).exact_region())
        material = MeshSnapshot.from_exact_region(region)
    claim = dict(frame=FRAME, order=list(ORDER), moments=[7/8, 0, 0, 0, 31/384, 31/384, 31/384, 0, 0, 0])
    return SealedSpinArtifact.create(problem, outer=outer, design_results=[('d', retained, cavity)],
        final_material=material, newly_carved=cavity, solver_claim=claim)


def test_owned_operation_roundtrip_recomputes_recipe_and_preserves_facts(tmp_path):
    original = artifact()
    with patch.object(SealedSpinArtifact, 'revalidate',
                      side_effect=AssertionError('save repeated an already completed seal')):
        loaded = save_spin_artifact(original, tmp_path/'bundle')
    assert loaded.facts == original.facts
    assert loaded.scientific_digest == original.scientific_digest
    assert np.signbit(loaded.problem.original.arrays()[0][0, 0])
    manifest = json.loads((tmp_path/'bundle/manifest.json').read_text())
    assert manifest['schema'] == 'sealed-spin-bundle-v3-owned-operations'
    recipe = manifest['final']['material']['construction']
    assert recipe['operation'] == 'STRICT_INTERIOR_DIFFERENCE_V1' and len(recipe['operands']) == 2
    # Binary operands and the retained mesh may be rebuilt; the actual final
    # rational boundary must come from the replayed owned operation.
    with patch.object(native(), 'make_exact_material_region_from_rationals',
                      side_effect=AssertionError('final boundary was independently rebuilt instead of replayed')):
        fresh = load_spin_artifact(tmp_path/'bundle', expected_problem_digest=original.problem.input_digest)
    assert fresh.scientific_digest == original.scientific_digest


def test_default_nonowned_bundle_keeps_v2_wire_format(tmp_path):
    save_spin_artifact(artifact(owned=False), tmp_path/'bundle')
    manifest = json.loads((tmp_path/'bundle/manifest.json').read_text())
    assert manifest['schema'] == 'sealed-spin-bundle-v2'
    assert 'construction' not in manifest['final']['material']


def test_reload_reuses_only_the_input_proof_it_just_constructed(tmp_path):
    original = artifact(owned=False)
    save_spin_artifact(original, tmp_path/'bundle')
    from make_it_stand.verification import spin_artifact
    with patch.object(spin_artifact, '_recheck_problem',
                      wraps=spin_artifact._recheck_problem) as recheck:
        loaded = load_spin_artifact(tmp_path/'bundle',
            expected_problem_digest=original.problem.input_digest)
    assert loaded.scientific_digest == original.scientific_digest
    assert any(call.kwargs.get('reuse_validation') is True for call in recheck.call_args_list)


@pytest.mark.parametrize('fault', ['unknown_operation', 'reversed_operands', 'self_reported_flag', 'removed_recipe'])
def test_modified_owned_recipe_is_rejected(tmp_path, fault):
    original = artifact()
    save_spin_artifact(original, tmp_path/'bundle')
    path = tmp_path/'bundle/manifest.json'
    manifest = json.loads(path.read_text())
    recipe = manifest['final']['material']['construction']
    if fault == 'unknown_operation':
        recipe['operation'] = 'TRUST_ME'
    elif fault == 'reversed_operands':
        recipe['operands'].reverse()
    elif fault == 'self_reported_flag':
        recipe['validated'] = True
    else:
        manifest['final']['material']['construction'] = None
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        load_spin_artifact(tmp_path/'bundle', expected_problem_digest=original.problem.input_digest)
