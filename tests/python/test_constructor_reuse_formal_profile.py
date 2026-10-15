"""Explicit input-proof reuse must retain real geometry and terminal checks."""
import json
from unittest.mock import patch

import pytest

from make_it_stand import spin
from make_it_stand.verification import faster_material as fm
from make_it_stand.verification import cherchi_inprocess as ci
from make_it_stand.verification import spin_problem as p
from make_it_stand.verification.spin_artifact import SealedSpinArtifact
from make_it_stand.verification.spin_common import FRAME, ORDER
from test_cherchi_rational_boundary_profile import CORE, actual_material
from test_cherchi_terminal_checkpoints import terminal_solver, settings
from test_exact_material_region import box, combine


def test_formal_entry_uses_constructor_certificate_and_records_explicit_profile(tmp_path, monkeypatch, actual_material):
    problem, _, q, _ = actual_material
    terminal_solver(monkeypatch, q)
    monkeypatch.setattr(ci.InProcessCherchiEvaluator, 'complement',
                        lambda *a, **k: (_ for _ in ()).throw(KeyboardInterrupt('test terminal checkpoint')))
    evaluator, solver = settings()
    with patch.object(p, 'partition', side_effect=AssertionError('duplicate input partition')):
        with pytest.raises(KeyboardInterrupt):
            ci.solve_cherchi_inprocess_and_seal_faster(problem,
                evaluator_settings=evaluator, solver_settings=solver, output_directory=tmp_path/'run',
                rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1',
                solver_boundary_profile='EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1',
                input_validation_profile='CONSTRUCTOR_PROOF_REUSE_V1')
    config = json.loads((tmp_path/'run/configuration.json').read_text())
    assert config['input_validation_profile'] == 'CONSTRUCTOR_PROOF_REUSE_V1'
    assert (tmp_path/'run/native-result.json').is_file()
    assert not (tmp_path/'run/complete.json').exists()


def test_terminal_input_reuse_keeps_complete_material_partition_and_moment_gate():
    outer, cavity = box(), box((.25,)*3, (.75,)*3)
    retained = combine(outer, (cavity[0], cavity[1][:, ::-1]))
    problem = p.SpinArtifactProblem(original=outer, fixed_regions=[], design_regions=[('d', outer, 1.)],
        length_unit='m', mass_unit='kg', minimum_outer_wall=.25,
        target=dict(mode='SPINNING_TOP', axis=[0., 0., 1.], origin=[0., 0., 0.], support_points=[]))
    claim = dict(frame=FRAME, order=list(ORDER), moments=[7/8, 0, 0, 0, 31/384, 31/384, 31/384, 0, 0, 0])
    def seal(material, solver_claim):
        return SealedSpinArtifact.create(problem, outer=outer, design_results=[('d', retained, cavity)],
            final_material=material, newly_carved=cavity, solver_claim=solver_claim,
            _reuse_problem_validation=True)
    with patch.object(p, 'partition', side_effect=AssertionError('duplicate input partition')):
        result = seal(retained, claim)
        assert result.facts['geometry_status'] == 'MESH_MATERIAL_VALIDATED'
        wrong = dict(claim, moments=[1., *claim['moments'][1:]])
        with pytest.raises(ValueError, match='moment claim mismatch'):
            seal(retained, wrong)
        with pytest.raises(ValueError):
            seal(outer, claim)  # Missing the actual cavity still fails.
