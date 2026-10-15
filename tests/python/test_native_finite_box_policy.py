import copy
import pytest
from test_native_tangent_policy import arguments, solve


def policy():
    return dict(direction_policy='FINITE_BOX_QP_V1', finite_box_maximum_faces=256,
                finite_box_maximum_iterations=1028, finite_box_residual_tolerance=1e-9)


def test_full_box_policy_roundtrip_and_default_is_unchanged():
    p, kw = arguments()
    baseline = solve(p, kw)
    kw['solver_settings']['slqp'].update(policy())
    try:
        result = solve(p, kw)
    except ValueError as exc:
        pytest.fail('native finite-box direction policy is missing: ' + str(exc))
    assert result['settings']['solver_settings']['slqp'] == kw['solver_settings']['slqp']
    assert result == solve(p, kw)
    assert result['solve']['input_hash'] != baseline['solve']['input_hash']
    events = [t['finite_box'] for p in result['solve']['passes'] for t in p['trace'] if 'finite_box' in t]
    assert events
    assert all(q['status'].startswith('BOX_QP_') for q in events)
    assert all('finite_box' not in t for p in baseline['solve']['passes'] for t in p['trace'])
    changed = copy.deepcopy(kw)
    changed['solver_settings']['slqp']['finite_box_maximum_faces'] = 128
    assert solve(p, changed)['solve']['input_hash'] != result['solve']['input_hash']


@pytest.mark.parametrize('change', [dict(finite_box_maximum_faces=0),
    dict(finite_box_maximum_iterations=0), dict(finite_box_residual_tolerance=float('nan'))])
def test_invalid_full_box_controls_rejected(change):
    p, kw = arguments()
    kw['solver_settings']['slqp'].update(policy())
    kw['solver_settings']['slqp'].update(change)
    with pytest.raises(ValueError):
        solve(p, kw)


def test_finite_box_disk_trace_roundtrip(tmp_path):
    from test_source_online_trace import run, reconstruct
    from dataclasses import asdict
    from test_source_spin_material_bridge import settings
    from make_it_stand.verification import source_spin_material as bridge
    controls = asdict(settings(bridge)['solver_settings'])['slqp']
    controls.update(policy())
    disk = run(tmp_path/'box', slqp=controls)
    memory = run(slqp=controls)
    assert disk.diagnostics == memory.diagnostics == reconstruct(disk.iter_diagnostics())
    assert any('finite_box' in t for p in disk.diagnostics['solve']['passes'] for t in p['trace'])
