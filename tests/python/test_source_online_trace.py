"""Exercise real source solves: spooling must preserve the scientific payload."""
import gc
from dataclasses import asdict
import pytest
import _make_it_stand as native
from make_it_stand.verification import source_spin_material as bridge
from test_source_spin_material_bridge import problem, settings
from test_source_trace_stream import reconstruct


def run(directory=None, max_bytes=1073741824, **changes):
    p = problem()
    kw = {key: asdict(value) for key, value in settings(bridge).items()}
    kw['solver_settings'].update(changes)
    if directory is not None:
        kw.update(iteration_store_directory=str(directory), iteration_store_max_bytes=max_bytes)
    return native.solve_source_spin2014_material(*p.original.arrays(),
        *p.design_regions[0].mesh.arrays(), **kw, density=1., minimum_outer_wall=.5)


@pytest.mark.parametrize('changes', [{}, {'maximum_leaf_count': 1},
    {'termination_policy': 'PAPER_POST_TRANSITION_V1'}, {'maximum_adaptive_passes': 1}])
def test_online_exact_science_and_bounded_history(tmp_path, changes):
    memory = run(**changes)
    disk = run(tmp_path / 'trace', **changes)
    expected = memory.diagnostics
    assert disk.diagnostics == expected
    assert reconstruct(disk.iter_diagnostics()) == expected
    stats = disk.trace_storage
    count = sum(len(p['trace']) for p in expected['solve']['passes'])
    assert stats['record_count'] == count
    assert stats['retained_iteration_records'] == 0
    assert stats['peak_retained_iteration_records'] <= 1
    if count:
        assert stats['bytes_written'] > 0
        assert expected['solve']['passes'][-1]['trace'][-1]['final_state']
    assert stats['resident_pass_count'] == len(expected['solve']['passes'])


def test_online_repeat_read_and_owner_lifetime(tmp_path):
    disk = run(tmp_path / 'trace')
    expected = disk.diagnostics
    first, second = disk.iter_diagnostics(), disk.iter_diagnostics()
    del disk
    gc.collect()
    assert reconstruct(first) == expected
    assert reconstruct(second) == expected
    assert len(expected['solve']['passes']) > 1


@pytest.mark.parametrize('mutation', ['truncate', 'payload', 'version', 'missing_seal', 'length', 'count'])
def test_online_corruption_fails_closed(tmp_path, mutation):
    disk = run(tmp_path / 'trace')
    path = next((tmp_path / 'trace').glob('*.slqptrace'))
    data = bytearray(path.read_bytes())
    if mutation == 'truncate': data = data[:len(data)//2]
    elif mutation == 'missing_seal': data = data[:-24]
    elif mutation == 'version': data[8] ^= 0x40
    elif mutation == 'length': data[24:32] = b'\xff' * 8
    elif mutation == 'count': data[-16:-8] = b'\xff' * 8
    else: data[100] ^= 0x40
    path.write_bytes(data)
    with pytest.raises(RuntimeError): list(disk.iter_diagnostics())


def test_existing_path_and_resource_limit_preserve_evidence(tmp_path):
    directory = tmp_path / 'existing'
    directory.mkdir()
    marker = directory / 'marker'
    marker.write_text('preserve')
    with pytest.raises(RuntimeError): run(directory)
    assert marker.read_text() == 'preserve'
    failed = tmp_path / 'limited'
    with pytest.raises(RuntimeError): run(failed, max_bytes=64)
    assert list(failed.glob('*.slqptrace'))
    with pytest.raises(RuntimeError): run(tmp_path / 'absent' / 'child')


def test_source_cumulative_budget_and_metric_tangent(tmp_path):
    controls = asdict(settings(bridge)['solver_settings'])['slqp']
    controls.update(direction_policy='METRIC_TANGENT_RESCUE_V1',
                    metric_tangent_maximum_faces=256,
                    metric_tangent_maximum_iterations=1028,
                    metric_tangent_residual_tolerance=1e-9)
    disk = run(tmp_path / 'metric', slqp=controls)
    assert disk.diagnostics == run(slqp=controls).diagnostics
    from pathlib import Path
    files = [Path(p) for p in disk.trace_storage['files']]
    assert len(files) > 1
    # A budget that fits the largest individual pass cannot fit the whole run.
    with pytest.raises(RuntimeError):
        run(tmp_path / 'cumulative', max_bytes=max(p.stat().st_size for p in files), slqp=controls)


def test_relative_store_survives_cwd_changes_before_first_and_repeat_reads(tmp_path, monkeypatch):
    from pathlib import Path
    expected = run().diagnostics
    monkeypatch.chdir(tmp_path)
    disk = run('relative-trace')
    first = disk.iter_diagnostics()
    monkeypatch.chdir(tmp_path.parent)
    assert reconstruct(first) == expected
    assert all(Path(path).is_absolute() for path in disk.trace_storage['files'])
    monkeypatch.chdir('/')
    assert reconstruct(disk.iter_diagnostics()) == expected
    assert disk.diagnostics == expected
