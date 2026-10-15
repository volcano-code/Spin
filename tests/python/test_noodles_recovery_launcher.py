import importlib.util
from pathlib import Path
from types import SimpleNamespace
import hashlib
import sys
import pytest


def test_offline_launcher_uses_independent_boundary_and_new_job(monkeypatch):
    root = Path(__file__).resolve().parents[2]
    path = root / 'experiments/faster_core_repair_20261012/launch_noodles_terminal_recovery_v2.py'
    spec = importlib.util.spec_from_file_location('noodles_recovery_launcher', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    calls = []
    def run(*args, **kwargs):
        calls.append((args, kwargs))
        return SimpleNamespace(status='PROCESS_EXITED', returncode=0)
    monkeypatch.setattr(module, 'run_serial_source_job', run)
    module.main()
    args, _ = calls[0]
    assert args[1].name == 'noodles-terminal-recovery-job-03'
    assert '--classified-boundary' not in args[2]
    assert '--reuse-input-proof' in args[2]
    assert '--owned-terminal' in args[2]


def test_recovery_preflight_rejects_backend_and_input_tampering(monkeypatch, tmp_path):
    root = Path(__file__).resolve().parents[2]
    path = root / 'experiments/faster_core_repair_20261012/run_noodles_terminal_recovery.py'
    spec = importlib.util.spec_from_file_location('noodles_recovery_worker', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    binary = tmp_path / 'bridge.so'
    binary.write_bytes(b'bridge')
    monkeypatch.setitem(sys.modules, '_cherchi_core_boolean', SimpleNamespace(__file__=str(binary)))
    identity = 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:' + hashlib.sha256(b'bridge').hexdigest()
    config = dict(boolean_backend_identity=identity)
    source = tmp_path / 'input.obj'
    source.write_bytes(b'original')
    frozen = dict(input_sources={str(source): hashlib.sha256(b'original').hexdigest()})
    module.preflight_saved_backend(config, frozen)
    with pytest.raises(ValueError, match='RECOVERY_NATIVE_BACKEND_CHANGED'):
        module.preflight_saved_backend(dict(boolean_backend_identity='wrong'), frozen)
    source.write_bytes(b'mutated')
    with pytest.raises(ValueError, match='FROZEN_INPUT_CHANGED'):
        module.preflight_saved_backend(config, frozen)
