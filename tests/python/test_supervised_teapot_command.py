"""Admission/command binding only; these files are not scientific mesh fixtures."""
import hashlib
import importlib.util
import importlib.machinery
import json
from pathlib import Path

import pytest


def adapter():
    path=Path(__file__).resolve().parents[2]/'experiments/teapot_spin2014/supervised.py'
    assert path.is_file(),'existing teapot runner has no whole-job adapter'
    spec=importlib.util.spec_from_file_location('supervised_teapot',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def inputs(tmp_path):
    data=tmp_path/'inputs';data.mkdir()
    hashes={}
    for name in ['outer','design','fixed']:
        p=data/(name+'.obj');p.write_text('admission identity fixture only')
        hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    (data/'manifest.json').write_text(json.dumps({'geometry_hashes':hashes,'source_files':{}}))
    binary=tmp_path/('_make_it_stand'+importlib.machinery.EXTENSION_SUFFIXES[0]);binary.write_bytes(b'not a loadable module: preflight only')
    return data,binary


def test_adapter_keeps_existing_scientific_runner_and_strategy(tmp_path):
    m=adapter();data,binary=inputs(tmp_path)
    request=m.prepare(data,'mass_lowering',binary)
    argv=request['argv']
    assert Path(argv[1]).name=='run.py'
    assert argv[argv.index('--strategy')+1]=='mass_lowering'
    assert argv[argv.index('--output')+1]=='{job}/payload'
    assert request['identity_files']['native']==binary.resolve()
    assert argv[argv.index('--expected-native-module')+1]==str(binary.resolve())
    assert 'maximum_iterations' not in request
    assert request['scope']=='RECONSTRUCTED_EXISTING_TEAPOT_PROTOCOL_NOT_PAPER_ORIGINAL'


def test_changed_or_missing_input_blocks_before_job(tmp_path):
    m=adapter();data,binary=inputs(tmp_path)
    (data/'design.obj').write_text('changed')
    with pytest.raises(ValueError,match='STALE'):m.prepare(data,'inertia',binary)


def test_incomplete_manifest_cannot_admit_inputs(tmp_path):
    m=adapter();data,binary=inputs(tmp_path)
    (data/'manifest.json').write_text(json.dumps({'geometry_hashes':{},'source_files':{}}))
    with pytest.raises(ValueError,match='GEOMETRY'):m.prepare(data,'inertia',binary)


def test_unknown_strategy_cannot_select_an_unreviewed_protocol(tmp_path):
    m=adapter();data,binary=inputs(tmp_path)
    with pytest.raises(ValueError,match='STRATEGY'):m.prepare(data,'faster',binary)


def test_invalid_native_suffix_is_not_a_module(tmp_path):
    m=adapter();data,binary=inputs(tmp_path)
    invalid=tmp_path/'_make_it_stand.not_importable.so';binary.rename(invalid)
    with pytest.raises(ValueError,match='NATIVE_MODULE'):m.prepare(data,'inertia',invalid)


def test_actual_loaded_module_must_be_pinned_file(tmp_path):
    m=adapter();data,binary=inputs(tmp_path)
    other=tmp_path/'different.so';other.write_bytes(b'other')
    with pytest.raises(ValueError,match='LOADED_NATIVE_MISMATCH'):
        m.verify_native_module(binary,other)
    m.verify_native_module(binary,binary)


def test_runner_checks_loaded_identity_before_geometry_or_output(tmp_path,monkeypatch):
    import runpy
    from types import SimpleNamespace
    from make_it_stand.verification import spin_common
    m=adapter();data,binary=inputs(tmp_path)
    other=tmp_path/'different.so';other.write_bytes(b'other')
    monkeypatch.syspath_prepend(str(Path(m.__file__).parent))
    monkeypatch.setattr(spin_common,'native',lambda:SimpleNamespace(__file__=str(other)))
    runner=runpy.run_path(str(Path(m.__file__).with_name('run.py')))
    output=tmp_path/'not-created'
    with pytest.raises(ValueError,match='LOADED_NATIVE_MISMATCH'):
        runner['execute'](data,output,'inertia',binary)
    assert not output.exists()


def test_relative_entrypoint_still_produces_absolute_runner(tmp_path,monkeypatch):
    import runpy
    m=adapter();data,binary=inputs(tmp_path)
    monkeypatch.chdir(Path(m.__file__).resolve().parents[2])
    entry=runpy.run_path('experiments/teapot_spin2014/supervised.py')
    request=entry['prepare'](data,'inertia',binary)
    assert Path(request['argv'][1]).is_absolute()
