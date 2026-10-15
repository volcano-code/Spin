"""Large provenance must not weaken transitive validation or durable publication."""
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import sys

import pytest

# This can exercise the legacy implementation to demonstrate the real regression.
import importlib
api = importlib.import_module(os.environ.get('SPIN_CHECKPOINT_TEST_MODULE',
    'make_it_stand.verification.source_checkpoint_store_v2'))


def open_store(tmp_path, files, **options):
    return api.StageCheckpointStore(tmp_path / 'store', files,
        maximum_wall_ns=100, maximum_trace_bytes=1000, maximum_payload_bytes=4096, **options)


def publish(store, index=0, **options):
    values = dict(payload=b'native-state\x00', cumulative_wall_ns=10+index,
                  cumulative_trace_bytes=20+index, artifacts={})
    values.update(options)
    return store.publish(index, **values)


@pytest.fixture
def files(tmp_path):
    path = tmp_path / 'source'; path.write_bytes(b'source bytes')
    return {'source': path}


def large_files(files):
    # Long labels reproduce real full-build provenance, without thousands of files.
    return {f'dependency-{i:05}-' + 'x'*220: files['source'] for i in range(5000)}


def test_large_identity_roundtrip_keeps_header_below_legacy_limit(tmp_path, files):
    inputs = large_files(files)
    store = open_store(tmp_path, inputs)
    first = publish(store)
    raw = first.path.read_bytes()
    _, header_size, _ = struct.unpack('<8sQQ', raw[:24])
    assert header_size < 1024*1024
    assert open_store(tmp_path, inputs).load_latest().payload == b'native-state\x00'
    files['source'].write_bytes(b'altered source')
    with pytest.raises(ValueError): store.load_latest()
    with pytest.raises(ValueError): open_store(tmp_path, inputs).load_latest()


def manifest_path(store):
    paths = list(store.directory.glob('identity-*.json'))
    assert len(paths) == 1
    return paths[0]


@pytest.mark.parametrize('damage', ['missing', 'changed', 'symlink', 'trailing'])
def test_manifest_damage_is_not_repaired_on_reopen(tmp_path, files, damage):
    store = open_store(tmp_path, files); publish(store)
    path = manifest_path(store); original = path.read_bytes()
    if damage == 'missing': path.unlink()
    elif damage == 'changed': path.write_bytes(original.replace(b'source bytes', b'other bytes') + b'x')
    elif damage == 'trailing': path.write_bytes(original + b'\n')
    else:
        outside = tmp_path / 'outside'; outside.write_bytes(original)
        path.unlink(); path.symlink_to(outside)
    with pytest.raises((ValueError, OSError)): open_store(tmp_path, files).load_latest()
    if damage == 'missing': assert not path.exists()


def test_missing_dependency_and_same_bytes_symlink_retarget_are_detected(tmp_path, files):
    alias = tmp_path / 'alias'; alias.symlink_to(files['source'])
    inputs = {'module': alias}; store = open_store(tmp_path, inputs); publish(store)
    alternate = tmp_path / 'alternate'; alternate.write_bytes(files['source'].read_bytes())
    alias.unlink(); alias.symlink_to(alternate)
    with pytest.raises(ValueError): store.load_latest()
    with pytest.raises(ValueError): open_store(tmp_path, inputs).load_latest()
    alternate.unlink()
    with pytest.raises(OSError): store.load_latest()


def test_oversize_identity_rejected_before_creating_store(tmp_path, files):
    with pytest.raises(ValueError, match='IDENTITY_MANIFEST_LIMIT'):
        open_store(tmp_path, large_files(files), maximum_identity_bytes=1024)
    assert not (tmp_path / 'store').exists()


def test_preflight_covers_future_artifacts_and_rejects_large_header(tmp_path, files):
    store = open_store(tmp_path, files)
    planned = {'trace': tmp_path / 'not-created-yet.trace'}
    info = store.preflight(artifact_paths=planned)
    assert info['header_upper_bound_bytes'] < 1024*1024
    assert info['identity_file_count'] == 1
    assert not planned['trace'].exists()
    assert store.load_latest() is None
    oversized = {'x'*1024*1024: tmp_path / 'future'}
    with pytest.raises(ValueError, match='CHECKPOINT_HEADER_LIMIT'):
        store.preflight(artifact_paths=oversized)
    assert not list(store.directory.glob('stage-*'))


def test_fresh_process_reads_same_payload_and_manifest_is_not_rewritten(tmp_path, files):
    store = open_store(tmp_path, files); first = publish(store)
    path = manifest_path(store); before = (path.read_bytes(), path.stat().st_ino, path.stat().st_mtime_ns)
    publish(store, 1)
    code = '''import json,sys
from make_it_stand.verification.source_checkpoint_store_v2 import StageCheckpointStore
s=StageCheckpointStore(sys.argv[1],json.loads(sys.argv[2]),maximum_wall_ns=100,maximum_trace_bytes=1000,maximum_payload_bytes=4096)
r=s.load_latest(); print(r.stage,r.payload.hex())
'''
    output = subprocess.check_output([sys.executable, '-c', code, str(store.directory),
        json.dumps({k:str(v) for k,v in files.items()})], text=True)
    assert output.strip() == '1 6e61746976652d737461746500'
    assert (path.read_bytes(), path.stat().st_ino, path.stat().st_mtime_ns) == before
    assert first.path.is_file()


@pytest.mark.parametrize('point', ['manifest-link', 'manifest-sync', 'stage-link', 'receipt-link'])
def test_interrupted_publication_never_yields_unbound_state(tmp_path, files, monkeypatch, point):
    store = open_store(tmp_path, files)
    link, sync = os.link, os.fsync
    manifest_linked = False
    def broken_link(src, dst):
        nonlocal manifest_linked
        name = Path(dst).name
        if point == 'manifest-link' and name.startswith('identity-'): raise OSError('injected')
        if point == 'stage-link' and name.startswith('stage-'): raise OSError('injected')
        if point == 'receipt-link' and name.startswith('commit-'): raise OSError('injected')
        link(src, dst)
        if name.startswith('identity-'): manifest_linked = True
    def broken_sync(fd):
        if point == 'manifest-sync' and manifest_linked: raise OSError('injected')
        sync(fd)
    with monkeypatch.context() as m:
        m.setattr(os, 'link', broken_link); m.setattr(os, 'fsync', broken_sync)
        with pytest.raises(OSError, match='injected'): publish(store)
    assert store.load_latest() is None
    assert list(store.directory.glob('.pending-*'))


def test_legacy_reader_remains_explicit_and_mixed_chain_is_rejected(tmp_path, files):
    from make_it_stand.verification.source_checkpoint_store import StageCheckpointStore as Legacy
    old = Legacy(tmp_path / 'old', files, maximum_wall_ns=100, maximum_trace_bytes=1000,
                 maximum_payload_bytes=4096)
    saved = publish(old); original = saved.path.read_bytes()
    modern = api.StageCheckpointStore(old.directory, files, maximum_wall_ns=100,
                 maximum_trace_bytes=1000, maximum_payload_bytes=4096)
    with pytest.raises(ValueError): modern.load_latest()
    with pytest.raises(ValueError): publish(modern, 1)
    assert saved.path.read_bytes() == original
    assert old.load_latest().sha256 == saved.sha256
    new = open_store(tmp_path, files); publish(new)
    reverse = Legacy(new.directory, files, maximum_wall_ns=100, maximum_trace_bytes=1000,
                     maximum_payload_bytes=4096)
    with pytest.raises(ValueError): reverse.load_latest()


def test_same_durability_contract_executes_against_v2(tmp_path):
    # Reuse the unchanged behavioral suite, including real subprocess interruption,
    # concurrent publishers, partial writes, counter rollback and corrupt frames.
    # Replace only its module import target, including imports in child programs.
    contract = Path(__file__).with_name('test_source_checkpoint_store.py').read_text()
    target = tmp_path / 'test_v2_durability_contract.py'
    target.write_text(contract.replace('make_it_stand.verification.source_checkpoint_store',
                                      'make_it_stand.verification.source_checkpoint_store_v2'))
    completed = subprocess.run([sys.executable, '-m', 'pytest', '-q', str(target),
        '-o', 'cache_dir='+str(tmp_path / 'cache')], capture_output=True, text=True, timeout=60)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert 'skipped' not in completed.stdout and 'deselected' not in completed.stdout


def test_identity_limit_boundary_has_no_implicit_cap_increase(tmp_path, files):
    store = open_store(tmp_path, files); publish(store)
    exact = len(manifest_path(store).read_bytes())
    equal = api.StageCheckpointStore(tmp_path / 'equal', files, maximum_wall_ns=100,
        maximum_trace_bytes=1000, maximum_payload_bytes=4096, maximum_identity_bytes=exact)
    assert publish(equal).payload == b'native-state\x00'
    with pytest.raises(ValueError, match='IDENTITY_MANIFEST_LIMIT'):
        api.StageCheckpointStore(tmp_path / 'less', files, maximum_wall_ns=100,
            maximum_trace_bytes=1000, maximum_identity_bytes=exact-1)
    assert not (tmp_path / 'less').exists()


@pytest.mark.parametrize('field,value', [('name','../outside'), ('size',1),
    ('file_count',0), ('sha256','0'*64)])
def test_rehashed_manifest_reference_forgery_is_rejected(tmp_path, files, field, value):
    store = open_store(tmp_path, files); saved = publish(store)
    raw = saved.path.read_bytes(); magic, size, payload_size = struct.unpack('<8sQQ', raw[:24])
    header = json.loads(raw[24:24+size]); header['identity_manifest'][field] = value
    encoded = json.dumps(header, sort_keys=True, separators=(',',':')).encode()
    body = struct.pack('<8sQQ', magic, len(encoded), payload_size)+encoded+raw[24+size:-64]
    digest = hashlib.sha256(body).hexdigest()
    saved.path.write_bytes(body+digest.encode())
    (store.directory/'commit-00000000.sha256').write_text(f'00000000 {digest}\n')
    with pytest.raises(ValueError): store.load_latest()
