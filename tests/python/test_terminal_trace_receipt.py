"""Tests for compact terminal receipts over native SLQP trace artifacts."""
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path

import pytest
from test_stage_audit import failed_pass_fixture


def receipt_api():
    assert importlib.util.find_spec('terminal_trace_receipt') is not None, 'native trace receipt module missing'
    return importlib.import_module('terminal_trace_receipt')


def fixture(tmp_path):
    header, stage, terminal = failed_pass_fixture()
    header_event = {'kind': 'header', 'trace_counts': [1], 'payload': header}
    events = [header_event, {'kind': 'pass', 'pass_index': 0, 'payload': stage},
              {'kind': 'iteration', 'pass_index': 0, 'trace_index': 0, 'payload': terminal}]
    encode = lambda value: (json.dumps(value, sort_keys=True, separators=(',', ':'),
                                        allow_nan=False) + '\n').encode()
    trace_root = tmp_path/'repeat-0'; trace = trace_root/'attempt-00000000/payload/trace/pass-000000.slqptrace'
    trace.parent.mkdir(parents=True); trace.write_bytes(b'synthetic sealed trace bytes')
    output = trace_root/'attempt-00000000/payload'; output.mkdir(parents=True, exist_ok=True)
    (output/'header.json').write_text(json.dumps(header_event, sort_keys=True, indent=2) + '\n')
    storage = {'mode': 'DISK_ITERATIONS_V1', 'resident_pass_count': 1,
               'record_count': 1, 'bytes_written': trace.stat().st_size,
               'files': [str(trace)]}
    return output, trace_root, trace, header_event, iter(encode(event) for event in events), storage


def test_receipt_binds_header_audit_and_each_native_trace(tmp_path):
    api = receipt_api()
    output, trace_root, trace, header_event, events, storage = fixture(tmp_path)
    receipt = api.write_terminal_trace_receipt(
        output, header_event, events, storage, trace_root=trace_root)
    checked = api.verify_terminal_trace_receipt(output/'trace-receipt.json')
    assert receipt['schema'] == 'spin2014-native-terminal-trace-receipt-v1'
    assert receipt['trace_counts'] == [1]
    assert receipt['trace_files'] == [dict(pass_index=0, path=str(trace.resolve()),
        size=trace.stat().st_size, sha256=hashlib.sha256(trace.read_bytes()).hexdigest())]
    assert checked['status'] == 'VERIFIED'
    assert checked['audit_sha256'] == receipt['audit_sha256']
    assert json.loads((output/'stage-audit.json').read_text())['scientific_status'] == \
        header_event['payload']['solve']['status']


def test_receipt_verifier_rejects_changed_or_missing_native_trace(tmp_path):
    api = receipt_api()
    output, trace_root, trace, header_event, events, storage = fixture(tmp_path)
    api.write_terminal_trace_receipt(output, header_event, events, storage, trace_root=trace_root)
    trace.write_bytes(b'changed after receipt')
    with pytest.raises(ValueError, match='TRACE_FILE_HASH_MISMATCH'):
        api.verify_terminal_trace_receipt(output/'trace-receipt.json')


def test_receipt_verifier_rejects_missing_native_trace(tmp_path):
    api = receipt_api()
    output, trace_root, trace, header_event, events, storage = fixture(tmp_path)
    api.write_terminal_trace_receipt(output, header_event, events, storage, trace_root=trace_root)
    trace.unlink()
    with pytest.raises(ValueError, match='TRACE_FILE_UNAVAILABLE'):
        api.verify_terminal_trace_receipt(output/'trace-receipt.json')


@pytest.mark.parametrize('name,reason', [
    ('header.json', 'TERMINAL_HEADER_HASH_MISMATCH'),
    ('stage-audit.json', 'STAGE_AUDIT_HASH_MISMATCH'),
])
def test_receipt_verifier_rejects_changed_header_or_audit(tmp_path, name, reason):
    api = receipt_api()
    output, trace_root, _, header_event, events, storage = fixture(tmp_path)
    api.write_terminal_trace_receipt(output, header_event, events, storage, trace_root=trace_root)
    path = output/name
    path.write_bytes(path.read_bytes() + b' ')
    with pytest.raises(ValueError, match=reason):
        api.verify_terminal_trace_receipt(output/'trace-receipt.json')


def test_receipt_writer_rejects_duplicate_trace_paths(tmp_path):
    api = receipt_api()
    output, trace_root, trace, header_event, events, storage = fixture(tmp_path)
    storage.update(resident_pass_count=2, files=[str(trace), str(trace)], bytes_written=2*trace.stat().st_size)
    header_event['trace_counts'] = [1, 0]
    (output/'header.json').write_text(json.dumps(header_event, sort_keys=True, indent=2) + '\n')
    with pytest.raises(ValueError, match='DUPLICATE_TRACE_FILE_PATH'):
        api.write_terminal_trace_receipt(output, header_event, iter(()), storage, trace_root=trace_root)


def test_receipt_writer_rejects_trace_count_mismatch_before_publishing(tmp_path):
    api = receipt_api()
    output, trace_root, _, header_event, events, storage = fixture(tmp_path)
    storage['files'] = []
    with pytest.raises(ValueError, match='TRACE_FILE_COUNT_MISMATCH'):
        api.write_terminal_trace_receipt(output, header_event, events,
                                         storage, trace_root=trace_root)
    assert not (output/'trace-receipt.json').exists()
    assert not (output/'stage-audit.json').exists()
