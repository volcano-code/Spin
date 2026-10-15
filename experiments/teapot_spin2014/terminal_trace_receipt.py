"""Durable receipt for direct audits of native SLQP trace files.

The receipt replaces the expanded JSONL transport only for the staged G7
campaign. It binds the audit and header to every native trace file; it does not
grant solver, geometry, or paper acceptance.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import stat
import uuid

from stage_audit_v2 import audit_native_diagnostics


SCHEMA = 'spin2014-native-terminal-trace-receipt-v1'


def _sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _json_unique(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('DUPLICATE_RECEIPT_KEY')
        result[key] = value
    return result


def _reject_constant(_):
    raise ValueError('NONFINITE_RECEIPT')


def _finite_float(text):
    value = float(text)
    if not math.isfinite(value):
        raise ValueError('NONFINITE_RECEIPT')
    return value


def _decode_event(data):
    if (not isinstance(data, bytes) or not data.endswith(b'\n') or
            data.count(b'\n') != 1):
        raise ValueError('TRUNCATED_NATIVE_EVENT')
    return json.loads(data, object_pairs_hook=_json_unique,
                      parse_constant=_reject_constant, parse_float=_finite_float)


def _write_durable_new(path, value):
    path = Path(path)
    if not path.parent.is_dir():
        raise ValueError('RECEIPT_DIRECTORY_MISSING')
    encoded = (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()
    temporary = path.parent / ('.pending-' + uuid.uuid4().hex)
    with temporary.open('xb') as stream:
        stream.write(encoded)
        stream.flush()
        os.fsync(stream.fileno())
    try:
        os.link(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        temporary.unlink(missing_ok=True)
    return hashlib.sha256(encoded).hexdigest()


def _regular_snapshot(path, root):
    path = Path(path)
    if not path.is_absolute():
        raise ValueError('TRACE_FILE_OUTSIDE_RUN_ROOT')
    try:
        resolved = path.resolve(strict=True)
    except OSError as exc:
        raise ValueError('TRACE_FILE_UNAVAILABLE') from exc
    if path != resolved or not path.is_relative_to(root):
        raise ValueError('TRACE_FILE_OUTSIDE_RUN_ROOT')
    try:
        info = path.lstat()
    except OSError as exc:
        raise ValueError('TRACE_FILE_UNAVAILABLE') from exc
    if not stat.S_ISREG(info.st_mode):
        raise ValueError('TRACE_FILE_NOT_REGULAR')
    return path, (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


def _file_receipts(paths, trace_root):
    root = Path(trace_root).resolve(strict=True)
    receipts = []
    total = 0
    snapshots = []
    seen = set()
    for index, raw_path in enumerate(paths):
        path, snapshot = _regular_snapshot(raw_path, root)
        if path in seen:
            raise ValueError('DUPLICATE_TRACE_FILE_PATH')
        seen.add(path)
        snapshots.append((path, snapshot))
        receipts.append({'pass_index': index, 'path': str(path), 'size': snapshot[2]})
        total += snapshot[2]
    return root, receipts, total, snapshots


def _seal_file_hashes(receipts, snapshots, trace_root):
    for receipt, (path, before) in zip(receipts, snapshots):
        hash_before = _regular_snapshot(path, trace_root)[1]
        actual = _sha256(path)
        after = _regular_snapshot(path, trace_root)[1]
        # The full root-containment check happened before audit; here compare
        # the no-follow metadata snapshot surrounding the native trace read.
        if before != hash_before or before != after:
            raise ValueError('TRACE_FILE_CHANGED_DURING_AUDIT')
        receipt['sha256'] = actual


def write_terminal_trace_receipt(output_directory, header_event, encoded_events,
                                 trace_storage, *, trace_root):
    output = Path(output_directory).resolve(strict=True)
    root = Path(trace_root).resolve(strict=True)
    if not output.is_relative_to(root):
        raise ValueError('OUTPUT_OUTSIDE_TRACE_ROOT')
    if not isinstance(header_event, dict) or header_event.get('kind') != 'header':
        raise ValueError('INVALID_TERMINAL_HEADER')
    header_path = output / 'header.json'
    audit_path = output / 'stage-audit.json'
    receipt_path = output / 'trace-receipt.json'
    try:
        header_path, _ = _regular_snapshot(header_path, root)
    except ValueError as exc:
        raise ValueError('TERMINAL_HEADER_FILE_REQUIRED') from exc
    try:
        saved_header = json.loads(header_path.read_text(), object_pairs_hook=_json_unique)
    except (OSError, ValueError) as exc:
        raise ValueError('TERMINAL_HEADER_INVALID') from exc
    if saved_header != header_event:
        raise ValueError('TERMINAL_HEADER_MISMATCH')

    counts = header_event.get('trace_counts')
    solve = header_event.get('payload', {}).get('solve', {})
    if (not isinstance(counts, list) or any(type(value) is not int or value < 0 for value in counts) or
            trace_storage.get('mode') != 'DISK_ITERATIONS_V1' or
            type(trace_storage.get('resident_pass_count')) is not int or
            trace_storage.get('resident_pass_count') != len(counts) or
            type(trace_storage.get('record_count')) is not int or
            trace_storage.get('record_count') != sum(counts)):
        raise ValueError('NATIVE_TRACE_COUNT_MISMATCH')
    paths = trace_storage.get('files')
    if not isinstance(paths, list) or len(paths) != len(counts):
        raise ValueError('TRACE_FILE_COUNT_MISMATCH')
    root, file_receipts, total_bytes, snapshots = _file_receipts(paths, root)
    if type(trace_storage.get('bytes_written')) is not int or trace_storage.get('bytes_written') != total_bytes:
        raise ValueError('TRACE_FILE_BYTES_MISMATCH')
    header_digest = _sha256(header_path)

    # The native iterator checks the sealed binary record frames as it yields
    # each line; this independent Python checker validates event JSON and order.
    source_events = iter(encoded_events)
    try:
        encoded_header = next(source_events)
    except StopIteration:
        raise ValueError('MISSING_HEADER') from None
    if _decode_event(encoded_header) != header_event:
        raise ValueError('NATIVE_HEADER_MISMATCH')

    def bound_events():
        yield encoded_header
        yield from source_events

    audit = audit_native_diagnostics(
        bound_events(), expected_pass_count=len(counts),
        expected_record_count=sum(counts))
    if audit.get('native_trace_counts') != counts:
        raise ValueError('AUDIT_TRACE_COUNT_MISMATCH')

    _seal_file_hashes(file_receipts, snapshots, root)
    if _sha256(header_path) != header_digest:
        raise ValueError('TERMINAL_HEADER_CHANGED')
    audit_sha = _write_durable_new(audit_path, audit)
    solve_hashes = {name: solve.get(name) for name in ('input_hash', 'trace_hash', 'output_hash')}
    if any(not isinstance(value, str) or len(value) != 64 or
           any(char not in '0123456789abcdef' for char in value) for value in solve_hashes.values()):
        raise ValueError('TERMINAL_SCIENTIFIC_HASHES_REQUIRED')
    receipt = {
        'schema': SCHEMA,
        'trace_root': str(root),
        'header_path': str(header_path),
        'header_sha256': header_digest,
        'trace_counts': counts,
        'pass_count': len(counts),
        'record_count': sum(counts),
        'native_trace_bytes': total_bytes,
        **solve_hashes,
        'trace_files': file_receipts,
        'audit_path': str(audit_path),
        'audit_sha256': audit_sha,
        'paper_acceptance': 'NOT_GRANTED',
        'sealed_artifact': 'NOT_GRANTED',
    }
    _write_durable_new(receipt_path, receipt)
    return receipt


def verify_terminal_trace_receipt(receipt_path):
    receipt_path = Path(receipt_path).absolute()
    receipt_root = receipt_path.parent.resolve(strict=True)
    receipt_path, _ = _regular_snapshot(receipt_path, receipt_root)
    try:
        receipt = json.loads(receipt_path.read_text(), object_pairs_hook=_json_unique,
                             parse_constant=_reject_constant, parse_float=_finite_float)
    except (OSError, ValueError) as exc:
        raise ValueError('TRACE_RECEIPT_INVALID') from exc
    if receipt.get('schema') != SCHEMA:
        raise ValueError('TRACE_RECEIPT_SCHEMA')
    root = Path(receipt.get('trace_root', '')).resolve(strict=True)
    header_path, _ = _regular_snapshot(receipt.get('header_path', ''), root)
    if _sha256(header_path) != receipt.get('header_sha256'):
        raise ValueError('TERMINAL_HEADER_HASH_MISMATCH')
    header = json.loads(header_path.read_text(), object_pairs_hook=_json_unique,
                        parse_constant=_reject_constant, parse_float=_finite_float)
    header_counts = header.get('trace_counts')
    if (header.get('kind') != 'header' or not isinstance(header_counts, list) or
            any(type(value) is not int or value < 0 for value in header_counts) or
            header_counts != receipt.get('trace_counts')):
        raise ValueError('TERMINAL_HEADER_MISMATCH')
    counts = receipt.get('trace_counts')
    if (not isinstance(counts, list) or any(type(value) is not int or value < 0 for value in counts) or
            receipt.get('pass_count') != len(counts) or
            receipt.get('record_count') != sum(counts)):
        raise ValueError('TRACE_RECEIPT_COUNT_MISMATCH')
    solve = header.get('payload', {}).get('solve', {})
    for name in ('input_hash', 'trace_hash', 'output_hash'):
        value = receipt.get(name)
        if (not isinstance(value, str) or len(value) != 64 or
                any(char not in '0123456789abcdef' for char in value) or solve.get(name) != value):
            raise ValueError('TERMINAL_SCIENTIFIC_HASH_MISMATCH')
    files = receipt.get('trace_files')
    if not isinstance(files, list) or len(files) != len(counts):
        raise ValueError('TRACE_FILE_COUNT_MISMATCH')
    total = 0
    for index, entry in enumerate(files):
        if not isinstance(entry, dict) or entry.get('pass_index') != index:
            raise ValueError('TRACE_FILE_SEQUENCE_MISMATCH')
        path, snapshot = _regular_snapshot(entry.get('path', ''), root)
        if (entry.get('size') != snapshot[2] or
                _sha256(path) != entry.get('sha256')):
            raise ValueError('TRACE_FILE_HASH_MISMATCH')
        total += snapshot[2]
    if total != receipt.get('native_trace_bytes'):
        raise ValueError('TRACE_FILE_BYTES_MISMATCH')
    audit_path, _ = _regular_snapshot(receipt.get('audit_path', ''), receipt_root)
    if _sha256(audit_path) != receipt.get('audit_sha256'):
        raise ValueError('STAGE_AUDIT_HASH_MISMATCH')
    audit = json.loads(audit_path.read_text(), object_pairs_hook=_json_unique,
                       parse_constant=_reject_constant, parse_float=_finite_float)
    if (audit.get('schema') != 'spin2014-stage-audit-v2' or
            audit.get('archive_integrity') != 'VERIFIED_NATIVE_TRACE_SEQUENCE_AND_RECORD_CHECKSUMS' or
            audit.get('native_trace_counts') != counts):
        raise ValueError('STAGE_AUDIT_RECEIPT_MISMATCH')
    return dict(status='VERIFIED', receipt_sha256=_sha256(receipt_path),
                audit_sha256=receipt['audit_sha256'], header_sha256=receipt['header_sha256'],
                trace_file_count=len(files), trace_bytes=total)
