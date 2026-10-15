"""Lossless sequential diagnostics, not paper/geometry acceptance evidence.

Export/read memory scales with the largest single record, not the whole trace.
The native result still retains all solver history. A reader is verified only
after exhaustion: yielded prefixes have NOT passed the terminal integrity check.
Checksums detect corruption; they do not authenticate an author or prove physics.
"""
import hashlib
import json
from pathlib import Path

DEFAULT_RECORD_LIMIT = 64 * 1024 * 1024


class _Sequence:
    def __init__(self):
        self.counts = None
        self.pass_index = -1
        self.trace_index = 0
        self.records = 0

    def accept(self, event):
        if not isinstance(event, dict) or not isinstance(event.get('payload'), dict):
            raise ValueError('INVALID_EVENT')
        kind = event.get('kind')
        if self.counts is None:
            counts = event.get('trace_counts')
            if (kind != 'header' or not isinstance(counts, list) or
                    any(type(n) is not int or n < 0 for n in counts)):
                raise ValueError('INVALID_HEADER')
            # The producer/reader owns the mutable event; admission state does not.
            self.counts = tuple(counts)
        elif kind == 'pass':
            if self.pass_index >= 0 and self.trace_index != self.counts[self.pass_index]:
                raise ValueError('INCOMPLETE_PASS')
            index = event.get('pass_index')
            if type(index) is not int or index != self.pass_index + 1 or index >= len(self.counts):
                raise ValueError('PASS_SEQUENCE')
            self.pass_index = index
            self.trace_index = 0
        elif kind == 'iteration':
            p, t = event.get('pass_index'), event.get('trace_index')
            if (type(p) is not int or type(t) is not int or self.pass_index < 0 or
                    p != self.pass_index or t != self.trace_index or t >= self.counts[p]):
                raise ValueError('ITERATION_SEQUENCE')
            self.trace_index += 1
        else:
            raise ValueError('UNEXPECTED_EVENT')
        self.records += 1

    def finish(self):
        if (self.counts is None or self.pass_index + 1 != len(self.counts) or
                (self.pass_index >= 0 and self.trace_index != self.counts[self.pass_index])):
            raise ValueError('INCOMPLETE_ARCHIVE')


def _encoded(event):
    return (json.dumps(event, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode('utf-8')


def _limit(value):
    if type(value) is not int or value <= 0:
        raise ValueError('INVALID_RECORD_LIMIT')
    return value


def write_archive(path, events, *, max_record_bytes=DEFAULT_RECORD_LIMIT):
    """Exclusive creation; failed writes remain partial, without a valid footer."""
    return _write_archive(path, events, max_record_bytes=max_record_bytes, encoded=False)


def write_encoded_archive(path, events, *, max_record_bytes=DEFAULT_RECORD_LIMIT):
    """Native JSON lines: validate each record, without reserializing its arrays.

    Uses the same v1 reader/sequence/footer contract. Bytes differ from the
    sorted-key dictionary writer; compare decoded scientific records, not the
    archive hash across encoders. No caller-supplied completion footer is used.
    """
    return _write_archive(path, events, max_record_bytes=max_record_bytes, encoded=True)


def write_trusted_native_archive(path, result, *, max_record_bytes=DEFAULT_RECORD_LIMIT):
    """Write the compiled native result's JSONL without reparsing every trace.

    This deliberately narrow fast path trusts only the byte iterator owned by
    ``NativeSourceSpinMaterialResult``.  It still checks the native event
    envelope/order/counts, enforces the same record limit, and writes the same
    v1 footer and SHA-256 digests.  The regular public writers remain strict
    JSON parsers for caller-provided records.  The native iterator itself reads
    and checksum-verifies each sealed SLQP binary record as it emits it.

    A complete independent JSON readback is still available through
    ``read_archive`` and is used in tests and offline audits; it is omitted from
    the time-critical worker because the supervisor subsequently hashes the
    durable file and the writer already computed its complete digest.
    """
    limit = _limit(max_record_bytes)
    try:
        from _make_it_stand import NativeSourceSpinMaterialResult
        if not isinstance(result, NativeSourceSpinMaterialResult):
            raise ValueError('NATIVE_RESULT_REQUIRED')
        storage = result.trace_storage
        pass_count = storage['resident_pass_count']
        iteration_count = storage['record_count']
        events = iter(result.iter_diagnostics_json())
    except ImportError as exc:
        raise ValueError('NATIVE_RESULT_REQUIRED') from exc
    except (AttributeError, KeyError, TypeError) as exc:
        raise ValueError('NATIVE_RESULT_REQUIRED') from exc
    if type(pass_count) is not int or pass_count < 0 or type(iteration_count) is not int or iteration_count < 0:
        raise ValueError('INVALID_NATIVE_TRACE_COUNTS')

    try:
        first = next(events)
    except StopIteration:
        raise ValueError('MISSING_HEADER') from None
    if (not isinstance(first, bytes) or len(first) > limit or
            first.count(b'\n') != 1 or not first.endswith(b'\n') or
            not first.startswith(b'{"kind":"header",') or not first.endswith(b'}\n')):
        raise ValueError('INVALID_NATIVE_HEADER_LINE')
    header = json.loads(first, object_pairs_hook=_pairs, parse_constant=_nonfinite, parse_float=_float)
    if not isinstance(header, dict) or header.get('kind') != 'header' or not isinstance(header.get('payload'), dict):
        raise ValueError('INVALID_NATIVE_HEADER')
    counts = header.get('trace_counts')
    if (not isinstance(counts, list) or len(counts) != pass_count or
            any(type(count) is not int or count < 0 for count in counts) or
            sum(counts) != iteration_count):
        raise ValueError('NATIVE_TRACE_COUNT_MISMATCH')

    def field(data, offset, suffix, error):
        end = data.find(suffix, offset)
        if end < 0:
            raise ValueError(error)
        digits = data[offset:end]
        if not digits or not digits.isdigit() or (len(digits) > 1 and digits.startswith(b'0')):
            raise ValueError(error)
        return int(digits), end + len(suffix)

    def validate_line(data, pass_index, trace_index):
        if not isinstance(data, bytes):
            raise TypeError('NATIVE_EVENT_REQUIRES_BYTES')
        if len(data) > limit:
            raise ValueError('RECORD_LIMIT')
        if data.count(b'\n') != 1 or not data.endswith(b'\n') or not data.endswith(b'}\n'):
            raise ValueError('INVALID_NATIVE_JSON_LINE')
        if data.startswith(b'{"kind":"pass","pass_index":'):
            if (pass_index >= 0 and trace_index + 1 != counts[pass_index]) or (
                    pass_index < 0 and trace_index != -1):
                raise ValueError('INCOMPLETE_NATIVE_PASS')
            actual, end = field(data, len(b'{"kind":"pass","pass_index":'), b',"payload":', 'PASS_SEQUENCE')
            if actual != pass_index + 1 or actual >= len(counts):
                raise ValueError('PASS_SEQUENCE')
            return actual, -1
        prefix = b'{"kind":"iteration","pass_index":'
        if data.startswith(prefix):
            actual_pass, offset = field(data, len(prefix), b',"trace_index":', 'ITERATION_SEQUENCE')
            actual_trace, _ = field(data, offset, b',"payload":', 'ITERATION_SEQUENCE')
            if (pass_index < 0 or actual_pass != pass_index or actual_trace != trace_index + 1 or
                    actual_trace >= counts[pass_index]):
                raise ValueError('ITERATION_SEQUENCE')
            return pass_index, actual_trace
        raise ValueError('UNEXPECTED_NATIVE_EVENT')

    sequence_pass = -1
    sequence_trace = -1
    records = 0
    digest = hashlib.sha256()
    whole = hashlib.sha256()
    with Path(path).open('xb') as out:
        def emit(data):
            nonlocal sequence_pass, sequence_trace, records
            if records == 0:
                # The first record was parsed completely above; bind its exact
                # bytes just as the ordinary writer does.
                pass
            else:
                sequence_pass, sequence_trace = validate_line(data, sequence_pass, sequence_trace)
            out.write(data)
            digest.update(data)
            whole.update(data)
            records += 1

        emit(first)
        for data in events:
            emit(data)
        if records != 1 + pass_count + iteration_count:
            raise ValueError('NATIVE_EVENT_COUNT_MISMATCH')
        if (pass_count == 0 and iteration_count != 0) or (
                pass_count > 0 and (sequence_pass != pass_count - 1 or
                                    sequence_trace != counts[-1] - 1)):
            raise ValueError('INCOMPLETE_NATIVE_ARCHIVE')
        footer = dict(kind='complete', schema='source-trace-archive-v1',
                      records=records, payload_sha256=digest.hexdigest())
        encoded_footer = _encoded(footer)
        if len(encoded_footer) > limit:
            raise ValueError('RECORD_LIMIT')
        out.write(encoded_footer)
        whole.update(encoded_footer)
    return header, dict(records=records, payload_sha256=digest.hexdigest(),
                        archive_sha256=whole.hexdigest())


def _write_archive(path, events, *, max_record_bytes, encoded):
    limit = _limit(max_record_bytes)
    sequence = _Sequence()
    digest = hashlib.sha256()
    whole = hashlib.sha256()
    with Path(path).open('xb') as out:
        for event in events:
            if encoded:
                data = event
                if not isinstance(data, bytes):
                    raise TypeError('ENCODED_EVENT_REQUIRES_BYTES')
                if len(data) > limit:
                    raise ValueError('RECORD_LIMIT')
                if not data.endswith(b'\n') or data.count(b'\n') != 1:
                    raise ValueError('INVALID_JSON_LINE')
                event = json.loads(data, object_pairs_hook=_pairs, parse_constant=_nonfinite, parse_float=_float)
            else:
                data = _encoded(event)
            sequence.accept(event)
            if len(data) > limit:
                raise ValueError('RECORD_LIMIT')
            out.write(data)
            digest.update(data)
            whole.update(data)
        sequence.finish()
        footer = dict(kind='complete', schema='source-trace-archive-v1',
                      records=sequence.records, payload_sha256=digest.hexdigest())
        data = _encoded(footer)
        if len(data) > limit:
            raise ValueError('RECORD_LIMIT')
        out.write(data)
        whole.update(data)
    return dict(records=sequence.records, payload_sha256=digest.hexdigest(),
                archive_sha256=whole.hexdigest())


def _pairs(pairs):
    out = {}
    for key, value in pairs:
        if key in out:
            raise ValueError('DUPLICATE_JSON_KEY')
        out[key] = value
    return out


def _nonfinite(value):
    raise ValueError('NONFINITE_JSON')


def _float(value):
    import math
    number = float(value)
    if not math.isfinite(number):
        raise ValueError('NONFINITE_JSON')
    return number


def read_archive(path, *, max_record_bytes=DEFAULT_RECORD_LIMIT):
    """Yield events; only full exhaustion verifies sequence, EOF and checksum."""
    limit = _limit(max_record_bytes)
    sequence = _Sequence()
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while True:
            line = stream.readline(limit + 1)
            if not line:
                raise ValueError('MISSING_COMPLETION_FOOTER')
            if len(line) > limit:
                raise ValueError('RECORD_LIMIT')
            if not line.endswith(b'\n'):
                raise ValueError('TRUNCATED_RECORD')
            event = json.loads(line, object_pairs_hook=_pairs, parse_constant=_nonfinite, parse_float=_float)
            if not isinstance(event, dict):
                raise ValueError('INVALID_EVENT')
            if event.get('kind') == 'complete':
                sequence.finish()
                if (event.get('schema') != 'source-trace-archive-v1' or
                        type(event.get('records')) is not int or event['records'] != sequence.records or
                        event.get('payload_sha256') != digest.hexdigest()):
                    raise ValueError('ARCHIVE_INTEGRITY_MISMATCH')
                if stream.read(1):
                    raise ValueError('TRAILING_RECORDS')
                return
            sequence.accept(event)
            digest.update(line)
            yield event


def read_archive_envelope(path, *, max_record_bytes=DEFAULT_RECORD_LIMIT):
    """Read the header and terminal footer only; does not validate payload bytes.

    Callers may use this only when they also bind the complete file SHA-256 to
    the digest returned by ``write_trusted_native_archive``.  For untrusted or
    caller-produced archives, exhaust ``read_archive`` instead.
    """
    limit = _limit(max_record_bytes)
    path = Path(path)
    with path.open('rb') as stream:
        first = stream.readline(limit + 1)
        if not first:
            raise ValueError('MISSING_HEADER')
        if len(first) > limit:
            raise ValueError('RECORD_LIMIT')
        if not first.endswith(b'\n'):
            raise ValueError('TRUNCATED_RECORD')
        header = json.loads(first, object_pairs_hook=_pairs,
                            parse_constant=_nonfinite, parse_float=_float)
        if not isinstance(header, dict) or header.get('kind') != 'header':
            raise ValueError('MISSING_HEADER')
        stream.seek(0, 2)
        size = stream.tell()
        if size < len(first) + 1:
            raise ValueError('MISSING_COMPLETION_FOOTER')
        tail_size = min(size, limit + 1)
        stream.seek(-tail_size, 2)
        tail = stream.read(tail_size)
    lines = tail.splitlines(keepends=True)
    if not lines or len(lines[-1]) > limit or not lines[-1].endswith(b'\n'):
        raise ValueError('MISSING_COMPLETION_FOOTER')
    footer = json.loads(lines[-1], object_pairs_hook=_pairs,
                        parse_constant=_nonfinite, parse_float=_float)
    if not isinstance(footer, dict) or footer.get('kind') != 'complete':
        raise ValueError('MISSING_COMPLETION_FOOTER')
    return header, footer
