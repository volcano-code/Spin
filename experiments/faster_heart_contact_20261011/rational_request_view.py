"""Streaming coordinate-only gate; never a geometry or solver certificate.

This attachment cannot make a lossy request eligible. It identifies where the
existing binary64 request interface needs a separately certified derivative
representation or a genuine rational-input implementation.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            result.update(chunk)
    return result.hexdigest()


class _Cursor:
    def __init__(self, stream, chunk_size):
        self.stream, self.chunk_size = stream, chunk_size
        self.buffer = ''; self.end = False
        self.decoder = json.JSONDecoder()

    def fill(self):
        part = self.stream.read(self.chunk_size)
        self.end = not part
        self.buffer += part

    def peek(self):
        self.buffer = self.buffer.lstrip()
        while not self.buffer and not self.end:
            self.fill(); self.buffer = self.buffer.lstrip()
        return self.buffer[:1]

    def take(self, char):
        if self.peek() != char: raise ValueError('MALFORMED_BOUNDARY_JSON')
        self.buffer = self.buffer[1:]

    def value(self):
        self.peek()
        while True:
            try:
                value, end = self.decoder.raw_decode(self.buffer)
                # Values here are complete strings/rows, never bare numbers.
                self.buffer = self.buffer[end:]
                return value
            except json.JSONDecodeError as exc:
                if self.end: raise ValueError('TRUNCATED_BOUNDARY_JSON') from exc
                if len(self.buffer) > (1 << 20):
                    raise ValueError('BOUNDARY_ROW_EXCEEDS_STREAMING_LIMIT') from exc
                self.fill()


def boundary_rows(path, chunk_size=65536):
    """Yield one vertex/face at a time; validate the complete JSON document."""
    if type(chunk_size) is not int or chunk_size <= 0:
        raise ValueError('invalid streaming chunk size')
    with Path(path).open() as stream:
        cursor = _Cursor(stream, chunk_size); cursor.take('{'); seen = set()
        while True:
            key = cursor.value()
            if not isinstance(key, str) or key in seen:
                raise ValueError('INVALID_OR_DUPLICATE_BOUNDARY_KEY')
            seen.add(key); cursor.take(':')
            if key in ('vertices', 'faces'):
                cursor.take('[')
                if cursor.peek() != ']':
                    while True:
                        yield key, cursor.value()
                        if cursor.peek() == ']': break
                        cursor.take(',')
                cursor.take(']')
            elif key == 'schema':
                if cursor.value() != 'authoritative-rational-input-boundary-v1':
                    raise ValueError('UNSUPPORTED_BOUNDARY_SCHEMA')
            elif key == 'obj_role':
                if not isinstance(cursor.value(), str):
                    raise ValueError('INVALID_DISPLAY_ROLE_METADATA')
            else:
                raise ValueError('UNSUPPORTED_BOUNDARY_FIELD: '+key)
            if cursor.peek() == '}': break
            cursor.take(',')
        cursor.take('}')
        if cursor.peek() or not {'schema', 'vertices', 'faces'} <= seen:
            raise ValueError('INCOMPLETE_OR_TRAILING_BOUNDARY_JSON')


def check_region(path, expected_sha256, *, chunk_size=65536):
    path = Path(path)
    if sha(path) != expected_sha256: raise ValueError('SOURCE_IDENTITY_MISMATCH')
    vertices = coordinates = faces = lossy = unavailable = 0
    first = None
    for kind, row in boundary_rows(path, chunk_size):
        if not isinstance(row, list) or len(row) != 3:
            raise ValueError('INVALID_VERTEX_OR_FACE_ROW')
        if kind == 'faces':
            if any(type(i) is not int or i < 0 for i in row):
                raise ValueError('INVALID_FACE_ROW')
            faces += 1; continue
        for axis, text in enumerate(row):
            if not isinstance(text, str): raise ValueError('NONCANONICAL_RATIONAL_COORDINATE')
            try: exact = Fraction(text)
            except (ValueError, ZeroDivisionError) as exc:
                raise ValueError('NONCANONICAL_RATIONAL_COORDINATE') from exc
            if text != f'{exact.numerator}/{exact.denominator}':
                raise ValueError('NONCANONICAL_RATIONAL_COORDINATE')
            try:
                floating = float(exact)
                represented = Fraction.from_float(floating) if math.isfinite(floating) else None
            except OverflowError:
                represented = None
            coordinates += 1
            if represented != exact:
                lossy += 1
                unavailable += represented is None
                if first is None:
                    first = dict(vertex=vertices, axis=axis, authority=text,
                        binary64_exact=None if represented is None else
                        f'{represented.numerator}/{represented.denominator}')
        vertices += 1
    if sha(path) != expected_sha256: raise ValueError('SOURCE_IDENTITY_CHANGED_DURING_CHECK')
    return dict(source_sha256=expected_sha256, vertex_count=vertices, face_count=faces,
        coordinate_count=coordinates, non_binary64_coordinates=lossy,
        unrepresentable_coordinates=unavailable, first_lossy_coordinate=first,
        status='NEEDS_BOUNDED_DERIVED_REQUEST_CERTIFICATION' if lossy else
               'LOSSLESS_BINARY64_COORDINATE_VIEW_ONLY',
        coordinate_view_lossless=lossy == 0,
        geometry_admission_granted=False, solver_admission_granted=False)


def inspect_admission(admission):
    admission = Path(admission)
    record_path = admission/'admission.json'; record_sha = sha(record_path)
    record = json.loads(record_path.read_text())
    regions = {name:check_region(admission/(name+'.rational.json'),
                               record['files'][name+'.rational.json'])
               for name in ('outer', 'fixed', 'design')}
    if sha(record_path) != record_sha: raise ValueError('ADMISSION_IDENTITY_CHANGED')
    lossy = any(not r['coordinate_view_lossless'] for r in regions.values())
    return dict(schema='rational-request-view-check-v1',
        checker_sha256=sha(__file__), admission_sha256=record_sha,
        problem_digest=record['problem_digest'], parent_admission_status=record['status'],
        regions=regions, status='NEEDS_BOUNDED_DERIVED_REQUEST_CERTIFICATION' if lossy else
        'LOSSLESS_BINARY64_COORDINATE_VIEW_ONLY',
        scope='COORDINATE_REPRESENTABILITY_ONLY_NOT_GEOMETRY_OR_SOLVER_ADMISSION',
        solver_admission_granted=False, geometry_admission_granted=False,
        optimization_invocations=0, simulation_invocations=0)


def require_lossless_admission_view(admission):
    """Block authority reuse; a return value still grants no solver admission."""
    report = inspect_admission(admission)
    if report['status'] == 'NEEDS_BOUNDED_DERIVED_REQUEST_CERTIFICATION':
        raise ValueError(report['status'])
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--admission', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = inspect_admission(args.admission)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.flush(); os.fsync(stream.fileno())
    directory = os.open(args.output.parent, os.O_RDONLY | os.O_DIRECTORY)
    try: os.fsync(directory)
    finally: os.close(directory)
    print(json.dumps(dict(status=result['status'], solver_admission_granted=False,
        counts={name:r['non_binary64_coordinates'] for name,r in result['regions'].items()})))


if __name__ == '__main__': main()
