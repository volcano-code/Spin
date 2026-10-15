"""Read-only receipt/derivative inventory; never grants a scientific gate.

Stream the large aggregation rows individually, and authenticate point packet
bytes without loading their geometry. No optimizer or Boolean calls occur.
"""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def aggregation_rows(path):
    """Read the existing sorted-key terminal format one JSON row at a time."""
    decoder = json.JSONDecoder()
    with Path(path).open() as stream:
        buffer = stream.read(65536)
        marker = '"rows": ['
        if marker not in buffer:
            raise ValueError('aggregation rows header absent')
        buffer = buffer.split(marker, 1)[1]
        while True:
            buffer = buffer.lstrip()
            if buffer.startswith(']'):
                return
            if buffer.startswith(','):
                buffer = buffer[1:].lstrip()
            while True:
                try:
                    row, end = decoder.raw_decode(buffer)
                    break
                except json.JSONDecodeError:
                    chunk = stream.read(65536)
                    if not chunk:
                        raise ValueError('truncated aggregation row')
                    buffer += chunk
            yield row
            buffer = buffer[end:]


def inventory(evidence):
    evidence = Path(evidence)
    manifest_path = evidence / 'split-production-manifest-04/manifest.json'
    manifest = json.loads(manifest_path.read_text())
    dependencies = {}
    for name, path in manifest['mapping_settings']['source_paths'].items():
        dependencies[name] = Path(path).is_file() and digest(path) == manifest['sources'][name]
    points = []
    for index in range(81):
        root = evidence / f'split-production-point-04-{index:03d}'
        record = dict(index=index, status='MISSING')
        try:
            terminal = json.loads((root / 'output/terminal.json').read_text())
            result = json.loads((root / 'result.json').read_text())
            complete = json.loads((root / 'complete.json').read_text())
            record.update(status=terminal['status'], elapsed_seconds=terminal['elapsed_seconds'],
                manifest_match=terminal['manifest_sha256'] == manifest['identity'],
                packet_authenticated=digest(root / 'output/point.json') == terminal['packet_sha256'],
                independent_authenticated=digest(root / 'output/independent.json') == terminal['independent_sha256'],
                process_success=result['status'] == 'PROCESS_EXITED' and result['returncode'] == 0,
                result_authenticated=digest(root / 'result.json') == complete['result_sha256'],
                resource_receipt=result['service']['MemoryMax'] == str(24 * 1024**3)
                    and result['service']['MemorySwapMax'] == '0'
                    and result['service']['MemorySwapPeak'] == '0'
                    and result['elapsed_ns'] < 3600_000_000_000)
        except (OSError, KeyError, TypeError, json.JSONDecodeError) as error:
            record.update(status='INCOMPLETE_OR_INVALID', error=str(error))
        points.append(record)
    aggregation = evidence / 'split-production-aggregation-04/output/terminal.json'
    rows = []
    if aggregation.is_file():
        for row in aggregation_rows(aggregation):
            rows.append({name: row.get(name) for name in
                ('column', 'factor', 'indices', 'fixed_branch', 'within_tolerance',
                 'maximum_component_error', 'passed', 'status')})
    review_root = Path(__file__).resolve().parents[2] / 'experiments/faster_heart_contact_20261011/rounded-asymmetric-v2/production-review-v1'
    review = dict(status='NOT_FOUND')
    review_terminal = review_root / 'output/terminal.json'
    if review_terminal.is_file():
        with review_terminal.open() as stream:
            prefix = stream.read(65536).split('"rows": [', 1)[0]
        header = json.loads(prefix.rstrip().rstrip(',') + '}')
        with review_terminal.open('rb') as stream:
            stream.seek(max(0, review_terminal.stat().st_size - 65536))
            suffix = stream.read().decode('utf-8')
        header.update(json.loads('{"schema":' + suffix.rsplit('], "schema":', 1)[1]))
        # The rows are large rational evidence; retain only their decisions.
        passed_rows = sum(row.get('passed') is True for row in aggregation_rows(review_terminal))
        parent = review_root.parents[1]
        result = json.loads((review_root/'result.json').read_text())
        complete = json.loads((review_root/'complete.json').read_text())
        review = dict(status='EXISTING_VERSIONED_REVIEW', manifest_match=header['manifest_sha256'] == manifest['identity'],
            protocol=header['protocol'], scope=header['scope'],
            recorded_derivative_gate=header['full_derivative_gate'], passed_rows=passed_rows,
            reviewer_authenticated=digest(parent/'review_saved_production.py') == header['reviewer_sha256'],
            witness_authenticated=digest(parent/'retriangulation_witness.py') == header['witness_sha256'],
            sources_match_manifest=header['sources'] == manifest['sources'],
            process_success=result['status'] == 'PROCESS_EXITED' and result['returncode'] == 0,
            result_authenticated=digest(review_root/'result.json') == complete['result_sha256'],
            terminal_sha256=digest(review_terminal), rounded_heart_admission_granted=header['rounded_heart_admission_granted'])
    return dict(schema='quadric-production-readonly-status-v1',
        manifest_identity=manifest['identity'], manifest_sha256=digest(manifest_path),
        dependencies=dependencies, points=points, rows=rows,
        complete_authenticated_points=sum(p['status'] == 'COMPLETE' and
            all(p.get(k) is True for k in ('manifest_match', 'packet_authenticated',
                'independent_authenticated', 'process_success', 'result_authenticated',
                'resource_receipt')) for p in points),
        derivative_pass_rows=sum(r['passed'] is True for r in rows),
        derivative_failed_rows=[r for r in rows if r['passed'] is not True],
        later_exact_retriangulation_review=review,
        aggregation_sha256=digest(aggregation) if aggregation.is_file() else None,
        native_build_receipt_present=(Path(manifest['mapping_settings']['source_paths']['native'])
            .parent.parent / 'build-receipt.json').is_file(),
        scientific_gate_granted=False, optimizer_invocations=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(inventory(args.evidence), sort_keys=True, allow_nan=False))
