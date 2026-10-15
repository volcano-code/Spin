"""Storage-only engineering comparison on unchanged reconstructed input.

Never grants geometry, physics or paper acceptance. Source observations are not
an attestation that arbitrary source files produced a previously frozen binary.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import subprocess
import time
import traceback
import uuid

from make_it_stand.verification import source_trace_archive as archive

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify_files(files):
    for path, digest in files.items():
        if not Path(path).is_file() or sha256(path) != digest:
            raise ValueError('STALE_RESULT: ' + str(path))


def write_new(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')


def compare_archives(left, right):
    left, right = Path(left).resolve(), Path(right).resolve()
    if left == right or left.samefile(right):
        raise ValueError('DISTINCT_RUN_ARCHIVES_REQUIRED')
    marker = object()
    count = 0
    header = None
    # Exhaust BOTH readers, including their terminal digests and trailing-byte
    # checks. zip() alone can stop before the longer archive is validated.
    for a, b in itertools.zip_longest(archive.read_archive(left),
                                    archive.read_archive(right), fillvalue=marker):
        # Python numeric equality aliases signed zero and int/float spellings.
        # Canonical JSON retains these observable scientific representation facts.
        if (a is marker or b is marker or
                json.dumps(a, sort_keys=True, separators=(',', ':'), allow_nan=False) !=
                json.dumps(b, sort_keys=True, separators=(',', ':'), allow_nan=False)):
            raise ValueError('RECORD_MISMATCH: ' + str(count))
        if count == 0:
            header = a
        count += 1
    if not count or header is None:
        raise ValueError('EMPTY_COMPARISON')
    return dict(archive_equivalence='VERIFIED', records=count,
                numerical_status=header['payload']['solve']['status'],
                paper_reproduction='NOT_EVALUATED',
                left_archive_sha256=sha256(left), right_archive_sha256=sha256(right))


def source_observation(source_root):
    """Observe effective files, including untracked sources, without build claims."""
    root = Path(source_root).resolve()
    selected = {root/'CMakeLists.txt'}
    for directory in ('include', 'src', 'bindings', 'cmake', 'python'):
        base = root/directory
        if base.is_dir():
            selected.update(p for p in base.rglob('*') if p.is_file() and
                            p.suffix in ('.h', '.hpp', '.hxx', '.c', '.cc', '.cpp', '.cxx', '.cmake', '.in', '.py'))
    selected.update((root/'experiments'/'teapot_spin2014').glob('*.py'))
    selected = sorted(p for p in selected if p.is_file())
    if not selected:
        raise ValueError('SOURCE_OBSERVATION_UNAVAILABLE')
    files = {str(p.resolve()): sha256(p) for p in selected}
    # Exact selected source paths avoid report/output circular fingerprints.
    paths = [str(p.relative_to(root)) for p in selected]
    diff = subprocess.run(['git', 'diff', '--binary', 'HEAD', '--', *paths],
                          cwd=root, capture_output=True, check=False)
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=root,
                            capture_output=True, text=True, check=False)
    available = diff.returncode == 0 and commit.returncode == 0
    return dict(source_root=str(root), files=files,
                commit_observed=commit.stdout.strip() if available else None,
                dirty_source_diff_sha256=hashlib.sha256(diff.stdout if available else b'UNAVAILABLE_NOT_GIT_CHECKOUT').hexdigest(),
                dirty_source_diff_scope='SELECTED_TRACKED_SOURCE_DIFF_FROM_HEAD' if available else 'UNAVAILABLE_NOT_GIT_CHECKOUT',
                untracked_source_coverage='INCLUDED_BY_CONTENT_IN_FILES',
                binding='CURRENT_SOURCE_OBSERVATION_NOT_COMPILATION_ATTESTATION')


def make_provenance(config, native_module, output, storage, store_max_bytes, *, source_root=ROOT):
    if storage not in ('memory', 'disk'):
        raise ValueError('INVALID_STORAGE_MODE')
    config, module, output = (Path(p).resolve() for p in (config, native_module, output))
    inputs = Path(json.loads(config.read_text())['input']).resolve()
    input_files = {str(p): sha256(p) for p in
                   (inputs/'outer.obj', inputs/'design.obj', inputs/'manifest.json')}
    checker_files = {str(p.resolve()): sha256(p) for p in (Path(__file__), Path(archive.__file__))}
    build = module.parent.parent
    cache = build/'CMakeCache.txt'
    if not cache.is_file():
        raise ValueError('NATIVE_BUILD_CACHE_REQUIRED: ' + str(cache))
    build_files = {str(cache): sha256(cache)}
    for relative in ('compile_commands.json', 'CMakeFiles/Makefile.cmake',
                     'CMakeFiles/_make_it_stand.dir/flags.make', 'CMakeFiles/_make_it_stand.dir/link.txt',
                     'CMakeFiles/make_it_stand_core.dir/flags.make', 'CMakeFiles/make_it_stand_core.dir/link.txt'):
        path = build/relative
        if path.is_file(): build_files[str(path)] = sha256(path)
    observed = source_observation(source_root)
    files = {**observed['files'], **input_files, **checker_files, **build_files,
             str(config): sha256(config), str(module): sha256(module)}
    return dict(schema='source-online-run-provenance-v2', run_id=str(uuid.uuid4()),
                output_directory=str(output), storage=storage, store_max_bytes=store_max_bytes,
                config_path=str(config), config_sha256=sha256(config), input_files=input_files,
                checker_files=checker_files, native_binary=str(module), native_binary_sha256=sha256(module),
                build_files=build_files, source_observation=observed, files=files,
                source_binding='CURRENT_SOURCE_AND_BUILD_FILES_OBSERVED_NOT_COMPILATION_ATTESTATION',
                frozen_binary_source_binding='HISTORICAL_COMPILE_SOURCES_NOT_ATTESTED',
                evidence='RECONSTRUCTED_STORAGE_REGRESSION', original_paper_case=False,
                authenticity='LOCAL_CONTENT_INTEGRITY_AND_ACCIDENTAL_COPY_PROTECTION_NOT_ADVERSARIAL_AUTHENTICATION')


def verify_provenance(provenance):
    verify_files(provenance['files'])
    observed = provenance['source_observation']
    if source_observation(observed['source_root']) != observed:
        raise ValueError('STALE_RESULT: effective source set or dirty-source fingerprint changed')


def output_hashes(output):
    """Only this run, recursively; bind durable native files as well as exports."""
    output = Path(output).resolve()
    files = {}
    for path in sorted(output.rglob('*')):
        if not path.is_file() or path == output/'outputs.json': continue
        resolved = path.resolve()
        if not resolved.is_relative_to(output):
            raise ValueError('OUTPUT_OUTSIDE_RUN_DIRECTORY')
        files[str(resolved)] = sha256(resolved)
    return files


def verified_run(directory):
    directory = Path(directory).resolve()
    try:
        p = json.loads((directory/'provenance.json').read_text())
        summary = json.loads((directory/'summary.json').read_text())
        recorded = json.loads((directory/'outputs.json').read_text())
        header = json.loads((directory/'header.json').read_text())
    except (OSError, ValueError) as exc:
        raise ValueError('RUN_MANIFEST_REQUIRED: ' + str(directory)) from exc
    if p.get('schema') != 'source-online-run-provenance-v2':
        raise ValueError('RUN_MANIFEST_VERSION')
    if p.get('output_directory') != str(directory) or summary.get('output_directory') != str(directory):
        raise ValueError('RUN_DIRECTORY_IDENTITY_MISMATCH')
    try: uuid.UUID(p['run_id'])
    except (KeyError, ValueError, TypeError) as exc: raise ValueError('RUN_ID_REQUIRED') from exc
    if (summary.get('run_id') != p['run_id'] or summary.get('storage') != p['storage'] or
            summary.get('provenance_sha256') != sha256(directory/'provenance.json')):
        raise ValueError('RUN_MANIFEST_IDENTITY_MISMATCH')
    verification = summary.get('storage_verification')
    if verification not in ('COMPLETE_ARCHIVE_READBACK', 'COMPLETE_NATIVE_STREAM_FOOTER_DIGEST') or (directory/'failure.json').exists():
        raise ValueError('INCOMPLETE_RUN')
    required = {str(directory/name) for name in ('provenance.json', 'summary.json', 'header.json', 'native-run.jsonl')}
    if not isinstance(recorded, dict) or not required.issubset(recorded):
        raise ValueError('INCOMPLETE_OUTPUT_MANIFEST')
    for path in recorded:
        resolved = Path(path).resolve()
        if str(resolved) != path or not resolved.is_relative_to(directory) or resolved == directory/'outputs.json':
            raise ValueError('INVALID_OUTPUT_MANIFEST_PATH')
    actual_outputs = output_hashes(directory)
    if actual_outputs != recorded:
        raise ValueError('OUTPUT_MANIFEST_FILE_SET_MISMATCH')
    verify_provenance(p)
    identities = {p['config_path']: p['config_sha256'], p['native_binary']: p['native_binary_sha256'],
                  **p['input_files'], **p['checker_files'], **p['build_files'], **p['source_observation']['files']}
    if identities != p['files']:
        raise ValueError('INCOMPLETE_PROVENANCE_FILE_BINDING')
    if p['storage'] == 'disk' and (summary.get('native_storage') or {}).get('mode') != 'DISK_ITERATIONS_V1':
        raise ValueError('NATIVE_STORAGE_MODE_MISMATCH')
    if p['storage'] == 'memory' and (summary.get('native_storage') or {}).get('mode', 'MEMORY') != 'MEMORY':
        raise ValueError('NATIVE_STORAGE_MODE_MISMATCH')
    # Bind summary and separately exported header to the archive. Native mode
    # was sequence-checked during generation and is bound to one complete file
    # digest here; its small envelope is read without reparsing every payload.
    archive_path = directory/'native-run.jsonl'
    if verification == 'COMPLETE_ARCHIVE_READBACK':
        first = None; count = 0
        for event in archive.read_archive(archive_path):
            if first is None: first = event
            count += 1
    else:
        first, footer = archive.read_archive_envelope(archive_path)
        counts = header.get('trace_counts')
        storage = summary.get('native_storage') or {}
        if (not isinstance(counts, list) or
                any(type(value) is not int or value < 0 for value in counts) or
                type(storage.get('resident_pass_count')) is not int or
                type(storage.get('record_count')) is not int or
                storage.get('resident_pass_count') != len(counts) or
                storage.get('record_count') != sum(counts)):
            raise ValueError('NATIVE_TRACE_COUNT_MISMATCH')
        count = 1 + len(counts) + sum(counts)
        if (footer.get('schema') != 'source-trace-archive-v1' or
                type(footer.get('records')) is not int or
                footer.get('records') != count or
                type(summary.get('archive', {}).get('records')) is not int or
                summary.get('archive', {}).get('records') != count):
            raise ValueError('ARCHIVE_EVENT_COUNT_MISMATCH')
        if summary.get('archive', {}).get('payload_sha256') != footer.get('payload_sha256'):
            raise ValueError('ARCHIVE_PAYLOAD_DIGEST_MISMATCH')
    if first != header or summary.get('trace_counts') != header['trace_counts']:
        raise ValueError('ARCHIVE_HEADER_BINDING_MISMATCH')
    solve = header['payload']['solve']
    for name in ('input_hash', 'trace_hash', 'output_hash'):
        if summary.get(name) != solve[name]: raise ValueError('SUMMARY_SCIENTIFIC_HASH_MISMATCH')
    if summary.get('numerical_status') != solve['status'] or summary['archive']['records'] != count:
        raise ValueError('SUMMARY_ARCHIVE_MISMATCH')
    if summary['archive'].get('archive_sha256') != recorded.get(str(archive_path)):
        raise ValueError('SUMMARY_ARCHIVE_DIGEST_MISMATCH')
    if verification == 'COMPLETE_ARCHIVE_READBACK':
        with archive_path.open('rb') as stream:
            stream.seek(max(0, archive_path.stat().st_size - 4096))
            footer = json.loads(stream.read(4096).splitlines()[-1])
        if summary['archive'].get('payload_sha256') != footer['payload_sha256']:
            raise ValueError('SUMMARY_ARCHIVE_PAYLOAD_DIGEST_MISMATCH')
    return p, summary


def compare_runs(left, right, *, kind):
    """Content integrity and accidental-copy checks, not adversarial authentication."""
    if kind not in ('cross-mode', 'legacy-binary', 'deterministic-repeat'):
        raise ValueError('EXPLICIT_COMPARISON_KIND_REQUIRED')
    left, right = Path(left).resolve(), Path(right).resolve()
    a, _ = verified_run(left); b, _ = verified_run(right)
    if left == right or a['run_id'] == b['run_id']:
        raise ValueError('DISTINCT_RUN_INVOCATIONS_REQUIRED')
    for name in ('config_sha256', 'input_files', 'checker_files'):
        if a[name] != b[name]: raise ValueError('INPUT_OR_CHECKER_IDENTITY_MISMATCH: ' + name)
    modes = (a['storage'], b['storage'])
    same_binary = a['native_binary_sha256'] == b['native_binary_sha256']
    if kind == 'cross-mode' and (set(modes) != {'memory', 'disk'} or not same_binary):
        raise ValueError('CROSS_MODE_REQUIRES_MEMORY_DISK_AND_SAME_BINARY_IDENTITY')
    if kind == 'legacy-binary' and (modes != ('memory', 'memory') or same_binary):
        raise ValueError('LEGACY_BINARY_REQUIRES_DISTINCT_BINARIES_IN_MEMORY_MODE')
    if kind == 'deterministic-repeat' and (modes != ('disk', 'disk') or not same_binary):
        raise ValueError('DETERMINISTIC_REPEAT_REQUIRES_DISK_MODE_AND_SAME_BINARY_IDENTITY')
    result = compare_archives(left/'native-run.jsonl', right/'native-run.jsonl')
    result.update(comparison_kind=kind, run_equivalence='VERIFIED',
                  left_run_id=a['run_id'], right_run_id=b['run_id'],
                  left_output_directory=str(left), right_output_directory=str(right),
                  authenticity='LOCAL_CONTENT_INTEGRITY_AND_ACCIDENTAL_COPY_PROTECTION_NOT_ADVERSARIAL_AUTHENTICATION')
    return result


def export_diagnostics(result, path, encoding='dictionary'):
    """Explicit v1 transport choice; never skips sequence/integrity checks."""
    if encoding == 'dictionary':
        events, writer = result.iter_diagnostics(), archive.write_archive
    elif encoding == 'native-json':
        # Native bytes are already serialized and are emitted from the sealed
        # binary SLQP reader. The specialized writer checks their event envelope
        # and sequence without parsing every multi-megabyte payload again.
        return archive.write_trusted_native_archive(path, result)
    else:
        raise ValueError('UNSUPPORTED_TRACE_ENCODING')
    first = next(events)
    header = json.loads(first) if encoding == 'native-json' else first
    facts = writer(path, itertools.chain([first], events))
    return header, facts


def run(config, output, storage, store_max_bytes, encoding='dictionary'):
    import numpy as np
    import trimesh
    import _make_it_stand as native

    config = Path(config).resolve()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    try:
        cfg = json.loads(config.read_text())
        if cfg.get('asset_evidence') != 'RECONSTRUCTED' or cfg.get('original_paper_case') is not False:
            raise ValueError('RECONSTRUCTED_ENGINEERING_INPUT_REQUIRED')
        inputs = Path(cfg['input']).resolve()
        geometry = [inputs/'outer.obj', inputs/'design.obj']
        manifest = json.loads((inputs/'manifest.json').read_text())
        verify_files(manifest['geometry_hashes'])
        provenance = make_provenance(config, native.__file__, output, storage, store_max_bytes)
        provenance.update(
            python_version=__import__('sys').version, numpy_version=np.__version__,
            trimesh_version=trimesh.__version__, trace_encoding=encoding)
        write_new(output/'provenance.json', provenance)
        verify_provenance(provenance)
        meshes = [trimesh.load(path, force='mesh', process=False) for path in geometry]
        args = [array for mesh in meshes for array in
                (np.asarray(mesh.vertices), np.asarray(mesh.faces))]
        extra = {} if storage == 'memory' else dict(
            iteration_store_directory=str(output/'iteration-store'),
            iteration_store_max_bytes=store_max_bytes)
        before = time.monotonic()
        result = native.solve_source_spin2014_material(*args,
            source_settings=cfg['geometry'], solver_settings=cfg['solver'],
            boundary_settings=cfg['boundary'], density=1.,
            minimum_outer_wall=cfg['minimum_wall'], **extra)
        solved = time.monotonic()
        solve_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print('NATIVE_DONE', 'seconds', solved-before, 'peak_RSS_KiB', solve_rss, flush=True)
        header, facts = export_diagnostics(result, output/'native-run.jsonl', encoding)
        exported = time.monotonic()
        if encoding == 'native-json':
            native_storage = result.trace_storage
            read_records = facts['records']
            expected_records = 1 + native_storage['resident_pass_count'] + native_storage['record_count']
            if read_records != expected_records:
                raise ValueError('NATIVE_ARCHIVE_RECORD_COUNT')
            storage_verification = 'COMPLETE_NATIVE_STREAM_FOOTER_DIGEST'
        else:
            read_records = sum(1 for _ in archive.read_archive(output/'native-run.jsonl'))
            if facts['records'] != read_records:
                raise ValueError('ARCHIVE_RECORD_COUNT')
            storage_verification = 'COMPLETE_ARCHIVE_READBACK'
        verify_provenance(provenance)
        write_new(output/'header.json', header)
        s = header['payload']['solve']
        write_new(output/'summary.json', dict(
            run_id=provenance['run_id'], output_directory=str(output), storage=storage, trace_encoding=encoding,
            provenance_sha256=sha256(output/'provenance.json'),
            storage_verification=storage_verification, numerical_status=s['status'],
            paper_reproduction='NOT_EVALUATED', physical_validation='NOT_RUN',
            input_hash=s['input_hash'], trace_hash=s['trace_hash'], output_hash=s['output_hash'],
            trace_counts=header['trace_counts'], archive=facts,
            native_seconds=solved-before, export_seconds=exported-solved,
            total_seconds=time.monotonic()-start, solve_peak_rss_kib=solve_rss,
            peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            native_storage=getattr(result, 'trace_storage', None),
            scope='ITERATION_STORAGE_ONLY_NOT_CONSTANT_TOTAL_ADAPTIVE_MEMORY'))
        write_new(output/'outputs.json', output_hashes(output))
        verify_provenance(provenance)
        print('STORAGE_ARCHIVE_COMPLETE', s['status'], 'records', read_records, flush=True)
    except Exception:
        write_new(output/'failure.json', dict(status='FAILED_OR_PARTIAL',
            elapsed_seconds=time.monotonic()-start, traceback=traceback.format_exc(),
            paper_reproduction='NOT_EVALUATED'))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    execute = sub.add_parser('run')
    execute.add_argument('--config', type=Path, required=True)
    execute.add_argument('--output', type=Path, required=True)
    execute.add_argument('--storage', choices=('memory', 'disk'), required=True)
    execute.add_argument('--store-max-bytes', type=int, default=8 * 1024**3)
    execute.add_argument('--trace-encoding', choices=('dictionary', 'native-json'), default='dictionary')
    compare = sub.add_parser('compare')
    compare.add_argument('left', type=Path)
    compare.add_argument('right', type=Path)
    compare.add_argument('--kind', choices=('cross-mode', 'legacy-binary', 'deterministic-repeat'), required=True)
    compare.add_argument('--output', type=Path, required=True)
    options = parser.parse_args()
    if options.command == 'run':
        resource.setrlimit(resource.RLIMIT_AS, (24 * 1024**3, 24 * 1024**3))
        run(options.config, options.output, options.storage, options.store_max_bytes, options.trace_encoding)
    else:
        destination = options.output.resolve()
        if any(destination.is_relative_to(p.resolve()) for p in (options.left, options.right)):
            raise ValueError('COMPARISON_OUTPUT_MUST_BE_OUTSIDE_RUNS')
        write_new(destination, compare_runs(options.left, options.right, kind=options.kind))
