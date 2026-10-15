"""Same native problem, lazy diagnostic export; no geometry/physics admission.

Each invocation is a fresh solve. This is an engineering replay, not a new
paper-original case. Never load the legacy result as an optimization warm start.
"""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import time
import traceback
import numpy as np
import trimesh
import _make_it_stand as native
from make_it_stand.verification import source_trace_archive as archive
from teapot import ROOT, write


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify(files):
    for path, expected in files.items():
        if not Path(path).is_file() or sha(path) != expected:
            raise ValueError('STALE_RESULT: '+path)


def execute(reference, output):
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    cfg = json.loads(reference.read_text())
    inputs = Path(cfg['input'])
    manifest = json.loads((inputs/'manifest.json').read_text())
    paths = [Path(native.__file__), Path(__file__), Path(__file__).with_name('teapot.py'),
             reference, inputs/'manifest.json', ROOT/'CMakeLists.txt']
    for directory in ('include','src','bindings','python/make_it_stand'):
        paths.extend(p for p in (ROOT/directory).rglob('*') if p.is_file() and
                     p.suffix in ('.h','.hpp','.cpp','.py'))
    build = Path(native.__file__).resolve().parent.parent
    for name in ('CMakeCache.txt','compile_commands.json','build.ninja','libmake_it_stand_core.a'):
        if (build/name).is_file(): paths.append(build/name)
    files = {str(p.resolve()):sha(p) for p in paths}
    files.update(manifest['geometry_hashes'])
    write(output/'source.json',dict(files=files,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope='SOURCE_AND_LOADED_BINARY_OBSERVATION_NOT_FROZEN_RELEASE_ATTESTATION'))
    solver = {**cfg['solver'],'termination_policy':'PAPER_POST_TRANSITION_V1'}
    write(output/'config.json',dict(reference_config=str(reference.resolve()),
        reference_sha256=sha(reference),reference_settings=cfg,
        source_settings=cfg['geometry'],solver_settings=solver,
        boundary_settings=cfg['boundary'],density=1.,minimum_outer_wall=cfg['minimum_wall'],
        scope='ENGINEERING_REPLAY_OF_EXPLICIT_PAPER_STOP_NOT_REFERENCE_DEFAULT_POLICY'))
    try:
        verify(files)
        meshes = [trimesh.load(inputs/(name+'.obj'),force='mesh',process=False)
                  for name in ('outer','design')]
        args = [a for m in meshes for a in (np.asarray(m.vertices),np.asarray(m.faces))]
        before = time.monotonic()
        result = native.solve_source_spin2014_material(*args,source_settings=cfg['geometry'],
            solver_settings=solver,boundary_settings=cfg['boundary'],density=1.,
            minimum_outer_wall=cfg['minimum_wall'])
        solved = time.monotonic()
        solve_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print('NATIVE_SOLVE_DONE', solved-before, 'peak_RSS_KiB', solve_rss, flush=True)
        # Header has terminal facts but no pass/iteration array. Never call diagnostics.
        header = None
        def events():
            nonlocal header
            for event in result.iter_diagnostics():
                if event['kind']=='header': header=event
                yield event
        archive_facts = archive.write_archive(output/'native-run.jsonl',events())
        exported = time.monotonic()
        export_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        del result
        gc.collect()
        # Full exhaustion checks EOF, digest, and missing/duplicated event sequence.
        read_records = sum(1 for _ in archive.read_archive(output/'native-run.jsonl'))
        if read_records != archive_facts['records']: raise ValueError('ARCHIVE_RECORD_COUNT')
        verify(files)
        s = header['payload']['solve']
        write(output/'header.json',header)
        summary = dict(status=s['status'],passes=len(header['trace_counts']),
            iterations=sum(header['trace_counts']),trace_hash=s['trace_hash'],output_hash=s['output_hash'],
            binary_check=s['binary_check'],archive=archive_facts,
            native_solve_seconds=solved-before,export_seconds=exported-solved,
            elapsed_seconds=time.monotonic()-start,solve_max_rss_kib=solve_rss,
            export_max_rss_kib=export_rss,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            original_paper_case=False,independent_python_seal=False,physical_validation='NOT_RUN',
            scope='POST_SOLVE_STREAMING_EXPORT_NOT_BOUNDED_NATIVE_SOLVER_MEMORY')
        write(output/'summary.json',summary)
        write(output/'outputs.json',{str(p.resolve()):sha(p) for p in output.iterdir() if p.is_file()})
        print('RESULT',s['status'],'archive_records',read_records,flush=True)
        return 0 if s['status']=='BINARY_MOMENT_VALIDATED' else 2
    except Exception:
        write(output/'failure.json',dict(traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-start))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--reference-config',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    raise SystemExit(execute(a.reference_config,a.output))
