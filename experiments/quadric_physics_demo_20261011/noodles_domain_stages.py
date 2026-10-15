"""Durable exact INPUT preprocessing; grants no optimization admission."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path[:0]=[str(ROOT/'python'),str(ROOT/'experiments/faster_heart_contact_20261011')]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def next_batch(names,excluded,count):
    if count<=0 or names[:len(excluded)]!=excluded or len(excluded)>len(names):
        raise ValueError('EXCLUSION_PREFIX_MISMATCH')
    return names[len(excluded):len(excluded)+count]


def worker(case,out,previous=None):
    import numpy as np
    import trimesh
    import _make_it_stand as native
    from noodles_admission import validate_manifest,publish
    from admit_rounded_input import exact_rows,construct,save_region
    from make_it_stand.verification.spin_problem import MeshSnapshot
    case,out=Path(case),Path(out)
    out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((case/'manifest.json').read_text())
    validate_manifest(manifest,'dense-v5')
    names=manifest['construction']['fixed_parts']
    identity=sha(case/'manifest.json');excluded=[]
    def load_source(name):
        entry=manifest['parts'][name];path=case/entry['model_obj']
        if sha(path)!=entry['sha256']['model']:raise ValueError('FROZEN_SOURCE_CHANGED:'+name)
        mesh=trimesh.load(path,force='mesh',process=False)
        MeshSnapshot((mesh.vertices,mesh.faces)).exact_region()
        return exact_rows(mesh.vertices),mesh.faces.tolist()
    if previous:
        previous=Path(previous)
        predecessor=json.loads((previous/'stage.json').read_text())
        if predecessor['manifest_sha256']!=identity:raise ValueError('PREDECESSOR_INPUT_MISMATCH')
        if predecessor['native_sha256']!=sha(native.__file__):raise ValueError('PREDECESSOR_NATIVE_MISMATCH')
        path=previous/'design.rational.json'
        if sha(path)!=predecessor['design_sha256']:raise ValueError('PREDECESSOR_BOUNDARY_CHANGED')
        boundary=json.loads(path.read_text())
        domain=boundary['vertices'],boundary['faces']
        MeshSnapshot.from_rationals(*domain).exact_region()
        excluded=predecessor['excluded_parts']
    else:domain=load_source('inner-ellipsoid-candidate')
    batch=next_batch(names,excluded,8)
    if not batch:raise ValueError('PREPROCESSING_ALREADY_COMPLETE')
    started=time.monotonic()
    with (out/'events.jsonl').open('x') as stream:
        for name in batch:
            def event(status):
                row=dict(status=status,part=name,elapsed_seconds=time.monotonic()-started)
                stream.write(json.dumps(row)+'\n');stream.flush();os.fsync(stream.fileno())
                print(json.dumps(row),flush=True)
            event('EXACT_EXCLUSION_STARTED')
            domain=construct(native,domain,load_source(name),'difference')
            event('EXACT_EXCLUSION_FINISHED')
    save_region(out/'design.obj',domain)
    region=MeshSnapshot.from_rationals(*domain).exact_region()
    done=excluded+batch
    publish(out/'stage.json',dict(schema='noodles-exact-domain-preprocessing-stage-v1',
        manifest_sha256=identity,native_sha256=sha(native.__file__),excluded_parts=done,part_count=len(names),
        design_sha256=sha(out/'design.rational.json'),region_facts=region.facts,
        predecessor_stage_sha256=sha(previous/'stage.json') if previous else None,
        preprocessing_complete=len(done)==len(names),geometry_admission=False,
        optimization_admission=False,simulation_admission=False))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--job',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--previous',type=Path)
    args=parser.parse_args()
    if args.output:
        worker(args.input,args.output,args.previous);return
    if args.job is None:raise ValueError('JOB_REQUIRED')
    from make_it_stand.verification.source_job import run_source_job
    native=ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    case=args.input.resolve()
    manifest=json.loads((case/'manifest.json').read_text())
    files=dict(worker=Path(__file__),native=native,manifest=case/'manifest.json',
        validation=HERE/'noodles_admission.py',generator=HERE/'noodles_dense_v5_input.py',
        base_generator=HERE/'noodles_dense_v3_input.py',tube_helper=HERE/'noodles_input.py',
        boolean_helper=ROOT/'experiments/faster_heart_contact_20261011/admit_rounded_input.py',
        independent_helper=ROOT/'experiments/faster_heart_contact_20261011/rounded_input.py')
    for name,entry in manifest['parts'].items():files[name]=case/entry['model_obj']
    command=[str(Path(__file__)),'--input',str(case),'--output','{job}/output']
    if args.previous:
        previous=args.previous.resolve()
        command+=['--previous',str(previous)]
        files.update(predecessor=previous/'stage.json',boundary=previous/'design.rational.json')
    result=run_source_job(args.job.resolve(),['/usr/bin/env',
        'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
        '/home/lin/anaconda3/bin/python3','-B',*command],files,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    print(json.dumps(dict(status=result.status,returncode=result.returncode)),flush=True)
    raise SystemExit(0 if result.status=='PROCESS_EXITED' and result.returncode==0 else 1)


if __name__=='__main__':main()
