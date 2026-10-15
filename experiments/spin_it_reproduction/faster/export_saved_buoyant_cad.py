#!/usr/bin/env python3
"""Export a saved, independently rechecked buoyant stage; never re-optimize."""
import argparse
from pathlib import Path
import time

import numpy as np
from make_it_stand import spin
from make_it_stand.verification.buoyant_cad import export_saved_buoyant_step
from make_it_stand.verification.buoyant_evidence import audit_buoyant_execution, JOB_FILES
from make_it_stand.verification.run_evidence import begin_run, finish_run, EXPLICIT_EFFECTIVE_SOURCE_POLICY
from make_it_stand.verification.spin_common import parse, canonical
from make_it_stand.verification.spin_bundle import load_spin_artifact
from make_it_stand.verification.cad_coordinates import export_original_frame_cad


def main():
    parser=argparse.ArgumentParser()
    for name in ('application','bundle-root','execution-evidence','source-job-root','occt-executable','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--stage',choices=('A','B'),default='B')
    parser.add_argument('--original-frame',action='store_true',
                        help='also restore model coordinates and declared STEP length unit')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[3]
    paths={key:getattr(args,key).resolve(strict=True) for key in
           ('application','bundle_root','execution_evidence','source_job_root')}
    occt=args.occt_executable.resolve(strict=True)
    assets={paths['application'],paths['execution_evidence'],occt,
            *[paths['source_job_root']/name for name in JOB_FILES],
            *[p.resolve(strict=True) for p in paths['bundle_root'].rglob('*') if p.is_file()]}
    args.output.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter()
    begin=begin_run(repository_root=root,runner_files=[Path(__file__)],
        checker_files=[*sorted((root/'python/make_it_stand/verification').glob('*.py')),
                       root/'python/make_it_stand/spin.py',root/'experiments/occt_step_roundtrip/main.cpp'],
        asset_files=sorted(assets),build_configuration_files=[root/'experiments/occt_step_roundtrip/CMakeLists.txt'],
        parameters={'operation':'SAVED_BUOYANT_STAGE_CAD','stage':args.stage,
                    'original_frame':args.original_frame,
                    **{key:str(value) for key,value in paths.items()},'occt_executable':str(occt)},
        seed_protocol={'seed':0,'protocol':'NO_OPTIMIZATION_SAVED_STAGE_CAD'},
        dependency_roots={'numpy':Path(np.__file__).resolve().parent},
        native_modules=[spin._native_module()],pure_python=False,
        source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
    (args.output/'begin.json').write_text(begin.to_json()+'\n')
    original=audit_buoyant_execution(**paths)
    case=parse(paths['application'].read_text())['cases'][0]
    cad=export_saved_buoyant_step(case,bundle_root=paths['bundle_root'],run_index=0,
        stage=args.stage,occt_executable=occt,output_directory=args.output/'cad')
    report={'schema':'saved-buoyant-stage-cad-evidence-v1','original_execution':original,
            'cad':cad,'physical_validation':'NOT_IN_SCOPE',
            'manufacturing_approval':'NOT_GRANTED','original_paper_approval':'NOT_GRANTED'}
    if args.original_frame:
        selected=next(s for s in case['runs'][0]['stages'] if s['stage']==args.stage)
        artifact=load_spin_artifact(paths['bundle_root']/selected['relative_bundle'],
                                   expected_problem_digest=selected['problem_digest'])
        report['original_frame_cad']=export_original_frame_cad(artifact,cad['cad'],
            normalized_brep=args.output/'cad/material.brep',occt_executable=occt,
            output_directory=args.output/'original-frame')
    scientific=args.output/'scientific.json'
    scientific.write_text(canonical(report)+'\n')
    evidence=finish_run(begin,events=[{'iteration':0,'status':'ok','data':{
        'stage':args.stage,'artifact_scientific_digest':cad['artifact_scientific_digest'],
        'cad_input_digest':cad['cad']['cad_input_digest'],
        'maximum_normalized_material_moment_error':cad['maximum_normalized_material_moment_error']}}],
        # Successful staged CAD paths are hard-link aliases of the published
        # files. Bind each actual scientific file once. Diagnostic logs can
        # legitimately be empty; their size/hash is bound by the CAD sidecar.
        output_files=[scientific,*sorted(p for sub in ('cad','original-frame')
                       for p in (args.output/sub).glob('*') if p.is_file() and p.suffix!='.log')],
        measurements={'elapsed_seconds':time.perf_counter()-started})
    (args.output/'evidence.json').write_text(evidence.to_json()+'\n')
    print(canonical({'stage':args.stage,'physical_mass':cad['physical_mass'],
                     'moment_error':cad['maximum_normalized_material_moment_error']}))


if __name__=='__main__':main()
