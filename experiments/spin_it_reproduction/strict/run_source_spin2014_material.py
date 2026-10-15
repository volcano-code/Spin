#!/usr/bin/env python3
"""SYNTHETIC_REGRESSION: fixed Task5e native fixture to independent sealed mesh.

Each invocation writes a new directory, including unsuccessful attempts. The
scientific payload is path independent. Evidence observes declared files on disk
before/after the actual solve, and cannot attest a build or authenticate execution.
The sealed subdirectory persists geometry validation only, not a solver-run seal.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import time
import traceback

import numpy as np

from make_it_stand.verification import source_spin_material as bridge
from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, finish_run, check_freshness,
)
from make_it_stand.verification.spin_bundle import save_spin_artifact, load_spin_artifact
from make_it_stand.verification.spin_common import canonical, native, parse
from make_it_stand.verification.spin_problem import SpinArtifactProblem

ROOT = Path(__file__).resolve().parents[3]


def _box(lo, hi):
    # Exact prior native fixture byte ordering; not an original paper model.
    v = np.array([[hi[a] if i & (1 << a) else lo[a] for a in range(3)] for i in range(8)],dtype='<f8')
    f = np.array([[0,2,1],[1,2,3],[4,5,6],[5,7,6],[0,1,4],[1,5,4],
                  [2,6,3],[3,6,7],[0,4,2],[2,4,6],[1,3,5],[3,7,5]],dtype='<i4')
    return v,f


def fixture():
    outer, design = _box((-1,-1,0),(1,1,2)), _box((-.5,-.5,1),(.5,.5,1.5))
    shell = (np.vstack((outer[0],design[0])),np.vstack((outer[1],design[1][:,::-1]+8)))
    problem = SpinArtifactProblem(original=outer,fixed_regions=[('shell',shell,1.)],
        design_regions=[('design',design,1.)],length_unit='m',mass_unit='kg',
        target={'mode':'STATIC_BALANCE','axis':[0,0,1],'origin':[0,0,0],'support_points':[[0,0,0]]},
        minimum_outer_wall=.5)
    settings = dict(source_settings=bridge.SourceGeometrySettings(root_origin=(-1,-1,0),root_side=2.,maximum_depth=4),
        solver_settings=bridge.SourceSolverSettings(maximum_optimization_level=4,laplacian_weight=1.,
            slqp=bridge.AuditedSlqpSettings(maximum_iterations=300,equality_tolerance=1e-10,projected_gradient_tolerance=1e-8)),
        boundary_settings=bridge.SourceBoundarySettings())
    return problem,settings


def _write(path, value):
    with path.open('x',encoding='ascii') as stream:
        stream.write(canonical(value)+'\n')


def run(output_directory, *, build_directory, dependency_roots, boundary_settings=None):
    output = Path(output_directory).absolute()
    output.mkdir(parents=False,exist_ok=False)
    started = datetime.now(timezone.utc).isoformat(); tick = time.perf_counter()
    begin = None
    events=[]
    try:
        p,settings = fixture()
        if boundary_settings is not None:
            if type(boundary_settings) is not bridge.SourceBoundarySettings:
                raise TypeError('expected immutable SourceBoundarySettings')
            settings['boundary_settings'] = boundary_settings
        inputs = {'schema':'source-spin2014-synthetic-inputs-v1','classification':'SYNTHETIC_REGRESSION',
                  'problem_identity':p.identity,'settings':{k:asdict(v) for k,v in settings.items()}}
        inputs = parse(canonical(inputs))
        assets = []
        for label,mesh in [('outer',p.original),('fixed_shell',p.fixed_regions[0].mesh),('design',p.design_regions[0].mesh)]:
            for kind,array in zip(('vertices','faces'),mesh.arrays()):
                path=output/(label+'_'+kind+'.npy')
                with path.open('xb') as stream: np.lib.format.write_array(stream,array,version=(1,0),allow_pickle=False)
                assets.append(path)
        input_path=output/'inputs.json';_write(input_path,inputs);assets.append(input_path)
        extension=native()
        build=Path(build_directory).resolve()
        if Path(extension.__file__).resolve().parents[1] != build:
            raise ValueError('declared build directory does not contain the actually loaded extension')
        begin=begin_run(repository_root=ROOT,runner_files=[Path(__file__)],
            checker_files=[Path(bridge.__file__),ROOT/'python/make_it_stand/verification/spin_artifact.py',
                ROOT/'python/make_it_stand/verification/spin_problem.py',ROOT/'python/make_it_stand/verification/mass_moments.py'],
            asset_files=assets,build_configuration_files=[ROOT/'CMakeLists.txt',build/'CMakeCache.txt'],
            parameters=inputs,seed_protocol={'seed':0,'protocol':'native audited HiGHS simplex: threads=1, parallel=off, random_seed=0'},
            dependency_roots=dependency_roots,native_modules=[extension],pure_python=False,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
        _write(output/'begin.json',begin.as_dict())
        result=bridge.solve_and_seal_source_spin2014(p,**settings)
        scientific={'schema':'source-spin2014-synthetic-scientific-v1','classification':'SYNTHETIC_REGRESSION',
                    'scientific_digest':result.scientific_digest,'identity':result.identity}
        _write(output/'scientific.json',scientific)
        # Persist the complete ordered trace independently for direct inspection.
        _write(output/'native-run.json',result.identity['native_run'])
        _write(output/'bridge-observations.json',result.evidence)
        for pass_index,pas in enumerate(result.identity['native_run']['solve']['passes']):
            for record in pas['trace']:
                events.append({'iteration':len(events),'status':'ok',
                    'data':{'kind':'SLQP_ITERATION_OBSERVATION','pass_index':pass_index,'record':record}})
        sealed=False
        if result.artifact is not None:
            saved=save_spin_artifact(result.artifact,output/'sealed')
            reloaded=load_spin_artifact(output/'sealed',expected_problem_digest=p.input_digest)
            if saved.scientific_digest != reloaded.scientific_digest or saved.scientific_digest != result.artifact.scientific_digest:
                raise ValueError('persisted independent seal identity mismatch')
            sealed=True
        events.append({'iteration':len(events),'status':'ok' if sealed else 'failed',
            'data':{'kind':'ACTUAL_SOLVE_AND_SEAL_OUTCOME','numerical_status':result.status,
                    'seal_error':result.seal_error,'scientific_digest':result.scientific_digest,
                    'artifact_scientific_digest':result.identity['artifact_scientific_digest']}})
        payloads=sorted(path for path in output.rglob('*') if path.is_file() and path.name!='begin.json')
        evidence=finish_run(begin,events=events,output_files=payloads,
            measurements={'started_at':started,'elapsed_seconds':time.perf_counter()-tick})
        _write(output/'evidence.json',evidence.as_dict())
        fresh=check_freshness(evidence)
        summary={'status':result.status,'sealed':sealed,'seal_error':result.seal_error,
            'scientific_digest':result.scientific_digest,'artifact_scientific_digest':result.identity['artifact_scientific_digest'],
            'freshness':asdict(fresh),'classification':'SYNTHETIC_REGRESSION','project_status':'PARTIAL_BLOCKED',
            'observation_scope':'OBSERVED_FILES_ON_DISK_ONLY','build_source_attestation':'NOT_VERIFIED',
            'physical_validation':'NOT_IN_SCOPE','original_evidence_admission':'BLOCKED_UNKNOWN',
            'solver_target_feasibility':'UNKNOWN','bundle_scope':'GEOMETRY_SEAL_ONLY_NOT_SOLVER_EXECUTION_AUTHENTICATION'}
        _write(output/'summary.json',summary)
        return summary
    except Exception as exc:
        failure={'schema':'source-spin2014-synthetic-failure-v1','status':'FAILED',
                 'error':str(exc),'exception_type':type(exc).__name__,'traceback':traceback.format_exc()}
        _write(output/'failure.json',failure)
        if begin is not None and not (output/'evidence.json').exists():
            events.append({'iteration':len(events),'status':'failed','data':failure})
            evidence=finish_run(begin,events=events,
                output_files=sorted(path for path in output.rglob('*') if path.is_file() and path.name!='begin.json'),
                measurements={'started_at':started,'elapsed_seconds':time.perf_counter()-tick})
            _write(output/'evidence.json',evidence.as_dict())
        raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-directory',type=Path,required=True)
    parser.add_argument('--build-directory',type=Path,required=True)
    parser.add_argument('--dependency-root',action='append',required=True,metavar='NAME=PATH',
                        help='actual declared local dependency tree to observe (repeatable)')
    parser.add_argument('--maximum-intermediate-complexity',type=int,default=4000000)
    args=parser.parse_args()
    roots={}
    for item in args.dependency_root:
        name,separator,path=item.partition('=')
        if not name or not separator or name in roots: parser.error('invalid or duplicate dependency root')
        roots[name]=Path(path).resolve()
    result=run(args.output_directory,build_directory=args.build_directory,dependency_roots=roots,
        boundary_settings=bridge.SourceBoundarySettings(maximum_intermediate_complexity=args.maximum_intermediate_complexity))
    print(canonical(result))
    return 0 if result['sealed'] and result['freshness']['status']=='CURRENT' else 1


if __name__=='__main__':
    raise SystemExit(main())
