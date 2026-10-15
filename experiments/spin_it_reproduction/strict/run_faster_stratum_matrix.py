#!/usr/bin/env python3
"""Capture the isolated Faster quadric-stratum matrix as strict evidence."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
from make_it_stand.verification.run_evidence import EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, finish_run, check_freshness
from make_it_stand.verification.spin_common import canonical

ROOT=Path(__file__).resolve().parents[3]
NAMES={"plane","elliptic-cylinder","elliptic-paraboloid","hyperbolic-paraboloid","parabolic-cylinder","cone","parallel-planes","intersecting-planes","double-plane","isolated-point","imaginary-ellipsoid","ellipsoid-control","one-sheet-hyperboloid-control","two-sheet-hyperboloid-control"}

def write(path,value):
    with path.open("x",encoding="ascii") as f:f.write(canonical(value)+"\n")

def run(out:Path, exe:Path, roots:dict[str,Path]):
    out=out.absolute(); out.mkdir(parents=False,exist_ok=False); exe=exe.resolve(strict=True)
    matrix=out/'stratum.json'; started=time.perf_counter()
    inputs={'schema':'faster-stratum-matrix-inputs-v1','executable':{'path':str(exe),'sha256':hashlib.sha256(exe.read_bytes()).hexdigest()},'expected_cases':sorted(NAMES)}
    ip=out/'inputs.json'; write(ip,inputs)
    begin=begin_run(repository_root=ROOT,runner_files=[Path(__file__)],checker_files=[ROOT/'tests/quadric_stratum_matrix.cpp',ROOT/'experiments/m1_stratum_matrix/CMakeLists.txt'],asset_files=[ip,exe],build_configuration_files=[ROOT/'CMakeLists.txt'],parameters=inputs,seed_protocol={'seed':0,'protocol':'deterministic isolated frozen stratum matrix'},dependency_roots=roots,native_modules=[],pure_python=True,source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
    write(out/'begin.json',begin.as_dict())
    p=subprocess.run([str(exe),str(matrix)],cwd=out,capture_output=True,text=True)
    log=out/'stratum.stdout'; log.write_text(p.stdout+p.stderr,encoding='utf-8')
    data=json.loads(matrix.read_text(encoding='utf-8')) if p.returncode==0 and matrix.is_file() else {}
    cases=data.get('cases',[])
    if p.returncode or data.get('schema')!='m1-stratum-matrix-result-v1' or not data.get('overall_pass') or {x.get('name') for x in cases}!=NAMES or any(x.get('violations') for x in cases): raise RuntimeError('isolated stratum matrix failed')
    result={'schema':'faster-stratum-matrix-scientific-v1','classification':'SYNTHETIC_REGRESSION','scope':'ISOLATED_FROZEN_STRATUM_MATRIX_ONLY_NOT_AUTHOR_CASE','matrix':str(matrix),'matrix_sha256':hashlib.sha256(matrix.read_bytes()).hexdigest(),'case_count':len(cases),'stdout_sha256':hashlib.sha256(log.read_bytes()).hexdigest()}
    write(out/'scientific.json',result)
    evidence=finish_run(begin,events=[{'iteration':0,'status':'ok','data':{'kind':'FASTER_STRATUM_MATRIX','matrix_sha256':result['matrix_sha256']}}],output_files=[matrix,log,out/'scientific.json'],measurements={'elapsed_seconds':time.perf_counter()-started})
    write(out/'evidence.json',evidence.as_dict()); fresh=check_freshness(evidence)
    summary={'status':'SYNTHETIC_STRATUM_MATRIX_VALIDATED','classification':'SYNTHETIC_REGRESSION','project_status':'PARTIAL_BLOCKED','original_evidence_admission':'BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS','freshness':{'status':fresh.status,'mismatches':list(fresh.mismatches)}};write(out/'summary.json',summary);return summary

def main():
 p=argparse.ArgumentParser();p.add_argument('--executable',type=Path,required=True);p.add_argument('--output-directory',type=Path,required=True);p.add_argument('--dependency-root',action='append',required=True);a=p.parse_args(); roots={}
 for x in a.dependency_root:
  n,s,v=x.partition('=')
  if not n or not s or n in roots:p.error('invalid dependency root')
  roots[n]=Path(v).resolve()
 try:r=run(a.output_directory,a.executable,roots)
 except Exception as e: print(str(e),file=sys.stderr);return 2
 print(canonical(r));return 0 if r['freshness']['status']=='CURRENT' else 1
if __name__=='__main__':raise SystemExit(main())
