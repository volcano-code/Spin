#!/usr/bin/env python3
"""Capture the native joint fill/cage-deformation verification suite.

Synthetic only: this records five complementary native checks, not a Dancing
Couple reproduction or an author-asset claim.
"""
from __future__ import annotations
import argparse, hashlib, subprocess, time, traceback
from datetime import datetime, timezone
from pathlib import Path
from make_it_stand.verification.run_evidence import EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run
from make_it_stand.verification.spin_common import canonical, native, parse

ROOT=Path(__file__).resolve().parents[3]
CHECKS=(
 ("functional","spin-it-2014-cage-functional-regression",("[PASS] scenario=spin-it-2014-cage-functional",)),
 ("solver","spin-it-2014-cage-solver-regression",("[PASS] scenario=spin-it-2014-cage-solver",)),
 ("shape","spin-it-2014-source-cage-shape-regression",("SYNTHETIC_REGRESSION source cage shape PASS",)),
 ("geometry","spin-it-2014-source-cage-geometry-regression",("actual-source geometry admission PASS",)),
 ("guarded","spin-it-2014-source-cage-guarded-slqp-regression",("actual-source guarded SLQP PASS",)),
)
def write(p,v):
 with p.open("x",encoding="ascii") as f:f.write(canonical(v)+"\n")
def run(out,*,build_directory,dependency_roots):
 out=Path(out).absolute();out.mkdir(parents=False,exist_ok=False); build=Path(build_directory).resolve()
 bins=[build/name for _,name,_ in CHECKS]
 if not all(x.is_file() and x.stat().st_mode&0o111 for x in bins):raise ValueError("joint-deformation native regressions unavailable")
 start=datetime.now(timezone.utc).isoformat(); tick=time.perf_counter(); begin=None; events=[]
 try:
  inputs=parse(canonical({"schema":"source-spin2014-deformation-inputs-v1","classification":"SYNTHETIC_REGRESSION","binaries":[{"name":x.name,"sha256":hashlib.sha256(x.read_bytes()).hexdigest()} for x in bins]})); inp=out/"inputs.json";write(inp,inputs)
  ext=native()
  if Path(ext.__file__).resolve().parents[1]!=build:raise ValueError("loaded extension differs from declared build")
  begin=begin_run(repository_root=ROOT,runner_files=[Path(__file__)],checker_files=[ROOT/"tests/spin_it_2014_cage_functional_regression.cpp",ROOT/"tests/spin_it_2014_cage_solver_regression.cpp",ROOT/"tests/spin_it_2014_source_cage_shape_regression.cpp",ROOT/"tests/spin_it_2014_source_cage_geometry_regression.cpp",ROOT/"tests/spin_it_2014_source_cage_guarded_slqp_regression.cpp"],asset_files=[inp,*bins],build_configuration_files=[ROOT/"CMakeLists.txt",build/"CMakeCache.txt"],parameters=inputs,seed_protocol={"seed":0,"protocol":"deterministic native joint deformation regressions"},dependency_roots=dependency_roots,native_modules=[ext],pure_python=False,source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY);write(out/"begin.json",begin.as_dict())
  records=[]
  for label,binary,markers in CHECKS:
   p=subprocess.run([str(build/binary)],cwd=build,capture_output=True,text=True); text=p.stdout+p.stderr; log=out/(label+".stdout");log.write_text(text,encoding="utf-8")
   if p.returncode or any(m not in text for m in markers):raise RuntimeError("joint deformation regression failed: "+binary)
   r={"label":label,"path":str(log),"sha256":hashlib.sha256(text.encode()).hexdigest(),"returncode":p.returncode,"required_markers":list(markers)};records.append(r);events.append({"iteration":len(events),"status":"ok","data":{"kind":"JOINT_DEFORMATION_NATIVE_CHECK",**r}})
  scientific={"schema":"source-spin2014-deformation-scientific-v1","classification":"SYNTHETIC_REGRESSION","checks":records};write(out/"scientific.json",scientific)
  ev=finish_run(begin,events=events,output_files=sorted(x for x in out.rglob("*") if x.is_file() and x.name!="begin.json"),measurements={"started_at":start,"elapsed_seconds":time.perf_counter()-tick});write(out/"evidence.json",ev.as_dict());fresh=check_freshness(ev)
  s={"status":"SYNTHETIC_JOINT_DEFORMATION_VALIDATED","classification":"SYNTHETIC_REGRESSION","project_status":"PARTIAL_BLOCKED","original_evidence_admission":"BLOCKED_UNKNOWN","freshness":{"status":fresh.status,"mismatches":list(fresh.mismatches)}};write(out/"summary.json",s);return s
 except Exception as e:
  if begin is not None and not (out/"evidence.json").exists():write(out/"evidence.json",finish_run(begin,events=events+[{"iteration":len(events),"status":"failed","data":{"error":str(e)}}],output_files=sorted(x for x in out.rglob("*") if x.is_file() and x.name!="begin.json"),measurements={"started_at":start,"elapsed_seconds":time.perf_counter()-tick}).as_dict())
  write(out/"failure.json",{"status":"FAILED","error":str(e),"traceback":traceback.format_exc()});raise
def main():
 p=argparse.ArgumentParser();p.add_argument("--output-directory",type=Path,required=True);p.add_argument("--build-directory",type=Path,required=True);p.add_argument("--dependency-root",action="append",required=True);a=p.parse_args();roots={}
 for x in a.dependency_root:
  n,sep,v=x.partition("=")
  if not n or not sep or n in roots:p.error("invalid or duplicate dependency root")
  roots[n]=Path(v).resolve()
 r=run(a.output_directory,build_directory=a.build_directory,dependency_roots=roots);print(canonical(r));return 0 if r["freshness"]["status"]=="CURRENT" else 1
if __name__=="__main__":raise SystemExit(main())
