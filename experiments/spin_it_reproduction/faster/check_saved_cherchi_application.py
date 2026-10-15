#!/usr/bin/env python3
"""Check persisted Cherchi runs and re-load exact materials; never optimize."""
import argparse
import json
from pathlib import Path
import time
import numpy as np

from make_it_stand import spin
from make_it_stand.verification.cherchi_application import audit_saved_cherchi_application,SCOPE
from make_it_stand.verification.run_evidence import begin_run,finish_run,EXPLICIT_EFFECTIVE_SOURCE_POLICY

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--application",type=Path,required=True)
    parser.add_argument("--bundle-root",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[3]
    application=args.application.resolve(strict=True)
    bundles=args.bundle_root.resolve(strict=True)
    args.output.mkdir(parents=True,exist_ok=False)
    checks=root/"python/make_it_stand/verification"
    started=time.perf_counter()
    begin=begin_run(repository_root=root,runner_files=[Path(__file__)],
        checker_files=[*checks.glob("*.py"),root/"experiments/spin_it_reproduction/faster/run_paper_table1.py",
                       root/"experiments/spin_it_reproduction/faster/run_paper_figures.py",
                       root/"experiments/spin_it_reproduction/strict/run_faster_cherchi_newton_integration.py"],
        asset_files=[application],build_configuration_files=[root/"CMakeLists.txt"],
        parameters={"scope":SCOPE,"operation":"READ_SAVED_CHERCHI_APPLICATION","bundle_root":str(bundles)},
        seed_protocol={"seed":0,"protocol":"NO_OPTIMIZATION_SAVED_RUN_AND_MATERIAL_RELOAD"},
        dependency_roots={"numpy":Path(np.__file__).resolve().parent},
        native_modules=[spin._native_module()],pure_python=False,
        source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
    (args.output/"begin.json").write_text(begin.to_json()+"\n")
    report=audit_saved_cherchi_application(application,bundle_root=bundles)
    scientific=args.output/"scientific.json"
    scientific.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    events=[{"iteration":i,"status":"ok" if row["material_verified"] else "failed",
             "data":{"case":case["id"],**row}}
            for case in report["cases"] for i,row in enumerate(case["runs"])]
    evidence=finish_run(begin,events=events,output_files=[scientific],
                        measurements={"elapsed_seconds":time.perf_counter()-started})
    (args.output/"evidence.json").write_text(evidence.to_json()+"\n")
    print(json.dumps({"scope":SCOPE,"qualified_material_runs":sum(c["qualified_material_runs"] for c in report["cases"])}))

if __name__=="__main__":main()
