import sys,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
if len(sys.argv)>1:
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
    for command in [['cmake','-S',str(HERE/'cherchi-traversal-diagnostic-binding-01'),'-B',str(out/'build'),'-DCMAKE_BUILD_TYPE=Release','-DPython_EXECUTABLE='+sys.executable],['cmake','--build',str(out/'build'),'--parallel','4']]:
        subprocess.run(command,check=True)
else:
    sys.path.insert(0,str(ROOT/'python'))
    from make_it_stand.verification.serial_source_job import run_serial_source_job
    files={'worker':Path(__file__)}
    for parent in ['cherchi-traversal-diagnostic-binding-01','cherchi-traversal-diagnostic-source-01/arrangements/code']:
        for p in (HERE/parent).glob('*'):
            if p.is_file():files[str(p.relative_to(HERE))]=p
    files['predicate_header']=HERE/'cherchi-traversal-diagnostic-source-01/arrangements/external/Indirect_Predicates/include/implicit_point.hpp'
    r=run_serial_source_job(ROOT,HERE/'cherchi-traversal-diagnostic-build-job-05',[sys.executable,'-B',str(Path(__file__)),'{job}/output/candidate'],files)
    print(r.status,r.returncode,flush=True)
