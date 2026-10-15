"""Same stored materials: exact boundary profile versus full Nef snapshot."""
import json
from pathlib import Path
import sys
import time
import numpy as np
from make_it_stand.verification.spin_problem import MeshSnapshot,rational_solver_candidate_snapshot
from make_it_stand.verification.faster_material import _raw_moments

HERE=Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    fixed=HERE/'projective-formal-pipeline-job-01/output/fixed'
    cases=[(name,fixed/(name+'-material.json'),2.) for name in ('ellipsoid','one_sheet','two_sheet')]
    run=HERE/'core-benchmark-faces-12288-seed-0-job-01/output/run'
    terminal=json.loads((run/'terminal-material-reference.json').read_text())['evaluation_file']
    density=json.loads((run/'configuration.json').read_text())['request']['material_density']
    cases.extend([('faces-12288-initial',run/'evaluation-000001.json',density),
                  ('faces-12288-terminal',run/terminal,density)])
    rows=[]
    for name,path,density in cases:
        material=json.loads(path.read_text());v=material['variable_rational_vertices'];f=material['variable_faces']
        start=time.monotonic();full=MeshSnapshot.from_rationals(v,f);full_seconds=time.monotonic()-start
        start=time.monotonic();candidate,mode=rational_solver_candidate_snapshot(v,f);candidate_seconds=time.monotonic()-start
        a,b=full.arrays(),candidate.arrays()
        row=dict(name=name,faces=len(f),full_snapshot_seconds=full_seconds,candidate_snapshot_seconds=candidate_seconds,
            candidate_mode=mode,arrays_equal=all(np.array_equal(x,y) for x,y in zip(a,b)),
            rational_coordinates_equal=full.rational_vertices==candidate.rational_vertices,
            independent_moments_equal=np.array_equal(_raw_moments(full,density),_raw_moments(candidate,density)),
            final_region_path_changed=False)
        rows.append(row);(output/'rows.json').write_text(json.dumps(rows))
    (output/'report.json').write_text(json.dumps(dict(rows=rows,core_admission=False,
        passed=all(r['arrays_equal'] and r['rational_coordinates_equal'] and r['independent_moments_equal'] for r in rows),
        scope='Saved material trial snapshot equivalence; full solve and final seal still required')))


if __name__=='__main__':main(Path(sys.argv[1]))
