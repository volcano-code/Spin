import json,time
from pathlib import Path
import numpy as np
import _make_it_stand as native
HERE=Path(__file__).resolve().parent
source=HERE/'noodles-terminal-recovery-job-03/output/material'
rows={}
for name,num in [('original',0),('fixed',1),('design',2),('retained',8),('cavity',9)]:
    from fractions import Fraction
    ep=source/f'mesh_{num:04d}_rational_vertices.json'
    v=np.array([[float(Fraction(x)) for x in p] for p in json.loads(ep.read_text())]) if ep.exists() else np.load(source/f'mesh_{num:04d}_vertices.npy')
    f=np.load(source/f'mesh_{num:04d}_faces.npy')
    t=time.monotonic();a=dict(native.audit_exact_triangle_boundary(v,f,'AnyNonzero'))
    rows[name]=dict(boundary_valid=a['boundary_valid'],intersecting_pairs=a['intersecting_pairs'],degenerate_faces=a['degenerate_faces'],seconds=time.monotonic()-t)
    print(name,rows[name],flush=True)
(HERE/'noodles-remaining-parts-audit-01.json').write_text(json.dumps(rows))
