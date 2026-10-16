"""Independent finite inverse-Gauss reference excludes cap as sole cause."""
import json,hashlib
from pathlib import Path
import numpy as np
import mpmath as m
from audit_vertex152_reference import normals
HERE=Path(__file__).resolve().parent
base=HERE/'heart-radial-source-derivatives-job-01/output/point-000.npz'
b=np.load(base);m.mp.dps=120;ns=normals();q=list(map(lambda x:m.mpf(float(x)),b['coefficients']))
rows=[]
for vertex,normal in [(27,92),(46,111)]:
    def value(t):
        p=q.copy();p[5]+=t
        a=m.matrix([[p[4],p[9],p[8]],[p[9],p[5],p[7]],[p[8],p[7],p[6]]]);ai=a**-1;linear=m.matrix(p[1:4]);n=ns[normal]
        den=(linear.T*ai*linear)[0]-4*p[0]
        return ai*(n/m.sqrt(4*(n.T*ai*n)[0]/den)-linear/2)
    exact=value(0);der=m.diff(value,m.mpf(0));prod=m.matrix(b['velocity'].reshape(-1,3,10)[vertex,:,5].tolist())
    checks=[]
    for h in [3e-6,1.5e-6,.75e-6]:
        h=m.mpf(h);fd=(value(-h)-8*value(-h/2)+8*value(h/2)-value(h))/(6*h)
        checks.append(dict(h=float(h),relative_error=float(m.norm(fd-der)/max(1,m.norm(der)))))
    rows.append(dict(vertex=vertex,normal=normal,source_type='FINITE_INVERSE_GAUSS_NO_CAP',coordinate_error=float(m.norm(exact-m.matrix(b['vertices'][vertex].tolist()))),analytic_relative_error=float(m.norm(der-prod)/max(1,m.norm(der))),checks=checks))
report=dict(decimal_precision=120,base_sha256=hashlib.sha256(base.read_bytes()).hexdigest(),records=rows,original_protocol_passed=all(c['relative_error']<=1e-6 for r in rows for c in r['checks']),admission=False)
(HERE/'heart-finite-reference-120digits-01.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
