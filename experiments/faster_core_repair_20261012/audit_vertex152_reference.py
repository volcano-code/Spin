"""Independent 90-digit inverse Gauss / ideal ray / radial cap reference."""
import json
from pathlib import Path
import numpy as np
import mpmath as m
import sys
import hashlib

def normals():
    phi=(1+np.sqrt(5.))/2
    v=[[-1,phi,0],[1,phi,0],[-1,-phi,0],[1,-phi,0],[0,-1,phi],[0,1,phi],[0,-1,-phi],[0,1,-phi],[phi,0,-1],[phi,0,1],[-phi,0,-1],[-phi,0,1]]
    v=[np.array(x)/np.linalg.norm(x) for x in v]
    f=[(0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),(3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1)]
    for _ in range(2):
        cache={};new=[]
        def mid(a,b):
            key=tuple(sorted((a,b)))
            if key not in cache:
                x=v[a]+v[b];cache[key]=len(v);v.append(x/np.linalg.norm(x))
            return cache[key]
        for a,b,c in f:
            ab,bc,ca=mid(a,b),mid(b,c),mid(c,a)
            new.extend([(a,ab,ca),(b,bc,ab),(c,ca,bc),(ab,bc,ca)])
        f=new
    return [m.matrix([m.mpf(float(t)) for t in x]) for x in v]

def run():
    m.mp.dps=int(sys.argv[1]) if len(sys.argv)>1 else 90
    root=Path(__file__).resolve().parent
    folder=root/'heart-radial-source-derivatives-job-01/output'
    base=np.load(folder/'point-000.npz');q=[m.mpf(float(x)) for x in base['coefficients']]
    ns=normals()
    # Ring point 90 lies on the registered truncation sphere.
    radius=m.mpf(float(np.linalg.norm(base['vertices'][90])))
    def evaluate(t):
        p=q.copy();p[5]+=t
        a=m.matrix([[p[4],p[9],p[8]],[p[9],p[5],p[7]],[p[8],p[7],p[6]]]);ai=a**-1;b=m.matrix(p[1:4]);c=-ai*b/2
        n=ns[111];den=(b.T*ai*b)[0]-4*p[0]
        origin=ai*(n/m.sqrt(4*(n.T*ai*n)[0]/den)-b/2)
        n=ns[143];delta=ns[160]-n
        aa=(delta.T*ai*delta)[0];bb=2*(n.T*ai*delta)[0];cc=(n.T*ai*n)[0]
        roots=sorted([(-bb-m.sqrt(bb*bb-4*aa*cc))/(2*aa),(-bb+m.sqrt(bb*bb-4*aa*cc))/(2*aa)])
        roots=[x for x in roots if 0<x<1];direction=ai*(n+roots[1]*delta)
        aa=(direction.T*direction)[0];bb=(origin.T*direction)[0];cc=(origin.T*origin)[0]-radius**2
        point=origin+direction*((-bb+m.sqrt(bb*bb-aa*cc))/aa)
        values,vectors=m.eigsy(a);axis=vectors[:,2]
        if (axis.T*(point-c))[0]<0:axis=-axis
        return c+(point-c)*(2*radius/(axis.T*(point-c))[0])
    reference=evaluate(m.mpf(0));der=m.diff(evaluate,m.mpf(0));prod=base['velocity'].reshape(-1,3,10)[152,:,5]
    report={'scope':'Independent formula reference; original three-scale gate unchanged', 'decimal_precision':m.mp.dps,
        'bindings':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__),folder/'point-000.npz',root/'heart-seed1-source-failure-registration-01.json')},
        'vertex':152,'source_vertex':90,'finite_normal':111,'ideal_edge':[143,160],
        'ideal_root_ordinal':1,'radius':float(radius),'original_protocol_passed':False,
        'base_coordinate_max_error':float(max(abs(reference[i]-float(base['vertices'][152,i])) for i in range(3))),
        'analytic_relative_error':float(m.norm(der-m.matrix(prod.tolist()))/max(1,m.norm(der))), 'checks':[]}
    for h in (3e-6,1.5e-6,.75e-6):
        h=m.mpf(h);fd=(evaluate(-h)-8*evaluate(-h/2)+8*evaluate(h/2)-evaluate(h))/(6*h)
        report['checks'].append({'h':float(h),'reference_fd_relative_error':float(m.norm(fd-der)/max(1,m.norm(der)))})
    print(json.dumps(report,indent=2))
    if len(sys.argv)>2:
        with Path(sys.argv[2]).open('x') as stream:json.dump(report,stream,indent=2)

if __name__=='__main__':run()
