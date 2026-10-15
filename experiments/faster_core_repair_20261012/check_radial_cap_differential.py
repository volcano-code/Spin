"""Independent differential of radial cap projection, literal coefficients."""
import json
from pathlib import Path
import numpy as np


def matrix(q):
    return np.array([[q[4],q[9],q[8]],[q[9],q[5],q[7]],[q[8],q[7],q[6]]])


def project(q, point, radius):
    a=matrix(q);center=-np.linalg.solve(a,q[1:4])/2
    eigenvalues,eigenvectors=np.linalg.eigh(a);axis=eigenvectors[:,2]
    if axis@point<0:axis=-axis
    offset=point-center;den=axis@offset
    if den<=0:raise ValueError('RADIAL_DENOMINATOR_NONPOSITIVE')
    return center+offset*(2*radius/den)


def differential(q, point, radius, dq, dp):
    a=matrix(q);da=matrix(dq);center=-np.linalg.solve(a,q[1:4])/2
    dc=-np.linalg.solve(a,da@center+dq[1:4]/2)
    values,vectors=np.linalg.eigh(a);axis=vectors[:,2]
    daxis=sum((vectors[:,k]@da@axis)/(values[2]-values[k])*vectors[:,k] for k in (0,1))
    if axis@point<0:axis=-axis;daxis=-daxis
    x=point-center;dx=dp-dc;den=axis@x;dden=daxis@x+axis@dx
    return dc+(2*radius/den)*dx-x*(2*radius*dden/den**2)


def main():
    root=Path(__file__).resolve().parent
    q=np.array(json.loads((root/'heart-seed1-source-failure-registration-01.json').read_text())['event']['quadric_coefficients'])
    error=json.loads((root/'heart-source-failure-replay-job-03/stdout.log').read_text())['error']
    fields=dict(a.split('=',1) for a in error.split(';')[1:] if '=' in a)
    vertices=np.array([[float(x) for x in a.split(',') if x] for a in fields['witness_vertices'].split('|') if a])
    point=vertices[89];radius=float(fields['truncation_radius']);rng=np.random.default_rng(0)
    maximum=0.
    radial=differential(q,point,radius,q,np.zeros(3))
    if np.linalg.norm(radial)>1e-8:
        raise ValueError('RADIAL_CAP_DEPENDS_ON_COEFFICIENT_SCALE')
    for _ in range(10):
        dq=rng.normal(size=10);dq/=np.linalg.norm(dq);dp=rng.normal(size=3)
        analytic=differential(q,point,radius,dq,dp)
        for h in (1e-4,1e-5,1e-6):
            def f(t):return project(q+t*dq,point+t*dp,radius)
            numerical=(-f(2*h)+8*f(h)-8*f(-h)+f(-2*h))/(12*h)
            error=np.linalg.norm(numerical-analytic)/max(1.,np.linalg.norm(analytic))
            maximum=max(maximum,error)
    if maximum>1e-6:raise ValueError(('RADIAL_DIFFERENTIAL_FAILED',maximum))
    print(json.dumps(dict(maximum_relative_error=maximum,
        coefficient_scale_direction_norm=float(np.linalg.norm(radial)),
        scope='Isolated radial projection formula; not production source or material derivative')))


if __name__=='__main__':main()
