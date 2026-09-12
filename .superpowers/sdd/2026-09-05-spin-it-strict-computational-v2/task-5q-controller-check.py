"""Independent fixed-rest BBW + signed tetrahedra / complex-step oracle.
Test-only dense arrays are small; no native/project imports or production formulas.
"""
from collections import defaultdict
import sys
import numpy as np

def indexed(rows, count=None, dtype=float):
    if count is None: count=len(rows)
    assert len(rows)==count and [int(r[0]) for r in rows]==list(range(count))
    a=np.array(rows,dtype=dtype)
    assert np.isfinite(a).all()
    return a[:,1:]

def tetra(points, faces):
    a,b,c=(points[faces[:,i]] for i in range(3))
    volume=np.einsum('ij,ij->i',a,np.cross(b,c))/6
    summed=a+b+c
    q=(np.einsum('ni,nj->nij',summed,summed)+np.einsum('ni,nj->nij',a,a)
       +np.einsum('ni,nj->nij',b,b)+np.einsum('ni,nj->nij',c,c))*volume[:,None,None]/20
    return np.column_stack((volume,volume[:,None]*summed/4,
        q[:,0,0],q[:,1,1],q[:,2,2],q[:,1,2],q[:,0,2],q[:,0,1]))

def verify(c):
    fixture=int(c['CASE'][0][0]);density=float(c['CASE'][0][1])
    nb,nh=map(int,c['CASE'][0][2:]);assert nh==26 and nb==8
    assert list(map(int,c['FLAGS'][0]))==[1,0,0,0]
    points=indexed(c['POINT']);moved=indexed(c['MOVED'],len(points))
    value_points=indexed(c['VALUEPOINT'],len(points))
    assert np.array_equal(moved,value_points)
    faces=indexed(c['OUTERFACE'],dtype=np.int64)
    triangles=indexed(c['TRI'],dtype=np.int64)
    assert faces.shape[1]==3 and triangles.shape[1]==3
    assert faces.min()>=0 and triangles.min()>=0
    assert faces.max()<len(points) and triangles.max()<len(points)
    beta=np.array(c['BETA'][0],float);params=np.array(c['PARAM'][0],float).reshape(nh,4)
    assert len(beta)==nb and np.isfinite(beta).all() and np.isfinite(params).all()
    nodes=indexed(c['NODE']);cells=indexed(c['CELL'],dtype=np.int64)
    handles=indexed(c['HANDLE'],nh)
    indices=cells[:,3:];lo=nodes[indices[:,0],:3];hi=nodes[indices[:,7],:3]
    weights=[]
    for p in points:
        choices=np.flatnonzero(np.all(p>=lo,axis=1)&np.all(p<=hi,axis=1))
        j=min(choices,key=lambda k:tuple(cells[k,:3]))
        u=(p-lo[j])/(hi[j]-lo[j])
        alpha=np.array([np.prod([u[a] if i&(1<<a) else 1-u[a] for a in range(3)]) for i in range(8)])
        weights.append(alpha@nodes[indices[j],3:])
    weights=np.array(weights);assert weights.shape==(len(points),nh)
    expected_points=np.sum(weights[:,:,None]*(handles[None,:,:]
        +params[None,:,0,None]*(points[:,None,:]-handles[None,:,:])+params[None,:,1:]),axis=1)
    motion_error=np.max(np.abs(expected_points-moved));assert motion_error<1e-12
    coefficients=np.zeros((len(triangles),nb))
    covered=set()
    for row in c['PATCH']:
        neg,pos,begin,count=map(int,row)
        assert -1<=neg<nb and -1<=pos<nb and count>=4 and begin>=0 and begin+count<=len(triangles)
        ids=range(begin,begin+count)
        assert not covered.intersection(ids);covered.update(ids)
        if neg>=0: coefficients[begin:begin+count,neg]-=1
        if pos>=0: coefficients[begin:begin+count,pos]+=1
    assert len(covered)==len(triangles)
    def moments(p,b):
        return density*(tetra(p,faces).sum(axis=0)+(coefficients@b)@tetra(p,triangles))
    expected=moments(expected_points,beta)
    native=np.array(c['MOMENTS'][0],float)
    assert np.isfinite(native).all() and np.array_equal(native,np.array(c['VALUEMOMENTS'][0],float))
    # Fixed pre-solve scales; never depend on optimized mass or measured errors.
    length=np.linalg.norm(np.ptp(points[:8],axis=0))
    scales=density*length**np.array([3,4,4,4,5,5,5,5,5,5])
    moment_error=np.max(np.abs(expected-native)/scales);assert moment_error<1e-10
    jb=indexed(c['JB'],nb).T;jh=indexed(c['JH'],4*nh).T
    expected_jb=np.empty_like(jb);expected_jh=np.empty_like(jh)
    epsilon=1e-30
    for k in range(nb):
        perturbed=beta.astype(complex);perturbed[k]+=1j*epsilon
        expected_jb[:,k]=moments(expected_points,perturbed).imag/epsilon
    for h in range(nh):
        for a in range(4):
            delta=weights[:,h,None]*(points-handles[h]) if a==0 else weights[:,h,None]*np.eye(3)[a-1]
            expected_jh[:,4*h+a]=moments(expected_points.astype(complex)+1j*epsilon*delta,beta).imag/epsilon
    expected_j=np.column_stack((expected_jb,expected_jh));native_j=np.column_stack((jb,jh))
    difference=np.abs(expected_j-native_j)
    assert np.all(difference<=1e-10*scales[:,None]+1e-7*np.abs(expected_j))
    scaled_j_error=np.max(difference/scales[:,None])
    # Beta derivatives have exact independent linear meaning, not only FD agreement.
    beta_direct=density*tetra(expected_points,triangles).T@coefficients
    assert np.max(np.abs(beta_direct-jb)/scales[:,None])<1e-10
    print(f'PASS fixture={fixture} variables={nb+4*nh} points={len(points)} '
          f'moment_scaled_error={moment_error:.3e} complexstep_jacobian_scaled_error={scaled_j_error:.3e} '
          f'bbw_motion_error={motion_error:.3e}')

case=None;count=0
with open(sys.argv[1],encoding='utf-8') as stream:
    for line in stream:
        parts=line.split();assert parts
        if parts[0]=='CASE':assert case is None;case=defaultdict(list)
        assert case is not None
        if parts[0]=='END':verify(case);count+=1;case=None
        else:case[parts[0]].append(parts[1:])
assert case is None and count==2
print('SYNTHETIC_REGRESSION only; no final geometry, joint optimizer or original-paper grant.')
