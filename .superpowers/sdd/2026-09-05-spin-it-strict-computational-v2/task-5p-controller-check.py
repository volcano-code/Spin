"""Independent Full-unit-face, topology, tetrahedron and BBW-motion checker.
No project Python or native solver/checker imports; finite fixtures only.
"""
from collections import defaultdict, Counter
import ctypes
import itertools
import sys
import numpy as np

fma=ctypes.CDLL('libm.so.6').fma
fma.argtypes=[ctypes.c_double]*3
fma.restype=ctypes.c_double

def indexed(rows, n=None, dtype=float):
    if n is None: n=len(rows)
    assert len(rows)==n and [int(r[0]) for r in rows]==list(range(n))
    a=np.array(rows,dtype=dtype)
    assert np.isfinite(a).all()
    return a

def box(lo,hi):
    vol=np.prod(hi-lo); mean=(hi+lo)/2
    return vol*np.r_[1,mean,(lo*lo+hi*lo+hi*hi)/3,mean[1]*mean[2],mean[0]*mean[2],mean[0]*mean[1]]

def triangle_moments(points,tri):
    a,b,c=(points[tri[:,i]] for i in range(3))
    v=np.einsum('ij,ij->i',a,np.cross(b,c))/6
    s=a+b+c
    q=(np.einsum('ni,nj->nij',s,s)+np.einsum('ni,nj->nij',a,a)+np.einsum('ni,nj->nij',b,b)+np.einsum('ni,nj->nij',c,c))*v[:,None,None]/20
    return np.column_stack((v,v[:,None]*s/4,q[:,0,0],q[:,1,1],q[:,2,2],q[:,1,2],q[:,0,2],q[:,0,1]))

def error(a,b):
    return np.max(np.abs(a-b))/max(1e-30,np.max(np.abs(b)))

def verify(c):
    fixture,d,dl,np_,npoints,npatch,nt,no,nh=map(int,c['CASE'][0])
    assert dl==d+1 and nh==(8 if fixture==2 else 26)
    points=indexed(c['POINT'],npoints)[:,1:]
    moved=indexed(c['MOVED'],npoints)[:,1:]
    outer=indexed(c['OUTER'],no)[:,1:]
    keys=indexed(c['KEY'],npoints-no,np.int64)[:,1:]
    assert points[:no].tobytes()==outer.tobytes()
    assert list(map(tuple,keys))==sorted(set(map(tuple,keys)))
    root=np.array(c['ROOT'][0],float); origin,side=root[:3],root[3]
    planes=lambda k,level: np.array([fma(float(k[a])/(1<<level),side,origin[a]) for a in range(3)])
    assert np.array_equal(points[no:],np.array([planes(k,dl) for k in keys]))
    prim=indexed(c['PRIM'],np_)
    patches=indexed(c['PATCH'],npatch,np.int64)
    tri=indexed(c['TRI'],nt,np.int64)[:,1:]
    assert (tri>=no).all() and (tri<npoints).all()
    occupied={}
    for p in prim:
        i=int(p[0]); level=int(p[2]); lo=[int(k)<<(d-level) for k in p[3:6]]; hi=[(int(k)+1)<<(d-level) for k in p[3:6]]
        assert np.array_equal(p[6:9],planes(lo,d)) and np.array_equal(p[9:12],planes(hi,d))
        for xyz in itertools.product(*(range(lo[a],hi[a]) for a in range(3))):
            assert xyz not in occupied
            occupied[xyz]=i
    lo=tuple(map(int,c['EXPECTED_LO'][0]));hi=tuple(map(int,c['EXPECTED_HI'][0]))
    assert set(occupied)==set(itertools.product(*(range(lo[a],hi[a]) for a in range(3))))
    expected_faces={}
    for xyz,owner in occupied.items():
        for axis in range(3):
            u,v=(axis+1)%3,(axis+2)%3
            for sign in (-1,1):
                neighbor=list(xyz);neighbor[axis]+=sign
                other=occupied.get(tuple(neighbor),-1)
                if owner==other:continue
                key=(axis,xyz[axis]+(sign==1),xyz[u],xyz[v])
                pair=(owner,other) if sign==1 else (other,owner)
                if key in expected_faces: assert expected_faces[key]==pair
                expected_faces[key]=pair
    got_faces={};expected_refs=set();all_triangle_ids=[];coarse_fine=0;inserted_edge_patches=0
    key_tri=keys[tri-no]
    cross=np.cross(key_tri[:,1]-key_tri[:,0],key_tri[:,2]-key_tri[:,0])
    for p in patches:
        _,axis,plane,level,ku,kv,negative,positive,start,count=map(int,p)
        assert negative>=0 or positive>=0
        assert negative<np_ and positive<np_ and negative!=positive
        width=1<<(d-level); uv=(ku*width,kv*width)
        for uu,vv in itertools.product(range(uv[0],uv[0]+width),range(uv[1],uv[1]+width)):
            key=(axis,plane,uu,vv);assert key not in got_faces
            got_faces[key]=(negative,positive)
        assert count>=4 and start>=0 and start+count<=nt
        ids=list(range(start,start+count));all_triangle_ids.extend(ids)
        assert (key_tri[ids,:,axis]==2*plane).all()
        assert (cross[ids,axis]>0).all()
        assert (cross[ids][:,[(axis+1)%3,(axis+2)%3]]==0).all()
        assert int(cross[ids,axis].sum())==2*(2*width)**2
        if count>4:inserted_edge_patches+=1
        if negative>=0 and positive>=0 and prim[negative,2]!=prim[positive,2]:coarse_fine+=1
        for t in ids:
            if negative>=0:expected_refs.add((negative,t,1))
            if positive>=0:expected_refs.add((positive,t,-1))
    assert got_faces==expected_faces and all_triangle_ids==list(range(nt))
    refs=[tuple(map(int,r)) for r in c['REF']]
    assert len(refs)==len(set(refs)) and set(refs)==expected_refs
    groups=defaultdict(list)
    for owner,t,sign in refs:groups[owner].append((t,sign))
    raw=triangle_moments(points,tri);warped=triangle_moments(moved,tri)
    native=indexed(c['MOMREST'],np_)[:,1:];native_moved=indexed(c['MOMMOVED'],np_)[:,1:]
    expected_moments=[];got_moved=[];max_error=0.
    for p in prim:
        owner=int(p[0]);edges=Counter();use=groups[owner]
        for t,sign in use:
            vertices=list(tri[t]);
            if sign<0:vertices[1],vertices[2]=vertices[2],vertices[1]
            for a,b in zip(vertices,vertices[1:]+vertices[:1]):edges[(a,b)]+=1
        assert all(n==1 and edges[(b,a)]==1 for (a,b),n in edges.items()), 'unmatched or nonmanifold primitive edge'
        ids=np.array([t for t,s in use]);sign=np.array([s for t,s in use])[:,None]
        actual=(raw[ids]*sign).sum(axis=0);deformed=(warped[ids]*sign).sum(axis=0)
        analytic=box(p[6:9],p[9:12]);expected_moments.append(analytic);got_moved.append(deformed)
        max_error=max(max_error,error(actual,analytic),error(native[owner],analytic),error(deformed,native_moved[owner]))
    assert max_error<1e-10
    # Native whole-exterior integral independently compared to primitive sums.
    union_expected=box(planes(lo,d),planes(hi,d))
    union_native=np.array(c['UNIONREST'][0],float)
    union_moved=np.array(c['UNIONMOVED'][0],float)
    union_error=max(error(union_native,union_expected),error(np.sum(got_moved,axis=0),union_moved))
    assert union_error<1e-10
    beta=.13+.017*prim[:,1]
    triangle_beta=np.zeros(nt)
    for owner,t,sign in refs:triangle_beta[t]+=sign*beta[owner]
    weighted=(warped*triangle_beta[:,None]).sum(axis=0)
    weighted_error=error(weighted,beta@native_moved)
    assert weighted_error<1e-10
    # Independently re-evaluate all actual BBW motion from exported nodal fields.
    nodes=indexed(c['NODE']);cells=indexed(c['GRIDCELL'],dtype=np.int64);handles=indexed(c['HANDLE'],nh)
    ni=cells[:,4:];low=nodes[ni[:,0],1:4];high=nodes[ni[:,7],1:4]
    move_error=0.
    for p,out in zip(points,moved):
        choices=np.flatnonzero(np.all(p>=low,axis=1)&np.all(p<=high,axis=1))
        j=min(choices,key=lambda j:tuple(cells[j,1:4]))
        unit=(p-low[j])/(high[j]-low[j])
        alpha=np.array([np.prod([unit[a] if b&(1<<a) else 1-unit[a] for a in range(3)]) for b in range(8)])
        w=alpha@nodes[ni[j],4:]
        value=np.sum(w[:,None]*(handles[:,1:4]+handles[:,4,None]*(p-handles[:,1:4])+handles[:,5:8]),axis=0)
        move_error=max(move_error,np.max(np.abs(value-out)))
    assert move_error<1e-12
    if fixture==0:assert np_==496 and coarse_fine>0 and inserted_edge_patches>0
    print(f'PASS fixture={fixture} sourceDepth={d} primitives={np_} points={npoints} patches={npatch} triangles={nt} coarseFine={coarse_fine} hangingEdgePatches={inserted_edge_patches} moment_relative_error={max_error:.3e} union_error={union_error:.3e} weighted_error={weighted_error:.3e} bbw_motion_error={move_error:.3e}')

case=None;count=0
with open(sys.argv[1],encoding='utf-8') as stream:
    for line in stream:
        parts=line.split();assert parts
        if parts[0]=='CASE':assert case is None;case=defaultdict(list)
        assert case is not None
        if parts[0]=='END':verify(case);count+=1;case=None
        else:case[parts[0]].append(parts[1:])
assert case is None and count==3
print('SYNTHETIC_REGRESSION: shared surfaces/moments under tested motions; no selfintersection, wall, solver or original-paper acceptance.')
