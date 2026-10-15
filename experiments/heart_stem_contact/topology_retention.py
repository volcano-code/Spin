"""Explicit conservative manufacturing extension, NOT paper binary rounding.

Only refill existing carved cells. Local tests identify critical cubical links;
they are NOT a replacement for exact final-boundary and material validation.
Bounded uniform expansion preserves existing cells exactly, never upsamples an
optimization or claims a new source resolution. Intended for small experiments.
"""
import itertools
import math
import numpy as np


def _cells(values,maximum_cells=200000):
    if len(values)>maximum_cells: raise ValueError('RESOURCE_LIMIT: topology cells')
    result=set()
    for value in values:
        if len(value)!=3 or any(isinstance(x,(bool,np.bool_)) or not isinstance(x,(int,np.integer)) for x in value):
            raise ValueError('integer lattice coordinates required')
        result.add(tuple(int(x) for x in value))
    return result


def _plane(index,axis,depth,origin,side):
    return math.fma(math.ldexp(float(index),-depth),side,float(origin[axis]))


def expand_cuboids(boxes,*,depth,origin,side,maximum_cells=200000):
    if not isinstance(depth,int) or isinstance(depth,bool) or not 0<=depth<=20:
        raise ValueError('invalid geometry depth')
    if len(origin)!=3 or not np.isfinite([*origin,side]).all() or side<=0:
        raise ValueError('invalid grid')
    count=0; checked=[]
    for box in boxes:
        level=box['cell']['level']; xyz=box['cell']['coordinates']
        if not isinstance(level,int) or isinstance(level,bool) or not 0<=level<=depth:
            raise ValueError('invalid cuboid level')
        key=next(iter(_cells([xyz])))
        if any(x<0 or x>=2**level for x in key): raise ValueError('cell outside root')
        factor=2**(depth-level); count+=factor**3
        if count>maximum_cells: raise ValueError('RESOURCE_LIMIT: exact uniform expansion')
        lo=tuple(x*factor for x in key)
        expected_lo=[_plane(x,a,depth,origin,side) for a,x in enumerate(lo)]
        expected_hi=[_plane(x+factor,a,depth,origin,side) for a,x in enumerate(lo)]
        if box['lower']!=expected_lo or box['upper']!=expected_hi:
            raise ValueError('GRID_GEOMETRY_MISMATCH: no snapping or resampling')
        checked.append((lo,factor))
    out=set()
    for lo,factor in checked:
        out.update(itertools.product(*(range(x,x+factor) for x in lo)))
    if len(out)!=count: raise ValueError('overlapping/duplicate source cuboids')
    return out


def _components(cells):
    remaining=set(cells); count=0
    while remaining:
        count+=1; stack=[remaining.pop()]
        while stack:
            c=stack.pop()
            for a in range(3):
                for sign in [-1,1]:
                    q=list(c); q[a]+=sign; q=tuple(q)
                    if q in remaining: remaining.remove(q); stack.append(q)
    return count


def critical_cells(values):
    cells=_cells(values); affected=set(); edges=set(); vertices=set()
    # A checkerboard cross-section around an edge has four boundary sheets.
    for c in sorted(cells):
        for a,b in [(0,1),(0,2),(1,2)]:
            for da,db in itertools.product([-1,1],repeat=2):
                d=list(c); d[a]+=da; d[b]+=db; d=tuple(d)
                x=list(c); x[a]+=da; y=list(c); y[b]+=db
                if c<d and d in cells and tuple(x) not in cells and tuple(y) not in cells:
                    edges.add((c,d)); affected.update([c,d])
    corners={tuple(c[a]+s[a] for a in range(3)) for c in cells
             for s in itertools.product([0,1],repeat=3)}
    for v in sorted(corners):
        octants={tuple(v[a]+s[a] for a in range(3)) for s in itertools.product([-1,0],repeat=3)}
        inside=octants&cells; outside=octants-cells
        if _components(inside)>1 or _components(outside)>1:
            vertices.add(v); affected.update(inside)
    return affected,dict(critical_edges=[[list(c),list(d)] for c,d in sorted(edges)],
                        critical_vertices=[list(v) for v in sorted(vertices)])


def retain_critical_cells(values,*,maximum_passes=100):
    if not isinstance(maximum_passes,int) or isinstance(maximum_passes,bool) or maximum_passes<=0:
        raise ValueError('invalid pass budget')
    cells=_cells(values); trace=[]
    for step in range(maximum_passes):
        if not cells: raise ValueError('EMPTY_CAVITY: not a successful hollowing result')
        affected,facts=critical_cells(cells)
        if not affected: return cells,trace
        before=len(cells); cells-=affected
        trace.append(dict(iteration=step,before=before,after=len(cells),
            retained_cells=[list(c) for c in sorted(affected)],**facts))
    if not cells: raise ValueError('EMPTY_CAVITY: not a successful hollowing result')
    if not critical_cells(cells)[0]: return cells,trace
    raise ValueError('RESOURCE_LIMIT: topology retention iterations')


def boundary_mesh(values,*,depth,origin,side):
    cells=_cells(values); points={}; triangles=[]
    def index(c):
        if c not in points: points[c]=len(points)
        return points[c]
    for c in sorted(cells):
        for a in range(3):
            b=(a+1)%3; d=(a+2)%3
            for sign in [-1,1]:
                neighbor=list(c); neighbor[a]+=sign
                if tuple(neighbor) in cells: continue
                corners=[]
                for u,v in [(0,0),(1,0),(1,1),(0,1)]:
                    q=list(c); q[a]+=int(sign>0); q[b]+=u; q[d]+=v; corners.append(index(tuple(q)))
                if sign<0: corners.reverse()
                triangles.extend([(corners[0],corners[1],corners[2]),(corners[0],corners[2],corners[3])])
    vertices=np.array([[_plane(c[a],a,depth,origin,side) for a in range(3)] for c in points],float).reshape(-1,3)
    return vertices,np.array(triangles,np.int32).reshape(-1,3)


def box_moments(values,*,depth,origin,side):
    total=np.zeros(10)
    for c in sorted(_cells(values)):
        lo=np.array([_plane(c[a],a,depth,origin,side) for a in range(3)])
        hi=np.array([_plane(c[a]+1,a,depth,origin,side) for a in range(3)])
        lengths=hi-lo; center=(hi+lo)/2; volume=float(np.prod(lengths))
        if volume<=0: raise ValueError('unrepresentable cell')
        second=np.outer(center,center)+np.diag(lengths**2/12)
        total+=volume*np.r_[1.,center,np.diag(second),second[1,2],second[0,2],second[0,1]]
    return total
