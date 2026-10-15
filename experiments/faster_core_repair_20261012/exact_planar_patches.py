"""Exact coplanar edge certificates for saved indexed material boundaries.

This does not certify embedding or correspondence of vertex provenance across
different parameters. It identifies whether connectivity changes disappear
when only exactly coplanar, consistently oriented internal edges are removed.
"""
from collections import defaultdict
from fractions import Fraction


def subtract(a,b):return tuple(x-y for x,y in zip(a,b))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b):return sum(x*y for x,y in zip(a,b))


def point_on_triangle(point,triangle):
    """Exact incidence of the closed triangle, not just its supporting plane."""
    p=tuple(Fraction(x) for x in point)
    a,b,c=[tuple(Fraction(x) for x in row) for row in triangle]
    normal=cross(subtract(b,a),subtract(c,a))
    if normal==(0,0,0):raise ValueError('EXACT_ZERO_AREA')
    if dot(normal,subtract(p,a))!=0:return False
    return all(dot(cross(subtract(v,u),subtract(p,u)),normal)>=0
               for u,v in ((a,b),(b,c),(c,a)))


def coalesce_collinear_cycle(points,cycle):
    """Certificate-view only: remove exact in-segment boundary subdivisions."""
    vertices=[tuple(Fraction(x) for x in row) for row in points]
    result=list(cycle)
    changed=True
    while changed and len(result)>3:
        changed=False
        for i,index in enumerate(result):
            a=vertices[result[i-1]];p=vertices[index];b=vertices[result[(i+1)%len(result)]]
            if a!=b and cross(subtract(p,a),subtract(b,a))==(0,0,0) and dot(subtract(p,a),subtract(p,b))<=0:
                result.pop(i);changed=True;break
    return result


def oriented_planar_chain_signature(vertices,faces):
    """Exact oriented surface-chain certificate, never a mesh admission check.

    Within each exact supporting plane a compact planar 2-chain is fixed by
    its oriented boundary 1-chain. Subdivide collinear boundary segments at
    exact endpoints, sum integer multiplicities and coalesce equal intervals.
    Exact-null triangles contribute the zero 2-chain ONLY in this proof view;
    their indices remain reported and the caller's geometry is not altered.
    """
    points=[tuple(Fraction(x) for x in p) for p in vertices]
    if any(len(p)!=3 for p in points):raise ValueError('BAD_COORDINATE_DIMENSION')
    planes={};null=[]
    for index,face in enumerate(faces):
        if len(face)!=3 or any(int(v)<0 or int(v)>=len(points) for v in face):raise ValueError('BAD_FACE')
        a,b,c=[points[int(v)] for v in face];normal=cross(subtract(b,a),subtract(c,a))
        if normal==(0,0,0):null.append(index);continue
        pivot=next(x for x in normal if x!=0)
        plane=tuple(x/pivot for x in normal)+(dot(normal,a)/pivot,)
        lines=planes.setdefault(plane,{})
        for u,v in ((a,b),(b,c),(c,a)):
            direction=subtract(v,u);axis=next(i for i,x in enumerate(direction) if x!=0)
            normalized=tuple(x/direction[axis] for x in direction)
            line=normalized+cross(u,normalized)
            events=lines.setdefault(line,{})
            low,high=sorted((u[axis],v[axis]));sign=1 if u[axis]<v[axis] else -1
            events[low]=events.get(low,0)+sign;events[high]=events.get(high,0)-sign
    signature=[]
    for plane,lines in sorted(planes.items()):
        line_signatures=[]
        for line,events in sorted(lines.items()):
            stops=sorted(events);active=0;intervals=[]
            for i,low in enumerate(stops[:-1]):
                active+=events[low];high=stops[i+1]
                if active:
                    if intervals and intervals[-1][1]==low and intervals[-1][2]==active:
                        intervals[-1]=(intervals[-1][0],high,active)
                    else:intervals.append((low,high,active))
            if stops and active+events[stops[-1]]!=0:raise ValueError('UNBALANCED_FINITE_EDGE_CHAIN')
            if intervals:line_signatures.append((line,tuple(intervals)))
        if line_signatures:signature.append((plane,tuple(line_signatures)))
    return dict(signature=tuple(signature),exact_null_triangle_indices=null,
                geometry_modified=False,mesh_admission=False)


def patch_boundary_signature(vertices,faces,*,closed=True):
    points=[tuple(Fraction(x) for x in row) for row in vertices]
    if any(len(p)!=3 for p in points):raise ValueError('BAD_COORDINATE_DIMENSION')
    edges=defaultdict(list);planes=[]
    for index,face in enumerate(faces):
        if len(face)!=3 or len(set(face))!=3 or any(v<0 or v>=len(points) for v in face):
            raise ValueError('BAD_FACE')
        a,b,c=[points[v] for v in face]
        n=cross(subtract(b,a),subtract(c,a))
        if n==(0,0,0):raise ValueError('EXACT_ZERO_AREA')
        pivot=next(x for x in n if x!=0)
        planes.append((tuple(x/pivot for x in n)+(dot(n,a)/pivot,),n))
        for j in range(3):
            u,v=face[j],face[(j+1)%3];edges[tuple(sorted((u,v)))].append((index,u,v))
    parent=list(range(len(faces)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    removed=[]
    for edge,uses in edges.items():
        if len(uses)>2 or (closed and len(uses)!=2):raise ValueError('NONMANIFOLD_OR_OPEN_EDGE')
        if len(uses)==1:continue
        left,right=uses
        if left[1:]!=right[1:][::-1]:raise ValueError('INCONSISTENT_WINDING')
        i,j=left[0],right[0]
        if planes[i][0]==planes[j][0]:
            if dot(planes[i][1],planes[j][1])<=0:raise ValueError('FOLDED_COPLANAR_FACES')
            parent[find(j)]=find(i)
            removed.append(dict(edge=list(edge),faces=[i,j],exact_coplanar=True))
    groups=defaultdict(list)
    for i in range(len(faces)):groups[find(i)].append(i)
    patches=[]
    for members in groups.values():
        member_set=set(members);boundary=[]
        for i in members:
            face=faces[i]
            for j in range(3):
                u,v=face[j],face[(j+1)%3]
                uses=edges[tuple(sorted((u,v)))]
                if not any(other!=i and other in member_set for other,_,_ in uses):boundary.append((u,v))
        successor={}
        for u,v in boundary:
            if u in successor:raise ValueError('PATCH_BOUNDARY_BRANCH')
            successor[u]=v
        if sorted(successor)!=sorted(successor.values()):raise ValueError('PATCH_BOUNDARY_OPEN')
        cycles=[]
        while successor:
            start=min(successor);loop=[];current=start
            while current in successor:
                loop.append(current);current=successor.pop(current)
            if current!=start:raise ValueError('PATCH_BOUNDARY_NOT_CYCLE')
            cycles.append(tuple(loop))
        if not cycles:raise ValueError('PATCH_WITHOUT_BOUNDARY')
        patches.append(tuple(sorted(cycles)))
    return dict(patches=sorted(patches),removed_edges=removed,
                exact_nonzero_face_count=len(faces),patch_count=len(patches))
