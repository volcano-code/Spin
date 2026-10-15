"""Exact closest-pair continua and coverage by a finite union of triangles.

For compact convex sets every minimum-distance pair has the same displacement
d (strict convexity of squared norm). Thus all source minimizers equal
A intersect (B-d). This may be a point, a segment or a convex polygon.
Coverage applies only to this global minimizer set, NOT all other pairs or
constrained local extrema that can arise after removing a source star.
"""
from fractions import Fraction as Q
from .near_feature_witness import add,sub,mul,dot,cross,triangle,point_triangle,disjoint_witness
from .segment_material_query import rational_point,triangle_interval


def dimension(poly):
    if not poly:return -1
    differences=[sub(p,poly[0]) for p in poly[1:] if p!=poly[0]]
    if not differences:return 0
    return 2 if any(any(cross(differences[0],d)) for d in differences[1:]) else 1


def clip(poly,score):
    if not poly:return []
    output=[]
    previous=poly[-1];ps=score(previous)
    for current in poly:
        cs=score(current)
        if (ps<0 and cs>=0) or (ps>=0 and cs<0):
            output.append(add(previous,mul(sub(current,previous),ps/(ps-cs))))
        if cs>=0:output.append(current)
        previous,ps=current,cs
    unique=[]
    for p in output:
        if not unique or unique[-1]!=p:unique.append(p)
    if len(unique)>1 and unique[0]==unique[-1]:unique.pop()
    return unique


def inside(point,t,n):return point_triangle(point,t,n)==point


def closest_set(a,b,*,expected_squared_exact):
    witness=disjoint_witness(a,b,expected_squared_exact=expected_squared_exact)
    p,q=map(rational_point,(witness['point_a_exact'],witness['point_b_exact']))
    displacement=sub(q,p)
    a,an=triangle(a);b,bn=triangle(b)
    shifted=tuple(sub(v,displacement) for v in b)
    signs=[dot(bn,sub(v,shifted[0])) for v in a]
    if all(s==0 for s in signs):
        poly=list(a)
        for i in range(3):
            start=shifted[i];edge=sub(shifted[(i+1)%3],start)
            poly=clip(poly,lambda v:dot(cross(edge,sub(v,start)),bn))
    else:
        section=[a[i] for i in range(3) if signs[i]==0]
        for i in range(3):
            j=(i+1)%3
            if signs[i]*signs[j]<0:
                section.append(add(a[i],mul(sub(a[j],a[i]),signs[i]/(signs[i]-signs[j]))))
        section=sorted(set(section))
        if len(section)==1:
            poly=section if inside(section[0],shifted,bn) else []
        elif len(section)==2:
            interval=triangle_interval(section[0],section[1],shifted,bn)
            if interval is None:poly=[]
            else:
                lo,hi,_=interval;d=sub(section[1],section[0])
                poly=list(dict.fromkeys((add(section[0],mul(d,lo)),add(section[0],mul(d,hi)))))
        else:raise ValueError('unexpected plane/triangle section')
    if not poly:raise ValueError('empty closest-pair set inconsistent with exact witness')
    dim=dimension(poly)
    if dim<2:poly=sorted(set(poly));poly=[poly[0]] if dim==0 else [poly[0],poly[-1]]
    for point in poly:
        if not inside(point,a,an) or not inside(add(point,displacement),b,bn):
            raise ValueError('closest-pair vertex outside triangle')
    return {'dimension':dim,'squared_exact':witness['squared_exact'],
            'displacement_exact':list(map(str,displacement)),
            'source_vertices_exact':[list(map(str,v)) for v in poly],
            'target_vertices_exact':[list(map(str,add(v,displacement))) for v in poly]}


def cover(poly,triangles):
    """Exact coverage of a convex point/segment/polygon by closed triangles.

    triangles are already exact (triangle, normal) pairs from the checked mesh.
    Uncovered points are recomputed and tested against every covering triangle.
    """
    poly=[rational_point(v) for v in poly];dim=dimension(poly)
    if dim<0:raise ValueError('nonempty convex minimizer set required')
    def uncovered(point):return not any(inside(point,t,n) for t,n in triangles)
    witness=None
    if dim==0:
        if uncovered(poly[0]):witness=poly[0]
    elif dim==1:
        p,q=min(poly),max(poly);intervals=[]
        for t,n in triangles:
            found=triangle_interval(p,q,t,n)
            if found is not None:intervals.append(found[:2])
        intervals.sort();reach=Q(0)
        for lo,hi in intervals:
            if lo>reach:
                witness=add(p,mul(sub(q,p),(reach+lo)/2));break
            reach=max(reach,hi)
        if witness is None and reach<1:witness=add(p,mul(sub(q,p),(reach+1)/2))
    else:
        normal=next(cross(sub(poly[i],poly[0]),sub(poly[j],poly[0]))
                    for i in range(1,len(poly)) for j in range(i+1,len(poly))
                    if any(cross(sub(poly[i],poly[0]),sub(poly[j],poly[0]))))
        if any(dot(normal,sub(p,poly[0]))!=0 for p in poly):
            raise ValueError('nonplanar minimizer polygon')
        for i in range(len(poly)):
            edge=sub(poly[(i+1)%len(poly)],poly[i])
            if any(dot(cross(edge,sub(p,poly[i])),normal)<0 for p in poly):
                raise ValueError('minimizer polygon must be convex and cyclically ordered')
        pieces=[poly]
        for t,n in triangles:
            # Noncoplanar triangles cover at most a line in this 2D set.
            if any(cross(normal,n)) or dot(n,sub(poly[0],t[0]))!=0:continue
            remaining_pieces=[]
            for piece in pieces:
                remainder=piece
                for i in range(3):
                    start=t[i];edge=sub(t[(i+1)%3],start)
                    def score(v):return dot(cross(edge,sub(v,start)),n)
                    outside=clip(remainder,lambda v:-score(v))
                    if dimension(outside)==2 and any(score(v)<0 for v in outside):
                        remaining_pieces.append(outside)
                    remainder=clip(remainder,score)
                    if dimension(remainder)<2:break
            pieces=remaining_pieces
            if not pieces:break
        if pieces:
            piece=pieces[0]
            center=tuple(sum(v[i] for v in piece)/len(piece) for i in range(3))
            if uncovered(center):witness=center
            else:
                a,b=next((a,b) for a in piece for b in piece if any(cross(sub(a,center),sub(b,center))))
                # A noncoplanar plane meets this nondegenerate quadratic curve
                # at <=2 points. Try more than 2 times the number of triangles.
                for k in range(2*len(triangles)+3):
                    t=Q(1,k+2)
                    candidate=add(add(mul(center,1-t-t*t),mul(a,t)),mul(b,t*t))
                    if uncovered(candidate):witness=candidate;break
                if witness is None:raise ValueError('uncovered polygon witness construction failed')
    if witness is not None and not uncovered(witness):raise ValueError('invalid uncovered witness')
    return {'minimizer_dimension':dim,'all_global_minimizers_covered':witness is None,
            'uncovered_target_exact':None if witness is None else list(map(str,witness)),
            'other_local_extrema_excluded':False,'complete_wall_thickness':'UNKNOWN'}
