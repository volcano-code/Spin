"""Independent rational closest points for authenticated disjoint triangles.

This deliberately does not admit intersecting triangles or certify that the
open closest-point segment is inside material. Local face halfspaces alone do
not describe the solid at edges, vertices or nonconvex neighborhoods.
"""
from fractions import Fraction as Q
import math


def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def mul(a,t): return tuple(x*t for x in a)
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def triangle(value):
    if len(value)!=3 or any(len(row)!=3 for row in value):
        raise ValueError('triangle required')
    for row in value:
        for x in row:
            if isinstance(x,bool) or not math.isfinite(x):
                raise ValueError('finite coordinates required')
    result=tuple(tuple(Q(float(x)) for x in row) for row in value)
    normal=cross(sub(result[1],result[0]),sub(result[2],result[0]))
    if dot(normal,normal)==0: raise ValueError('degenerate triangle')
    return result,normal


def point_segment(p,a,b):
    u=sub(b,a)
    t=max(Q(0),min(Q(1),dot(sub(p,a),u)/dot(u,u)))
    return add(a,mul(u,t))


def point_triangle(p,t,n):
    projected=sub(p,mul(n,dot(sub(p,t[0]),n)/dot(n,n)))
    if all(dot(cross(sub(t[(i+1)%3],t[i]),sub(projected,t[i])),n)>=0 for i in range(3)):
        return projected
    candidates=[point_segment(p,t[i],t[(i+1)%3]) for i in range(3)]
    return min(candidates,key=lambda q:(dot(sub(p,q),sub(p,q)),q))


def segment_pairs(p,p1,q,q1):
    candidates=[(p,point_segment(p,q,q1)),(p1,point_segment(p1,q,q1)),
                (point_segment(q,p,p1),q),(point_segment(q1,p,p1),q1)]
    u,v,w=sub(p1,p),sub(q1,q),sub(p,q)
    a,b,c,d,e=dot(u,u),dot(u,v),dot(v,v),dot(u,w),dot(v,w)
    denominator=a*c-b*b
    if denominator:
        s,t=(b*e-c*d)/denominator,(a*e-b*d)/denominator
        if 0<=s<=1 and 0<=t<=1:
            candidates.append((add(p,mul(u,s)),add(q,mul(v,t))))
    return candidates


def disjoint_witness(a,b,*,expected_squared_exact):
    """Must match a separately authenticated positive triangle distance.

    A zero oracle value, including edge-through-face intersection, is rejected:
    point/face and edge/edge candidates alone cannot detect that case.
    expected_squared_exact is a cross-check, not a caller-supplied pass flag.
    The caller must bind it to the exact native enumeration of the same input.
    """
    expected=Q(expected_squared_exact)
    if expected<=0: raise ValueError('positive independently checked separation required')
    a,an=triangle(a); b,bn=triangle(b)
    candidates=[(p,point_triangle(p,b,bn)) for p in a]
    candidates.extend((point_triangle(p,a,an),p) for p in b)
    for i in range(3):
        for j in range(3):
            candidates.extend(segment_pairs(a[i],a[(i+1)%3],b[j],b[(j+1)%3]))
    p,q=min(candidates,key=lambda pq:(dot(sub(pq[0],pq[1]),sub(pq[0],pq[1])),pq))
    distance=dot(sub(p,q),sub(p,q))
    if distance!=expected: raise ValueError('independent nearest-point distance mismatch')
    da,db=dot(an,sub(q,p)),dot(bn,sub(p,q))
    # Signs refer only to supplied triangle winding, not a global solid.
    side=('BOTH_NEGATIVE' if da<0 and db<0 else 'BOTH_POSITIVE' if da>0 and db>0 else
          'TANGENT_AT_LEAST_ONE' if da==0 or db==0 else 'MIXED_SIGNS')
    return {'squared_exact':str(distance),'point_a_exact':list(map(str,p)),
            'point_b_exact':list(map(str,q)),'face_dot_a_exact':str(da),'face_dot_b_exact':str(db),
            'local_halfspaces':side,'open_segment_material_membership':'UNKNOWN',
            'wall_thickness_admission':'NOT_GRANTED'}
