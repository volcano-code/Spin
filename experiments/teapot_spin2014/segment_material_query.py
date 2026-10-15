"""Exact segment intersections and parity queries on a closed triangle boundary.

Topology is checked here; geometric validity/self-intersection must be bound to
the independent exact mesh audit by the caller. Results do not certify thickness.
Broad-phase bounds round OUTWARD; all final predicates use rational arithmetic.
"""
from collections import defaultdict
from fractions import Fraction as Q
import math
import numpy as np

from .near_feature_witness import add, sub, mul, dot, cross


def rational_point(values):
    if len(values)!=3: raise ValueError('three coordinates required')
    if any(isinstance(x,bool) for x in values): raise ValueError('boolean coordinate')
    result=tuple(Q(x) if isinstance(x,(str,Q,int)) else Q(float(x)) for x in values)
    if not all(math.isfinite(float(x)) for x in result): raise ValueError('nonfinite point')
    return result


def triangle_interval(p,q,triangle,normal):
    """Closed segment parameter interval meeting a closed triangle, or None."""
    delta=sub(q,p)
    offset=dot(normal,sub(p,triangle[0])); slope=dot(normal,delta)
    low,high=Q(0),Q(1)
    if slope:
        low=high=-offset/slope
        if low<0 or low>1: return None
    elif offset:
        return None
    for i in range(3):
        edge=sub(triangle[(i+1)%3],triangle[i])
        base=dot(cross(edge,sub(p,triangle[i])),normal)
        change=dot(cross(edge,delta),normal)
        if change>0: low=max(low,-base/change)
        elif change<0: high=min(high,-base/change)
        elif base<0: return None
        if low>high: return None
    point=add(p,mul(delta,low))
    strict=all(dot(cross(sub(triangle[(i+1)%3],triangle[i]),sub(point,triangle[i])),normal)>0
               for i in range(3))
    return low,high,strict


class ClosedBoundaryQueries:
    def __init__(self,vertices,faces):
        self.vertices=np.asarray(vertices)
        self.faces=np.asarray(faces)
        if (self.vertices.ndim!=2 or self.vertices.shape[1]!=3 or len(self.vertices)==0 or
            self.vertices.dtype!=np.dtype('float64') or not np.isfinite(self.vertices).all() or
            self.faces.ndim!=2 or self.faces.shape[1]!=3 or len(self.faces)==0 or
            self.faces.dtype.kind not in 'iu' or self.faces.min()<0 or self.faces.max()>=len(self.vertices)):
            raise ValueError('finite float64 indexed triangle mesh required')
        self.vertices=self.vertices.copy(); self.faces=self.faces.copy()
        edges=defaultdict(list)
        for face in self.faces:
            if len(set(map(int,face)))!=3: raise ValueError('repeated face index')
            for i in range(3):
                a,b=int(face[i]),int(face[(i+1)%3])
                edges[tuple(sorted((a,b)))].append(1 if a<b else -1)
        if any(len(values)!=2 or sum(values)!=0 for values in edges.values()):
            raise ValueError('closed consistently wound edge-manifold required')
        triangles=self.vertices[self.faces]
        self.low=triangles.min(axis=1); self.high=triangles.max(axis=1)
        self.cache={}
        self.xmax=Q(float(self.vertices[:,0].max()))
        self.margin=max(Q(1),self.xmax-Q(float(self.vertices[:,0].min())))

    def _triangle(self,index):
        if index not in self.cache:
            t=tuple(rational_point(self.vertices[i]) for i in self.faces[index])
            n=cross(sub(t[1],t[0]),sub(t[2],t[0]))
            if dot(n,n)==0: raise ValueError('degenerate triangle')
            self.cache[index]=(t,n)
        return self.cache[index]

    def intersections(self,p,q):
        p,q=rational_point(p),rational_point(q)
        if p==q: raise ValueError('nonzero segment required')
        # Conversion of a rational to double is nearest-rounded. One outward
        # nextafter encloses it, also for exact zero and subnormal values.
        low=np.array([math.nextafter(float(min(a,b)),-math.inf) for a,b in zip(p,q)])
        high=np.array([math.nextafter(float(max(a,b)),math.inf) for a,b in zip(p,q)])
        candidates=np.flatnonzero(np.all((self.high>=low)&(self.low<=high),axis=1))
        result=[]
        for raw in candidates:
            i=int(raw); t,n=self._triangle(i)
            interval=triangle_interval(p,q,t,n)
            if interval is not None: result.append((i,*interval))
        return result

    def parity(self,point):
        point=rational_point(point)
        # End of ray is strictly beyond the mesh x maximum. Boundary hits on
        # triangle edges/vertices trigger a new direction, never a guessed count.
        extent=max(self.margin,self.xmax-point[0]+self.margin)
        for attempt in range(33):
            direction=(Q(1),Q(0),Q(0)) if attempt==0 else (Q(1),Q(attempt),Q(attempt*attempt+1))
            endpoint=add(point,mul(direction,extent))
            hits=self.intersections(point,endpoint)
            if any(low==0 for _,low,high,strict in hits): return 'BOUNDARY',attempt+1
            if any(low!=high or not strict for _,low,high,strict in hits): continue
            count=len(hits)
            return ('PARITY_INSIDE' if count%2 else 'PARITY_OUTSIDE'),attempt+1
        return 'UNKNOWN_DEGENERATE_RAYS',33

    def audit_segment(self,p,q):
        p,q=rational_point(p),rational_point(q)
        hits=self.intersections(p,q)
        endpoints=[any(low<=t<=high for _,low,high,_ in hits) for t in (Q(0),Q(1))]
        overlap=any(low<high for _,low,high,_ in hits)
        internal=[row for row in hits if row[2]>0 and row[1]<1]
        midpoint=mul(add(p,q),Q(1,2))
        location,attempts=self.parity(midpoint)
        if not all(endpoints): state='ENDPOINT_NOT_ON_BOUNDARY'
        elif overlap: state='BOUNDARY_OVERLAP'
        elif internal: state='INTERIOR_BOUNDARY_CONTACT'
        elif location=='PARITY_INSIDE': state='OPEN_SEGMENT_PARITY_INSIDE'
        elif location=='PARITY_OUTSIDE': state='OPEN_SEGMENT_PARITY_OUTSIDE'
        else: state='UNKNOWN'
        return {'state':state,'endpoint_on_boundary':endpoints,
                'midpoint_parity':location,'ray_attempts':attempts,
                'intersections':[{'face':i,'t0':str(lo),'t1':str(hi),'face_interior':strict}
                                 for i,lo,hi,strict in hits],
                'complete_wall_thickness':'UNKNOWN','wall_thickness_admission':'NOT_GRANTED'}
