"""Exact initial-star exclusion and moving-center orthogonality diagnostics.

Implements only necessary local checks from Epsilon-shapes Sections 5.2/5.3.
A single excluded minimizer does NOT exclude its whole triangle pair: a nearest
point continuum can include both excluded and non-excluded points. No global
thickness or region-growth completion status can be obtained from this module.
"""
import math
import numpy as np
from .near_feature_witness import sub,dot,cross
from .segment_material_query import rational_point


class SourceStarQuery:
    def __init__(self,boundary):
        self.boundary=boundary
        self.point_cache={}

    def incident_faces(self,point):
        point=rational_point(point)
        if point not in self.point_cache:
            lo=np.array([math.nextafter(float(x),-math.inf) for x in point])
            hi=np.array([math.nextafter(float(x),math.inf) for x in point])
            b=self.boundary
            candidates=np.flatnonzero(np.all((b.high>=lo)&(b.low<=hi),axis=1))
            faces=[]
            for raw in candidates:
                i=int(raw);t,n=b._triangle(i)
                if dot(n,sub(point,t[0]))!=0:continue
                if all(dot(cross(sub(t[(j+1)%3],t[j]),sub(point,t[j])),n)>=0 for j in range(3)):
                    faces.append(i)
            self.point_cache[point]=tuple(faces)
        return self.point_cache[point]

    def audit(self,source_face,center,target,*,source_edge=None):
        b=self.boundary
        if type(source_face) is not int or not 0<=source_face<len(b.faces):
            raise ValueError('invalid source triangle')
        p,q=rational_point(center),rational_point(target)
        if p==q:raise ValueError('nonzero center-target displacement required')
        if source_face not in self.incident_faces(p):raise ValueError('center not on source triangle')
        seeds=set(map(int,b.faces[source_face]))
        mode='TRIANGLE'
        if source_edge is not None:
            if (len(source_edge)!=2 or any(type(i) is not int for i in source_edge)
                or len(set(source_edge))!=2 or not set(source_edge)<=seeds):
                raise ValueError('source edge must be an edge of the source triangle')
            seeds=set(source_edge);a,c=[rational_point(b.vertices[i]) for i in source_edge]
            if any(cross(sub(p,a),sub(c,a))) or dot(sub(p,a),sub(p,c))>0:
                raise ValueError('center not on source edge')
            mode='EDGE'
        target_faces=self.incident_faces(q)
        if not target_faces:raise ValueError('target not on mesh boundary')
        star_hits=[{'face':i,'shared_source_vertices':sorted(seeds&set(map(int,b.faces[i])))}
                   for i in target_faces if seeds&set(map(int,b.faces[i]))]
        displacement=sub(q,p)
        orthogonal_edges=set();orthogonal_faces=[]
        for i in target_faces:
            t,n=b._triangle(i)
            if not any(cross(displacement,n)):orthogonal_faces.append(i)
            for j in range(3):
                a,c=t[j],t[(j+1)%3]
                if not any(cross(sub(q,a),sub(c,a))) and dot(sub(q,a),sub(q,c))<=0:
                    if dot(displacement,sub(c,a))==0:
                        orthogonal_edges.add(tuple(sorted((int(b.faces[i,j]),int(b.faces[i,(j+1)%3])))))
        return {'source_mode':mode,'source_face':source_face,'source_vertices':sorted(seeds),
                'target_incident_faces':list(target_faces),'initial_star_witnesses':star_hits,
                'target_in_initial_star':bool(star_hits),
                'orthogonal_target_edges':[list(e) for e in sorted(orthogonal_edges)],
                'orthogonal_target_faces':orthogonal_faces,
                'moving_event_algebraic_condition':bool(orthogonal_edges or orthogonal_faces),
                'candidate_admission':('EXCLUDED_BY_INITIAL_STAR' if star_hits else 'NOT_EXCLUDED_NOT_VALIDATED'),
                'whole_triangle_pair_excluded':False,'region_growth_complete':False,
                'complete_wall_thickness':'UNKNOWN','wall_thickness_admission':'NOT_GRANTED'}
