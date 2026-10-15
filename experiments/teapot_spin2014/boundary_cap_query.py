"""Exact fixed-center open-ball boundary topology filtration.

Requires an independently authenticated embedded closed orientable PL surface.
Each triangle/ball and shared-edge/ball intersection is convex. The finite good
cover changes only at the exact minimum squared distances to faces, edges and
vertices. Connected components and chi = F-E+V characterize each cap component;
a connected orientable finite-type open surface with chi=1 is a disk.
This checks ALL positive radii through the threshold at ONE center, not all
surface centers, nor material-side epsilon modifications or global thickness.
"""
from collections import defaultdict
from fractions import Fraction as Q
import math
import numpy as np
from .near_feature_witness import sub,dot,point_triangle,point_segment
from .segment_material_query import rational_point


class BoundaryCapQuery:
    def __init__(self,boundary):
        self.boundary=boundary
        self.points=[rational_point(v) for v in boundary.vertices]
        self.edges=defaultdict(list);self.vertex_faces=defaultdict(list)
        for i,face in enumerate(boundary.faces):
            boundary._triangle(i)  # exact degeneracy check, in addition to edge topology
            for vertex in face:self.vertex_faces[int(vertex)].append(i)
            for j in range(3):self.edges[tuple(sorted((int(face[j]),int(face[(j+1)%3]))))].append(i)
        # The underlying query checks closed oriented edge links. Also reject
        # pinched vertices: each vertex link must be a single cycle.
        for vertex,faces in self.vertex_faces.items():
            link=defaultdict(set)
            for i in faces:
                a,b=[int(v) for v in boundary.faces[i] if int(v)!=vertex]
                link[a].add(b);link[b].add(a)
            if any(len(adj)!=2 for adj in link.values()):raise ValueError('nonmanifold vertex link')
            reached=set();pending=[next(iter(link))]
            while pending:
                x=pending.pop()
                if x in reached:continue
                reached.add(x);pending.extend(link[x]-reached)
            if len(reached)!=len(link):raise ValueError('disconnected vertex link')

    def audit(self,center,threshold_squared):
        if isinstance(threshold_squared,bool):raise ValueError('boolean squared radius')
        center=rational_point(center);limit=Q(threshold_squared)
        if limit<=0:raise ValueError('positive squared radius required')
        # Exact dyadic upper bound; no square-root rounding or underflow risk.
        k=limit.numerator.bit_length()-limit.denominator.bit_length()
        radius=Q(2)**((k+2)//2)
        while radius*radius<limit:radius*=2
        while (radius/2)**2>=limit:radius/=2
        lo=np.array([math.nextafter(float(x-radius),-math.inf) for x in center])
        hi=np.array([math.nextafter(float(x+radius),math.inf) for x in center])
        b=self.boundary
        candidates=np.flatnonzero(np.all((b.high>=lo)&(b.low<=hi),axis=1))
        events=defaultdict(lambda:{'faces':[],'edges':[],'vertices':[]})
        face_set=set();edge_set=set();vertex_set=set();on_boundary=False
        for raw in candidates:
            i=int(raw);t,n=b._triangle(i);nearest=point_triangle(center,t,n)
            squared=dot(sub(center,nearest),sub(center,nearest))
            if squared==0:on_boundary=True
            if squared<limit:
                events[squared]['faces'].append(i);face_set.add(i)
                face=b.faces[i]
                vertex_set.update(map(int,face))
                edge_set.update(tuple(sorted((int(face[j]),int(face[(j+1)%3])))) for j in range(3))
        if not on_boundary:raise ValueError('center must lie on authenticated boundary')
        for edge in sorted(edge_set):
            nearest=point_segment(center,self.points[edge[0]],self.points[edge[1]])
            squared=dot(sub(center,nearest),sub(center,nearest))
            if squared<limit:events[squared]['edges'].append(edge)
        for vertex in sorted(vertex_set):
            delta=sub(center,self.points[vertex]);squared=dot(delta,delta)
            if squared<limit:events[squared]['vertices'].append(vertex)
        keys=sorted(events);parent={};active_edges=[];active_vertices=[];trace=[]
        def root(i):
            while parent[i]!=i:
                parent[i]=parent[parent[i]];i=parent[i]
            return i
        def union(a,c):
            a,c=root(a),root(c)
            if a!=c:parent[max(a,c)]=min(a,c)
        first_bad=None
        for index,value in enumerate(keys):
            event=events[value]
            for face in sorted(event['faces']):parent[face]=face
            for edge in event['edges']:
                a,c=self.edges[edge]
                if a not in parent or c not in parent:raise ValueError('edge filtration inconsistency')
                union(a,c);active_edges.append(edge)
            for vertex in event['vertices']:
                faces=self.vertex_faces[vertex]
                if not all(i in parent for i in faces):raise ValueError('vertex filtration inconsistency')
                for i in faces[1:]:union(faces[0],i)
                active_vertices.append(vertex)
            counts=defaultdict(lambda:[0,0,0])
            for face in parent:counts[root(face)][0]+=1
            for edge in active_edges:counts[root(self.edges[edge][0])][1]+=1
            for vertex in active_vertices:counts[root(self.vertex_faces[vertex][0])][2]+=1
            components=[{'representative_face':i,'faces':f,'edges':e,'vertices':v,'euler':f-e+v}
                        for i,(f,e,v) in sorted(counts.items())]
            disk=len(components)==1 and components[0]['euler']==1
            upper=keys[index+1] if index+1<len(keys) else limit
            record={'after_squared_radius_exact':str(value),'next_squared_radius_exact':str(upper),
                    'components':components,'single_disk':disk}
            trace.append(record)
            if not disk and first_bad is None:
                first_bad={'after_squared_radius_exact':str(value),'before_squared_radius_exact':str(upper),
                           'witness_squared_radius_exact':str((value+upper)/2),'components':components}
        if not trace or keys[0]!=0:raise ValueError('missing boundary-center initial cap')
        return {'schema':'fixed-center-boundary-cap-filtration-v1',
                'center_exact':list(map(str,center)),'threshold_squared_exact':str(limit),
                'state':('SINGLE_DISK_ALL_POSITIVE_RADII_TO_THRESHOLD' if first_bad is None else 'NON_DISK_INTERVAL_FOUND'),
                'first_non_disk_interval':first_bad,'threshold_cap_single_disk':trace[-1]['single_disk'],
                'critical_intervals':trace,'independent_embedded_surface_audit_required':True,
                'all_surface_centers_checked':False,'material_side_classified':False,
                'complete_wall_thickness':'UNKNOWN','wall_thickness_admission':'NOT_GRANTED'}
