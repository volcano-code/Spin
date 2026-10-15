"""Exact endpoint-link diagnostics, not a complete epsilon-map algorithm.

Reports Section 5.1's local vertex/edge event test from Cabiddu & Attene,
https://arxiv.org/html/1704.08049 . No tetrahedral region growth, subdivision or
global coverage is claimed. Sections 5.2/5.3 have additional event conditions.
"""
from fractions import Fraction as Q
from collections import defaultdict
from .near_feature_witness import sub,dot,cross,point_triangle
from .segment_material_query import rational_point


def endpoint_link(query,center,endpoint,incident_faces):
    p,q=rational_point(center),rational_point(endpoint)
    radius2=dot(sub(p,q),sub(p,q))
    if radius2==0 or not incident_faces: raise ValueError('separated endpoints and full incidence required')
    ids=sorted(set(incident_faces))
    for i in ids:
        t,n=query._triangle(i)
        closest=point_triangle(q,t,n)
        if closest!=q: raise ValueError('endpoint not on listed triangle')
    vertices={int(v) for i in ids for v in query.faces[i]
              if rational_point(query.vertices[v])==q}
    if len(vertices)>1: raise ValueError('ambiguous coincident vertices')
    details={}
    if vertices:
        vertex=vertices.pop(); links=defaultdict(set)
        for i in ids:
            others=[int(v) for v in query.faces[i] if int(v)!=vertex]
            if len(others)!=2: raise ValueError('incomplete vertex incidence')
            a,b=others; links[a].add(b); links[b].add(a)
        if any(len(values)!=2 for values in links.values()): raise ValueError('vertex link not a cycle')
        start=min(links); order=[start]; previous=None; current=start
        while True:
            choices=sorted(links[current]-({previous} if previous is not None else set()))
            following=choices[0]
            if following==start: break
            if following in order: raise ValueError('nonmanifold vertex link')
            order.append(following); previous,current=current,following
        if len(order)!=len(links): raise ValueError('disconnected vertex link')
        # An edge from q enters the open ball iff its initial distance derivative
        # is negative, equivalently (neighbor-q) dot (p-q) > 0.
        internal=[dot(sub(rational_point(query.vertices[v]),q),sub(p,q))>0 for v in order]
        switches=sum(a!=b for a,b in zip(internal,internal[1:]+internal[:1]))
        event=(switches==0 or switches>2)
        kind='VERTEX'
        details={'vertex':vertex,'ordered_neighbor_vertices':order,
                 'incident_edge_internal':internal,'cyclic_switches':switches}
    elif len(ids)==2:
        common=set(map(int,query.faces[ids[0]]))&set(map(int,query.faces[ids[1]]))
        if len(common)!=2: raise ValueError('edge incidence mismatch')
        internal=[]
        for i in ids:
            t,n=query._triangle(i); nearest=point_triangle(p,t,n)
            distance=dot(sub(p,nearest),sub(p,nearest))
            if distance>radius2: raise ValueError('invalid incident face distance')
            internal.append(distance<radius2)
        event=internal[0]==internal[1]; kind='EDGE'
        details={'edge_vertices':sorted(common),'incident_face_internal':internal}
    elif len(ids)==1:
        t,n=query._triangle(ids[0])
        if not all(dot(cross(sub(t[(i+1)%3],t[i]),sub(q,t[i])),n)>0 for i in range(3)):
            raise ValueError('missing adjacent triangle at boundary edge')
        kind='FACE_INTERIOR'; event=True
    else: raise ValueError('unsupported endpoint feature')
    return {'feature':kind,'local_section_5_1_event':event,'details':details,
            'full_epsilon_map':'NOT_IMPLEMENTED','wall_thickness_admission':'NOT_GRANTED'}
