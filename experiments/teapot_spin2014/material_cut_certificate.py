"""Conditional local ball-removal topology witness on an audited material mesh.

For a regular ball with exactly two disjoint properly embedded boundary disks,
the disks divide the ball into two end regions and a middle region. An interior
material chord connecting the disks identifies the middle region as material.
Its intersection with the sphere is an annulus. Additivity gives
chi(material minus open ball) = chi(material)-1. This is a topology-change
witness, not a minimum wall thickness value or a general printability metric.
Embedding/material validity remain prerequisites of the independent mesh audit.
"""
from collections import defaultdict
from fractions import Fraction as Q
import math
import numpy as np
from .boundary_cap_query import BoundaryCapQuery
from .segment_material_query import rational_point
from .near_feature_witness import sub,dot,point_triangle,point_segment


def material_cut_witness(caps,center,threshold_squared):
    if not isinstance(caps,BoundaryCapQuery):raise TypeError('computed boundary cap query required')
    query=caps.boundary;p=rational_point(center)
    filtration=caps.audit(p,threshold_squared)
    base={'schema':'conditional-material-cut-topology-witness-v1',
          'center_exact':list(map(str,p)),'threshold_squared_exact':str(Q(threshold_squared)),
          'independent_embedded_material_audit_required':True,'wall_metric_adopted':False,
          'complete_wall_thickness':'UNKNOWN','wall_thickness_admission':'NOT_GRANTED'}
    first=filtration['first_non_disk_interval']
    if first is None:return {**base,'state':'NO_FIXED_CENTER_CAP_EVENT'}
    if len(first['components'])!=2 or any(c['euler']!=1 for c in first['components']):
        return {**base,'state':'UNSUPPORTED_CAP_TOPOLOGY','cap_event':first}
    squared=Q(first['witness_squared_radius_exact'])
    radius=Q(2)**((squared.numerator.bit_length()-squared.denominator.bit_length()+2)//2)
    while radius*radius<squared:radius*=2
    lo=np.array([math.nextafter(float(x-radius),-math.inf) for x in p])
    hi=np.array([math.nextafter(float(x+radius),math.inf) for x in p])
    candidates=np.flatnonzero(np.all((query.high>=lo)&(query.low<=hi),axis=1))
    active={};edge_candidates=set()
    for raw in candidates:
        i=int(raw);t,n=query._triangle(i);nearest=point_triangle(p,t,n)
        d2=dot(sub(nearest,p),sub(nearest,p))
        if d2<squared:
            active[i]=(d2,nearest)
            face=query.faces[i]
            edge_candidates.update(tuple(sorted((int(face[j]),int(face[(j+1)%3])))) for j in range(3))
    graph={i:set() for i in active};active_edges=[]
    for edge in edge_candidates:
        nearest=point_segment(p,caps.points[edge[0]],caps.points[edge[1]])
        if dot(sub(nearest,p),sub(nearest,p))<squared:
            a,b=caps.edges[edge]
            if a not in active or b not in active:raise ValueError('inconsistent active boundary')
            graph[a].add(b);graph[b].add(a);active_edges.append(edge)
    components=[];unseen=set(active)
    while unseen:
        pending=[min(unseen)];component=set()
        while pending:
            i=pending.pop()
            if i in component:continue
            component.add(i);pending.extend(graph[i]-component)
        unseen-=component;components.append(component)
    observations=[];center_components=[]
    for index,component in enumerate(components):
        vertices={int(v) for i in component for v in query.faces[i]}
        vc=sum(dot(sub(caps.points[v],p),sub(caps.points[v],p))<squared for v in vertices)
        ec=sum(caps.edges[e][0] in component for e in active_edges)
        observations.append({'representative_face':min(component),'faces':len(component),'edges':ec,'vertices':vc,
                             'euler':len(component)-ec+vc})
        if any(active[i][0]==0 for i in component):center_components.append(index)
    if observations!=first['components'] or len(center_components)!=1:
        raise ValueError('independent cap reconstruction mismatch')
    other=components[1-center_components[0]]
    target_face=min(other,key=lambda i:(active[i][0],i));_,target=active[target_face]
    chord=query.audit_segment(p,target)
    details={'cap_event':first,'ball_squared_radius_exact':str(squared),
             'boundary_disk_faces':[sorted(c) for c in components],
             'target_face':target_face,'target_exact':list(map(str,target)),'material_chord':chord}
    if chord['state']!='OPEN_SEGMENT_PARITY_INSIDE':
        return {**base,**details,'state':'NO_MATERIAL_CUT_WITNESS','cut_euler_change':None}
    return {**base,**details,'state':'MATERIAL_CUT_EULER_CHANGE_WITNESS','cut_euler_change':-1,
            'material_in_ball_euler':1,'material_on_sphere_euler':0,
            'interpretation':'Conditional on the bound embedded material audit; not a global wall metric.'}
