"""One rounded planar arrangement; exact intersections and global error bounds.

The unchanged analytic circle remains authority. This module defines a distinct
bounded polyhedral derivative, not an exact representation of that circle.
"""
from fractions import Fraction as F


def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def segment_crossing(a,b,c,d):
    u=sub(b,a);v=sub(d,c);den=cross(u,v)
    if den==0:return None
    t=cross(sub(c,a),v)/den;s=cross(sub(c,a),u)/den
    if 0<=t<=1 and 0<=s<=1:return (a[0]+t*u[0],a[1]+t*u[1])
    return None


def circle_chord_bound(points,center,radius):
    """Radial Hausdorff bound via exact quadratic extrema on every chord.

With positive angular ordering and a center-containing convex polygon, every
radial ray hits one chord. |sqrt(q)-R| <= |q-R²|/R applies along its whole length.
    """
    if radius<=0:raise ValueError('INVALID_CIRCLE_RADIUS')
    bound=F(0)
    for p,q in zip(points,points[1:]+points[:1]):
        a=sub(p,center);d=sub(q,p);dd=d[0]**2+d[1]**2
        if dd==0:raise ValueError('COLLAPSED_CIRCLE_EDGE')
        if cross(a,sub(q,center))<=0:raise ValueError('NONPOSITIVE_CIRCLE_ANGULAR_ORDER')
        t=max(F(0),min(F(1),-(a[0]*d[0]+a[1]*d[1])/dd))
        values=[sum(v*v for v in a),sum(v*v for v in sub(q,center)),
                (a[0]+t*d[0])**2+(a[1]+t*d[1])**2]
        bound=max(bound,max(abs(v-radius**2)/radius for v in values))
    return bound


def _loops(geometry):
    from paper_flower_fixture import boundary_loops
    return [[tuple(F(float(x)) for x in p) for p in loop] for loop in boundary_loops(geometry)]


def _inside(point,loops):
    # Signed winding over all CCW outer / CW hole loops, in exact arithmetic.
    winding=0;x,y=point
    for loop in loops:
        for a,b in zip(loop,loop[1:]+loop[:1]):
            if a[1]<=y<b[1]:
                if cross(sub(b,a),sub(point,a))>0:winding+=1
            elif b[1]<=y<a[1]:
                if cross(sub(b,a),sub(point,a))<0:winding-=1
    return winding!=0


def shared_partition(original,design,cutter,length):
    from shapely.geometry import LineString,Polygon
    from shapely.ops import polygonize,unary_union
    from shapely.strtree import STRtree
    from paper_flower_output import profile_geometry
    originals=_loops(original);designs=_loops(design);cuts=_loops(cutter)
    if len(cuts)!=1:raise ValueError('CUTTER_MUST_BE_SINGLE_RING')
    de=[(a,b) for loop in designs for a,b in zip(loop,loop[1:]+loop[:1])]
    oe=[(a,b) for loop in originals for a,b in zip(loop,loop[1:]+loop[:1])]
    same_source_design=set(oe)==set(de)
    ce=[(a,b) for loop in cuts for a,b in zip(loop,loop[1:]+loop[:1])]
    dn=[{a,b} for a,b in de];cn=[{a,b} for a,b in ce]
    tree=STRtree([LineString([tuple(map(float,a)),tuple(map(float,b))]) for a,b in ce])
    intersections=[]
    for i,(a,b) in enumerate(de):
        for j in tree.query(LineString([tuple(map(float,a)),tuple(map(float,b))])):
            p=segment_crossing(a,b,*ce[int(j)])
            if p is not None:
                dn[i].add(p);cn[int(j)].add(p)
                intersections.append(dict(design_edge=i,circle_edge=int(j),point=[str(x) for x in p]))
    def ordered(edge,nodes):
        a,b=edge;axis=0 if a[0]!=b[0] else 1
        return sorted(nodes,key=lambda p:(p[axis]-a[axis])/(b[axis]-a[axis]))
    dchains=[ordered(e,n) for e,n in zip(de,dn)]
    cchains=[ordered(e,n) for e,n in zip(ce,cn)]
    scale=F(float(length));mapping={};rounding=F(0)
    def convert(p):
        nonlocal rounding
        if p not in mapping:
            exact=tuple(v*scale for v in p);stored=tuple(float(v) for v in exact)
            rounding=max(rounding,sum(abs(F(v)-x) for v,x in zip(stored,exact)))
            mapping[p]=stored
        return mapping[p]
    # Original and design/rim share their transformed source vertices; each
    # new cut vertex is converted once and reused by all adjacent cells.
    source_segments=[]
    if not same_source_design:
        source_segments.extend(oe)
    for chain in dchains:source_segments.extend(zip(chain,chain[1:]))
    for chain in cchains:
        for a,b in zip(chain,chain[1:]):
            if _inside(((a[0]+b[0])/2,(a[1]+b[1])/2),designs):source_segments.append((a,b))
    edges={tuple(sorted((convert(a),convert(b)))) for a,b in source_segments}
    for loop in cuts:
        for p in loop:convert(p)
    if len(set(mapping.values()))!=len(mapping):raise ValueError('DISTINCT_AUTHORITY_VERTICES_COLLAPSE')
    # Confirm rounding did not introduce crossings that would require a
    # second noding operation and therefore invalidate the shared graph.
    lines=[LineString(edge) for edge in sorted(edges)];tree2=STRtree(lines)
    for i,line in enumerate(lines):
        a,b=[tuple(F(float(x)) for x in p) for p in line.coords]
        for j in tree2.query(line):
            if int(j)<=i:continue
            c,d=[tuple(F(float(x)) for x in p) for p in lines[int(j)].coords]
            p=segment_crossing(a,b,c,d)
            if p is not None and not (p in (a,b) and p in (c,d)):
                raise ValueError('ROUNDED_ARRANGEMENT_INTRODUCED_CROSSING')
    nodal_design=[];offset=0
    for loop in designs:
        ring=[]
        for _ in loop:ring.extend(convert(p) for p in dchains[offset][:-1]);offset+=1
        nodal_design.append(ring)
    derived_design=profile_geometry(nodal_design)
    derived_original=(derived_design if same_source_design else
                      profile_geometry([[convert(p) for p in loop] for loop in originals]))
    derived_cutter=Polygon([convert(p) for p in cuts[0]])
    buckets={key:[] for key in ['rim','retained','cavity']}
    vertices=set(mapping.values())
    for face in polygonize(lines):
        for ring in [face.exterior,*face.interiors]:
            if any(tuple(p) not in vertices for p in ring.coords):raise ValueError('POLYGONIZE_ADDED_VERTEX')
        p=face.representative_point()
        if not derived_original.contains(p):continue
        if not derived_design.contains(p):key='rim'
        else:key='cavity' if derived_cutter.contains(p) else 'retained'
        buckets[key].append(face)
    result={k:unary_union(v) for k,v in buckets.items()}
    result.update(design=derived_design,original=derived_original)
    # Producer diagnostics are not a final native solid certification.
    if result['design'].symmetric_difference(result['retained'].union(result['cavity'])).area!=0:
        raise ValueError('SHARED_EDITABLE_PLANAR_PARTITION_FAILED')
    if result['original'].symmetric_difference(result['design'].union(result['rim'])).area!=0:
        raise ValueError('SHARED_SOURCE_PLANAR_PARTITION_FAILED')
    authority_edges=de+[(a,b) for loop in originals for a,b in zip(loop,loop[1:]+loop[:1])]
    perimeter=sum(abs(a[0]-b[0])+abs(a[1]-b[1]) for a,b in authority_edges)*scale
    # Swept edge tubes bound the polygon interface perturbation. Every edge
    # uses a fixed endpoint homotopy, with displacement bounded by rounding.
    tube=2*rounding*perimeter+F(355,113)*rounding**2*len(authority_edges)
    result['certificate']=dict(schema='shared-rational-arrangement-v2',
        source_and_design_share_noded_boundary=same_source_design,
        authority_design_loops_model=[[[str(x) for x in p] for p in loop] for loop in designs],
        authority_original_loops_model=[[[str(x) for x in p] for p in loop] for loop in originals],
        edge_index_scope='AUTHENTICATED_ORDERED_AUTHORITY_LOOPS_NOT_MANIFEST_ARRAY_OFFSETS',
        intersections=intersections,maximum_coordinate_rounding_l1_m=str(rounding),
        design_boundary_tube_area_bound_m2=str(tube),
        boundary_tube_scope='SOURCE_AND_DESIGN_ALL_EDGES_L1_LENGTH_UPPER_BOUND',
        noded_design_loops_si=nodal_design,
        vertices=[dict(authority_model=[str(x) for x in p],stored_si=list(v)) for p,v in sorted(mapping.items())],
        curve_bound_scope='EXACT_CHORD_QUADRATIC_EXTREMA_AND_ENDPOINT_HOMOTOPY_NOT_SAMPLING',
        native_partition_status='PENDING')
    return result
