"""Multiprecision value formulas, no production/native derivatives imported.

Development reference only. Consumers must bind source recipes and validate
reconstructed geometry before comparing any differentiated material moments.
"""
import mpmath as mp


def ideal_edge(q, first, second):
    q=list(map(mp.mpf,q)); a=mp.matrix(first); delta=mp.matrix(second)-a
    matrix=mp.matrix([[q[4],q[9],q[8]],[q[9],q[5],q[7]],[q[8],q[7],q[6]]])
    inverse=matrix**-1
    quadratic=(delta.T*inverse*delta)[0]
    linear=2*(a.T*inverse*delta)[0]
    constant=(a.T*inverse*a)[0]
    if constant==0 or quadratic+linear+constant==0:
        raise ValueError('IDEAL_ENDPOINT_BRANCH')
    if quadratic==0:
        roots=[] if linear==0 else [-constant/linear]
    else:
        discriminant=linear**2-4*quadratic*constant
        if discriminant==0 and 0 < -linear/(2*quadratic) < 1:
            raise ValueError('TANGENT_IDEAL_BRANCH')
        roots=[] if discriminant<=0 else [(-linear+s*mp.sqrt(discriminant))/(2*quadratic) for s in (-1,1)]
    return [(t,list(inverse*(a+t*delta))) for t in sorted(roots) if 0<t<1]


def truncate_ray(origin, direction, radius):
    origin=mp.matrix(origin); direction=mp.matrix(direction); radius=mp.mpf(radius)
    aa=(direction.T*direction)[0]; bb=(origin.T*direction)[0]
    cc=(origin.T*origin)[0]-radius**2
    if radius<=0 or aa<=0 or cc>=0:
        raise ValueError('RAY_REQUIRES_INTERIOR_ORIGIN_AND_NONZERO_DIRECTION')
    # Solve the positive ray/sphere intersection directly at working precision.
    parameter=(-bb+mp.sqrt(bb**2-aa*cc))/aa
    return list(origin+parameter*direction)


def remote_cap(loop, radius):
    """Project a registered ordered loop onto its outward remote plane.

    This is a value formula, not a source topology constructor. A caller must
    independently establish the loop and that its cap is outside the domain.
    """
    vertices=[list(map(mp.mpf,point)) for point in loop]
    radius=mp.mpf(radius)
    if len(vertices)<3 or radius<=0 or any(len(point)!=3 for point in vertices):
        raise ValueError('INVALID_CAP_LOOP')
    area=[mp.mpf(0) for _ in range(3)]
    for a,b in zip(vertices,vertices[1:]+vertices[:1]):
        for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
            area[i]+=a[j]*b[k]-a[k]*b[j]
    norm=mp.sqrt(sum(x*x for x in area))
    if not norm:raise ValueError('DEGENERATE_CAP_LOOP')
    axis=[x/norm for x in area]
    side=sum(axis[i]*sum(point[i] for point in vertices) for i in range(3))
    if not side:raise ValueError('AMBIGUOUS_CAP_SIDE')
    if side<0:axis=[-x for x in axis]
    return [[point[i]+axis[i]*(2*radius-sum(a*x for a,x in zip(axis,point)))
             for i in range(3)] for point in vertices]


def projective_points(q, fixed_normals, normal_ids, rays, interface_count, faces, radius):
    """Evaluate a frozen construction recipe, not discover a new triangulation.

    Moving interval normals are recomputed from both ideal roots. This value
    function can be differentiated independently with multiprecision arithmetic.
    The consumer must check its output against the actual saved source mesh.
    """
    from independent_directional_reference import finite_point
    if len(normal_ids)!=len(rays) or not 0<=interface_count<=len(rays):
        raise ValueError('INCOMPLETE_SOURCE_RECIPES')
    moving={}
    for i in range(interface_count):
        if rays[i][0]==-3:
            nid=normal_ids[i]
            if nid<0 or (nid in moving and moving[nid]!=rays[i]):
                raise ValueError('CONFLICTING_MOVING_NORMAL')
            moving[nid]=rays[i]
    edges={};finite={}
    def edge(a,b):
        if not 0<=a<b<len(fixed_normals):raise ValueError('UNKNOWN_IDEAL_EDGE')
        if (a,b) not in edges:edges[a,b]=ideal_edge(q,fixed_normals[a],fixed_normals[b])
        return edges[a,b]
    def anchor(nid):
        if nid not in finite:
            if nid in moving:
                _,a,b,_=moving[nid]; roots=edge(a,b)
                if len(roots)!=2:raise ValueError('MOVING_INTERVAL_BRANCH_CHANGED')
                t=(roots[0][0]+roots[1][0])/2
                normal=[mp.mpf(x)+t*(mp.mpf(y)-mp.mpf(x)) for x,y in zip(fixed_normals[a],fixed_normals[b])]
            elif 0<=nid<len(fixed_normals):normal=fixed_normals[nid]
            else:raise ValueError('UNKNOWN_FINITE_NORMAL')
            finite[nid]=finite_point(list(map(mp.mpf,q)),normal)
        return finite[nid]
    points=[]
    for nid,recipe in zip(normal_ids[:interface_count],rays[:interface_count]):
        if nid>=0:
            if list(recipe)!=[-1]*4 and recipe[0]!=-3:raise ValueError('FINITE_RAY_CONFLICT')
            points.append(anchor(nid))
        else:
            nid,a,b,ordinal=recipe; roots=edge(a,b)
            if nid<0 or not 0<=ordinal<len(roots):raise ValueError('UNKNOWN_RAY_BRANCH')
            points.append(truncate_ray(anchor(nid),roots[ordinal][1],radius))
    loops={}
    for recipe in rays[interface_count:]:
        if recipe[0]==-2:
            _,source,loop,ordinal=recipe
            if not 0<=source<interface_count or loop<0 or ordinal<0:
                raise ValueError('INVALID_CAP_RECIPE')
            entries=loops.setdefault(loop,{})
            if ordinal in entries:raise ValueError('DUPLICATE_CAP_ORDINAL')
            entries[ordinal]=source
    caps={}
    for loop,entries in loops.items():
        if sorted(entries)!=list(range(len(entries))):raise ValueError('INCOMPLETE_CAP_LOOP')
        caps[loop]=remote_cap([points[entries[i]] for i in range(len(entries))],radius)
    for i in range(interface_count,len(rays)):
        recipe=rays[i]
        if normal_ids[i]!=-1:raise ValueError('CAP_WITH_FINITE_NORMAL')
        if recipe[0]==-2:points.append(caps[recipe[2]][recipe[3]])
        elif list(recipe)==[-1]*4:
            neighbors=sorted({v for face in faces if i in face for v in face if v!=i})
            if len(neighbors)<3 or any(not 0<=v<interface_count for v in neighbors):
                raise ValueError('UNSUPPORTED_CENTROID_CAP')
            points.append([sum(points[v][a] for v in neighbors)/len(neighbors) for a in range(3)])
        else:raise ValueError('UNKNOWN_CAP_RECIPE')
    return points
