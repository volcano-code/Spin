"""Conservative exact full-boundary cover certificate, not finite sampling.

Each convex covered cell lies in a source triangle and maps affinely into
one target triangle along the source normal. The squared norm of that affine
displacement is convex, hence its maximum on a convex cell occurs at a vertex.
Exact half-plane subtraction proves coverage of the complete source triangle.
Positive-area coverage plus closed target triangles covers the limiting edges.
Both directions are required. Missing coverage is UNKNOWN, never small error.
"""
import gmpy2
from admit_bounded_request import cross,text


def _dot(a,b):return sum(x*y for x,y in zip(a,b))
def _minus(a,b):return tuple(x-y for x,y in zip(a,b))
def _normal(points):return cross(_minus(points[1],points[0]),_minus(points[2],points[0]))
def _area(poly):
    return len(poly)>=3 and any(any(cross(_minus(poly[i],poly[0]),_minus(poly[i+1],poly[0])))
                               for i in range(1,len(poly)-1))
def _edge_value(p,a,b,drop):
    axes=[k for k in range(3) if k!=drop];i,j=axes
    return (b[i]-a[i])*(p[j]-a[j])-(b[j]-a[j])*(p[i]-a[i])


def _clip(poly,a,b,drop,sign):
    result=[]
    if not poly:return result
    last=poly[-1];last_value=sign*_edge_value(last,a,b,drop)
    for p in poly:
        value=sign*_edge_value(p,a,b,drop)
        if (value>=0)!=(last_value>=0):
            t=last_value/(last_value-value)
            result.append(tuple(last[k]+t*(p[k]-last[k]) for k in range(3)))
        if value>=0 and (not result or result[-1]!=p):result.append(p)
        last,last_value=p,value
    if len(result)>1 and result[0]==result[-1]:result.pop()
    return result


def _subtract(poly,triangle,drop):
    inside=poly;outside=[]
    for a,b in zip(triangle,triangle[1:]+triangle[:1]):
        piece=_clip(inside,a,b,drop,-1)
        if _area(piece):outside.append(piece)
        inside=_clip(inside,a,b,drop,1)
        if not _area(inside):return [],outside
    return inside,outside


def _triangles(packet):
    points=[tuple(gmpy2.mpq(x) for x in row) for row in packet['vertices']]
    triangles=[]
    for face in packet['faces']:
        if len(face)!=3 or any(type(i)is not int or not 0<=i<len(points) for i in face):
            raise ValueError('INVALID_CERTIFICATE_FACE')
        triangle=[points[i] for i in face];normal=_normal(triangle)
        if not any(normal):raise ValueError('DEGENERATE_CERTIFICATE_TRIANGLE')
        triangles.append((triangle,normal,[(min(p[k] for p in triangle),max(p[k] for p in triangle)) for k in range(3)]))
    return triangles


def _directed(source,target,epsilon,maximum_cover_cells):
    covered=[];missing=[];maximum=gmpy2.mpq(0);cells=0
    for source_index,(triangle,normal,box) in enumerate(source):
        remaining=[triangle];drop=max(range(3),key=lambda k:abs(normal[k]));nn=_dot(normal,normal)
        for target_index,(target_triangle,target_normal,target_box) in enumerate(target):
            if not remaining:break
            if any(a[1]+epsilon<b[0] or b[1]+epsilon<a[0] for a,b in zip(box,target_box)):continue
            denominator=_dot(target_normal,normal)
            if not denominator:continue
            projected=[tuple(p[k]-normal[k]*_dot(normal,_minus(p,triangle[0]))/nn for k in range(3))
                       for p in target_triangle]
            if not _area(projected):continue
            if _edge_value(projected[2],projected[0],projected[1],drop)<0:projected.reverse()
            next_remaining=[]
            for poly in remaining:
                intersection,outside=_subtract(poly,projected,drop)
                if not intersection:
                    next_remaining.append(poly);continue
                values=[(_dot(target_normal,_minus(target_triangle[0],p))/denominator)**2*nn for p in intersection]
                bound=max(values)
                if bound>epsilon**2:
                    next_remaining.append(poly);continue
                cells+=1
                if cells>maximum_cover_cells:
                    return dict(status='UNKNOWN',reason='EXACT_COVER_CELL_BUDGET',covered_cells=covered,
                                uncovered_source_faces=[source_index],maximum_squared_bound_exact=text(maximum))
                maximum=max(maximum,bound)
                covered.append(dict(source_face=source_index,target_face=target_index,
                    source_cell_vertices_exact=[[text(x) for x in p] for p in intersection],
                    distance_squared_upper_bound_exact=text(bound)))
                next_remaining.extend(outside)
            remaining=next_remaining
        if remaining:
            missing.append(dict(source_face=source_index,
                uncovered_convex_cells_exact=[[[text(x) for x in p] for p in poly] for poly in remaining]))
    return dict(status='CERTIFIED' if source and not missing else 'UNKNOWN',covered_cells=covered,
                uncovered_source_faces=missing,maximum_squared_bound_exact=text(maximum))


def boundary_distance_certificate(authority,derived,maximum_distance,*,maximum_cover_cells=200000):
    epsilon=gmpy2.mpq(maximum_distance)
    if epsilon<=0:raise ValueError('NONPOSITIVE_BOUNDARY_BOUND')
    a,b=_triangles(authority),_triangles(derived)
    forward=_directed(a,b,epsilon,maximum_cover_cells)
    reverse=_directed(b,a,epsilon,maximum_cover_cells)
    valid=forward['status']==reverse['status']=='CERTIFIED'
    bound=max(gmpy2.mpq(forward['maximum_squared_bound_exact']),gmpy2.mpq(reverse['maximum_squared_bound_exact']))
    return dict(schema='bounded-pin-full-boundary-cover-v1',status='CERTIFIED' if valid else 'UNKNOWN',
        coverage_scope='ALL_TRIANGLE_POINTS_VIA_EXACT_CONVEX_COVER',forward=forward,reverse=reverse,
        bidirectional_distance_squared_upper_bound_exact=text(bound) if valid else None,
        requested_distance_upper_bound_exact=text(epsilon),geometry_admission_granted=False,
        scope='BOUNDARY_DISTANCE_ONLY_REQUIRES_INDEPENDENT_NATIVE_SOLID_AND_SUBSET_CHECKS')


def numerical_pin_gate(boundary,authority_moments,derived_moments,difference_volume,length):
    """Necessary numeric gates only; this helper never grants native admission."""
    length=gmpy2.mpq(length);volume=gmpy2.mpq(difference_volume);tol=gmpy2.mpq(1,10**6)
    if length<=0 or len(authority_moments)!=10 or len(derived_moments)!=10:
        raise ValueError('INVALID_PIN_MOMENT_OR_LENGTH')
    a,b=[[gmpy2.mpq(x) for x in row] for row in (authority_moments,derived_moments)]
    if a[0]<=0 or b[0]<=0 or volume<0 or a[0]-b[0]!=volume:
        raise ValueError('PIN_DIFFERENCE_VOLUME_INDEPENDENT_MOMENT_MISMATCH')
    errors=[abs(x-y)/length**(3 if i==0 else 4 if i<4 else 5) for i,(x,y) in enumerate(zip(a,b))]
    numeric=volume<=tol*length**3 and all(x<=tol for x in errors)
    full=boundary.get('status')=='CERTIFIED'
    distance=gmpy2.mpq(boundary['bidirectional_distance_squared_upper_bound_exact']) if full else None
    if full and distance>tol**2*length**2:numeric=False
    return dict(schema='bounded-pin-numerical-gate-v1',
        status='CERTIFIED' if full and numeric else 'REJECTED' if not numeric else 'UNKNOWN',
        normalized_moment_errors_exact=[text(x) for x in errors],
        normalized_removed_volume_exact=text(volume/length**3),
        authority_normalization_length_exact=text(length),tolerance_exact=text(tol),
        geometry_admission_granted=False,solver_admission_granted=False)
