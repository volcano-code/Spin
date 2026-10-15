"""Reuse independently integrated fixed moments for the existing z-axis bound."""
import gmpy2
from rounded_input import _q,_half_surface_first

def capacity_from_fixed_moments(fixed,design_vertices,design_faces,origin):
    if len(fixed)!=10 or len(origin)!=3:raise ValueError('CAPACITY_FRAME_DIMENSIONS')
    fixed=[_q(v) for v in fixed];origin=[_q(v) for v in origin]
    if fixed[0]<=0:raise ValueError('CAPACITY_REQUIRES_POSITIVE_FIXED_MASS')
    points=[tuple(_q(v)-origin[i] for i,v in enumerate(row)) for row in design_vertices]
    result=dict(schema='saved-fixed-moments-first-moment-capacity-v1',
        scope='NECESSARY_CONTINUOUS_SUBSET_BOUND_GIVEN_VALID_INPUT_PARTITION',
        coordinate_frame='TARGET_ORIGIN_TRANSLATED_MODEL',target_origin=list(map(float,origin)),
        target_axis=[0.,0.,1.],density_model=1.,feasible_witness=False,sealed_certificate='NOT_ISSUED')
    excluded=False
    for axis,label in enumerate(('x','y')):
        first=fixed[axis+1]-fixed[0]*origin[axis]
        low=_half_surface_first(points,design_faces,axis,-1)
        high=_half_surface_first(points,design_faces,axis,1)
        total_low,total_high=first+low,first+high
        rejected=total_low>0 or total_high<0;excluded|=rejected
        result[label]=dict(fixed_first_exact=str(first),design_min_exact=str(low),design_max_exact=str(high),
            total_min_exact=str(total_low),total_max_exact=str(total_high),
            total_min=float(total_low),total_max=float(total_high),excludes_zero=rejected)
    result['status']='EXCLUDED_BY_FIRST_MOMENT_BOUND' if excluded else 'UNKNOWN_NOT_EXCLUDED_BY_FIRST_MOMENT_BOUNDS'
    return result
