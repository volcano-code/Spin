"""Independent necessary COM capacity bounds for frozen Flower input."""
from pathlib import Path
import json
import hashlib
import sys
from fractions import Fraction as F
from shapely.geometry import Polygon, box
from shapely.ops import unary_union
from paper_flower_fixture import exact_profile_moments, boundary_loops


def check(folder):
    record=json.loads((folder/'input.json').read_text())
    problem=json.loads((folder/'flower-cases.json').read_text())['cases'][0]['planar_problem']
    loops=problem['design_region_loops']
    signed=[]
    for loop in loops:
        p=Polygon(loop)
        signed.append((exact_profile_moments([loop])[0]>0,p))
    positive=unary_union([p for sign,p in signed if sign])
    negative=unary_union([p for sign,p in signed if not sign])
    design=positive.difference(negative)
    fixed=[F(x) for x in record['fixed_exact_rational_moments']]
    density=F(problem['design_areal_density'])
    xmin,ymin,xmax,ymax=design.bounds
    radius=2*max(abs(xmin),abs(ymin),abs(xmax),abs(ymax),1.)
    rows=[]
    for axis,minimum_box,maximum_box in [
            ('x',box(-radius,-radius,0,radius),box(0,-radius,radius,radius)),
            ('y',box(-radius,-radius,radius,0),box(-radius,0,radius,radius))]:
        i=1 if axis=='x' else 2
        low=fixed[i]+density*exact_profile_moments(boundary_loops(design.intersection(minimum_box)))[i]
        high=fixed[i]+density*exact_profile_moments(boundary_loops(design.intersection(maximum_box)))[i]
        rows.append(dict(axis=axis,lower_exact=str(low),upper_exact=str(high),
                         lower=float(low),upper=float(high),zero_within_capacity=low<=0<=high))
    return dict(schema='paper-flower-necessary-capacity-v1',
                input_sha256=hashlib.sha256((folder/'input.json').read_bytes()).hexdigest(),
                request_sha256=hashlib.sha256((folder/'flower-cases.json').read_bytes()).hexdigest(),
                directions=rows,feasibility='UNKNOWN' if all(r['zero_within_capacity'] for r in rows) else 'EXCLUDED_BY_PLANAR_CAPACITY',
                limitation='Necessary independent coordinate bounds only; not a joint feasibility witness. Shapely clipping vertices use binary64.',
                algorithm='retain every design point of negative/positive coordinate to bound arbitrary material subsets')


if __name__=='__main__':
    result=json.dumps(check(Path(sys.argv[1])),indent=2)+'\n'
    if len(sys.argv)>2:
        with Path(sys.argv[2]).open('x') as out:out.write(result)
    print(result)
