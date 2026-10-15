"""Read-only exact signed surface clipping; no optimizer or certificate grant.

For x<=0, omitted closing facets on x=0 have zero origin-based tetrahedral
volume and zero x first moment. Thus clipped oriented surface triangles suffice
for these two integrals. Input material containment remains a separate gate.
"""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import numpy as np
import _make_it_stand as native
import sys

CASE = Path(__file__).parent / 'asymmetric-candidate-v1'


def validate_diagnostic_assumptions(identity):
    target = identity['target']
    if (target['mode'] != 'SPINNING_TOP' or target['origin'] != [0., 0., 0.]
            or target['axis'] != [0., 0., 1.]):
        raise ValueError('capacity diagnostic requires frozen zero origin and z axis')
    for kind in ('fixed_regions', 'design_regions'):
        if len(identity[kind]) != 1 or identity[kind][0]['density'] != 1.:
            raise ValueError('capacity diagnostic requires one unit-density region per kind')


def require_partition_identity(actual, frozen):
    if actual != frozen:
        raise ValueError('capacity diagnostic partition differs from frozen admission')


def load(name):
    vertices, faces = [], []
    for line in (CASE / name).read_text().splitlines():
        fields = line.split()
        if fields and fields[0] == 'v':
            vertices.append(tuple(Q(float(x)) for x in fields[1:4]))
        if fields and fields[0] == 'f':
            faces.append(tuple(int(x.split('/')[0])-1 for x in fields[1:4]))
    # Fixed material legitimately includes an inward cavity component. This
    # orientation audit alone is not a material containment/partition proof.
    audit = native.audit_exact_triangle_boundary(
        np.array([[float(x) for x in p] for p in vertices]),
        np.array(faces, dtype=np.int32),
        'AnyNonzero' if name == 'fixed.obj' else 'OutwardPositive')
    if not audit['boundary_valid']:
        raise ValueError(f'{name}: exact boundary invalid')
    return vertices, faces


def clip_negative_x(polygon):
    result = []
    for a,b in zip(polygon, polygon[1:]+polygon[:1]):
        inside, next_inside = a[0]<=0, b[0]<=0
        if inside:
            result.append(a)
        if inside != next_inside:
            s = -a[0]/(b[0]-a[0])
            result.append(tuple(a[k]+s*(b[k]-a[k]) for k in range(3)))
    return result


def integrate(vertices, faces, negative=False):
    mass, lx = Q(0), Q(0)
    for face in faces:
        polygon = [vertices[i] for i in face]
        if negative:
            polygon = clip_negative_x(polygon)
        for j in range(1,len(polygon)-1):
            a,b,c = polygon[0],polygon[j],polygon[j+1]
            dv = (a[0]*(b[1]*c[2]-b[2]*c[1])+
                  a[1]*(b[2]*c[0]-b[0]*c[2])+
                  a[2]*(b[0]*c[1]-b[1]*c[0]))/6
            mass += dv
            lx += dv*(a[0]+b[0]+c[0])/4
    return mass,lx


def main():
    box = [tuple(map(Q,p)) for p in [(-1,-1,-1),(1,-1,-1),
        (1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]]
    faces = [(0,2,1),(0,3,2),(4,5,6),(4,6,7),(0,1,5),(0,5,4),
        (3,7,6),(3,6,2),(0,4,7),(0,7,3),(1,2,6),(1,6,5)]
    # Literal analytical calibrations, not expectations from the implementation.
    if integrate(box,faces,True) != (Q(4),Q(-2)):
        raise ValueError('negative half-box analytical calibration failed')
    shifted = [tuple(p[k]+(Q(1,2) if k==0 else 0) for k in range(3)) for p in box]
    if integrate(shifted,faces,True) != (Q(2),Q(-1,2)):
        raise ValueError('shifted half-box analytical calibration failed')
    positive = [(p[0]+2,p[1],p[2]) for p in box]
    if integrate(positive,faces,True) != (Q(0),Q(0)):
        raise ValueError('positive box analytical calibration failed')
    admitted = json.loads((CASE / 'input-admission.json').read_text())
    identity = admitted['problem_identity']
    validate_diagnostic_assumptions(identity)
    fixed = integrate(*load('fixed.obj'))
    negative = integrate(*load('design.obj'),negative=True)
    lower = fixed[1]+negative[1]
    partition_digest = None
    if '--verify-partition' in sys.argv:
        from make_it_stand.verification.spin_problem import SpinArtifactProblem
        def mesh(name):
            vertices,faces = load(name)
            return (np.array([[float(x) for x in p] for p in vertices]),
                    np.array(faces,dtype=np.int32))
        problem = SpinArtifactProblem(original=mesh('outer.obj'),
            fixed_regions=[('fixed-shell-and-pin',mesh('fixed.obj'),1.)],
            design_regions=[('heart-inner-excluding-pin',mesh('design.obj'),1.)],
            length_unit=identity['length_unit'], mass_unit=identity['mass_unit'],
            target=identity['target'], minimum_outer_wall=identity['minimum_outer_wall'])
        require_partition_identity(problem.input_digest, admitted['problem_digest'])
        partition_digest = problem.input_digest
    record = dict(schema='heart-lateral-capacity-diagnostic-v2',
        scope='EXACT_BINARY64_MODEL_COORDINATE_SURFACE_INTEGRALS',
        target_origin_model=[0,0,0], density_model=1,
        fixed_mass=float(fixed[0]), fixed_Lx=float(fixed[1]),
        negative_design_mass=float(negative[0]), negative_design_Lx=float(negative[1]),
        total_Lx_lower_bound=float(lower), exact_lower_bound=str(lower),
        strict_positive_lower_bound=lower>0,
        sealed_certificate='NOT_YET_ISSUED',
        rechecked_partition_digest=partition_digest,
        input_sources={str(CASE/name):hashlib.sha256((CASE/name).read_bytes()).hexdigest()
            for name in ('fixed.obj','design.obj')})
    print(json.dumps(record, indent=2, allow_nan=False), flush=True)
    if '--output' in sys.argv:
        dest = Path(sys.argv[sys.argv.index('--output')+1])
        with dest.open('x') as out:
            json.dump(record,out,indent=2,allow_nan=False)


if __name__ == '__main__':
    main()
