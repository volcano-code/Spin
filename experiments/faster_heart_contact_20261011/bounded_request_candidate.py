"""Explicit, globally consistent binary64 derivative; no solver admission.

Geometric bound: T + [-r,r]^3 is obtained by three orthogonal segment
extrusions. Cavalieri gives volume increment 2r times projected area; a planar
convex extrusion adds 2r times projected width. Thus its EXACT volume is
r*sum(abs(edge_cross)) + 4*r²*sum(axis_widths) + 8*r³.
The straight vertex homotopy stays in this tube for every matching triangle.
For two certified closed oriented solids, winding numbers agree outside the
union of tubes, so symmetric difference volume is at most the sum. This file
checks local indexed triangles only; global solid validity and partition must
be certified independently before applying the conditional entity bound.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
import json
import math
import os
from pathlib import Path
import gmpy2

from rational_request_view import boundary_rows, sha


def text(q):
    return str(gmpy2.mpz(q.numerator))+'/'+str(gmpy2.mpz(q.denominator))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def triangle_tube_volume(vertices, radius):
    if radius < 0: raise ValueError('NEGATIVE_RADIUS')
    a,b,c = vertices
    normal = cross(tuple(b[k]-a[k] for k in range(3)),tuple(c[k]-a[k] for k in range(3)))
    widths = sum(max(v[k] for v in vertices)-min(v[k] for v in vertices) for k in range(3))
    return radius*sum(abs(x) for x in normal)+4*radius**2*widths+8*radius**3


class GlobalVertexMap:
    def __init__(self):
        self.mapping = {}; self.reverse = {}; self.radius = Q(0)

    @property
    def distinct_authority_vertices(self): return len(self.mapping)

    def convert(self, row):
        if len(row) != 3: raise ValueError('INVALID_VERTEX')
        original = tuple(Q(x) for x in row)
        if tuple(text(x) for x in original) != tuple(row):
            raise ValueError('NONCANONICAL_COORDINATE')
        if original in self.mapping: return self.mapping[original]
        try:
            floating = tuple(float(x) for x in original)
            if not all(math.isfinite(x) for x in floating): raise OverflowError()
            derived = tuple(Q.from_float(x) for x in floating)
        except OverflowError as exc: raise ValueError('UNREPRESENTABLE_BINARY64_VERTEX') from exc
        if derived in self.reverse and self.reverse[derived] != original:
            raise ValueError('GLOBAL_VERTEX_COLLISION: '+str(tuple(row)))
        self.radius = max(self.radius, *(abs(a-b) for a,b in zip(original,derived)))
        self.reverse[derived] = original; self.mapping[original] = derived
        return derived


def check_local_triangles(original, derived, faces):
    if len(original) != len(derived): raise ValueError('CONNECTIVITY_VERTEX_MISMATCH')
    edges = defaultdict(list)
    for index, face in enumerate(faces):
        if len(face) != 3 or any(type(i) is not int or not 0 <= i < len(original) for i in face):
            raise ValueError('INVALID_FACE: '+str(index))
        normals = []
        for vertices in (original,derived):
            a,b,c = [vertices[i] for i in face]
            normal = cross(tuple(b[k]-a[k] for k in range(3)),tuple(c[k]-a[k] for k in range(3)))
            if not any(normal): raise ValueError('ZERO_AREA_FACE: '+str(index))
            normals.append(normal)
        if sum(a*b for a,b in zip(*normals)) <= 0:
            raise ValueError('FACE_ORIENTATION_REVERSED: '+str(index))
        for a,b in zip(face,face[1:]+face[:1]):
            edges[min(a,b),max(a,b)].append((a,b))
    closed = all(len(e)==2 and e[0]==e[1][::-1] for e in edges.values()) and bool(faces)
    if not closed: raise ValueError('NOT_CLOSED_ORIENTED_INDEXED_SURFACE')
    return dict(closed_oriented_indexed_surface=True, unchanged_connectivity=True,
        zero_area_faces=0, reversed_faces=0,
        status='PENDING_NATIVE_ENTITY_AND_PARTITION_CERTIFICATION',
        self_intersection='NOT_CHECKED_NATIVE_REQUIRED', solver_admission_granted=False)


def save(path, value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False)
        stream.flush(); os.fsync(stream.fileno())
    fd=os.open(path.parent,os.O_RDONLY | os.O_DIRECTORY)
    try: os.fsync(fd)
    finally: os.close(fd)


def build_candidate(admission, output):
    from rounded_input import exact_moments
    admission,output=Path(admission),Path(output)
    output.mkdir(parents=True,exist_ok=False)
    record_sha=sha(admission/'admission.json')
    record=json.loads((admission/'admission.json').read_text())
    mapping=GlobalVertexMap(); reports={}
    try:
        for name in ('outer','fixed','design'):
            source=admission/(name+'.rational.json'); expected=record['files'][source.name]
            if sha(source)!=expected: raise ValueError('STALE_AUTHORITY: '+name)
            original=[]; derived=[]; faces=[]
            for kind,row in boundary_rows(source):
                if kind=='vertices':
                    derived.append(mapping.convert(row)); original.append(tuple(Q(x) for x in row))
                else: faces.append(row)
            local=check_local_triangles(original,derived,faces)
            radius=max((abs(a-b) for p,q in zip(original,derived) for a,b in zip(p,q)),default=Q(0))
            bound=sum((triangle_tube_volume([original[i] for i in face],radius) for face in faces),Q(0))
            a=exact_moments([[text(x) for x in row] for row in original],faces)
            b=exact_moments([[text(x) for x in row] for row in derived],faces)
            reports[name]=dict(local, source_sha256=expected,
                vertex_count=len(original),face_count=len(faces),
                maximum_linf_vertex_shift_exact=text(radius),
                bidirectional_boundary_hausdorff_linf_upper_bound_exact=text(radius),
                boundary_hausdorff_l2_upper_bound_squared_exact=text(3*radius**2),
                conditional_symmetric_difference_volume_upper_bound_exact=text(bound),
                entity_bound_applied=False,
                independent_authority_moments_exact=[str(x) for x in a],
                independent_derived_moments_exact=[str(x) for x in b],
                independent_moment_delta_exact=[str(y-x) for x,y in zip(a,b)])
            save(output/(name+'.rational.json'),dict(schema='authoritative-rational-input-boundary-v1',
                obj_role='DERIVED_BINARY64_REQUEST_NOT_ORIGINAL_AUTHORITY',
                vertices=[[text(x) for x in row] for row in derived],faces=faces))
            save(output/(name+'.binary64.json'),dict(vertices=[[float(x) for x in row] for row in derived],faces=faces))
            save(output/(name+'-local-check.json'),reports[name])
            if sha(source)!=expected: raise ValueError('AUTHORITY_CHANGED_DURING_CHECK: '+name)
        if sha(admission/'admission.json')!=record_sha: raise ValueError('ADMISSION_CHANGED')
        report=dict(schema='globally-mapped-bounded-request-candidate-v1',
            source_admission_sha256=record_sha,source_problem_digest=record['problem_digest'],
            creator_sha256=sha(__file__),global_map_distinct_vertices=mapping.distinct_authority_vertices,
            global_maximum_linf_shift_exact=text(mapping.radius),regions=reports,
            status='PENDING_NATIVE_ENTITY_AND_PARTITION_CERTIFICATION',
            original_authority_unchanged=True,solver_admission_granted=False,
            geometry_admission_granted=False,optimization_invocations=0,simulation_invocations=0,
            required_checks=['NATIVE_VALID_SOLIDS','DERIVED_COMPLETE_PARTITION',
                'DESIGN_CONTAINMENT','FIXED_PIN_ROLE','ACTUAL_CLEARANCE',
                'BOUND_AND_MOMENT_TOLERANCE_ACCEPTANCE','DERIVED_REQUEST_IDENTITY'],
            files={p.name:sha(p) for p in output.iterdir() if p.is_file()})
        save(output/'candidate.json',report)
        return report
    except Exception as exc:
        save(output/'failure.json',dict(status='DERIVED_REQUEST_REJECTED',
            type=type(exc).__name__,message=str(exc),source_admission_sha256=record_sha,
            solver_admission_granted=False,optimization_invocations=0,simulation_invocations=0))
        raise


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--admission',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=build_candidate(args.admission,args.output)
    print(json.dumps(dict(status=result['status'],solver_admission_granted=False)))


if __name__=='__main__':main()
