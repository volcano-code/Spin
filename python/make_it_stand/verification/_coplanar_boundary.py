"""Propose exact-surface-preserving local changes to a rounded mesh.

Only strictly convex, consistently oriented coplanar two-triangle patches
within one labelled region qualify. Vertices are immutable; no welding,
deletion, approximate coplanarity or Newton-map replacement is performed.
The caller MUST independently audit the whole derived material afterwards.
Unsupported patches stay unchanged with explicit reasons.
"""
from fractions import Fraction as Q
import math


def _sub(a,b): return tuple(x-y for x,y in zip(a,b))
def _dot(a,b): return sum(x*y for x,y in zip(a,b))
def _cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def _normal(points, face):
    a,b,c=(points[i] for i in face)
    return _cross(_sub(b,a),_sub(c,a))
def _side(points,a,b,c,d):
    edge=_sub(points[b],points[a])
    return _dot(_cross(edge,_sub(points[c],points[a])),
                _cross(edge,_sub(points[d],points[a])))
def _edges(face):
    return tuple(zip(face,face[1:]+face[:1]))
def _boundary(faces):
    counts={}
    for face in faces:
        for a,b in _edges(face):
            edge=tuple(sorted((a,b)))
            counts[edge]=counts.get(edge,0)+(1 if a<b else -1)
    return {edge:n for edge,n in counts.items() if n}


def propose_coplanar_diagonal_repairs(exact_vertices, display_vertices, faces,
                                      rejected_pairs, *, face_regions,
                                      allow_positive_display_patch=False):
    """Return proposals only, never a passed boundary or material status.

    Exact coordinates must already be authenticated Fractions; displays must
    be their binary64 rounding. This low-level API cannot attest provenance.
    Rejected pairs are hints, not trusted geometric proof. Every local decision
    uses exact rational predicates on both authority and represented display.
    """
    def require(condition, message):
        if not condition: raise ValueError(message)
    sequence=lambda x:isinstance(x,(list,tuple))
    require(sequence(exact_vertices) and sequence(display_vertices) and
            len(exact_vertices)==len(display_vertices)>0,'invalid coordinate arrays')
    exact,display=[],[]
    for e,v in zip(exact_vertices,display_vertices):
        require(sequence(e) and len(e)==3 and all(type(x) is Q for x in e),
                'exact coordinates require authenticated Fractions')
        require(sequence(v) and len(v)==3 and all(type(x) is float and math.isfinite(x) for x in v),
                'display coordinates require finite binary64 values')
        try: rounded=tuple(float(x) for x in e)
        except OverflowError as error: raise ValueError('unrepresentable exact coordinate') from error
        require(rounded==tuple(v),'display differs from canonical rounding')
        exact.append(tuple(e));display.append(tuple(Q(x) for x in v))
    require(sequence(faces) and len(faces)>0,'missing triangle faces')
    triangles=[]
    for face in faces:
        require(sequence(face) and len(face)==3 and
                all(type(i) is int and 0<=i<len(exact) for i in face) and len(set(face))==3,
                'invalid triangle indices')
        triangles.append(tuple(face))
    require(sequence(face_regions) and len(face_regions)==len(triangles) and
            all(type(label) is str and label for label in face_regions),'invalid face region labels')
    require(sequence(rejected_pairs),'invalid rejected pair collection')
    pairs=set()
    for pair in rejected_pairs:
        require(sequence(pair) and len(pair)==2 and
                all(type(i) is int and 0<=i<len(triangles) for i in pair) and pair[0]!=pair[1],
                'invalid rejected pair indices')
        pairs.add(tuple(sorted(pair)))
    adjacency={}
    for index,face in enumerate(triangles):
        for edge in _edges(face):adjacency.setdefault(tuple(sorted(edge)),set()).add(index)
    output=list(triangles);patches=[];skipped=[];touched=set()
    for first,second in sorted(pairs):
        def skip(reason): skipped.append({'face_indices':[first,second],'reason':reason})
        if {first,second}&touched:
            skip('OVERLAPPING_PROPOSAL_REQUIRES_NEW_WHOLE_AUDIT');continue
        old=[output[first],output[second]]
        shared=set(old[0])&set(old[1])
        if len(shared)!=2:
            skip('NOT_A_SHARED_EDGE_PAIR');continue
        if face_regions[first]!=face_regions[second]:
            skip('REGION_INTERFACE');continue
        a,b=next((a,b) for a,b in _edges(old[0]) if a in shared and b in shared)
        c=next(i for i in old[0] if i not in shared)
        d=next(i for i in old[1] if i not in shared)
        if ((b,a) not in _edges(old[1]) or
                adjacency[tuple(sorted((a,b)))]!={first,second}):
            skip('INVALID_OR_NONMANIFOLD_SHARED_EDGE');continue
        if tuple(sorted((c,d))) in adjacency:
            skip('REPLACEMENT_DIAGONAL_ALREADY_EXISTS');continue
        n=_normal(exact,old[0])
        if (not any(n) or _dot(n,_sub(exact[d],exact[a]))!=0 or
                _side(exact,a,b,c,d)>=0 or _side(exact,c,d,a,b)>=0):
            skip('NOT_A_STRICTLY_CONVEX_EXACT_COPLANAR_PATCH');continue
        new=[(c,d,b),(c,a,d)]
        if (_boundary(old)!=_boundary(new) or
                any(_dot(n,_normal(exact,f))<=0 for f in old+new)):
            skip('ORIENTED_PATCH_BOUNDARY_MISMATCH');continue
        if not allow_positive_display_patch and all(_dot(n,_normal(display,f))>0 for f in old):
            skip('NO_LOCAL_ORIENTATION_LOSS');continue
        if any(_dot(n,_normal(display,f))<=0 for f in new):
            skip('NEW_DISPLAY_TRIANGLES_NOT_POSITIVELY_ORIENTED');continue
        for index,face in zip((first,second),old):
            for edge in _edges(face):
                key=tuple(sorted(edge));adjacency[key].remove(index)
                if not adjacency[key]:del adjacency[key]
        for index,face in zip((first,second),new):
            output[index]=face
            for edge in _edges(face):adjacency.setdefault(tuple(sorted(edge)),set()).add(index)
        touched.update((first,second))
        patches.append({'face_indices':[first,second],'old_faces':[list(f) for f in old],
            'new_faces':[list(f) for f in new], 'region':face_regions[first],
            'exact_relation':'EQUAL_ORIENTED_CONVEX_COPLANAR_QUAD',
            'old_exact_side_dot':str(_side(exact,a,b,c,d)),
            'new_exact_side_dot':str(_side(exact,c,d,a,b))})
    return {'status':'PROPOSED_REQUIRES_FULL_MATERIAL_AUDIT',
            'scope':'EXACT_LOCAL_PATCH_EQUALITY_NOT_FULL_BOUNDARY_OR_MATERIAL_APPROVAL',
            'faces':output,'patches':patches,'skipped_pairs':skipped,
            'vertex_positions_changed':False}
