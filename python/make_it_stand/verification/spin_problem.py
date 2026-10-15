"""Immutable indexed geometry and locked input for hollow-only material checking."""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
import math
import re
import weakref

import numpy as np

from ..spin import SPIN_MODES
from .spin_common import (SCOPE, TOLERANCE, byte_digest, canonical, digest, keys,
                          native, parse, require, scalar, text_label, vector)
from .rational_text import fraction_from_rational_text, fraction_to_rational_text


def _mesh_input_array(value, *, indices=False):
    """Freeze sequence scalars before common dtype inference can erase types.

    Existing arrays retain their supplied dtype for native lossless validation.
    Sequences need scalar-level checks first: np.array on a mixed sequence can
    erase bools, round integers, or promote signed/unsigned indices to floats.
    The staging arrays below contain only values already checked losslessly.
    """
    if isinstance(value, np.ndarray):
        return np.array(value, copy=True, order='C')
    label = 'triangle indices' if indices else 'vertices'
    require(isinstance(value, (list, tuple)), label + ': expected array, list or tuple')
    checked = []
    for row in tuple(value):
        require(isinstance(row, (list, tuple, np.ndarray)), label + ': expected rows of three scalars')
        if isinstance(row, np.ndarray):
            kinds = ('i', 'u') if indices else ('i', 'u', 'f')
            require(row.ndim == 1 and row.dtype.kind in kinds and row.dtype.itemsize <= 8,
                    label + ': unsupported row array shape or dtype')
        values = tuple(row)
        require(len(values) == 3, label + ': expected rows of three scalars')
        if indices:
            require(all(isinstance(x, (int, np.integer)) and not isinstance(x, (bool, np.bool_))
                        for x in values), label + ': expected integer scalars; bool is invalid')
            converted = tuple(int(x) for x in values)
            require(all(-(2**31) <= x < 2**31 for x in converted), label + ': exceeds signed32 range')
        else:
            converted = tuple(scalar(x, 'vertex coordinate') for x in values)
        checked.append(converted)
    return np.array(checked, dtype=np.int64 if indices else np.float64, order='C')


@dataclass(frozen=True, init=False)
class MeshSnapshot:
    """Canonical little-endian binary64/int32 bytes, without repair or reordering.

    None means an explicit empty set. All other inputs must be valid, nonempty,
    embedded signed regions. Native lossless parsing precedes conversion.
    """
    _vertices: bytes
    _faces: bytes
    vertex_count: int
    face_count: int
    explicit_empty: bool
    _exact_region: object = field(compare=False, repr=False)
    _rational_vertices_json: bytes | None = field(repr=False)

    def __init__(self, mesh, *, rational_vertices=None):
        authenticated_region = None
        authenticated_rationals = None
        if isinstance(mesh, MeshSnapshot):
            authenticated_region = mesh._exact_region
            authenticated_rationals = mesh.rational_vertices
            mesh = None if mesh.explicit_empty else mesh.arrays()
        if mesh is None:
            region = (authenticated_region if authenticated_region is not None
                      else native().empty_exact_material_region())
            v, f = np.empty((0, 3), dtype='<f8'), np.empty((0, 3), dtype='<i4')
            exact_json = None
        else:
            require(isinstance(mesh, (tuple, list)) and len(mesh) == 2,
                    'mesh must be a (vertices, triangles) pair or explicit None')
            # Freeze and validate supplied scalar types before sequence inference;
            # then validate both actual staging arrays with the native checker.
            v = _mesh_input_array(mesh[0])
            f = _mesh_input_array(mesh[1], indices=True)
            rational_vertices = (authenticated_rationals if authenticated_rationals is not None
                                 else rational_vertices)
            exact_json = None
            if rational_vertices is not None:
                exact_rows = _canonical_rational_vertices(rational_vertices)
                require(len(exact_rows) == len(v), 'rational vertex count mismatch')
                exact_float = np.asarray([[float(fraction_from_rational_text(x)) for x in row] for row in exact_rows], dtype='<f8')
                require(exact_float.tobytes(order='C') == np.asarray(v, dtype='<f8', order='C').tobytes(order='C'),
                        'rational coordinates do not match their binary64 derived mesh')
                region = (authenticated_region if authenticated_region is not None else
                          native().make_exact_material_region_from_rationals(exact_rows, f))
                exact_json = canonical([list(row) for row in exact_rows]).encode('ascii')
            else:
                region = (authenticated_region if authenticated_region is not None
                          else native().make_exact_material_region(v, f))
            v, f = np.asarray(v, dtype='<f8', order='C'), np.asarray(f, dtype='<i4', order='C')
        require(region.facts['completed'] and region.facts['valid'],
                'invalid material boundary: ' + region.facts['error'])
        object.__setattr__(self, '_vertices', v.tobytes(order='C'))
        object.__setattr__(self, '_faces', f.tobytes(order='C'))
        object.__setattr__(self, 'vertex_count', len(v))
        object.__setattr__(self, 'face_count', len(f))
        object.__setattr__(self, 'explicit_empty', mesh is None)
        object.__setattr__(self, '_exact_region', region)

        object.__setattr__(self, '_rational_vertices_json', exact_json)
    @classmethod
    def from_rationals(cls, vertices, faces):
        """Create a snapshot with authoritative exact rational coordinates."""
        exact_rows = _canonical_rational_vertices(vertices)
        face_array = _mesh_input_array(faces, indices=True)
        display_vertices = np.asarray([[float(fraction_from_rational_text(x)) for x in row] for row in exact_rows],
                                      dtype='<f8')
        return cls((display_vertices, face_array), rational_vertices=exact_rows)

    @classmethod
    def from_exact_region(cls, region):
        """Export the actual immutable owned exact region without a second Nef.

        No caller geometry or passed flags are accepted. Boundary coordinates
        and faces come from the same native object that owns the region proof.
        Final partition, mass and saved-file checks remain independent.
        """
        require(type(region) is native().ExactMaterialRegion,
                'expected an owned native ExactMaterialRegion')
        facts = region.facts
        require(facts['completed'] and facts['valid'], 'invalid owned exact material region')
        if facts['explicit_empty']:
            return cls(None)
        boundary = region.rational_boundary
        rows = _canonical_rational_vertices(boundary['vertices'])
        faces = np.asarray(_mesh_input_array(boundary['faces'], indices=True), dtype='<i4', order='C')
        policy = getattr(region, 'boundary_coordinate_policy', 'RATIONAL_INDEXED')
        require(policy in ('BINARY64_INDEXED', 'RATIONAL_INDEXED'), 'unknown owned coordinate policy')
        vertices = (np.asarray(region.binary64_boundary_vertices, dtype='<f8', order='C')
                    if policy == 'BINARY64_INDEXED' else
                    np.asarray([[float(fraction_from_rational_text(x)) for x in row] for row in rows], dtype='<f8'))
        require(np.isfinite(vertices).all(), 'owned exact region display is nonfinite')
        require(region.matches_rational_boundary(rows, faces), 'owned boundary export mismatch')
        obj = object.__new__(cls)
        attrs = dict(_vertices=vertices.tobytes(order='C'), _faces=faces.tobytes(order='C'),
                     vertex_count=len(vertices), face_count=len(faces), explicit_empty=False,
                     _exact_region=region,
                     _rational_vertices_json=(canonical([list(row) for row in rows]).encode('ascii')
                         if policy == 'RATIONAL_INDEXED' else None))
        for name, value in attrs.items():
            object.__setattr__(obj, name, value)
        return obj

    def arrays(self):
        return (np.frombuffer(self._vertices, dtype='<f8').reshape((-1, 3)).copy(),
                np.frombuffer(self._faces, dtype='<i4').reshape((-1, 3)).copy())

    @property
    def rational_vertices(self):
        if self._rational_vertices_json is None:
            return None
        return tuple(tuple(row) for row in parse(self._rational_vertices_json.decode('ascii')))

    def exact_region(self):
        out = self._exact_region
        require(out.facts['completed'] and out.facts['valid'], 'invalid authenticated region')
        return out

    @property
    def identity(self):
        result = {'explicit_empty': self.explicit_empty,
                  'vertices': {'dtype': '<f8', 'shape': [self.vertex_count, 3],
                               'sha256': byte_digest(self._vertices)},
                  'faces': {'dtype': '<i4', 'shape': [self.face_count, 3],
                            'sha256': byte_digest(self._faces)}}
        if self._rational_vertices_json is not None:
            result['rational_vertices'] = {
                'format': 'CANONICAL_FRACTION_V1',
                'sha256': byte_digest(self._rational_vertices_json)}
        return result


def _canonical_rational_vertices(vertices):
    require(isinstance(vertices, (list, tuple)) and len(vertices) > 0,
            'rational vertices must be a nonempty sequence')
    rows = []
    token = re.compile(r'-?(?:0|[1-9][0-9]*)/(?:[1-9][0-9]*)\Z')
    for row in vertices:
        require(isinstance(row, (list, tuple)) and len(row) == 3,
                'rational vertex must have three coordinates')
        converted = []
        for value in row:
            require(isinstance(value, str) and token.fullmatch(value),
                    'rational coordinate must be numerator/positive-denominator text')
            number = fraction_from_rational_text(value)
            converted.append(fraction_to_rational_text(number))
        rows.append(tuple(converted))
    return tuple(rows)
class SolverBoundarySnapshot:
    """Exact single-shell boundary used for intermediate Faster trial maps.

    This proves a closed, consistently oriented, exactly self-intersection-free
    represented boundary. It intentionally does not construct a Nef solid and
    is not accepted by sealed-artifact code. Multi-shell candidates must use a
    full :class:`MeshSnapshot` so nesting and material-region semantics are
    checked. Final outputs always use ``MeshSnapshot`` independently.
    """
    _vertices: bytes
    _faces: bytes
    vertex_count: int
    face_count: int
    boundary_digest: str

    def __init__(self, mesh):
        require(isinstance(mesh, (tuple, list)) and len(mesh) == 2,
                'solver boundary must be a (vertices, triangles) pair')
        vertices = _mesh_input_array(mesh[0])
        faces = _mesh_input_array(mesh[1], indices=True)
        facts = native().audit_exact_triangle_boundary(vertices, faces, 'AnyNonzero')
        require(facts['boundary_valid'],
                'invalid exact solver candidate boundary: ' + canonical(facts))
        signs = tuple(int(value) for value in facts['component_volume_signs'])
        require(signs == (1,),
                'solver boundary snapshot requires one outward-oriented shell')
        self._store(mesh[0], mesh[1], facts['input_digest'])

    def _store(self, vertices, faces, boundary_digest):
        vertices = _mesh_input_array(vertices)
        faces = _mesh_input_array(faces, indices=True)
        vertices = np.asarray(vertices, dtype='<f8', order='C')
        faces = np.asarray(faces, dtype='<i4', order='C')
        object.__setattr__(self, '_vertices', vertices.tobytes(order='C'))
        object.__setattr__(self, '_faces', faces.tobytes(order='C'))
        object.__setattr__(self, 'vertex_count', len(vertices))
        object.__setattr__(self, 'face_count', len(faces))
        object.__setattr__(self, 'boundary_digest', boundary_digest)

    @classmethod
    def _from_audit(cls, vertices, faces, facts):
        snapshot = cls.__new__(cls)
        snapshot._store(vertices, faces, facts['input_digest'])
        return snapshot

    def arrays(self):
        return (np.frombuffer(self._vertices, dtype='<f8').reshape((-1, 3)).copy(),
                np.frombuffer(self._faces, dtype='<i4').reshape((-1, 3)).copy())

    @property
    def identity(self):
        return {'scope': 'EXACT_SINGLE_OUTWARD_BOUNDARY',
                'boundary_digest': self.boundary_digest,
                'vertices': {'dtype': '<f8', 'shape': [self.vertex_count, 3],
                             'sha256': byte_digest(self._vertices)},
                'faces': {'dtype': '<i4', 'shape': [self.face_count, 3],
                          'sha256': byte_digest(self._faces)}}


def solver_candidate_snapshot(mesh):
    """Use a Nef region check only when a candidate has multiple shells.

    For one outward shell, the exact boundary audit proves a compact embedded
    oriented 2-manifold, so Jordan--Brouwer gives an unambiguous bounded solid.
    Multiple boundary components need nesting/alternating-orientation checks
    and retain the full ``MeshSnapshot`` path.
    """
    require(isinstance(mesh, (tuple, list)) and len(mesh) == 2,
            'solver candidate must be a (vertices, triangles) pair')
    vertices = _mesh_input_array(mesh[0])
    faces = _mesh_input_array(mesh[1], indices=True)
    facts = native().audit_exact_triangle_boundary(vertices, faces, 'AnyNonzero')
    require(facts['boundary_valid'],
            'invalid exact solver candidate boundary: ' + canonical(facts))
    signs = tuple(int(value) for value in facts['component_volume_signs'])
    if signs == (1,):
        return (SolverBoundarySnapshot._from_audit(vertices, faces, facts),
                'EXACT_SINGLE_OUTWARD_BOUNDARY')
    return MeshSnapshot((vertices, faces)), 'EXACT_NEF_SIGNED_REGION'


@dataclass(frozen=True, init=False)
class RationalSolverBoundarySnapshot(SolverBoundarySnapshot):
    """Checked single rational shell for trial maps; final material uses Nef."""
    _rational_vertices_json: bytes = field(repr=False)

    @classmethod
    def _from_rational_audit(cls, vertices, faces, rational_vertices, facts):
        snapshot=cls.__new__(cls)
        snapshot._store(vertices,faces,facts['input_digest'])
        object.__setattr__(snapshot,'_rational_vertices_json',
            canonical([list(row) for row in rational_vertices]).encode('ascii'))
        return snapshot

    @property
    def rational_vertices(self):
        return tuple(tuple(row) for row in parse(self._rational_vertices_json.decode('ascii')))

    @property
    def identity(self):
        return {**super().identity,'scope':'EXACT_RATIONAL_SINGLE_OUTWARD_BOUNDARY',
            'rational_vertices':{'format':'CANONICAL_FRACTION_V1',
                'sha256':byte_digest(self._rational_vertices_json)}}


@dataclass(frozen=True, init=False)
class RationalSeparatedSolverBoundarySnapshot(RationalSolverBoundarySnapshot):
    """Embedded positive root shells with exact mutual exterior witnesses.

    Trial-only proof; never accepted as a final MeshSnapshot/Nef region.
    """
    _separation_json: str = field(repr=False)

    @property
    def separation_proof(self):
        return parse(self._separation_json)

    @property
    def identity(self):
        return {**super().identity, 'scope': 'EXACT_RATIONAL_SEPARATED_OUTWARD_BOUNDARIES',
                'separation_digest': byte_digest(self._separation_json.encode('ascii'))}


@dataclass(frozen=True, init=False)
class RationalClassifiedSolverBoundarySnapshot(RationalSolverBoundarySnapshot):
    """Trial-only rational boundary with native exact signed nesting facts."""
    _classification_json: str = field(repr=False)

    @property
    def classification_proof(self):
        return parse(self._classification_json)

    @property
    def identity(self):
        return {**super().identity, 'scope': 'EXACT_RATIONAL_SIGNED_BOUNDARY_CLASSIFICATION',
                'classification_digest': byte_digest(self._classification_json.encode('ascii'))}


def _positive_root_exterior_witnesses(rows, faces, facts):
    """Called only after the actual complete exact embedded-boundary audit.

    Disjoint connected embedded closed boundaries cannot change inside/outside
    without intersecting. A vertex strictly outside the other's exact bounding
    box proves that entire boundary is exterior. Requiring this in BOTH ordered
    directions rules out nesting, without assuming disjoint bounding boxes.
    Inconclusive cases continue through the full signed-region constructor.
    """
    signs = facts['component_volume_signs']
    if len(signs) < 2 or any(sign != 1 for sign in signs):
        return None
    points = [tuple(fraction_from_rational_text(value) for value in row) for row in rows]
    members = [set() for _ in signs]
    for face, component in zip(faces, facts['face_components'], strict=True):
        members[component].update(map(int, face))
    extrema = []
    for indices in members:
        extrema.append([(min(indices, key=lambda i: (points[i][axis], i)),
                         max(indices, key=lambda i: (points[i][axis], -i))) for axis in range(3)])
    witnesses = []
    for i in range(len(signs)):
        for j in range(len(signs)):
            if i == j:
                continue
            found = None
            for axis in range(3):
                low, high = extrema[i][axis]
                other_low, other_high = extrema[j][axis]
                if points[low][axis] < points[other_low][axis]:
                    found = dict(component=i, exterior_to=j, vertex=low, axis=axis, side='BELOW',
                                 bound=rows[other_low][axis])
                elif points[high][axis] > points[other_high][axis]:
                    found = dict(component=i, exterior_to=j, vertex=high, axis=axis, side='ABOVE',
                                 bound=rows[other_high][axis])
                if found is not None:
                    break
            if found is None:
                return None
            witnesses.append(found)
    return dict(schema='exact-mutual-exterior-bbox-witnesses-v1',
                boundary_digest=facts['input_digest'], parents=[-1]*len(signs),
                outside_witnesses=witnesses)


def rational_solver_candidate_snapshot(rational_vertices, faces, *, fixed_boundary_proof=None,
                                      allow_separated_roots=False,
                                      allow_signed_classification=False):
    """Use exact rational boundary checking for one embedded outward shell.

    The same topology, vertex manifold, winding and exact self-intersection
    predicates used by the region constructor run here. Multiple shells and
    cavity nesting still require the full signed-region constructor. This
    intermediate snapshot is not a MeshSnapshot or a material seal.
    """
    require(type(allow_separated_roots) is bool, 'allow_separated_roots must be boolean')
    require(type(allow_signed_classification) is bool, 'allow_signed_classification must be boolean')
    require(not (allow_separated_roots and allow_signed_classification),
            'select only one exact multi-component classification profile')
    rows=_canonical_rational_vertices(rational_vertices)
    face_array=_mesh_input_array(faces,indices=True)
    n=native()
    if allow_signed_classification:
        operation = getattr(n, 'audit_exact_rational_material_boundary', None)
        require(callable(operation), 'native exact material classifier unavailable')
        classification = operation(rows, face_array, fixed_boundary_proof)
        facts = classification.get('boundary', {})
        require(classification.get('completed') and classification.get('classification_valid') and
                facts.get('boundary_valid'),
                'invalid exact rational material classification: ' + canonical(classification))
        require(facts.get('input_digest') and
                classification['original_signs'] == facts['component_volume_signs'],
                'exact material classification does not bind the audited boundary')
        display=np.asarray([[float(fraction_from_rational_text(x)) for x in row] for row in rows],dtype='<f8')
        require(np.isfinite(display).all(),'rational boundary display is nonfinite')
        snapshot=RationalClassifiedSolverBoundarySnapshot._from_rational_audit(
            display,face_array,rows,facts)
        object.__setattr__(snapshot,'_classification_json',canonical({
            key: classification[key] for key in ('algorithm','parents','depths','original_signs',
                'material_components','bounded_void_components','fixed_proof_identity')}))
        return snapshot, 'EXACT_RATIONAL_SIGNED_BOUNDARY_CLASSIFICATION'
    if not hasattr(n,'audit_exact_rational_triangle_boundary'):
        raise RuntimeError('native rational boundary audit unavailable')
    facts=(n.audit_exact_rational_triangle_boundary(rows,face_array,'AnyNonzero')
        if fixed_boundary_proof is None else
        n.audit_rational_boundary_with_fixed_proof(rows,face_array,'AnyNonzero',fixed_boundary_proof))
    require(facts['boundary_valid'],
        'invalid exact rational solver candidate boundary: '+canonical(facts))
    if tuple(int(value) for value in facts['component_volume_signs'])==(1,):
        display=np.asarray([[float(fraction_from_rational_text(x)) for x in row] for row in rows],dtype='<f8')
        require(np.isfinite(display).all(),'rational boundary display is nonfinite')
        return (RationalSolverBoundarySnapshot._from_rational_audit(display,face_array,rows,facts),
            'EXACT_RATIONAL_SINGLE_OUTWARD_BOUNDARY')
    if allow_separated_roots:
        proof = _positive_root_exterior_witnesses(rows, face_array, facts)
        if proof is not None:
            display = np.asarray([[float(fraction_from_rational_text(x)) for x in row] for row in rows], dtype='<f8')
            require(np.isfinite(display).all(), 'rational boundary display is nonfinite')
            snapshot = RationalSeparatedSolverBoundarySnapshot._from_rational_audit(display, face_array, rows, facts)
            object.__setattr__(snapshot, '_separation_json', canonical(proof))
            return snapshot, 'EXACT_RATIONAL_SEPARATED_OUTWARD_BOUNDARIES'
    return MeshSnapshot.from_rationals(rows,face_array),'EXACT_NEF_SIGNED_REGION'


@dataclass(frozen=True)
class MaterialSnapshot:
    region_id: str
    mesh: MeshSnapshot
    density: float

    @property
    def identity(self):
        return {'id': self.region_id, 'mesh': self.mesh.identity, 'density': self.density}


def partition(whole, parts):
    result = native().check_exact_material_partition(whole, parts)
    require(result['completed'] and result['complete_partition'],
            'incomplete or overlapping exact material partition: ' + canonical(result))
    return result


_VALIDATED_PROBLEMS = {}


def _live_problem_validation_signature(problem):
    def atom(value):
        return type(value), value
    def mesh(value):
        require(type(value) is MeshSnapshot, 'problem validation requires an actual MeshSnapshot')
        return tuple(atom(getattr(value, key)) for key in
            ('_vertices', '_faces', 'vertex_count', 'face_count', 'explicit_empty', '_rational_vertices_json')) + (value._exact_region,)
    def regions(values):
        require(type(values) is tuple, 'problem validation requires immutable region tuples')
        result = []
        for region in values:
            require(type(region) is MaterialSnapshot, 'problem validation requires actual material snapshots')
            result.append((atom(region.region_id), atom(region.density), mesh(region.mesh)))
        return tuple(result)
    return (mesh(problem.original), regions(problem.fixed_regions), regions(problem.design_regions),
            tuple(atom(getattr(problem, key)) for key in
                ('length_unit', 'mass_unit', 'minimum_outer_wall', 'validation_scope',
                 '_target', '_normalization', '_identity', 'input_digest')))


def _record_constructed_problem_validation(problem):
    if type(problem) is not SpinArtifactProblem:
        return  # Old subclass behavior remains available through full recheck.
    from .spin_common import implementation_observation
    key = id(problem)
    def remove(reference):
        entry = _VALIDATED_PROBLEMS.get(key)
        if entry is not None and entry[0] is reference:
            del _VALIDATED_PROBLEMS[key]
    _VALIDATED_PROBLEMS[key] = (weakref.ref(problem, remove),
        _live_problem_validation_signature(problem), native(), implementation_observation())


def _reuse_constructed_problem_validation(problem):
    """Run-local constructor certificate, bound to actual immutable content.

    Native handles and coordinate bytes are part of the witness, not only a
    caller-readable digest. Unregistered copies and any live dependency change
    are rejected. Persisted inputs always enter through a fresh constructor.
    """
    from .spin_common import implementation_observation
    require(type(problem) is SpinArtifactProblem, 'expected constructor-certified SpinArtifactProblem')
    entry = _VALIDATED_PROBLEMS.get(id(problem))
    require(entry is not None and entry[0]() is problem, 'problem was not constructor certified')
    require(native() is entry[2], 'problem native validation context changed')
    require(_live_problem_validation_signature(problem) == entry[1], 'problem validation content changed')
    require(implementation_observation() == entry[3], 'problem validation implementation changed')
    return problem


@dataclass(frozen=True, init=False)
class SpinArtifactProblem:
    """Input locked before a final geometry/claim is presented.

    Regions are ordered triples (unique ID, actual mesh, positive density).
    Unit labels define model units; density is mass_unit / length_unit**3.
    Target data is context only. Unknown original evidence always blocks admission.
    """
    original: MeshSnapshot
    fixed_regions: tuple[MaterialSnapshot, ...]
    design_regions: tuple[MaterialSnapshot, ...]
    length_unit: str
    mass_unit: str
    minimum_outer_wall: float
    validation_scope: str
    _target: str
    _normalization: str
    _identity: str
    input_digest: str

    def __init__(self, *, original, fixed_regions, design_regions, length_unit,
                 mass_unit, target, minimum_outer_wall, validation_scope=SCOPE):
        require(validation_scope == SCOPE, 'unsupported validation profile: ' + str(validation_scope))
        outer = MeshSnapshot(original)
        require(not outer.explicit_empty, 'original domain must be nonempty')
        seen = set()
        def regions(values):
            require(isinstance(values, (tuple, list)), 'regions must be an explicit list')
            out = []
            for row in values:
                require(isinstance(row, (tuple, list)) and len(row) == 3,
                        'region must contain id, mesh, density')
                region_id = text_label(row[0], 'region ID')
                require(region_id not in seen, 'duplicate material ID')
                seen.add(region_id)
                mesh = MeshSnapshot(row[1])
                require(not mesh.explicit_empty, 'original density region must be nonempty')
                out.append(MaterialSnapshot(region_id, mesh, scalar(row[2], 'density', positive=True)))
            return tuple(out)
        fixed, design = regions(fixed_regions), regions(design_regions)
        require(bool(design), 'at least one design region is required')
        initial_partition = partition(outer.exact_region(), [r.mesh.exact_region() for r in fixed + design])
        keys(target, ('mode', 'axis', 'origin', 'support_points'), 'target')
        require(isinstance(target['mode'], str) and target['mode'] in SPIN_MODES, 'unknown target mode')
        axis, origin = vector(target['axis'], 'axis'), vector(target['origin'], 'origin')
        require(any(x != 0. for x in axis), 'axis must be nonzero')
        require(isinstance(target['support_points'], (tuple, list)), 'support points must be explicit list')
        points = [vector(p, 'support point') for p in target['support_points']]
        require(target['mode'] != 'STATIC_BALANCE' or bool(points), 'static balance needs a support point')
        target_data = {'mode': target['mode'], 'axis': axis, 'origin': origin, 'support_points': points}
        v, _ = outer.arrays()
        lo, hi = v.min(axis=0), v.max(axis=0)
        with np.errstate(over='ignore', invalid='ignore'):
            spans = hi - lo
            center = lo + spans * .5
        length = scalar(float(max(spans)), 'normalization length', positive=True)
        density = max(r.density for r in fixed + design)
        try:
            scales = [length**power * density for power in (3, 4, 5)]
        except OverflowError as exc:
            raise ValueError('normalization units are unrepresentable') from exc
        require(all(math.isfinite(x) and x > 0 for x in scales), 'normalization units are unrepresentable')
        require(all(r.density / density > 0 for r in fixed + design), 'normalized density underflow')
        normalization = {'origin': vector(center, 'normalization origin'), 'L': length, 'R': density,
                         'axes': [[1.,0.,0.], [0.,1.,0.], [0.,0.,1.]],
                         'authority': 'IMPLEMENTATION_DEFAULT', 'frame': 'INPUT_NORMALIZED_V1'}
        attrs = {'original': outer, 'fixed_regions': fixed, 'design_regions': design,
                 'length_unit': text_label(length_unit, 'length unit'),
                 'mass_unit': text_label(mass_unit, 'mass unit'),
                 'minimum_outer_wall': scalar(minimum_outer_wall, 'minimum outer wall', positive=True),
                 'validation_scope': validation_scope, '_target': canonical(target_data),
                 '_normalization': canonical(normalization)}
        for name, value in attrs.items():
            object.__setattr__(self, name, value)
        identity = {'schema': 'spin-artifact-problem-v1', 'validation_scope': SCOPE,
                    'original': outer.identity, 'fixed_regions': [r.identity for r in fixed],
                    'design_regions': [r.identity for r in design],
                    'length_unit': self.length_unit, 'mass_unit': self.mass_unit,
                    'target': target_data, 'minimum_outer_wall': self.minimum_outer_wall,
                    'wall_definition': 'NEWLY_CARVED_BOUNDARY_TO_ORIGINAL_BOUNDARY',
                    'moment_absolute_tolerance': TOLERANCE, 'normalization': normalization,
                    'original_evidence_admission': 'BLOCKED_UNKNOWN',
                    'partition': initial_partition}
        object.__setattr__(self, '_identity', canonical(identity))
        object.__setattr__(self, 'input_digest', digest(identity))
        _record_constructed_problem_validation(self)

    @property
    def normalization(self):
        return parse(self._normalization)

    @property
    def target(self):
        return parse(self._target)

    @property
    def identity(self):
        return parse(self._identity)
