"""Immutable indexed geometry and locked input for hollow-only material checking."""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from ..spin import SPIN_MODES
from .spin_common import (SCOPE, TOLERANCE, byte_digest, canonical, digest, keys,
                          native, parse, require, scalar, text_label, vector)


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

    def __init__(self, mesh):
        if isinstance(mesh, MeshSnapshot):
            mesh = None if mesh.explicit_empty else mesh.arrays()
        if mesh is None:
            native().empty_exact_material_region()
            v, f = np.empty((0, 3), dtype='<f8'), np.empty((0, 3), dtype='<i4')
        else:
            require(isinstance(mesh, (tuple, list)) and len(mesh) == 2,
                    'mesh must be a (vertices, triangles) pair or explicit None')
            # Freeze and validate supplied scalar types before sequence inference;
            # then validate both actual staging arrays with the native checker.
            v = _mesh_input_array(mesh[0])
            f = _mesh_input_array(mesh[1], indices=True)
            region = native().make_exact_material_region(v, f)
            require(region.facts['completed'] and region.facts['valid'],
                    'invalid material boundary: ' + region.facts['error'])
            v, f = np.asarray(v, dtype='<f8', order='C'), np.asarray(f, dtype='<i4', order='C')
        object.__setattr__(self, '_vertices', v.tobytes(order='C'))
        object.__setattr__(self, '_faces', f.tobytes(order='C'))
        object.__setattr__(self, 'vertex_count', len(v))
        object.__setattr__(self, 'face_count', len(f))
        object.__setattr__(self, 'explicit_empty', mesh is None)

    def arrays(self):
        return (np.frombuffer(self._vertices, dtype='<f8').reshape((-1, 3)).copy(),
                np.frombuffer(self._faces, dtype='<i4').reshape((-1, 3)).copy())

    def exact_region(self):
        n = native()
        out = n.empty_exact_material_region() if self.explicit_empty else n.make_exact_material_region(*self.arrays())
        require(out.facts['completed'] and out.facts['valid'], 'invalid authenticated region')
        return out

    @property
    def identity(self):
        return {'explicit_empty': self.explicit_empty,
                'vertices': {'dtype': '<f8', 'shape': [self.vertex_count, 3],
                             'sha256': byte_digest(self._vertices)},
                'faces': {'dtype': '<i4', 'shape': [self.face_count, 3],
                          'sha256': byte_digest(self._faces)}}


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

    @property
    def normalization(self):
        return parse(self._normalization)

    @property
    def target(self):
        return parse(self._target)

    @property
    def identity(self):
        return parse(self._identity)
