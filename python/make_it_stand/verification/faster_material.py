"""Locked hollow-only Faster candidates/runs and independently sealed material.

Native coordinates are input-normalized model coordinates; quadric coefficients
are native p (matrix-entry convention) in the derived input-normalized target
frame. This module never accepts caller-supplied moments or repaired geometry.
The immutable result binds the run to its material seal. Saving just its artifact
with spin_bundle persists only the material seal, not this run binding.
"""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, replace
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from ..spin import (SpinEvaluatorConfig, SpinItFasterRequest, SpinSolverConfig,
                    _finite_tuple, _native_module, compute_hydrostatic_geometry,
                    find_hydrostatic_equilibrium,
                    evaluate_spin_it_faster_plane_halfspaces,
                    evaluate_spin_it_faster_quadric,
                    evaluate_spin_it_faster_fullrank_quadric,
                    integrate_raw_mass_moments, solve_spin_it_faster,
                    solve_spin_it_faster_fullrank_quadric,
                    solve_spin_it_faster_with_external_moment_evaluator,
                    solve_spin_it_faster_buoyant_plane,
                    solve_spin_it_faster_static_plane_family)
from .mass_moments import TriangleMesh, compute_mass_properties
from .spin_artifact import SealedSpinArtifact, _recheck_problem
from .spin_common import FRAME, ORDER, SCOPE, canonical, digest, native, parse, require, scalar, text_label, vector
from .spin_problem import MeshSnapshot, SpinArtifactProblem, solver_candidate_snapshot, _mesh_input_array
from .rational_text import (fraction_from_rational_text, fraction_to_rational_text,
                            integer_from_decimal_text)


@dataclass(frozen=True, init=False)
class FasterMaterialResult:
    """Run binding kept separate from the independently reloadable material seal."""
    status: str
    artifact: SealedSpinArtifact | None
    seal_error: str | None
    scientific_digest: str
    _identity: str

    def __init__(self):
        raise TypeError('use a public Faster material candidate or solver bridge')

    @property
    def identity(self):
        return parse(self._identity)


def _parse_boolean_exact_rational(token: str) -> Fraction:
    if not isinstance(token, str) or token.count('/') != 1:
        raise ValueError('exact Boolean coordinate must be one rational token')
    numerator_text, denominator_text = token.split('/')

    def integer(text: str) -> int:
        sign = -1 if text.startswith('-') else 1
        unsigned = text[1:] if sign < 0 else text
        if unsigned.startswith('0b'):
            return sign * int(unsigned[2:], 2)
        return sign * integer_from_decimal_text(unsigned)

    numerator = integer(numerator_text)
    denominator = integer(denominator_text)
    if denominator <= 0:
        raise ValueError('exact Boolean coordinate denominator must be positive')
    return Fraction(numerator, denominator)


def read_exact_boolean_sidecar_coordinates(vertices, sidecar_path):
    """Bind Cherchi rational coordinates to an OBJ, without material authority.

    The OBJ must exactly match the sidecar's in-memory and serialized floating
    coordinates. Exact rationals remain authoritative; their floating display
    coordinates may differ from the OBJ only within a fixed roundoff envelope
    of 64 eps times the locked coordinate scale. Incomplete, reordered,
    unmatched, or ambiguous sidecars fail closed. This coordinate-only step
    may precede orientation normalization; closure, embedding and material
    validity must still be checked by MeshSnapshot or the exact region API.
    """
    display = np.asarray(vertices, dtype=np.float64)
    if display.ndim != 2 or display.shape[1:] != (3,) or not np.isfinite(display).all():
        raise ValueError('Boolean display vertices must be finite N-by-3 coordinates')
    coordinate_scale = float(max(
        1.0,
        float(np.max(np.abs(display))) if display.size else 0.0,
        float(np.linalg.norm(np.ptp(display, axis=0))) if len(display) else 0.0,
    ))
    absolute_error_bound = float(64.0 * np.finfo(np.float64).eps * coordinate_scale)
    expected_header = ('vertex_index', 'status', 'exact_x', 'exact_y', 'exact_z',
                       'memory_x', 'memory_y', 'memory_z', 'obj_x', 'obj_y', 'obj_z')
    exact_vertices = []
    maximum_absolute_error = 0.0
    approximate_vertex_count = 0
    sidecar_path = Path(sidecar_path)
    try:
        with sidecar_path.open('r', encoding='ascii', newline='') as stream:
            reader = csv.DictReader(stream, delimiter=chr(9))
            if tuple(reader.fieldnames or ()) != expected_header:
                raise ValueError('exact Boolean sidecar has an unsupported header')
            for index, row in enumerate(reader):
                if index >= len(display):
                    raise ValueError('exact Boolean sidecar has extra vertex rows')
                status = row.get('status')
                if row.get('vertex_index') != str(index) or status not in (
                        'MATCHED', 'APPROXIMATION_DISAGREES'):
                    raise ValueError('exact Boolean sidecar is incomplete, reordered, or unmatched')
                if status == 'APPROXIMATION_DISAGREES':
                    approximate_vertex_count += 1
                exact_row = tuple(_parse_boolean_exact_rational(row[f'exact_{axis}'])
                                  for axis in 'xyz')
                for axis_index, axis in enumerate('xyz'):
                    try:
                        memory_value = float(row[f'memory_{axis}'])
                        obj_value = float(row[f'obj_{axis}'])
                    except (TypeError, ValueError, OverflowError) as error:
                        raise ValueError('exact Boolean sidecar has invalid display coordinates') from error
                    exact_value = float(exact_row[axis_index])
                    display_value = float(display[index, axis_index])
                    if (not all(math.isfinite(value) for value in
                                (memory_value, obj_value, exact_value, display_value)) or
                            memory_value != obj_value or obj_value != display_value):
                        raise ValueError('exact Boolean sidecar does not bind the supplied OBJ vertex')
                    error = abs(exact_value - display_value)
                    maximum_absolute_error = max(maximum_absolute_error, error)
                    if error > absolute_error_bound:
                        raise ValueError(
                            'exact Boolean rational does not bind the supplied OBJ vertex within the bounded derivation error'
                        )
                exact_vertices.append(tuple(
                    fraction_to_rational_text(value) for value in exact_row))
    except OSError as error:
        raise ValueError('exact Boolean sidecar cannot be read') from error
    if len(exact_vertices) != len(display):
        raise ValueError('exact Boolean sidecar vertex count does not match the OBJ')
    audit = {
        'schema': 'exact-boolean-float-derivation-audit-v1',
        'status': 'BOUNDED_FLOAT_DERIVATION' if maximum_absolute_error else 'EXACT_FLOAT_MATCH',
        'vertex_count': int(len(display)),
        'approximation_status_vertex_count': approximate_vertex_count,
        'coordinate_scale': coordinate_scale,
        'maximum_absolute_error': maximum_absolute_error,
        'absolute_error_bound': absolute_error_bound,
        'sidecar_sha256': hashlib.sha256(sidecar_path.read_bytes()).hexdigest(),
    }
    return tuple(exact_vertices), audit


def mesh_snapshot_from_exact_boolean_sidecar(vertices, faces, sidecar_path, *, return_audit=False):
    """Authenticate a rational boundary after binding its original sidecar.

    A coordinate-only read cannot pass closure, winding, embedding or material
    gates. Preserve the original complete audit shape for legacy consumers.
    """
    exact_vertices, audit = read_exact_boolean_sidecar_coordinates(vertices, sidecar_path)
    snapshot = MeshSnapshot.from_rationals(exact_vertices, faces)
    audit['rational_geometry_identity'] = snapshot.identity
    return (snapshot, audit) if return_audit else snapshot

def rational_target_affine_binding(request):
    """Bind the actual binary64-input affine map and its exact rational inverse.

    A numerically orthonormal input need not be exactly orthonormal as a
    rational matrix. This is coordinate evidence, not material acceptance.
    """
    axes = np.asarray(request.target_axes, dtype=np.float64)
    origin = np.asarray(request.target_origin, dtype=np.float64)
    require(axes.shape == (3, 3) and origin.shape == (3,) and
            np.isfinite(axes).all() and np.isfinite(origin).all(),
            'exact target-frame transform requires a finite locked frame')
    matrix = [[Fraction.from_float(float(axes[j, i])) for j in range(3)] for i in range(3)]
    offset = [Fraction.from_float(float(x)) for x in origin]
    augmented = [row + [Fraction(int(i == j)) for j in range(3)]
                 for i, row in enumerate(matrix)]
    for column in range(3):
        pivot = next((row for row in range(column, 3) if augmented[row][column]), None)
        require(pivot is not None, 'singular actual target-frame affine map')
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        denominator = augmented[column][column]
        augmented[column] = [x/denominator for x in augmented[column]]
        for row in range(3):
            if row != column:
                multiplier = augmented[row][column]
                augmented[row] = [a-multiplier*b for a, b in zip(augmented[row], augmented[column])]
    inverse = [row[3:] for row in augmented]
    require(all(sum(matrix[i][k]*inverse[k][j] for k in range(3)) == int(i == j)
                for i in range(3) for j in range(3)), 'actual affine inverse identity failed')
    encode = fraction_to_rational_text
    body = dict(schema='rational-target-affine-v1',
        coordinate_semantics='ACTUAL_BINARY64_FRAME_ELEMENTS_INTERPRETED_RATIONALLY',
        model_to_target_matrix=[[encode(x) for x in row] for row in matrix],
        model_to_target_translation=[encode(-sum(matrix[i][k]*offset[k] for k in range(3))) for i in range(3)],
        target_to_model_matrix=[[encode(x) for x in row] for row in inverse],
        target_to_model_translation=[encode(x) for x in offset],
        original_units_per_normalized_unit=encode(Fraction.from_float(float(
            scalar(request.original_units_per_normalized_unit, 'original unit scale', positive=True)))))
    return dict(body, affine_digest=digest(body))


def verify_rational_target_affine_binding(binding, request):
    """Recompute from the locked request; a copied digest cannot bless mutation."""
    require(isinstance(binding, dict) and binding == rational_target_affine_binding(request),
            'saved actual target-frame affine binding mismatch')


def _target_rationals_to_normalized_model(vertices, request, *, exact_inverse=False,
                                         return_binding=False):
    """Apply the locked binary64 target frame as an exact rational affine map."""
    axes = np.asarray(request.target_axes, dtype=np.float64)
    origin = np.asarray(request.target_origin, dtype=np.float64)
    require(axes.shape == (3, 3) and origin.shape == (3,) and
            np.isfinite(axes).all() and np.isfinite(origin).all(),
            'exact target-frame transform requires a finite locked frame')
    require(type(exact_inverse) is bool and type(return_binding) is bool,
            'target-frame inverse policy must be explicit bool')
    require(not return_binding or exact_inverse, 'affine inverse binding requires exact inverse policy')
    binding = rational_target_affine_binding(request) if exact_inverse else None
    inverse = ([[fraction_from_rational_text(x) for x in row]
                for row in binding['target_to_model_matrix']] if binding is not None else None)
    transformed = []
    for row in vertices:
        target = tuple(fraction_from_rational_text(value) for value in row)
        model = tuple(sum((target[i] * (inverse[j][i] if inverse is not None
                           else Fraction.from_float(float(axes[j, i])))
                           for i in range(3)), Fraction.from_float(float(origin[j])))
                      for j in range(3))
        transformed.append(tuple(fraction_to_rational_text(value) for value in model))
    return (tuple(transformed), binding) if return_binding else tuple(transformed)


@dataclass(frozen=True, init=False)
class FrozenFasterMaterialComposition:
    """A prior sealed material transferred as fixed geometry into a new request.

    This records a geometry dependency, not proof that the source object was
    optimized to a paper target.  The upstream run evidence must remain bound
    separately by its caller.
    """
    request: SpinItFasterRequest
    _binding: str

    def __init__(self):
        raise TypeError('use compose_sealed_faster_material_as_fixed_region')

    @property
    def binding(self):
        return parse(self._binding)


@dataclass(frozen=True, init=False)
class FloatingIslandChainResult:
    """B solve bound to the frozen A material and final independent audits."""
    status: str
    composition: FrozenFasterMaterialComposition
    _solver_result: str
    _identity: str

    def __init__(self):
        raise TypeError('use solve_floating_island_b_from_sealed_a')

    @property
    def solver_result(self):
        return parse(self._solver_result)

    @property
    def identity(self):
        return parse(self._identity)


def _raw_moments(mesh, density):
    """Independent signed-tetrahedron mass moments in Faster's ten-term order."""
    v, f = mesh.arrays()
    rational_vertices = getattr(mesh, 'rational_vertices', None)
    return _raw_moments_from_arrays(v, f, density, rational_vertices=rational_vertices)


def _raw_moments_from_arrays(vertices, faces, density, *, rational_vertices=None):
    """Integrate an already-authenticated boundary without rebuilding its Nef.

    Callers may use this for an invertible affine image of a locked material:
    the source MeshSnapshot owns the topology/embedding proof, while the
    transform changes only the coordinates used by the signed-tetrahedron
    integral. New or solver-produced boundaries must still pass MeshSnapshot
    or the exact rational candidate audit before reaching this helper.
    """
    v = _mesh_input_array(vertices)
    f = _mesh_input_array(faces, indices=True)
    if rational_vertices is not None:
        # The authenticated exact boundary is authoritative. Its display mesh
        # can contain collapsed slivers even when every exact face is valid;
        # integrating the display would silently change the represented solid.
        from .mass_moments import integrate_exact_rational_raw_moments
        weight = Fraction(scalar(density, 'material density', positive=True))
        geometric = integrate_exact_rational_raw_moments(rational_vertices, f.tolist())
        moments = tuple(float(value * weight) for value in geometric)
        require(all(math.isfinite(value) for value in moments) and moments[0] > 0.,
                'exact material moments are nonfinite or have nonpositive mass')
        return moments
    result = compute_mass_properties(
        TriangleMesh(tuple(map(tuple, v.tolist())), tuple(map(tuple, f.tolist()))), density)
    q = result.second_moment
    return (result.mass, *result.first_moment,
            q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1])


def _add_moments(*rows):
    return tuple(math.fsum(float(row[i]) for row in rows) for i in range(10))


@dataclass(frozen=True, init=False)
class Spin2014StaticMaterialOutput:
    """Geometry attachment, not a SealedSpinArtifact or manufacturing grant.

    The saved candidate is never reoptimized or reevaluated at an opposite
    plane. Its MODEL-frame retained boundary and an exact rational complement
    are checked by the existing material-region and independent moment paths.
    """
    outer: MeshSnapshot
    fixed: MeshSnapshot
    design: MeshSnapshot
    retained: MeshSnapshot
    cavity: MeshSnapshot
    material: MeshSnapshot
    scientific_digest: str
    _identity: str
    _audit: str

    def __init__(self):
        raise TypeError('use prepare_spin2014_static_material_output')

    @property
    def identity(self):
        return parse(self._identity)

    @property
    def audit(self):
        return parse(self._audit)


def _nested_material_boundary(outer, removed):
    """Retain the original outer indexing; reverse the actual cavity boundary.

    This is a representation operation, not a geometry repair. Touching or
    overlapping boundaries fail the existing exact signed-region audit.
    """
    if removed.explicit_empty:
        return outer
    if outer.identity == removed.identity:
        return MeshSnapshot(None)
    vertices, faces = outer.arrays()
    inner_vertices, inner_faces = removed.arrays()
    def rational_rows(mesh, rows):
        return mesh.rational_vertices or tuple(tuple(
            fraction_to_rational_text(Fraction.from_float(float(value)))
            for value in row) for row in rows)
    exact = tuple(rational_rows(outer, vertices)) + tuple(rational_rows(removed, inner_vertices))
    return MeshSnapshot.from_rationals(exact,
        np.vstack((faces, inner_faces[:, ::-1] + len(vertices))).tolist())


def _canonical_exact_boolean_boundary(vertices, faces):
    """Canonicalize indexing only, then authenticate the exact boundary once.

    Do not merge vertices, remove faces, flip orientation or retriangulate a
    coplanar patch. Alternative triangulations still differ and fail reload;
    this only removes vertex/face traversal order from the scientific digest.
    """
    exact_rows = tuple(tuple(fraction_to_rational_text(_parse_boolean_exact_rational(token))
                             for token in row) for row in vertices)
    # Validate indices losslessly before permutation. Authentication of the
    # canonically reindexed boundary below still performs all exact geometric
    # checks; a bijective reindexing and oriented cyclic face rotation neither
    # repair nor alter the represented region. No cached region is trusted.
    old_faces = _mesh_input_array(faces, indices=True)
    if old_faces.dtype.kind not in ('i', 'u') or old_faces.dtype.itemsize > 8:
        raise TypeError('faces must use an integer dtype')
    require(old_faces.ndim == 2 and old_faces.shape[1] == 3 and len(old_faces) > 0,
            'triangle indices: expected nonempty rows of three scalars')
    require(all(0 <= int(index) < len(exact_rows) and int(index) < 2**31
                for index in old_faces.flat), 'triangle index out of range')
    order = sorted(range(len(exact_rows)),
                   key=lambda i: tuple(fraction_from_rational_text(x) for x in exact_rows[i]))
    inverse = {old: new for new, old in enumerate(order)}
    canonical_faces = []
    for face in old_faces:
        mapped = tuple(inverse[int(x)] for x in face)
        canonical_faces.append(min(mapped, mapped[1:]+mapped[:1], mapped[2:]+mapped[:2]))
    return MeshSnapshot.from_rationals(tuple(exact_rows[i] for i in order),
                                      sorted(canonical_faces))


def prepare_spin2014_static_material_output(run):
    """Authenticate the material represented by one actual static candidate.

    A qualified squared-height attachment is necessary but is not trusted as
    geometry evidence: regions, partitions, density-weighted moments, target
    frame and objective are freshly checked below. The result remains an
    algorithm geometry attachment; the existing seal's wall/manufacturing
    requirements are neither bypassed nor replaced by this entrypoint.
    """
    from ..spin import Spin2014StaticQuadricMaterialRun
    require(isinstance(run, Spin2014StaticQuadricMaterialRun), 'expected a static quadric material run')
    request, candidate, checked = run.measured_request, run.candidate, run.squared_height
    require(run.objective_profile == 'SPIN2014_OBJECTIVES_QUADRIC_V1' and
            request.mode == 'STATIC_BALANCE' and run.objective_settings.objective == 'STATIC_LOWERING',
            'static material output requires the original squared-height profile')
    require(candidate.get('status') == 'CONVERGED', 'static candidate is not a converged terminal result')
    require(checked.get('squared_height_conditions_verified') is True,
            'static squared-height terminal conditions are not qualified')
    fixed_density = scalar(run.fixed_material_density, 'fixed density', positive=True)
    variable_density = scalar(request.material_density, 'design density', positive=True)
    origin = np.asarray(request.target_origin, dtype=np.float64)
    axes = np.asarray(request.target_axes, dtype=np.float64)
    require(origin.shape == (3,) and axes.shape == (3,3) and np.isfinite(origin).all() and
            np.isfinite(axes).all() and np.max(abs(axes.T @ axes - np.eye(3))) <= 1e-12 and
            np.linalg.det(axes) > 0., 'invalid target coordinate frame')
    outer = MeshSnapshot((request.outer_vertices, request.outer_faces))
    design = MeshSnapshot((request.design_vertices, request.design_faces))
    n = native()
    def subset(one, two, label):
        audit = n.check_exact_material_subset(one.exact_region(), two.exact_region())
        require(audit['completed'] and audit['a_subset_b'], label + ': exact containment failed')
        return audit
    subset(design, outer, 'design')
    target_vertices = np.asarray(candidate['variable_material_vertices'], dtype=np.float64)
    target_faces = np.asarray(candidate['variable_material_faces'])
    require(bool(len(target_vertices)) == bool(len(target_faces)), 'missing retained material boundary')
    coordinate_map = 'SAVED_TARGET_TO_MODEL_BINARY64_V1'
    full_source_audit = None
    if len(target_vertices) and candidate.get('material_region') == 'STRICT_FULL':
        # A full-domain terminal mesh may have reordered vertices/triangles.
        # Prove entity equality in the actual TARGET map before reusing the
        # original MODEL design; a binary64 inverse roundtrip is not identity.
        dv, df = design.arrays()
        target_design = MeshSnapshot(((dv-origin) @ axes, df))
        target_retained = MeshSnapshot((target_vertices, target_faces))
        full_source_audit = n.check_exact_material_partition(target_design.exact_region(),
            [target_retained.exact_region()])
        require(full_source_audit['completed'] and full_source_audit['complete_partition'],
                'STRICT_FULL retained boundary does not represent the actual whole target design')
        retained = design
        coordinate_map = 'EXACT_WHOLE_DESIGN_SOURCE_IDENTITY_V1'
    elif len(target_vertices):
        # This is exactly the saved candidate's declared affine export. Never
        # move a vertex to a guessed plane or weld two rounded evaluations.
        retained = MeshSnapshot((target_vertices @ axes.T + origin, target_faces))
    else:
        retained = MeshSnapshot(None)
    subset(retained, design, 'retained')
    if retained.explicit_empty:
        cavity = design
        construction = 'EMPTY_RETAINED_DESIGN_COMPLEMENT_V1'
    else:
        equality = n.check_exact_material_partition(design.exact_region(), [retained.exact_region()])
        if equality['completed'] and equality['complete_partition']:
            cavity = MeshSnapshot(None)
            construction = 'FULL_RETAINED_EMPTY_COMPLEMENT_V1'
        else:
            dv, df = design.arrays()
            rv, rf = retained.arrays()
            exact = n.libigl_exact_boolean_construction(dv, df, rv, rf, 'difference')
            require(exact['succeeded'], 'exact same-retained complement construction failed')
            cavity = _canonical_exact_boolean_boundary(exact['vertices'], exact['faces'])
            construction = 'DESIGN_MINUS_SAVED_RETAINED_EXACT_V1'
    fixed = _nested_material_boundary(outer, design)
    material = _nested_material_boundary(outer, cavity)
    def partition(whole, parts, label):
        audit = n.check_exact_material_partition(whole.exact_region(), [p.exact_region() for p in parts])
        require(audit['completed'] and audit['complete_partition'] and
                audit['no_positive_volume_overlap'] and all(audit['part_subset_whole']),
                label + ': exact partition failed')
        return audit
    partitions = {'design_partition': partition(design, [retained, cavity], 'design'),
                  'whole_partition': partition(outer, [fixed, retained, cavity], 'whole'),
                  'material_partition': partition(material, [fixed, retained], 'material')}
    def moments(mesh, density):
        return np.zeros(10) if mesh.explicit_empty else np.asarray(_raw_moments(mesh, density))
    fixed_model = moments(fixed, fixed_density)
    retained_model = moments(retained, variable_density)
    total_model = fixed_model + retained_model
    model_to_target = axes.T
    translation = -model_to_target @ origin
    total_target = np.asarray(_transform_moments(total_model, model_to_target, translation))
    variable_target = np.asarray(_transform_moments(retained_model, model_to_target, translation))
    outer_vertices = outer.arrays()[0]
    outer_target = (outer_vertices-origin) @ axes
    length = float(np.max(np.ptp(outer_target, axis=0)))
    require(length > 0. and abs(length-run.target_length) <= 16*np.finfo(float).eps*length and
            np.min(outer_target[:,2]) >= 0., 'invalid input-derived target scale or support plane')
    scales = max(fixed_density, variable_density)*length**np.array([3,4,4,4,5,5,5,5,5,5])
    require(np.isfinite(scales).all() and (scales > 0.).all(), 'unrepresentable moment normalization')
    errors = []
    for actual, recorded in ((fixed_model, request.fixed_shell_moments),
                             (variable_target, candidate['variable_raw_moments']),
                             (total_target, candidate['total_raw_moments'])):
        recorded = np.asarray(recorded, dtype=np.float64)
        require(recorded.shape == (10,) and np.isfinite(recorded).all(), 'invalid recorded material moments')
        errors.append(float(np.max(abs(actual-recorded)/scales)))
    require(max(errors) <= 1e-6, 'geometry/density does not reproduce the saved material moments')
    mass = scalar(float(total_target[0]), 'material mass', positive=True)
    height = float(total_target[3]/mass)
    length_units = scalar(run.objective_settings.solver_length_unit_in_declared_units, 'declared length', positive=True)
    factor = scalar(run.objective_settings.center_weight, 'center weight', positive=True)*length_units**2
    objective = factor*height**2
    require(height > 0. and math.isfinite(objective) and
            abs(objective-scalar(checked['objective'], 'recorded squared objective')) <= 1e-10*max(1.,objective),
            'independent squared objective mismatch')
    # Preserve independent moments of each actual material, including all
    # off-diagonal inertia terms; caller can recheck the numerical gate/support.
    second = np.array([[total_model[4],total_model[9],total_model[8]],
                       [total_model[9],total_model[5],total_model[7]],
                       [total_model[8],total_model[7],total_model[6]]])
    centered = second - np.outer(total_model[1:4],total_model[1:4])/mass
    inertia = np.trace(centered)*np.eye(3)-centered
    audit = dict(schema='spin2014-static-material-output-audit-v1',
        scope='EXACT_REGION_PARTITION_AND_INDEPENDENT_MOMENTS_NOT_MANUFACTURING',
        complement_construction=construction, coordinate_map=coordinate_map,
        full_design_source_identity_audit=full_source_audit,
        fixed_density=fixed_density, design_density=variable_density,
        fixed_moments_model=fixed_model.tolist(), retained_moments_model=retained_model.tolist(),
        total_moments_model=total_model.tolist(), total_moments_target=total_target.tolist(),
        cavity_geometric_moments=moments(cavity, 1.).tolist(),
        center_of_mass_model=(total_model[1:4]/mass).tolist(), center_inertia_model=inertia.tolist(),
        squared_objective=objective, normalized_moment_error=max(errors), independent_moments_verified=True,
        manufacturing_approved=False, sealed_artifact=False, **partitions)
    identity = {'schema': 'spin2014-static-material-output-v1',
        'objective_profile': run.objective_profile, 'objective_identity': run.objective_identity,
        'source_run_sha256': digest(_json_value(asdict(run))),
        'meshes': {name: mesh.identity for name,mesh in (('outer',outer),('fixed',fixed),
            ('design',design),('retained',retained),('cavity',cavity),('material',material))},
        'audit': audit}
    result = object.__new__(Spin2014StaticMaterialOutput)
    for name, value in (('outer',outer),('fixed',fixed),('design',design),('retained',retained),
                        ('cavity',cavity),('material',material),('_audit',canonical(_json_value(audit))),
                        ('_identity',canonical(_json_value(identity))),('scientific_digest',digest(_json_value(identity)))):
        object.__setattr__(result, name, value)
    return result


def _json_value(value):
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        token = 'nan' if math.isnan(value) else '+inf' if value > 0. else '-inf'
        return {'$nonfinite_float': token}
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


@dataclass(frozen=True, init=False)
class Spin2014QuadricMaterialOutput:
    """Existing material-system attachment; never a manufacturing seal."""
    outer: MeshSnapshot
    fixed: MeshSnapshot
    design: MeshSnapshot
    retained: MeshSnapshot
    cavity: MeshSnapshot
    material: MeshSnapshot
    scientific_digest: str
    _identity: str
    _audit: str

    def __init__(self):
        raise TypeError('use prepare_spin2014_quadric_material_output')

    @property
    def identity(self):
        return parse(self._identity)

    @property
    def audit(self):
        return parse(self._audit)


def _independent_spin2014_rotational_physics(moments, objective):
    """Original equations checked independently of the native AD evaluator.

    Supports complex inputs for independent complex-step differentiation.
    In particular mass-weighted lowering is L_z squared in declared units.
    """
    m=objective.solver_mass_unit_in_declared_units
    l=objective.solver_length_unit_in_declared_units
    d=np.array([m,m*l,m*l,m*l,*([m*l*l]*6)])
    s=np.asarray(moments)*d
    mass=s[0];first=s[1:4]
    q=np.array([[s[4],s[9],s[8]],[s[9],s[5],s[7]],[s[8],s[7],s[6]]])
    centered=q-np.outer(first,first)/mass
    inertia=np.trace(centered)*np.eye(3)-centered
    ratio=(inertia[0,0]**2+inertia[1,1]**2+2*inertia[0,1]**2)/inertia[2,2]**2
    f=objective.inertia_weight*ratio
    if objective.objective=='TOP_LOWERING_ONLY':
        f+=objective.center_weight*(first[2]/mass)**2
    elif objective.objective=='TOP_MASS_AND_LOWERING':
        f+=objective.center_weight*first[2]**2
    if objective.objective=='YOYO':
        constraints=np.array([*first,inertia[0,2],inertia[1,2]])
    else:
        constraints=np.array([first[0],first[1],inertia[0,2],inertia[1,2]])
    return f,constraints


def _oriented_coordinate_faces(vertices, faces):
    """Index-independent, exact-binary64 oriented triangle multiset."""
    v=np.asarray(vertices,dtype=float);f=np.asarray(faces)
    require(v.ndim==2 and v.shape[1:]==(3,) and np.isfinite(v).all() and
            f.ndim==2 and f.shape[1:]==(3,) and np.issubdtype(f.dtype,np.integer) and
            (f>=0).all() and (f<len(v)).all(),'invalid comparison boundary')
    rows=[]
    for face in f:
        row=tuple(tuple(float(x) for x in v[int(i)]) for i in face)
        rows.append(min(row,row[1:]+row[:1],row[2:]+row[:2]))
    return sorted(rows)


def _exact_spin2014_saved_retained(run, evaluator_settings):
    """Recover rationals only when they round to the actual saved endpoint.

    This is geometry replay, not another optimizer or a repair. The normal
    settings, Boolean operation and both operand boundaries are explicit.
    """
    require(isinstance(evaluator_settings,SpinEvaluatorConfig),
            'saved evaluator settings are required for exact material replay')
    if evaluator_settings.rational_source_coordinates:
        return _rational_spin2014_saved_retained(run, evaluator_settings)
    r,s=run.measured_request,run.solve
    replay=evaluate_spin_it_faster_quadric(r,s.quadric_coefficients,
        evaluator_settings=evaluator_settings,coefficients_are_normalized=True)
    require(replay.geometry_hash==s.final_geometry_hash,
            'actual replay geometry identity differs from the saved endpoint')
    origin=np.asarray(r.target_origin);axes=np.asarray(r.target_axes)
    design_target=(np.asarray(r.design_vertices)-origin)@axes
    operation='intersection' if replay.quadric_solid_is_positive_region else 'difference'
    exact=native().libigl_exact_boolean_construction(design_target,r.design_faces,
        replay.interface_vertices,replay.interface_faces,operation)
    require(exact['succeeded'],'exact saved-operand Boolean replay failed')
    rationals=tuple(tuple(fraction_to_rational_text(_parse_boolean_exact_rational(x)) for x in row)
                    for row in exact['vertices'])
    rounded=np.asarray([[float(fraction_from_rational_text(x)) for x in row] for row in rationals])
    require(_oriented_coordinate_faces(rounded,np.asarray(exact['faces'],dtype=np.int32))==
            _oriented_coordinate_faces(s.material_vertices,np.asarray(s.material_faces,dtype=np.int32)),
            'exact replay does not round to the same saved oriented material boundary')
    retained=_canonical_exact_boolean_boundary(
        _target_rationals_to_normalized_model(rationals,r),exact['faces'])
    # Independently audit the floating derivative, never substitute it for the
    # authority. Its containment failure is preserved even when rationals pass.
    display=MeshSnapshot((np.asarray(s.material_vertices)@axes.T+origin,s.material_faces))
    design=MeshSnapshot((r.design_vertices,r.design_faces))
    display_subset=native().check_exact_material_subset(display.exact_region(),design.exact_region())
    binding=dict(schema='spin2014-saved-exact-boolean-replay-v2',
        coefficients_are_normalized=True,
        evaluator_settings=evaluator_settings.identity_dict(),saved_geometry_hash=s.final_geometry_hash,
        replay_geometry_hash=replay.geometry_hash,boolean_backend=replay.boolean_backend_identity,
        operation=operation,rounded_oriented_boundary_identical=True,
        target_operand_digest=digest({'vertices':design_target.tolist(),'faces':r.design_faces}),
        source_operand_digest=digest({'vertices':replay.interface_vertices,'faces':replay.interface_faces}),
        floating_derived_containment=display_subset,
        floating_derived_containment_verified=bool(display_subset['completed'] and display_subset['a_subset_b']))
    return retained,binding


def _rational_spin2014_saved_retained(run, evaluator_settings):
    """Consume the saved same-call authority, never a floating reconstruction.

    Re-evaluation authenticates its complete source recipe and provenance. The
    retained region itself comes from the saved rational tape. This local
    geometry check grants neither solver convergence nor manufacturing status.
    """
    from .mass_moments import integrate_exact_rational_raw_moments
    request, solved = run.measured_request, run.solve
    authority = getattr(solved, 'rational_authority', None)
    require(isinstance(authority, dict) and
            authority.get('schema') == 'rational-quadric-production-authority-v1' and
            authority.get('coordinate_space') == 'TARGET_FRAME' and
            authority.get('analytic_coordinates_exact') is False and
            authority.get('variable_moment_semantics') == 'DENSITY_WEIGHTED' and
            authority.get('geometry_hash') == solved.final_geometry_hash,
            'missing or mismatched saved rational production authority')
    faces = solved.material_faces
    raw = integrate_exact_rational_raw_moments(authority['material_vertices'], faces)
    density = Fraction.from_float(scalar(request.material_density, 'material density', positive=True))
    require(tuple(density*x for x in raw) == tuple(
                fraction_from_rational_text(x) for x in authority['variable_moments']),
            'saved authority independent exact moments or density mismatch')
    replay = evaluate_spin_it_faster_quadric(request, solved.quadric_coefficients,
        evaluator_settings=evaluator_settings, coefficients_are_normalized=True)
    require(replay.geometry_hash == solved.final_geometry_hash and
            canonical(replay.rational_authority) == canonical(authority) and
            tuple(tuple(f) for f in replay.material_faces) == tuple(tuple(f) for f in faces),
            'saved rational authority differs from actual source, geometry or provenance')
    rounded = np.asarray([[float(fraction_from_rational_text(x)) for x in row]
                          for row in authority['material_vertices']])
    require(np.array_equal(rounded, np.asarray(solved.material_vertices)),
            'saved display vertices do not derive from authority')
    model, affine = _target_rationals_to_normalized_model(
        authority['material_vertices'], request, exact_inverse=True, return_binding=True)
    design_model = _target_rationals_to_normalized_model(
        authority['design_vertices'], request, exact_inverse=True)
    expected_design = tuple(tuple(fraction_to_rational_text(Fraction.from_float(float(x)))
                                  for x in row) for row in request.design_vertices)
    require(design_model == expected_design,
            'authority design does not invert to actual input coordinates')
    retained = _canonical_exact_boolean_boundary(model, faces)
    # A rejected display is a separate diagnostic, not a reason to replace or
    # weld the exact material. Only ordinary geometry rejection is caught.
    try:
        display_exact = _target_rationals_to_normalized_model(tuple(tuple(
            fraction_to_rational_text(Fraction.from_float(float(x))) for x in row)
            for row in solved.material_vertices), request, exact_inverse=True)
        display = MeshSnapshot((np.asarray([[float(fraction_from_rational_text(x))
                    for x in row] for row in display_exact]), faces))
        design = MeshSnapshot((request.design_vertices, request.design_faces))
        subset = native().check_exact_material_subset(display.exact_region(), design.exact_region())
        floating = dict(status='CHECKED', containment=subset,
            verified=bool(subset['completed'] and subset['a_subset_b']))
    except ValueError as error:
        floating = dict(status='REJECTED', verified=False, diagnostic=str(error))
    binding = dict(schema='spin2014-saved-rational-authority-v1',
        coefficients_are_normalized=True, evaluator_settings=evaluator_settings.identity_dict(),
        saved_geometry_hash=solved.final_geometry_hash, replay_geometry_hash=replay.geometry_hash,
        boolean_backend=replay.boolean_backend_identity,
        source_construction_identity=authority['source_construction_identity'],
        authority_digest=digest(authority), affine=affine,
        independent_exact_moments_identical=True,
        floating_derived=floating, manufacturing_approved=False)
    return retained, binding


def _check_spin2014_rotational_ledger(solved, settings):
    """Reject truncated attempt/evaluation evidence before geometry replay.

    This checks accounting, not stationarity, geometry or initialization
    independence. It never substitutes a success flag for those other gates.
    """
    def check(condition, detail):
        require(condition, 'rotational ledger: '+detail)
    counts=solved.restart_function_evaluations
    check(type(solved.restart_count) is int and
          1<=solved.restart_count<=settings.maximum_restarts and
          len(counts)==solved.restart_count,'missing or over-budget attempts')
    check(all(type(x) is int and 1<=x<=settings.maximum_function_evaluations_per_restart
              for x in counts),'invalid per-attempt evaluation count')
    check(type(solved.function_evaluations) is int and
          sum(counts)==solved.function_evaluations,'total evaluation count differs')
    check(bool(solved.events),'missing events')
    coefficient_hash=_native_module().spin_it_faster_quadric_hash
    total=0;current=-1;local=0;active=False;terminals=[];evaluations=[]
    for event in solved.events:
        check(event.equation_profile==settings.equation_profile,
              'event equation profile differs')
        if event.type=='attempt_start':
            check(not active and event.restart==current+1 and event.restart<len(counts),
                  'duplicate, missing or overlapping attempt start')
            current=event.restart;local=0;active=True
            check(event.outcome=='started' and event.evaluation_index==total and
                  event.restart_evaluation_index==0,'invalid start counters')
        else:
            check(active and event.restart==current,'event outside its attempt')
            if event.type=='evaluation':
                total+=1;local+=1
                check(local<=counts[current] and event.outcome in ('accepted','rejected','failed'),
                      'extra evaluation or unclassified outcome')
                check(bool(event.reason),'missing evaluation diagnostic')
                if settings.capture_evaluation_payload and event.outcome in ('accepted','rejected'):
                    check(event.evaluation_payload is not None,'missing returned evaluation payload')
                evaluations.append((current,total,local))
            elif event.type=='attempt_terminal':
                check(event.outcome=='terminated' and bool(event.reason) and
                      local==counts[current],'missing evaluation or terminal diagnostic')
                terminals.append(event.reason);active=False
            elif event.type not in ('newton_direction','newton_linear_solve'):
                check(False,'unrecognized event type')
            check(event.evaluation_index==total and event.restart_evaluation_index==local,
                  'noncontiguous evaluation counters')
        if event.type in ('attempt_start','evaluation'):
            q=event.quadric_coefficients
            check(q is not None and len(q)==10 and all(math.isfinite(float(x)) for x in q)
                  and bool(event.quadric_hash),'missing finite attempted coefficients')
            try:
                actual_hash=coefficient_hash(np.asarray(q,dtype=float))
            except (TypeError,ValueError) as error:
                raise ValueError('rotational ledger: invalid coefficient identity') from error
            check(event.quadric_hash==actual_hash,'coefficient identity differs')
    check(not active and len(terminals)==len(counts) and total==solved.function_evaluations,
          'truncated final attempt')
    check(solved.status=='CONVERGED' and terminals[-1]=='converged' and
          all(reason!='converged' for reason in terminals[:-1]),
          'terminal status does not match attempts')
    timings=solved.evaluation_timings
    check(len(timings)==len(evaluations),'missing evaluation timings')
    check([(row.restart,row.evaluation_index,row.restart_evaluation_index) for row in timings]
          ==evaluations,'timing ledger counters differ')


def prepare_spin2014_quadric_material_output(run, *, evaluator_settings=None):
    """Recheck one saved rotational endpoint without optimizing or repairing.

    The authority is the exact construction bound to the actual saved retained
    boundary, not a second evaluation of -q. Build DESIGN-minus-that-boundary
    and use existing exact partitions and independent tetrahedron integration.
    A failed floating derivation remains a failure; no vertex is moved/welded.
    """
    from ..spin import Spin2014QuadricMaterialRun
    require(isinstance(run,Spin2014QuadricMaterialRun),'expected a rotational quadric run')
    request,solved,objective=run.measured_request,run.solve,run.objective_settings
    require(solved.status=='CONVERGED' and solved.final_geometry_verified,
            'rotational material requires a verified converged terminal result')
    modes={'YOYO':'YO_YO','TOP_INERTIA_ONLY':'SPINNING_TOP',
           'TOP_LOWERING_ONLY':'SPINNING_TOP','TOP_MASS_AND_LOWERING':'SPINNING_TOP'}
    require(run.objective_profile=='SPIN2014_OBJECTIVES_QUADRIC_V1' and
            modes.get(objective.objective)==request.mode,'wrong rotational objective profile/mode')
    # A saved result does not inherit live-native admission. Recheck the same
    # protocol before reading its terminal label/moments; in particular a
    # caller-modified tolerance cannot turn this attachment into a weaker gate.
    require(request.settings.equation_profile=='DUAL_MONOMIAL_RAW_V1' and
            solved.equation_profile=='DUAL_MONOMIAL_RAW_V1' and
            request.settings.equation_tolerance<=1e-12 and
            request.settings.maximum_restarts<=8 and
            request.settings.maximum_function_evaluations_per_restart<=60,
            'saved rotational source violates the strict raw protocol/tolerance')
    _check_spin2014_rotational_ledger(solved,request.settings)
    fixed_density=scalar(run.fixed_material_density,'fixed density',positive=True)
    design_density=scalar(request.material_density,'design density',positive=True)
    origin=np.asarray(request.target_origin,dtype=float);axes=np.asarray(request.target_axes,dtype=float)
    require(origin.shape==(3,) and axes.shape==(3,3) and np.isfinite(origin).all() and
            np.isfinite(axes).all() and np.max(abs(axes.T@axes-np.eye(3)))<=1e-12 and
            np.linalg.det(axes)>0.,'invalid target coordinate frame')
    outer=MeshSnapshot((request.outer_vertices,request.outer_faces))
    design=MeshSnapshot((request.design_vertices,request.design_faces))
    require(solved.material_vertices is not None and solved.material_faces is not None,
            'missing actual retained material')
    retained,exact_binding=_exact_spin2014_saved_retained(run,evaluator_settings)
    n=native()
    def subset(one,two,label):
        checked=n.check_exact_material_subset(one.exact_region(),two.exact_region())
        require(checked['completed'] and checked['a_subset_b'],label+': exact containment failed')
    subset(design,outer,'design');subset(retained,design,'retained')
    from fractions import Fraction
    dv,df=design.arrays();_,rf=retained.arrays()
    design_exact=design.rational_vertices or tuple(tuple(
        fraction_to_rational_text(Fraction.from_float(float(x))) for x in row) for row in dv)
    require(retained.rational_vertices is not None,'missing authoritative retained coordinates')
    exact=n.libigl_exact_boolean_construction_from_rationals(
        design_exact,df,retained.rational_vertices,rf,'difference')
    require(exact['succeeded'],'same-retained exact complement construction failed')
    cavity=(_canonical_exact_boolean_boundary(exact['vertices'],exact['faces'])
            if len(exact['faces']) else MeshSnapshot(None))
    fixed=_nested_material_boundary(outer,design)
    material=_nested_material_boundary(outer,cavity)
    def partition(whole,parts,label):
        checked=n.check_exact_material_partition(whole.exact_region(),[x.exact_region() for x in parts])
        require(checked['completed'] and checked['complete_partition'] and
                checked['no_positive_volume_overlap'] and all(checked['part_subset_whole']),
                label+': exact partition failed')
        return checked
    partitions={'design_partition':partition(design,[retained,cavity],'design'),
                'whole_partition':partition(outer,[fixed,retained,cavity],'whole'),
                'material_partition':partition(material,[fixed,retained],'material')}
    def moments(mesh,density):
        return np.zeros(10) if mesh.explicit_empty else np.asarray(_raw_moments(mesh,density))
    fixed_model=moments(fixed,fixed_density);retained_model=moments(retained,design_density)
    total_model=fixed_model+retained_model
    rotation=axes.T;translation=-rotation@origin
    total_target=np.asarray(_transform_moments(total_model,rotation,translation))
    retained_target=np.asarray(_transform_moments(retained_model,rotation,translation))
    length=float(np.max(np.ptp((outer.arrays()[0]-origin)@axes,axis=0)))
    require(math.isfinite(length) and length>0.,'invalid fixed input normalization')
    scales=max(fixed_density,design_density)*length**np.array([3,4,4,4,5,5,5,5,5,5])
    require(np.isfinite(scales).all() and (scales>0).all(),'invalid moment normalization')
    errors=[]
    for actual,saved in ((fixed_model,request.fixed_shell_moments),
                         (retained_target,solved.final_variable_raw_moments),
                         (total_target,solved.final_total_raw_moments),
                         (total_target,solved.independent_total_raw_moments)):
        recorded=np.asarray(saved,dtype=float)
        require(recorded.shape==(10,) and np.isfinite(recorded).all(),'missing finite saved mass moments')
        errors.append(float(np.max(abs(actual-recorded)/scales)))
    require(max(errors)<=1e-6,'geometry/density does not reproduce recorded material moments')
    require(total_target[0]>0. and np.isfinite(total_target).all(),'nonpositive or nonfinite material mass')
    second=np.array([[total_target[4],total_target[9],total_target[8]],
                     [total_target[9],total_target[5],total_target[7]],
                     [total_target[8],total_target[7],total_target[6]]])
    centered=second-np.outer(total_target[1:4],total_target[1:4])/total_target[0]
    target_inertia=np.trace(centered)*np.eye(3)-centered
    require(target_inertia[2,2]>1e-14,'axial inertia is zero or unqualified')
    f,constraints=_independent_spin2014_rotational_physics(total_target,objective)
    # This checker never calls the native objective/AD implementation. Complex
    # step avoids a subtractive finite-difference noise floor at the raw gate.
    gradient=np.zeros(10);jacobian=np.zeros((len(constraints),10))
    for i in range(10):
        perturbed=total_target.astype(complex);perturbed[i]+=1e-24j
        fi,ci=_independent_spin2014_rotational_physics(perturbed,objective)
        gradient[i]=np.imag(fi)/1e-24;jacobian[:,i]=np.imag(ci)/1e-24
    multipliers=np.asarray(solved.multipliers,dtype=float)
    q=np.asarray(solved.quadric_coefficients,dtype=float)
    require(multipliers.shape==constraints.shape and np.isfinite(multipliers).all() and
            q.shape==(10,) and np.isfinite(q).all(),'missing terminal coefficients/multipliers')
    mu=scalar(solved.multiplier_mu,'mu',positive=True)
    dual=q.copy();dual[7:]*=2;dual/=np.linalg.norm(dual)
    raw=float(np.linalg.norm(np.r_[gradient-jacobian.T@multipliers+mu*dual,constraints]))
    require(math.isfinite(raw) and raw<=request.settings.equation_tolerance,
            'independent original-objective raw residual fails terminal tolerance; witness '+canonical(dict(
                independent_raw_residual=raw,terminal_tolerance=request.settings.equation_tolerance,
                independently_integrated_total_moments=total_target.tolist(),
                saved_total_moments=list(solved.final_total_raw_moments),
                normalized_moment_errors=errors,
                independent_objective_gradient=gradient.tolist(),constraint_residuals=constraints.tolist())))
    require(solved.trace and math.isfinite(float(f)) and
            abs(float(f)-solved.trace[-1].objective)<=1e-10*max(1.,abs(float(f))),
            'independent original objective differs from the saved terminal')
    unit_mass=objective.solver_mass_unit_in_declared_units
    unit_length=objective.solver_length_unit_in_declared_units
    constraint_scales=np.array(([scales[1]*unit_mass*unit_length]*3 if objective.objective=='YOYO'
                               else [scales[1]*unit_mass*unit_length]*2)+
                              [scales[4]*unit_mass*unit_length**2]*2)
    normalized_constraints=float(np.max(abs(constraints)/constraint_scales))
    require(normalized_constraints<=1e-6,'independent physical constraints fail')
    cavity_moments=moments(cavity,1.)
    audit=dict(schema='spin2014-rotational-material-output-audit-v1',
        scope='EXACT_REGION_PARTITION_AND_INDEPENDENT_MOMENTS_NOT_MANUFACTURING',
        complement_construction='DESIGN_MINUS_SAVED_RETAINED_EXACT_V1',
        coordinate_map='SAVED_EXACT_TARGET_TO_MODEL_RATIONAL_AFFINE_V1',
        exact_replay=exact_binding,fixed_density=fixed_density,design_density=design_density,
        fixed_moments_model=fixed_model.tolist(),retained_moments_model=retained_model.tolist(),
        total_moments_model=total_model.tolist(),total_moments_target=total_target.tolist(),
        cavity_geometric_moments=cavity_moments.tolist(),center_of_mass_model=(total_model[1:4]/total_model[0]).tolist(),
        center_inertia_model=(axes@target_inertia@axes.T).tolist(),original_objective=float(f),
        independent_objective_gradient=gradient.tolist(),independent_constraint_jacobian=jacobian.tolist(),
        independent_raw_residual=raw,normalized_constraint_residual=normalized_constraints,
        normalized_moment_error=max(errors),independent_moments_verified=True,
        accepted_updates=sum(row.accepted_geodesic_step>0 for row in solved.trace),
        actual_removed_volume=float(cavity_moments[0]),manufacturing_approved=False,sealed_artifact=False,**partitions)
    identity=dict(schema='spin2014-rotational-material-output-v1',objective_profile=run.objective_profile,
        objective_identity=run.objective_identity,source_run_sha256=digest(_json_value(asdict(run))),
        meshes={name:mesh.identity for name,mesh in (('outer',outer),('fixed',fixed),('design',design),
            ('retained',retained),('cavity',cavity),('material',material))},audit=audit)
    # The full source binding retains every timing and diagnostic. It is not
    # itself a scientific fact: two independent processes can measure different
    # costs while producing the same material. Do not discard the source record
    # or weaken its separate integrity check when comparing material science.
    scientific_identity={key:value for key,value in identity.items() if key!='source_run_sha256'}
    result=object.__new__(Spin2014QuadricMaterialOutput)
    for name,value in (('outer',outer),('fixed',fixed),('design',design),('retained',retained),
                       ('cavity',cavity),('material',material),('_audit',canonical(_json_value(audit))),
                       ('_identity',canonical(_json_value(identity))),
                       ('scientific_digest',digest(_json_value(scientific_identity)))):
        object.__setattr__(result,name,value)
    return result


def _exact_partition_audit(request, regions):
    n = native()
    outer = MeshSnapshot((request.outer_vertices, request.outer_faces)).exact_region()
    parts = [region.exact_region() for _, region, _ in regions]
    audit = n.check_exact_material_partition(outer, parts)
    require(audit['completed'] and audit['complete_partition'] and
            audit['no_positive_volume_overlap'] and
            all(audit['part_subset_whole']),
            'exact geometry check failed: ' + canonical(audit))
    return audit


def compose_sealed_faster_material_as_fixed_region(
        request_b, source_artifact, *, region_prefix,
        source_to_request_rotation, source_to_request_translation,
        source_length_to_request_unit, reference_density,
        source_length_unit, target_length_unit,
        source_mass_unit, target_mass_unit, existing_fixed_regions=()):
    """Freeze sealed object A as actual fixed material in a Faster request for B.

    ``existing_fixed_regions`` contains ``(region_id, vertices, faces,
    density_ratio)`` tuples already expressed in B request coordinates.  Their
    independently integrated moments must match the request's prior fixed
    moment claim.  The native exact partition audit then proves the represented
    regions (A, other fixed regions and B's design domain) exactly tile B's
    outer domain without positive-volume overlap.
    """
    require(isinstance(request_b, SpinItFasterRequest), 'expected SpinItFasterRequest for target B')
    require(isinstance(source_artifact, SealedSpinArtifact), 'expected sealed source material artifact A')
    source_artifact.revalidate()
    region_prefix = text_label(region_prefix, 'region_prefix')
    source_length_unit = text_label(source_length_unit, 'source_length_unit')
    target_length_unit = text_label(target_length_unit, 'target_length_unit')
    source_mass_unit = text_label(source_mass_unit, 'source_mass_unit')
    target_mass_unit = text_label(target_mass_unit, 'target_mass_unit')
    require(source_length_unit == target_length_unit and source_mass_unit == target_mass_unit,
            'source and target units must match')
    scale = scalar(source_length_to_request_unit, 'source_length_to_request_unit', positive=True)
    require(scale == 1.0, 'matching unit labels require a unit-scale rigid transform')
    rho_ref = scalar(reference_density, 'reference_density', positive=True)
    rotation = np.asarray(source_to_request_rotation, dtype=float)
    require(rotation.shape == (3, 3) and np.isfinite(rotation).all(),
            'source-to-request transform must be a proper rigid rotation')
    require(np.allclose(rotation.T @ rotation, np.eye(3), rtol=0., atol=1e-12) and
            abs(float(np.linalg.det(rotation)) - 1.) <= 1e-12,
            'source-to-request transform must be a proper rigid rotation')
    translation = np.asarray(vector(source_to_request_translation, 'source_to_request_translation'))

    source_problem = source_artifact.problem
    require(source_problem.length_unit == source_length_unit and
            source_problem.mass_unit == source_mass_unit,
            'sealed source artifact units must match declared source units')
    design_density = {r.region_id: r.density for r in source_problem.design_regions}
    source_regions = []
    for region in source_problem.fixed_regions:
        if not region.mesh.explicit_empty:
            source_regions.append((region.region_id, region.mesh, region.density))
    seen_design = set()
    for region_id, retained, _ in source_artifact.design_results:
        require(region_id in design_density and region_id not in seen_design,
                'sealed source design-region identity mismatch')
        seen_design.add(region_id)
        if not retained.explicit_empty:
            source_regions.append((region_id, retained, design_density[region_id]))
    require(seen_design == set(design_density), 'sealed source material omits a design region')
    require(bool(source_regions), 'sealed source has no nonempty material regions')

    imported = []
    imported_moments = []
    imported_identities = []
    for region_id, mesh, physical_density in source_regions:
        vertices, faces = mesh.arrays()
        mapped = scale * (vertices @ rotation.T) + translation
        mapped_mesh = MeshSnapshot((mapped, faces))
        density_ratio = float(physical_density / rho_ref)
        require(math.isfinite(density_ratio) and density_ratio > 0.,
                'transferred material density ratio must be finite and positive')
        new_id = region_prefix + ':' + region_id
        imported.append((new_id, mapped_mesh, density_ratio))
        imported_moments.append(_raw_moments(mapped_mesh, density_ratio))
        imported_identities.append({'region_id': new_id, 'source_region_id': region_id,
                                    'mesh': mapped_mesh.identity, 'density_ratio': density_ratio})

    existing = []
    existing_moments = []
    existing_identities = []
    identifiers = {entry[0] for entry in imported}
    for index, item in enumerate(tuple(existing_fixed_regions)):
        require(isinstance(item, (tuple, list)) and len(item) == 4,
                'existing fixed region must be (id, vertices, faces, density_ratio)')
        region_id = text_label(item[0], f'existing_fixed_regions[{index}].id')
        require(region_id not in identifiers, 'duplicate fixed material region id')
        identifiers.add(region_id)
        mesh = MeshSnapshot((item[1], item[2]))
        rho = scalar(item[3], f'existing_fixed_regions[{index}].density_ratio', positive=True)
        existing.append((region_id, mesh, rho))
        existing_moments.append(_raw_moments(mesh, rho))
        existing_identities.append({'region_id': region_id, 'mesh': mesh.identity, 'density_ratio': rho})

    expected_existing = _add_moments(*existing_moments) if existing_moments else (0.,) * 10
    # Geometry and fixed moments are in B's model coordinates. Native evaluation
    # transforms the fixed moments to the target frame exactly once.
    # Compare with the explicitly supplied geometry; caller values alone are
    # never accepted as evidence.
    root_side = max(float(np.ptp(np.asarray(request_b.outer_vertices), axis=0).max()), 1e-300)
    moment_scales = (root_side**3, root_side**4, root_side**4, root_side**4,
                     root_side**5, root_side**5, root_side**5, root_side**5,
                     root_side**5, root_side**5)
    fixed_residual = max(abs(expected_existing[i] - request_b.fixed_shell_moments[i]) /
                         moment_scales[i] for i in range(10))
    require(fixed_residual <= 1e-10,
            'existing fixed geometry does not match request fixed-shell moments')

    design = MeshSnapshot((request_b.design_vertices, request_b.design_faces))
    partition_regions = [(rid, mesh, density) for rid, mesh, density in existing + imported]
    partition_regions.append(('design:B', design, request_b.material_density))
    geometry_audit = _exact_partition_audit(request_b, partition_regions)
    total_fixed = _add_moments(expected_existing, _add_moments(*imported_moments))
    new_request = replace(request_b, fixed_shell_moments=total_fixed)
    binding = {
        'schema': 'faster-floating-island-material-transfer-v1',
        'status': 'FROZEN_AS_FIXED_GEOMETRY',
        'optimization_claim': 'NOT_ASSERTED_BY_GEOMETRY_TRANSFER',
        'source_artifact_scientific_digest': source_artifact.scientific_digest,
        'source_artifact_evidence_digest': source_artifact.evidence_digest,
        'source_final_material': source_artifact.final_material.identity,
        'source_material_regions': imported_identities,
        'existing_fixed_regions': existing_identities,
        'existing_fixed_moment_residual_normalized': fixed_residual,
        'rigid_transform': {'rotation': rotation.tolist(), 'translation': translation.tolist(),
                            'scale': scale, 'source_length_unit': source_length_unit,
                            'target_length_unit': target_length_unit},
        'density_reference': {'reference_density': rho_ref,
                              'source_mass_unit': source_mass_unit,
                              'target_mass_unit': target_mass_unit},
        'exact_geometry_check': geometry_audit,
        'request_digest': digest(asdict(new_request)),
    }
    result = object.__new__(FrozenFasterMaterialComposition)
    object.__setattr__(result, 'request', new_request)
    object.__setattr__(result, '_binding', canonical(binding))
    return result


def seal_buoyant_plane_result_as_artifact(
        problem, request, solved, *, displacement_vertices, displacement_faces,
        waterline_origin, waterline_normal):
    """Seal a converged prescribed-waterline plane solve from actual geometry.

    The request must describe exactly the locked outer and design meshes in
    ``problem``. Fixed material moments, displaced-fluid geometry, both sides
    of the optimized plane, and the final solver claim are independently
    recomputed. Only a converged, equilibrium-checked, positive-stability
    result can produce a seal. This records prescribed-waterline numerical
    stability only; it is not a free-waterline or physical validation.
    """
    require(isinstance(problem, SpinArtifactProblem), 'expected locked SpinArtifactProblem')
    require(isinstance(request, SpinItFasterRequest), 'expected SpinItFasterRequest')
    require(request.mode == 'BUOYANT_STABILITY', 'A request must use BUOYANT_STABILITY mode')
    problem = _recheck_problem(problem)
    require(problem.target['mode'] == 'BUOYANT_STABILITY',
            'locked A problem target must be BUOYANT_STABILITY')
    require(len(problem.design_regions) == 1,
            'buoyant plane sealing currently requires exactly one design region')

    def same_mesh(snapshot, vertices, faces, label):
        actual_v, actual_f = snapshot.arrays()
        expected_v = np.asarray(vertices, dtype='<f8')
        expected_f = np.asarray(faces, dtype='<i4')
        require(actual_v.shape == expected_v.shape and actual_f.shape == expected_f.shape and
                actual_v.tobytes() == expected_v.tobytes() and
                actual_f.tobytes() == expected_f.tobytes(),
                f'{label} differs from locked A problem geometry')

    same_mesh(problem.original, request.outer_vertices, request.outer_faces, 'request outer')
    same_mesh(problem.design_regions[0].mesh, request.design_vertices,
              request.design_faces, 'request design')
    require(problem.design_regions[0].density == request.material_density,
            'request design density differs from locked A problem')
    require(request.original_units_per_normalized_unit == 1.0,
            'buoyant seal requires request coordinates in the locked problem units')

    axes = np.asarray(request.target_axes, dtype=np.float64)
    origin = np.asarray(request.target_origin, dtype=np.float64)
    require(np.allclose(axes.T @ axes, np.eye(3), rtol=0.0, atol=1e-10) and
            abs(float(np.linalg.det(axes)) - 1.0) <= 1e-10,
            'request target frame must be a proper rigid frame')

    fixed_target = []
    fixed_volume = 0.0
    for region in problem.fixed_regions:
        vertices, faces = region.mesh.arrays()
        target_vertices = (vertices - origin) @ axes
        fixed_target.append(_raw_moments(
            MeshSnapshot((target_vertices, faces)), region.density))
        fixed_volume += _raw_moments(region.mesh, 1.0)[0]
    fixed_moments = _add_moments(*fixed_target) if fixed_target else (0.0,) * 10
    side = max(float(np.ptp(np.asarray(request.outer_vertices) @ axes, axis=0).max()), 1e-300)
    moment_scales = (side**3, side**4, side**4, side**4,
                     side**5, side**5, side**5, side**5, side**5, side**5)
    request_fixed_target = _transform_moments(request.fixed_shell_moments,
        axes.T, -axes.T @ origin)
    fixed_error = max(abs(fixed_moments[i] - request_fixed_target[i]) /
                      moment_scales[i] for i in range(10))
    require(fixed_error <= 1e-10,
            'request fixed moments do not match independently integrated A materials')
    require(abs(fixed_volume - request.fixed_shell_volume) <=
            1e-10 * max(1.0, fixed_volume),
            'request fixed-shell volume does not match actual A fixed material')

    displacement = MeshSnapshot((displacement_vertices, displacement_faces))
    disp_v, disp_f = displacement.arrays()
    displacement_target = (disp_v - origin) @ axes
    waterline_world_origin = np.asarray(_finite_tuple(
        waterline_origin, 3, 'waterline_origin'), dtype=np.float64)
    waterline_world_normal = np.asarray(_finite_tuple(
        waterline_normal, 3, 'waterline_normal'), dtype=np.float64)
    normal_length = float(np.linalg.norm(waterline_world_normal))
    require(math.isfinite(normal_length) and normal_length > 0.0,
            'waterline normal must be nonzero and finite')
    normal_target = axes.T @ (waterline_world_normal / normal_length)
    require(np.allclose(normal_target, (0.0, 0.0, 1.0), rtol=0.0, atol=1e-10),
            'waterline normal must align with target-frame vertical')
    waterline_target_origin = (waterline_world_origin - origin) @ axes
    measured = compute_hydrostatic_geometry(
        displacement_target.tolist(), disp_f.tolist(),
        waterline_origin=waterline_target_origin,
        waterline_normal=(0.0, 0.0, 1.0))
    hydro_scale = max(float(request.displaced_liquid_volume), side**3, 1e-300)
    require(abs(measured.displaced_volume - request.displaced_liquid_volume) <=
            1e-10 * hydro_scale and
            abs(measured.waterline_area - request.waterline_area) <= 1e-10 * max(side**2, 1e-300) and
            np.allclose(measured.center_of_buoyancy,
                        request.center_of_buoyancy_target, rtol=0.0,
                        atol=1e-10 * max(side, 1e-300)) and
            np.allclose(measured.waterline_central_second_moment,
                        request.waterline_central_second_moment, rtol=0.0,
                        atol=1e-10 * max(side**4, 1e-300)),
            'request hydrostatics differ from the supplied displacement mesh and waterline')

    require(isinstance(solved, dict) and solved.get('status') == 'CONVERGED',
            'A buoyant solve must converge before sealing')
    require(solved.get('equilibrium_constraints_verified') is True and
            solved.get('eq7_numerical_conditions_verified') is True,
            'A buoyant solve lacks equilibrium or finite-precision Eq. 7 verification')
    require(solved.get('stable_at_prescribed_waterline') is True and
            float(solved.get('stability_margin', float('-inf'))) > 0.0,
            'A solve must have positive prescribed-waterline stability')
    coefficients = np.asarray(solved.get('plane_coefficients'), dtype=np.float64)
    require(coefficients.shape == (4,) and np.isfinite(coefficients).all(),
            'converged A result lacks finite plane coefficients')
    split = evaluate_spin_it_faster_plane_halfspaces(
        request, coefficients=coefficients)
    retained_target = split['positive']
    carved_target = split['negative']
    solver_boundary = np.asarray(solved.get('variable_material_vertices'), dtype=np.float64)
    solver_faces = np.asarray(solved.get('variable_material_faces'), dtype=np.int32)
    require(solver_boundary.shape == retained_target['vertices'].shape and
            solver_faces.shape == retained_target['faces'].shape and
            solver_boundary.tobytes() == retained_target['vertices'].tobytes() and
            solver_faces.tobytes() == retained_target['faces'].tobytes(),
            'A solver boundary differs from the independently reevaluated plane material')

    def model_mesh(target_mesh):
        return MeshSnapshot((target_mesh['vertices'] @ axes.T + origin,
                             target_mesh['faces']))

    retained = model_mesh(retained_target)
    carved = model_mesh(carved_target)
    outer_vertices, outer_faces = problem.original.arrays()
    cavity_vertices, cavity_faces = carved.arrays()
    final_vertices = np.vstack((outer_vertices, cavity_vertices))
    final_faces = np.vstack((outer_faces,
                             cavity_faces[:, ::-1] + len(outer_vertices)))

    fixed_plus_design = _add_moments(
        fixed_moments,
        _raw_moments(MeshSnapshot((retained_target['vertices'], retained_target['faces'])),
                     problem.design_regions[0].density))
    solver_total = tuple(float(x) for x in solved.get('total_raw_moments', ()))
    require(len(solver_total) == 10 and all(math.isfinite(x) for x in solver_total),
            'A solver result lacks ten finite total mass moments')
    total_error = max(abs(fixed_plus_design[i] - solver_total[i]) /
                      moment_scales[i] for i in range(10))
    require(total_error <= 1e-10,
            'A solver moments differ from independently integrated fixed and retained materials')

    frame = problem.normalization
    normalized_claim = _transform_moments(
        _transform_moments(solver_total, axes, origin),
        np.eye(3) / frame['L'], -np.asarray(frame['origin']) / frame['L'])
    normalized_claim = [x / (frame['L']**3 * frame['R']) for x in normalized_claim]
    artifact = SealedSpinArtifact.create(
        problem, outer=problem.original,
        design_results=[(problem.design_regions[0].region_id, retained, carved)],
        final_material=(final_vertices, final_faces), newly_carved=carved,
        solver_claim={'frame': FRAME, 'order': list(ORDER),
                      'moments': normalized_claim})
    identity = {
        'schema': 'faster-buoyant-plane-a-seal-v1',
        'operation': 'SEAL_BUOYANT_PLANE_A',
        'scope': 'SYNTHETIC_METHOD_VALIDATION',
        'problem_digest': problem.input_digest,
        'request_digest': request.request_digest,
        'request': _json_value(asdict(request)),
        'plane_coefficients_target': coefficients.tolist(),
        'target_axes': axes.tolist(),
        'target_origin': origin.tolist(),
        'displacement_mesh': displacement.identity,
        'waterline': {'origin_world': waterline_world_origin.tolist(),
                      'normal_world': (waterline_world_normal / normal_length).tolist(),
                      'origin_target': waterline_target_origin.tolist()},
        'recomputed_hydrostatics': {
            'displaced_volume': measured.displaced_volume,
            'center_of_buoyancy_target': list(measured.center_of_buoyancy),
            'waterline_area': measured.waterline_area,
            'waterline_central_second_moment': [list(row) for row in
                                                 measured.waterline_central_second_moment],
        },
        'solver_result_digest': digest(_json_value(solved)),
        'fixed_moment_error': fixed_error,
        'independent_moment_error': total_error,
        'prescribed_waterline_stable': True,
        'stability_margin': float(solved['stability_margin']),
        'equilibrium_constraints_verified': True,
        'eq7_numerical_conditions_verified': True,
        'plane_partition': split['exact_partition_audit'],
        'plane_interface_stitch': split['interface_stitch'],
        'plane_mesh_independent_moment_error': split['independent_mesh_moment_error'],
        'plane_partition_independent_moment_error': split['independent_partition_moment_error'],
        'retained_region': retained.identity,
        'carved_region': carved.identity,
        'artifact_scientific_digest': artifact.scientific_digest,
        'artifact_evidence_digest': artifact.evidence_digest,
        'free_waterline_validation': 'NOT_RUN',
    }
    return _result('CONVERGED', identity, artifact)


def solve_floating_island_b_from_sealed_a(
        request_b, source_artifact, *, region_prefix,
        source_to_request_rotation, source_to_request_translation,
        source_length_to_request_unit, reference_density,
        source_length_unit, target_length_unit,
        source_mass_unit, target_mass_unit, existing_fixed_regions=(),
        initial_coefficients, maximum_iterations=60,
        maximum_line_search_steps=24, constraint_tolerance=1e-10):
    """Freeze sealed A into B, solve B's static plane family, and audit its output.

    This is a computational method adapter, not evidence for an author
    Floating-Island asset or experiment. Every fixed-region moment is bound
    through :func:`compose_sealed_faster_material_as_fixed_region`. The initial
    fixed/A/B-design domains must exactly partition B's outer domain; after
    solving, B's retained mesh is independently integrated and checked as an
    exact subset of its design domain. Carved-away void is not misreported as
    material. A chain reports ``CONVERGED`` only after these checks and the
    solver's optimality flag pass.
    """
    existing_fixed_regions = tuple(existing_fixed_regions)
    composition = compose_sealed_faster_material_as_fixed_region(
        request_b, source_artifact, region_prefix=region_prefix,
        source_to_request_rotation=source_to_request_rotation,
        source_to_request_translation=source_to_request_translation,
        source_length_to_request_unit=source_length_to_request_unit,
        reference_density=reference_density,
        source_length_unit=source_length_unit, target_length_unit=target_length_unit,
        source_mass_unit=source_mass_unit, target_mass_unit=target_mass_unit,
        existing_fixed_regions=existing_fixed_regions)
    solved = solve_spin_it_faster_static_plane_family(
        composition.request, initial_coefficients=initial_coefficients,
        maximum_iterations=maximum_iterations,
        maximum_line_search_steps=maximum_line_search_steps,
        constraint_tolerance=constraint_tolerance)

    identity = {
        'schema': 'faster-floating-island-a-to-b-solve-v1',
        'scope': 'SYNTHETIC_METHOD_VALIDATION',
        'source_artifact_scientific_digest': source_artifact.scientific_digest,
        'source_artifact_evidence_digest': source_artifact.evidence_digest,
        'composition_binding': composition.binding,
        'request_digest': composition.request.request_digest,
        'solver_status': solved.get('status'),
        'solver_diagnostic': solved.get('diagnostic'),
        'solver_iterations': solved.get('iterations'),
        'solver_accepted_steps': solved.get('accepted_steps'),
        'normalized_constraints_norm': solved.get('normalized_constraints_norm'),
        'independent_total_moment_residual': None,
        'independent_normalized_constraints_norm': None,
        'final_material_domain_audit': {
            'status': 'NOT_RUN_SOLVER_NOT_CONVERGED',
            'b_material_subset_design': False,
        },
    }
    status = str(solved.get('status', 'UNKNOWN'))
    if status == 'CONVERGED' and not solved.get('paper_optimality_conditions_verified', False):
        status = 'OPTIMALITY_NOT_VERIFIED'
    if status == 'CONVERGED':
        vertices = np.asarray(solved.get('variable_material_vertices'), dtype=float)
        faces = np.asarray(solved.get('variable_material_faces'), dtype=np.int64)
        require(vertices.ndim == 2 and vertices.shape[1] == 3 and len(vertices) >= 4 and
                faces.ndim == 2 and faces.shape[1] == 3 and len(faces) >= 4,
                'converged B solve did not return a closed material boundary')
        b_mesh = MeshSnapshot((vertices, faces))
        b_moments = _raw_moments(b_mesh, composition.request.material_density)
        target_axes = np.asarray(composition.request.target_axes, dtype=float)
        model_to_target = target_axes.T
        target_translation = -model_to_target @ np.asarray(
            composition.request.target_origin, dtype=float)
        fixed_target_moments = _transform_moments(
            composition.request.fixed_shell_moments,
            model_to_target, target_translation)
        independent_total = _add_moments(fixed_target_moments, b_moments)
        reported_total = tuple(float(value) for value in solved['total_raw_moments'])
        side = max(float(np.ptp(np.asarray(composition.request.outer_vertices), axis=0).max()),
                   1e-300)
        scales = (side**3, side**4, side**4, side**4,
                  side**5, side**5, side**5, side**5, side**5, side**5)
        moment_residual = max(abs(independent_total[i] - reported_total[i]) / scales[i]
                              for i in range(10))
        require(moment_residual <= 1e-10,
                'independent chained material moments differ from solver total')
        outer_vertices = np.asarray(composition.request.outer_vertices, dtype=float)
        outer_target = (outer_vertices - np.asarray(
            composition.request.target_origin, dtype=float)) @ model_to_target.T
        target_length = float(np.ptp(outer_target, axis=0).max())
        constraint_scale = composition.request.material_density * target_length**4
        independent_constraint_norm = math.hypot(independent_total[1], independent_total[2]) / constraint_scale
        require(independent_constraint_norm <= constraint_tolerance,
                'independent chained material fails static horizontal balance constraints')

        design_b = MeshSnapshot((composition.request.design_vertices,
                                 composition.request.design_faces)).exact_region()
        model_vertices = vertices @ target_axes.T + np.asarray(
            composition.request.target_origin, dtype=float)
        model_b_mesh = MeshSnapshot((model_vertices, faces))
        retained_b = model_b_mesh.exact_region()
        containment = native().relate_exact_material_regions(retained_b, design_b)
        require(containment['completed'] and containment['a_subset_b'],
                'solved B material is not an exactly certified subset of its design domain')
        identity['independent_total_moment_residual'] = moment_residual
        identity['independent_normalized_constraints_norm'] = independent_constraint_norm
        identity['independent_total_raw_moments'] = list(independent_total)
        identity['final_material_domain_audit'] = {
            'status': 'PASS',
            'b_material_subset_design': True,
            'relation': containment,
            'base_partition': composition.binding['exact_geometry_check'],
            'target_to_model_transform': {
                'rotation': model_to_target.T.tolist(),
                'origin': list(composition.request.target_origin),
            },
            'cavity_note': 'carved design-domain volume is intentionally not part of the material partition',
        }
        identity['final_b_material_identity'] = model_b_mesh.identity
    identity['status'] = status
    solver_json = _json_value(solved)
    identity['solver_result_digest'] = digest(solver_json)
    result = object.__new__(FloatingIslandChainResult)
    object.__setattr__(result, 'status', status)
    object.__setattr__(result, 'composition', composition)
    object.__setattr__(result, '_solver_result', canonical(solver_json))
    object.__setattr__(result, '_identity', canonical(identity))
    return result


def solve_floating_island_b_buoyant_from_sealed_a(
        request_b, source_artifact, *, region_prefix,
        source_to_request_rotation, source_to_request_translation,
        source_length_to_request_unit, reference_density,
        source_length_unit, target_length_unit,
        source_mass_unit, target_mass_unit, displacement_vertices,
        displacement_faces, waterline_origin, waterline_normal, liquid_density,
        existing_fixed_regions=(), initial_parameters=(0., 0., 0.),
        maximum_iterations=60, maximum_line_search_steps=24,
        constraint_tolerance=1e-10, capture_seal_failure=False):
    """Optimize and seal buoyant B with sealed A as actual fixed material.

    Unlike :func:`solve_floating_island_b_from_sealed_a`, this path recomputes
    displacement and waterline geometry from a closed mesh, derives B's
    heterogeneous fixed-material mass from A plus any other fixed regions,
    runs the prescribed-waterline plane solver, and seals the resulting B
    material. Free-waterline or perturbed-pose stability remains a separate
    validation and is not inferred here.

    ``capture_seal_failure`` preserves the completed native ledger when the
    post-solve material check rejects it. The default keeps legacy exceptions.
    """
    require(isinstance(request_b, SpinItFasterRequest),
            'expected SpinItFasterRequest for target B')
    require(request_b.mode == 'STATIC_BALANCE',
            'buoyant B adapter expects a geometry-only static base request')
    liquid_density = scalar(liquid_density, 'liquid_density', positive=True)
    axes = np.asarray(request_b.target_axes, dtype=np.float64)
    origin = np.asarray(request_b.target_origin, dtype=np.float64)
    require(np.allclose(axes.T @ axes, np.eye(3), rtol=0., atol=1e-10) and
            abs(float(np.linalg.det(axes)) - 1.) <= 1e-10,
            'B target axes must be a right-handed orthonormal frame')
    world_normal = np.asarray(_finite_tuple(waterline_normal, 3, 'waterline_normal'))
    normal_length = float(np.linalg.norm(world_normal))
    require(math.isfinite(normal_length) and normal_length > 0.,
            'waterline_normal must be finite and nonzero')
    world_normal = world_normal / normal_length
    require(np.allclose(axes.T @ world_normal, (0., 0., 1.), rtol=0., atol=1e-10),
            'waterline normal must align with B target-frame vertical')
    world_waterline_origin = np.asarray(_finite_tuple(
        waterline_origin, 3, 'waterline_origin'), dtype=np.float64)

    composition = compose_sealed_faster_material_as_fixed_region(
        request_b, source_artifact, region_prefix=region_prefix,
        source_to_request_rotation=source_to_request_rotation,
        source_to_request_translation=source_to_request_translation,
        source_length_to_request_unit=source_length_to_request_unit,
        reference_density=reference_density,
        source_length_unit=source_length_unit, target_length_unit=target_length_unit,
        source_mass_unit=source_mass_unit, target_mass_unit=target_mass_unit,
        existing_fixed_regions=existing_fixed_regions)

    displacement = MeshSnapshot((displacement_vertices, displacement_faces))
    displacement_model, displacement_triangles = displacement.arrays()
    displacement_target = (displacement_model - origin) @ axes
    target_waterline_origin = (world_waterline_origin - origin) @ axes
    hydrostatic = compute_hydrostatic_geometry(
        displacement_target, displacement_triangles.tolist(),
        waterline_origin=target_waterline_origin,
        waterline_normal=(0., 0., 1.))

    rotation = np.asarray(source_to_request_rotation, dtype=np.float64)
    translation = np.asarray(source_to_request_translation, dtype=np.float64)
    scale = scalar(source_length_to_request_unit, 'source_length_to_request_unit', positive=True)
    reference = scalar(reference_density, 'reference_density', positive=True)
    source_problem = source_artifact.problem
    source_density = {region.region_id: region.density
                      for region in source_problem.design_regions}
    fixed_problem_regions = []
    fixed_volumes = []
    for region in source_problem.fixed_regions:
        if region.mesh.explicit_empty:
            continue
        vertices, faces = region.mesh.arrays()
        mapped = scale * (vertices @ rotation.T) + translation
        ratio = float(region.density / reference)
        fixed_problem_regions.append((region_prefix + ':' + region.region_id,
                                      MeshSnapshot((mapped, faces)), ratio))
        fixed_volumes.append(float(_raw_moments(MeshSnapshot((mapped, faces)), 1.)[0]))
    for region_id, retained, _ in source_artifact.design_results:
        if retained.explicit_empty:
            continue
        vertices, faces = retained.arrays()
        mapped = scale * (vertices @ rotation.T) + translation
        ratio = float(source_density[region_id] / reference)
        fixed_problem_regions.append((region_prefix + ':' + region_id,
                                      MeshSnapshot((mapped, faces)), ratio))
        fixed_volumes.append(float(_raw_moments(MeshSnapshot((mapped, faces)), 1.)[0]))
    for index, (region_id, vertices, faces, density_ratio) in enumerate(tuple(existing_fixed_regions)):
        mesh = MeshSnapshot((vertices, faces))
        fixed_problem_regions.append((str(region_id), mesh, float(density_ratio)))
        fixed_volumes.append(float(_raw_moments(mesh, 1.)[0]))

    fixed_volume = math.fsum(fixed_volumes)
    fixed_mass = float(composition.request.fixed_shell_moments[0])
    require(math.isfinite(fixed_volume) and fixed_volume > 0. and
            math.isfinite(fixed_mass) and fixed_mass > 0.,
            'A plus fixed B material must have positive independently measured volume and mass')
    effective_density = fixed_mass / fixed_volume
    buoyant_request = replace(
        composition.request, mode='BUOYANT_STABILITY',
        liquid_density=liquid_density,
        displaced_liquid_volume=float(hydrostatic.displaced_volume),
        center_of_buoyancy_target=tuple(hydrostatic.center_of_buoyancy),
        waterline_area=float(hydrostatic.waterline_area),
        waterline_central_second_moment=tuple(tuple(float(x) for x in row)
                                              for row in hydrostatic.waterline_central_second_moment),
        fixed_shell_volume=fixed_volume,
        fixed_shell_density=effective_density)
    solved = solve_spin_it_faster_buoyant_plane(
        buoyant_request, initial_parameters=initial_parameters,
        maximum_iterations=maximum_iterations,
        maximum_line_search_steps=maximum_line_search_steps,
        constraint_tolerance=constraint_tolerance)

    source_to_problem_regions = [(region_id, mesh, density)
                                 for region_id, mesh, density in fixed_problem_regions]
    outer_mesh = MeshSnapshot((request_b.outer_vertices, request_b.outer_faces))
    design_mesh = MeshSnapshot((request_b.design_vertices, request_b.design_faces))
    problem_b = SpinArtifactProblem(
        original=outer_mesh,
        fixed_regions=source_to_problem_regions,
        design_regions=[('design:B', design_mesh, request_b.material_density)],
        length_unit=target_length_unit, mass_unit=target_mass_unit,
        target={'mode': 'BUOYANT_STABILITY', 'axis': world_normal.tolist(),
                'origin': world_waterline_origin.tolist(), 'support_points': []},
        minimum_outer_wall=request_b.minimum_wall_thickness)

    identity = {
        'schema': 'faster-floating-island-b-buoyant-chain-v1',
        'operation': 'FLOATING_ISLAND_BUOYANT_B_FROM_SEALED_A',
        'scope': 'SYNTHETIC_METHOD_VALIDATION',
        'source_a_scientific_digest': source_artifact.scientific_digest,
        'source_a_evidence_digest': source_artifact.evidence_digest,
        'composition_binding': composition.binding,
        'request_digest': buoyant_request.request_digest,
        'hydrostatic_measurement': {
            'displacement_mesh': displacement.identity,
            'displaced_volume': hydrostatic.displaced_volume,
            'center_of_buoyancy_target': list(hydrostatic.center_of_buoyancy),
            'waterline_area': hydrostatic.waterline_area,
            'waterline_central_second_moment': [list(row) for row in
                                                 hydrostatic.waterline_central_second_moment],
            'world_waterline_origin': world_waterline_origin.tolist(),
            'world_waterline_normal': world_normal.tolist(),
        },
        'fixed_material_measurement': {
            'volume': fixed_volume, 'mass': fixed_mass,
            'effective_density': effective_density,
            'source_region_count': len(fixed_problem_regions),
        },
        'solver_status': solved.get('status'),
        'solver_diagnostic': solved.get('diagnostic'),
        'solver_iterations': solved.get('iterations'),
        'solver_accepted_steps': solved.get('accepted_steps'),
        'plane_coefficients': _json_value(solved.get('plane_coefficients')),
        'normalized_constraints': _json_value(solved.get('normalized_constraints')),
        'solver_trace': _json_value(solved.get('trace')),
        'solver_result': _json_value(solved),
        'solver_result_digest': digest(_json_value(solved)),
        'prescribed_waterline_stable': solved.get('stable_at_prescribed_waterline'),
        'stability_margin': solved.get('stability_margin'),
        'free_waterline_validation': 'NOT_RUN',
    }
    if solved.get('status') != 'CONVERGED':
        return _result(str(solved.get('status', 'UNKNOWN')), identity)

    try:
        sealed = seal_buoyant_plane_result_as_artifact(
            problem_b, buoyant_request, solved,
            displacement_vertices=displacement_model,
            displacement_faces=displacement_triangles,
            waterline_origin=world_waterline_origin,
            waterline_normal=world_normal)
    except (ValueError, RuntimeError) as error:
        if not capture_seal_failure:
            raise
        identity['seal_failure'] = {'type': type(error).__name__, 'message': str(error)}
        return _result('MATERIAL_VALIDATION_FAILED', identity, seal_error=str(error))
    require(sealed.status == 'CONVERGED' and sealed.artifact is not None,
            'B buoyancy solve did not produce an independently sealed material')
    identity['b_seal'] = sealed.identity
    identity['b_artifact_scientific_digest'] = sealed.artifact.scientific_digest
    identity['b_artifact_evidence_digest'] = sealed.artifact.evidence_digest
    return _result('CONVERGED', identity, sealed.artifact)


def recompute_sealed_buoyant_stability_at_pose(
        artifact, *, displacement_vertices, displacement_faces, rotation,
        translation, reference_origin, waterline_normal, liquid_density,
        relative_volume_tolerance=1e-10, maximum_iterations=80):
    """Recompute free-waterline equilibrium and small-angle stability at a pose.

    The displacement hull must be the exact outer boundary sealed in the
    artifact. A rigid pose is applied to both the hull and its independently
    integrated material center of mass. The returned metacentric matrix is
    the linearized hydrostatic stability result at the volume-matched
    waterline; it is not a fluid-dynamics simulation or a substitute for
    author-specific flooding/displacement semantics.
    """
    require(isinstance(artifact, SealedSpinArtifact),
            'expected an independently sealed Spin artifact')
    outer = MeshSnapshot((displacement_vertices, displacement_faces))
    require(outer.identity == artifact.outer.identity,
            'displacement hull must match the sealed artifact outer boundary')
    hull_vertices, hull_faces = outer.arrays()

    rotation = np.asarray(rotation, dtype=np.float64)
    require(rotation.shape == (3, 3) and np.isfinite(rotation).all(),
            'pose rotation must be a finite 3x3 matrix')
    require(np.allclose(rotation.T @ rotation, np.eye(3), rtol=0., atol=1e-10),
            'pose rotation must be orthonormal')
    require(abs(float(np.linalg.det(rotation)) - 1.) <= 1e-10,
            'pose rotation must be proper, not reflected')
    translation = np.asarray(_finite_tuple(translation, 3, 'pose translation'))
    reference_origin = np.asarray(_finite_tuple(reference_origin, 3, 'reference origin'))
    normal = np.asarray(_finite_tuple(waterline_normal, 3, 'waterline normal'))
    normal_norm = float(np.linalg.norm(normal))
    require(normal_norm > 0., 'waterline normal must be nonzero')
    normal /= normal_norm
    liquid_density = scalar(liquid_density, 'liquid density', positive=True)
    require(type(maximum_iterations) is int and maximum_iterations > 0,
            'maximum iterations must be a positive integer')
    relative_volume_tolerance = scalar(
        relative_volume_tolerance, 'relative volume tolerance', positive=True)

    facts = artifact.facts
    physical = facts['physical']
    mass = scalar(physical['mass'], 'sealed material mass', positive=True)
    center = np.asarray(_finite_tuple(physical['center_of_mass'], 3,
                                      'sealed material center of mass'))
    posed_vertices = hull_vertices @ rotation.T + translation
    posed_center = center @ rotation.T + translation
    equilibrium = find_hydrostatic_equilibrium(
        posed_vertices.tolist(), hull_faces.tolist(),
        reference_origin=reference_origin.tolist(),
        waterline_normal=normal.tolist(),
        target_displaced_volume=mass / liquid_density,
        relative_volume_tolerance=relative_volume_tolerance,
        maximum_iterations=maximum_iterations)

    volume = scalar(equilibrium.geometry.displaced_volume,
                    'recomputed displaced volume', positive=True)
    buoyancy_center = np.asarray(equilibrium.geometry.center_of_buoyancy,
                                 dtype=np.float64)
    waterplane = np.asarray(
        equilibrium.geometry.waterline_central_second_moment,
        dtype=np.float64)
    require(waterplane.shape == (2, 2) and np.isfinite(waterplane).all() and
            np.allclose(waterplane, waterplane.T, rtol=0., atol=1e-12),
            'recomputed waterplane moment must be finite and symmetric')
    eigenvalues = np.linalg.eigvalsh(waterplane)
    require(np.isfinite(eigenvalues).all() and eigenvalues[0] > 0.,
            'recomputed waterplane moment must be positive definite')

    vertical_cg_above_cb = float(np.dot(posed_center - buoyancy_center, normal))
    metacentric = waterplane / volume - vertical_cg_above_cb * np.eye(2)
    margins = np.linalg.eigvalsh(metacentric)
    require(np.isfinite(margins).all(), 'metacentric stability margins are non-finite')
    residual = abs(float(equilibrium.volume_residual))
    return {
        'schema': 'sealed-buoyant-pose-stability-v1',
        'scope': 'SEALED_OUTER_HULL_SMALL_ANGLE_HYDROSTATIC_LINEARIZATION',
        'artifact_scientific_digest': artifact.scientific_digest,
        'artifact_evidence_digest': artifact.evidence_digest,
        'displacement_hull_digest': outer.identity,
        'pose': {'rotation': rotation.tolist(), 'translation': translation.tolist()},
        'mass': mass,
        'liquid_density': liquid_density,
        'target_displaced_volume': equilibrium.target_displaced_volume,
        'recomputed_displaced_volume': volume,
        'volume_residual': residual,
        'waterline_signed_offset': equilibrium.waterline_signed_offset,
        'waterline_area': equilibrium.geometry.waterline_area,
        'center_of_gravity': posed_center.tolist(),
        'center_of_buoyancy': buoyancy_center.tolist(),
        'waterplane_central_second_moment': waterplane.tolist(),
        'metacentric_stability_matrix': metacentric.tolist(),
        'metacentric_eigenvalue_margins': margins.tolist(),
        'minimum_metacentric_margin': float(margins[0]),
        'stable_at_pose': bool(margins[0] > 0.),
        'iterations': equilibrium.iterations,
    }


def material_run_scientific_identity(identity):
    """Scientific projection of a saved supported-material run, not acceptance.

    Legacy identities stay unchanged. The versioned directional path retains
    runtime timing and original Python stacks in its full identity, but these
    named diagnostics do not enter its scientific digest. Native evaluation
    errors, types, events and all numerical arrays remain in the projection.
    """
    if identity.get('scientific_identity_profile') != 'FASTER_SUPPORTED_MATERIAL_V1':
        return identity
    projected = parse(canonical(identity))
    native = projected.get('native_run', {})
    native.pop('evaluation_timings', None)
    native.get('supported_space', {}).pop('python_callback_diagnostics', None)
    return projected


def _result(status, identity, artifact=None, seal_error=None):
    identity = {**identity, 'status': status, 'seal_error': seal_error,
                'artifact_scientific_digest': None if artifact is None else artifact.scientific_digest}
    if identity.get('directional_mapping', {}).get('schema') == 'faster-material-supported-directions-v1':
        identity['scientific_identity_profile'] = 'FASTER_SUPPORTED_MATERIAL_V1'
    result = object.__new__(FasterMaterialResult)
    for name, value in dict(status=status, artifact=artifact, seal_error=seal_error,
                            scientific_digest=digest(material_run_scientific_identity(identity)),
                            _identity=canonical(identity)).items():
        object.__setattr__(result, name, value)
    return result


def quadric_coefficients_to_cad_frame(
        coefficients, *, target_axes, target_origin,
        original_units_per_normalized_unit):
    """Pull back Faster target-frame coefficients into scaled CAD coordinates.

    The request frame maps model coordinates ``x`` to target coordinates
    ``y = target_axes.T @ (x - target_origin)``. The OCCT mesh import applies
    ``x_cad = scale * x``; this function therefore composes that exact affine
    map with the homogeneous quadric matrix, then returns the ten coefficients
    in the solver's canonical order with unit Euclidean norm.
    """
    values = np.asarray(_finite_tuple(coefficients, 10, 'quadric coefficients'),
                        dtype=np.float64)
    norm = float(np.linalg.norm(values))
    require(abs(norm - 1.) <= 1e-10,
            'solver quadric coefficients must have unit Euclidean norm')
    axes = np.asarray(target_axes, dtype=np.float64)
    require(axes.shape == (3, 3) and np.isfinite(axes).all(),
            'target_axes must be a finite 3x3 frame')
    require(np.allclose(axes.T @ axes, np.eye(3), rtol=0., atol=1e-10) and
            abs(float(np.linalg.det(axes)) - 1.) <= 1e-10,
            'target_axes must be right-handed and orthonormal')
    origin = np.asarray(_finite_tuple(target_origin, 3, 'target_origin'),
                        dtype=np.float64)
    scale = scalar(original_units_per_normalized_unit,
                   'original_units_per_normalized_unit', positive=True)

    quadratic = np.asarray(((values[4], values[9], values[8]),
                            (values[9], values[5], values[7]),
                            (values[8], values[7], values[6])), dtype=np.float64)
    linear = values[1:4]
    homogeneous = np.zeros((4, 4), dtype=np.float64)
    homogeneous[:3, :3] = quadratic
    homogeneous[:3, 3] = 0.5 * linear
    homogeneous[3, :3] = 0.5 * linear
    homogeneous[3, 3] = values[0]

    target_from_cad = np.eye(4, dtype=np.float64)
    target_from_cad[:3, :3] = axes.T / scale
    target_from_cad[:3, 3] = -axes.T @ origin
    cad_homogeneous = target_from_cad.T @ homogeneous @ target_from_cad
    cad_quadratic = cad_homogeneous[:3, :3]
    cad_linear = 2. * cad_homogeneous[:3, 3]
    transformed = np.asarray((
        cad_homogeneous[3, 3], cad_linear[0], cad_linear[1], cad_linear[2],
        cad_quadratic[0, 0], cad_quadratic[1, 1], cad_quadratic[2, 2],
        cad_quadratic[1, 2], cad_quadratic[0, 2], cad_quadratic[0, 1],
    ), dtype=np.float64)
    transformed_norm = float(np.linalg.norm(transformed))
    require(np.isfinite(transformed).all() and math.isfinite(transformed_norm) and
            transformed_norm > 0., 'transformed CAD quadric is non-finite or zero')
    return tuple(float(value) for value in transformed / transformed_norm)


def _transform_moments(moments, axes, translation):
    """Push x to A*x+t; mass unchanged, full symmetric second moments retained."""
    m, lx, ly, lz, xx, yy, zz, yz, xz, xy = moments
    a, t = np.asarray(axes), np.asarray(translation)
    first = a @ np.array([lx, ly, lz])
    second = a @ np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]) @ a.T
    second += np.outer(first, t) + np.outer(t, first) + m*np.outer(t, t)
    first += m*t
    return [float(m), *first.tolist(), *np.diag(second).tolist(),
            float(second[1,2]), float(second[0,2]), float(second[0,1])]


def _target_axes(axis):
    # Scaling first avoids overflow for finite very large direction components.
    z = np.asarray(axis, dtype=float)
    z = z / np.max(np.abs(z))
    z = z / np.linalg.norm(z)
    seed = np.eye(3)[int(np.argmin(np.abs(z)))]
    x = seed - z*np.dot(seed, z)
    x /= np.linalg.norm(x)
    return np.column_stack((x, np.cross(z, x), z))


def _cad_quadric_type(evaluation, coefficients=None):
    """Map native broad type/subtype to the exact typed CAD protocol token."""
    interface_type = evaluation.get('interface_type')
    if interface_type == 'HYPERBOLOID':
        return evaluation.get('interface_subtype')
    if interface_type == 'PARABOLOID':
        require(coefficients is not None,
                'paraboloid CAD classification requires its locked coefficients')
        q = np.asarray(_finite_tuple(coefficients, 10, 'paraboloid coefficients'),
                       dtype=np.float64)
        quadratic = np.asarray(((q[4], q[9], q[8]),
                                (q[9], q[5], q[7]),
                                (q[8], q[7], q[6])), dtype=np.float64)
        eigenvalues = np.linalg.eigvalsh(quadratic)
        threshold = 1e-12 * max(1.0, float(np.linalg.norm(quadratic)))
        transverse = eigenvalues[np.abs(eigenvalues) > threshold]
        require(transverse.shape == (2,),
                'paraboloid CAD classification requires exactly two transverse eigenvalues')
        if transverse[0] * transverse[1] < 0.0:
            return 'HYPERBOLIC_PARABOLOID'
        return 'PARABOLOID'
    return interface_type


def _prepare(problem, evaluator_settings, solver_settings, *, allow_reduced_stratum=False,
             reuse_input_validation=False):
    n = _native_module()
    if not all(hasattr(n, name) for name in ('SpinEvaluatorSettings', 'evaluate_spin_it_faster_quadric',
                                            'solve_spin_it_faster', 'integrate_raw_mass_moments')):
        raise RuntimeError('strict Faster material bridge requires native Faster evaluator support')
    require(isinstance(problem, SpinArtifactProblem), 'expected locked SpinArtifactProblem')
    problem = _recheck_problem(problem, reuse_validation=reuse_input_validation)
    require(problem.validation_scope == SCOPE, 'NOT_SUPPORTED: material scope')
    require(len(problem.design_regions) == 1, 'NOT_SUPPORTED: exactly one design region required')
    supported_methods = {'DISCRETE_FULL_RANK_QUADRIC', 'DISCRETE_ELLIPSOID'}
    if allow_reduced_stratum:
        supported_methods.add('DISCRETE_REDUCED_STRATUM')
    require(isinstance(evaluator_settings, SpinEvaluatorConfig) and
            evaluator_settings.derivative_method in supported_methods,
            'strict bridge requires explicit DISCRETE_FULL_RANK_QUADRIC settings')
    require(isinstance(solver_settings, SpinSolverConfig) and
            solver_settings.equation_profile == 'DUAL_MONOMIAL_RAW_V1',
            'strict bridge requires explicit DUAL_MONOMIAL_RAW_V1 settings')
    target = problem.target
    require(target['mode'] != 'BUOYANT_STABILITY', 'NOT_SUPPORTED: buoyancy requires actual hydrostatics')
    if target['mode'] == 'STATIC_BALANCE':
        require(target['support_points'] == [target['origin']],
                'NOT_SUPPORTED: static target needs exactly one support point equal to target origin')
    else:
        require(not target['support_points'], 'NOT_SUPPORTED: rotational target with support points')
    norm = problem.normalization
    origin, length, density = np.array(norm['origin']), norm['L'], norm['R']
    def normalized_mesh(mesh):
        v, f = mesh.arrays()
        return ((v-origin)/length).tolist(), f.tolist()
    fixed = []
    for region in problem.fixed_regions:
        v, f = normalized_mesh(region.mesh)
        fixed.append(integrate_raw_mass_moments(v, f, region.density/density))
    moments = tuple(math.fsum(row[i] for row in fixed) for i in range(10))
    outer_v, outer_f = normalized_mesh(problem.original)
    design_v, design_f = normalized_mesh(problem.design_regions[0].mesh)
    axes = _target_axes(target['axis'])
    target_origin = (np.array(target['origin']) - origin)/length
    request = SpinItFasterRequest(outer_v, outer_f, design_v, design_f,
        fixed_shell_moments=moments, target_origin=target_origin.tolist(), target_axes=axes.tolist(),
        mode=target['mode'], material_density=problem.design_regions[0].density/density,
        original_units_per_normalized_unit=length, minimum_wall_thickness=problem.minimum_outer_wall/length,
        settings=solver_settings)
    identity = {'schema': 'faster-locked-material-run-v1', 'problem_digest': problem.input_digest,
                'scope': SCOPE, 'project_status': 'PARTIAL_BLOCKED',
                'solver_target_feasibility': 'UNKNOWN', 'physical_validation': 'NOT_IN_SCOPE',
                'evaluator_settings': evaluator_settings.identity_dict(), 'solver_settings': asdict(solver_settings),
                'native_request_digest': digest(asdict(request)),
                'coefficient_frame': 'INPUT_NORMALIZED_TARGET_V1',
                'coefficient_convention': 'NATIVE_MATRIX_ENTRY_P',
                'normalization': norm, 'target_axes': axes.tolist(), 'target_origin': target_origin.tolist(),
                'fixed_moments_input_normalized': list(moments)}
    # The opt-in map must not alter historical default configuration identities.
    if not evaluator_settings.fixed_normal_paraboloid:
        identity['evaluator_settings'].pop('fixed_normal_paraboloid')
    return problem, request, identity


@dataclass(frozen=True, init=False)
class _PreparedFasterProblem:
    """Run-local capability for one validated immutable problem and settings."""
    problem: SpinArtifactProblem
    request: SpinItFasterRequest
    identity_json: str
    evaluator_settings: SpinEvaluatorConfig
    solver_settings: SpinSolverConfig
    native_module: object
    input_validation_profile: str

    def __init__(self):
        raise TypeError('prepared problems must be created by the validating entry point')


def _prepare_faster_problem(problem, evaluator_settings, solver_settings, *,
                            input_validation_profile='RECHECK_INPUT_V1'):
    require(input_validation_profile in ('RECHECK_INPUT_V1', 'CONSTRUCTOR_PROOF_REUSE_V1'),
            'unknown input validation profile')
    problem, request, identity = _prepare(problem, evaluator_settings, solver_settings,
        reuse_input_validation=input_validation_profile=='CONSTRUCTOR_PROOF_REUSE_V1')
    if input_validation_profile != 'RECHECK_INPUT_V1':
        identity['input_validation_profile'] = input_validation_profile
    prepared = object.__new__(_PreparedFasterProblem)
    for name, value in dict(problem=problem, request=request,
            identity_json=json.dumps(identity, sort_keys=True, allow_nan=False),
            evaluator_settings=evaluator_settings, solver_settings=solver_settings,
            native_module=_native_module(), input_validation_profile=input_validation_profile).items():
        object.__setattr__(prepared, name, value)
    return prepared


def _frozen_support_basis(value):
    if value is None:
        return None
    basis = np.array(value, dtype=np.float64, copy=True)
    require(basis.ndim == 2 and basis.shape[0] == 10 and
            2 <= basis.shape[1] <= 10 and np.isfinite(basis).all(),
            'support basis must be finite 10-by-k, 2 <= k <= 10')
    require(np.linalg.norm(basis.T @ basis - np.eye(basis.shape[1])) <= 1e-12,
            'support basis must be orthonormal')
    basis.flags.writeable = False
    return basis


@dataclass(frozen=True)
class _PreparedExternalFixedMaterial:
    problem: SpinArtifactProblem
    request: SpinItFasterRequest
    native_module: object
    moments: tuple
    identity: str


def _require_same_call_producer(material_validator):
    from .cherchi_inprocess import InProcessCherchiEvaluator
    producer = getattr(material_validator, '__self__', None)
    require(type(producer) is InProcessCherchiEvaluator and
            getattr(material_validator, '__func__', None) is
            InProcessCherchiEvaluator.validate_same_call_material,
            'layered validation requires a bound production material validator')


def _prepare_external_fixed_material(problem, request, *, stage_observer=None):
    """Independently integrate immutable fixed inputs once for one run."""
    axes = np.asarray(request.target_axes, dtype=np.float64)
    origin = np.asarray(request.target_origin, dtype=np.float64)
    normalization = problem.normalization
    fixed_model = []
    for index, region in enumerate(problem.fixed_regions):
        if region.mesh.explicit_empty:
            continue
        vertices, faces = region.mesh.arrays()
        target_vertices = (((vertices - np.asarray(normalization['origin'])) /
                            normalization['L']) - origin) @ axes
        def integrate(target_vertices=target_vertices, faces=faces, region=region):
            return _raw_moments_from_arrays(target_vertices, faces,
                                            region.density / normalization['R'])
        fixed_model.append(integrate() if stage_observer is None else stage_observer(
            f'map_fixed_material_integration_{index}', integrate))
    moments = tuple(_add_moments(*fixed_model) if fixed_model else (0.,) * 10)
    side = max(float(np.max(np.ptp(np.asarray(request.design_vertices), axis=0))), 1e-300)
    scales = np.asarray((side**3, side**4, side**4, side**4,
                         side**5, side**5, side**5, side**5, side**5, side**5))
    expected = np.asarray(_transform_moments(request.fixed_shell_moments,
                                             axes.T, -axes.T @ origin))
    error = float(np.max(np.abs(expected - moments) / scales))
    require(math.isfinite(error) and error <= 1e-10,
            'prepared fixed moments differ from independently integrated fixed material')
    module = _native_module()
    identity = digest(dict(problem=problem.input_digest, request=asdict(request),
        native_sha256=hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest(),
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        moments=moments))
    return _PreparedExternalFixedMaterial(problem, request, module, moments, identity)


def _external_geometry_moment_map(
        problem, request, evaluation, *, require_derivative=True,
        expected_coefficients=None, stage_observer=None, support_basis=None,
        solver_boundary_profile='EXACT_NEF_REGION_V1', prepared_fixed_material=None,
        validation_profile='INDEPENDENT_EACH_EVALUATION_V1', material_validator=None):
    """Validate one backend mesh and bind its variable moments to the request.

    External Boolean backends return the *variable* material and its Jacobian.
    Fixed material is independently integrated here and added exactly once in
    the solver's target frame.  ``variable_raw_moments`` is required as a
    cross-check, never trusted as the source of the solver residual.
    """
    require(stage_observer is None or callable(stage_observer),
            'stage_observer must be callable or None')
    require(solver_boundary_profile in ('EXACT_NEF_REGION_V1','EXACT_RATIONAL_BOUNDARY_V1','EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1'),
            'unknown solver boundary profile')
    require(validation_profile in ('INDEPENDENT_EACH_EVALUATION_V1', 'LAYERED_SAME_CALL_V1'),
            'unknown material validation profile')
    fixed_boundary_proof = None
    if (validation_profile == 'LAYERED_SAME_CALL_V1' or
            solver_boundary_profile in ('EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1')):
        require(callable(material_validator) and expected_coefficients is not None,
                'layered validation requires a same-call material validator')
        _require_same_call_producer(material_validator)
        if stage_observer is None:
            material_validator(evaluation, expected_coefficients, request)
        else:
            stage_observer('map_same_call_material_validation',
                           lambda: material_validator(evaluation, expected_coefficients, request))
        if solver_boundary_profile in ('EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1'):
            fixed_boundary_proof = material_validator.__self__.fixed_boundary_proof
            require(fixed_boundary_proof is not None,
                    'fixed subset profile requires a prepared native fixed boundary proof')
    if prepared_fixed_material is not None:
        require(isinstance(prepared_fixed_material, _PreparedExternalFixedMaterial) and
                prepared_fixed_material.problem is problem and
                prepared_fixed_material.request is request and
                prepared_fixed_material.native_module is _native_module(),
                'fixed material preparation does not bind this immutable problem/request/build')

    def run_stage(name, operation):
        return operation() if stage_observer is None else stage_observer(name, operation)

    require(isinstance(evaluation, dict), 'external geometry evaluator must return a mapping')
    support = _frozen_support_basis(support_basis)
    required = ('variable_vertices', 'variable_faces', 'variable_raw_moments',
                'derivative_defined', 'geometry_hash',
                'topology_hash', 'evaluation_map_identity',
                'boolean_backend_identity', 'derivative_identity')
    if support is None:
        required += ('moment_jacobian',)
    else:
        require('moment_jacobian' not in evaluation,
                'supported material map must not provide an ambient Jacobian')
        if require_derivative:
            require(expected_coefficients is not None,
                    'directional material map requires requested coefficients')
            required += ('native_tangent_basis', 'native_directional_jacobian')
    for key in required:
        require(key in evaluation, f'external geometry evaluator omitted {key}')
    if expected_coefficients is not None:
        require('quadric_coefficients' in evaluation,
                'external geometry evaluator omitted quadric_coefficients')
        expected = np.asarray(_finite_tuple(
            expected_coefficients, 10, 'expected quadric coefficients'), dtype=np.float64)
        actual = np.asarray(_finite_tuple(
            evaluation['quadric_coefficients'], 10,
            'evaluated quadric coefficients'), dtype=np.float64)
        require(actual.tobytes(order='C') == expected.tobytes(order='C'),
                'external geometry quadric coefficients do not match the requested coefficients')
    raw_vertices = evaluation['variable_vertices']
    raw_faces = evaluation['variable_faces']
    vertices = np.asarray(raw_vertices, dtype=np.float64)
    faces = np.asarray(raw_faces)
    require(vertices.ndim == 2 and vertices.shape[1:] == (3,) and
            faces.ndim == 2 and faces.shape[1:] == (3,),
            'external variable material must be a triangle mesh')
    rational_vertices = evaluation.get('variable_rational_vertices')
    exact_snapshot = evaluation.get('variable_exact_snapshot')
    exact_geometry_identity = evaluation.get('exact_geometry_identity')
    if len(vertices):
        def build_variable_snapshot():
            if exact_snapshot is not None:
                require(isinstance(exact_snapshot, MeshSnapshot),
                        'external exact snapshot must be a MeshSnapshot')
                snapshot = exact_snapshot
                snapshot_vertices, snapshot_faces = snapshot.arrays()
                display_vertices = np.asarray(vertices, dtype='<f8', order='C')
                display_faces = np.asarray(faces, dtype='<i4', order='C')
                require(snapshot_vertices.shape == display_vertices.shape and
                        snapshot_faces.shape == display_faces.shape and
                        snapshot_faces.tobytes(order='C') == display_faces.tobytes(order='C'),
                        'external exact snapshot does not bind callback topology')
                require(snapshot.rational_vertices is not None,
                        'external exact snapshot omitted rational coordinates')
                float_audit = evaluation.get('exact_float_derivation_audit')
                if float_audit is None:
                    require(snapshot_vertices.tobytes(order='C') ==
                            display_vertices.tobytes(order='C'),
                            'external exact snapshot does not bind callback mesh bytes')
                else:
                    require(isinstance(float_audit, dict) and
                            float_audit.get('schema') == 'exact-boolean-float-derivation-audit-v1' and
                            float_audit.get('rational_geometry_identity') == snapshot.identity,
                            'external float derivation audit is not bound to the rational snapshot')
                    scale = float(max(
                        1.0,
                        float(np.max(np.abs(display_vertices))) if display_vertices.size else 0.0,
                        float(np.linalg.norm(np.ptp(display_vertices, axis=0)))
                        if len(display_vertices) else 0.0,
                    ))
                    bound = float(64.0 * np.finfo(np.float64).eps * scale)
                    maximum_error = float(np.max(np.abs(snapshot_vertices - display_vertices)))                         if display_vertices.size else 0.0
                    expected_status = ('BOUNDED_FLOAT_DERIVATION' if maximum_error
                                       else 'EXACT_FLOAT_MATCH')
                    try:
                        reported_error = float(float_audit['maximum_absolute_error'])
                        reported_bound = float(float_audit['absolute_error_bound'])
                        reported_scale = float(float_audit['coordinate_scale'])
                        approximate_count = float_audit['approximation_status_vertex_count']
                        vertex_count = float_audit['vertex_count']
                        sidecar_sha = float_audit['sidecar_sha256']
                    except (KeyError, TypeError, ValueError, OverflowError) as error:
                        raise ValueError('external float derivation audit is incomplete') from error
                    require(float_audit.get('status') == expected_status and
                            reported_error == maximum_error and reported_bound == bound and
                            reported_scale == scale and 0 <= approximate_count <= len(display_vertices) and
                            vertex_count == len(display_vertices) and
                            isinstance(sidecar_sha, str) and len(sidecar_sha) == 64 and
                            all(character in '0123456789abcdef' for character in sidecar_sha) and
                            maximum_error <= bound,
                            'external float derivation audit does not match independently measured geometry error')
                require(snapshot.rational_vertices is not None,
                        'external exact snapshot omitted rational coordinates')
                return snapshot, 'EXACT_RATIONAL_BOOLEAN_SIDECAR'
            if rational_vertices is None:
                return solver_candidate_snapshot((raw_vertices, raw_faces))
            if solver_boundary_profile in ('EXACT_RATIONAL_BOUNDARY_V1','EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1'):
                from .spin_problem import rational_solver_candidate_snapshot,_mesh_input_array
                snapshot,mode=rational_solver_candidate_snapshot(rational_vertices,raw_faces,
                    fixed_boundary_proof=fixed_boundary_proof,
                    allow_separated_roots=solver_boundary_profile=='EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1',
                    allow_signed_classification=solver_boundary_profile=='EXACT_RATIONAL_FIXED_CLASSIFIED_V1')
                supplied=_mesh_input_array(raw_vertices)
                canonical_vertices,_=snapshot.arrays()
                require(np.asarray(supplied,dtype='<f8',order='C').tobytes(order='C')==
                        canonical_vertices.tobytes(order='C'),
                        'external rational boundary does not bind callback mesh bytes')
                return snapshot,mode
            snapshot = MeshSnapshot((raw_vertices, raw_faces),
                                    rational_vertices=rational_vertices)
            return snapshot, 'EXACT_RATIONAL_BOOLEAN_SIDECAR'

        snapshot, variable_geometry_validation = run_stage(
            'map_variable_mesh_snapshot', build_variable_snapshot)
        if exact_geometry_identity is not None:
            float_audit = evaluation.get('exact_float_derivation_audit')
            expected_identity = (digest({
                'rational_geometry': snapshot.identity,
                'float_derivation_audit': float_audit,
            }) if float_audit is not None else digest(snapshot.identity))
            require(text_label(exact_geometry_identity, 'exact_geometry_identity') ==
                    expected_identity,
                    'external exact geometry identity does not match its rational snapshot and derivation audit')
        rational_vertices = getattr(snapshot, 'rational_vertices', None)
        vertices, faces = snapshot.arrays()
        independently_integrated = np.asarray(
            _finite_tuple(evaluation['variable_raw_moments'], 10, 'same-call variable moments')
            if validation_profile == 'LAYERED_SAME_CALL_V1' else run_stage(
                'map_variable_independent_integration',
                lambda: _raw_moments(snapshot, request.material_density)), dtype=np.float64)
    else:
        require(len(faces) == 0 and faces.dtype.kind in ('i', 'u'),
                'empty external material must not contain invalid faces')
        vertices = np.empty((0, 3), dtype=np.float64)
        faces = np.empty((0, 3), dtype=np.int32)
        independently_integrated = np.zeros(10, dtype=np.float64)
        variable_geometry_validation = 'EXACT_EMPTY_MATERIAL'
    reported_variable = np.asarray(
        _finite_tuple(evaluation['variable_raw_moments'], 10,
                      'external variable_raw_moments'), dtype=np.float64)
    design_vertices = np.asarray(request.design_vertices, dtype=np.float64)
    side = max(float(np.max(np.ptp(design_vertices, axis=0))), 1e-300)
    scales = np.asarray((side**3, side**4, side**4, side**4,
                         side**5, side**5, side**5, side**5, side**5, side**5))
    variable_error = float(np.max(np.abs(reported_variable - independently_integrated) / scales))
    require(math.isfinite(variable_error) and variable_error <= 1e-10,
            'external variable moments differ from independently integrated variable material')

    axes = np.asarray(request.target_axes, dtype=np.float64)
    origin = np.asarray(request.target_origin, dtype=np.float64)
    prepared_fixed_material = (prepared_fixed_material or
        _prepare_external_fixed_material(problem, request, stage_observer=stage_observer))
    fixed_target = np.asarray(prepared_fixed_material.moments, dtype=np.float64)
    # Cross-check the separately prepared fixed-moment input after the same
    # model-to-target rigid transform.  It catches frame/scale drift before
    # the first Newton residual is formed.
    expected_fixed = np.asarray(_transform_moments(
        request.fixed_shell_moments, axes.T, -axes.T @ origin), dtype=np.float64)
    fixed_error = float(np.max(np.abs(expected_fixed - fixed_target) / scales))
    require(math.isfinite(fixed_error) and fixed_error <= 1e-10,
            'prepared fixed moments differ from independently integrated fixed material')

    model_vertices = None
    if 'variable_vertices_model' in evaluation:
        raw_model_vertices = evaluation['variable_vertices_model']
        model_vertices = np.asarray(raw_model_vertices, dtype=np.float64)
        require(model_vertices.shape == vertices.shape,
                'external target/model material meshes must have matching vertex counts')
        if len(model_vertices):
            # The Cherchi boundary has already passed exact represented-region
            # validation in its actual solver frame above.  This second array
            # is only that same indexed boundary expressed through the locked
            # rigid target transform; rebuilding MeshSnapshot here repeats an
            # expensive exact self-intersection/region check on a derived copy.
            # Verify the full indexed coordinate correspondence instead.  A
            # converged endpoint is independently validated again when sealed
            # in model coordinates.
            def verify_model_frame_identity():
                mapped_vertices = (model_vertices - origin) @ axes
                coordinate_error = float(np.max(np.abs(mapped_vertices - vertices)))
                require(math.isfinite(coordinate_error) and coordinate_error <= 1e-12 * side,
                        'external target/model mesh coordinates do not match the locked target transform')

            run_stage('map_target_model_frame_identity', verify_model_frame_identity)
        else:
            require(len(faces) == 0,
                    'empty target material must have an empty model-frame counterpart')
            model_vertices = np.empty((0, 3), dtype=np.float64)
    else:
        require(np.array_equal(axes, np.eye(3)) and np.array_equal(origin, np.zeros(3)),
                'non-identity target frames require source-preserving model-frame Boolean vertices')

    differential = {}
    if support is None:
        jacobian = np.asarray(evaluation['moment_jacobian'], dtype=np.float64)
        require(jacobian.shape == (10, 10) and np.isfinite(jacobian).all(),
                'external variable-material Jacobian must be finite 10-by-10')
        differential['moment_jacobian'] = jacobian
    elif require_derivative:
        basis = np.array(evaluation['native_tangent_basis'], dtype=np.float64, copy=True)
        jacobian = np.array(evaluation['native_directional_jacobian'], dtype=np.float64, copy=True)
        shape = (10, support.shape[1]-1)
        require(basis.shape == shape and jacobian.shape == shape and
                np.isfinite(basis).all() and np.isfinite(jacobian).all(),
                'external directional derivative dimensions must be finite 10-by-(k-1)')
        require(np.linalg.norm(basis.T @ basis - np.eye(shape[1])) <= 1e-12 and
                np.linalg.norm(basis.T @ expected) <= 1e-12,
                'external native basis must be orthonormal and tangent to requested coefficients')
        # Dual/native coordinate span and omitted gradient equations are checked
        # by the existing native assembly, not inferred by zero-padding here.
        differential.update(native_tangent_basis=basis, native_directional_jacobian=jacobian)
    derivative_defined = evaluation['derivative_defined']
    require(type(derivative_defined) is bool,
            'external derivative_defined must be boolean')
    if require_derivative:
        require(derivative_defined,
                'external Boolean geometry derivative is unavailable')
    elif support is not None:
        require(not derivative_defined,
                'complement-only directional map must explicitly mark derivative unavailable')
    labels = {}
    for key in ('geometry_hash', 'topology_hash', 'evaluation_map_identity',
                'boolean_backend_identity', 'derivative_identity'):
        labels[key] = text_label(evaluation[key], key)
    require(labels['boolean_backend_identity'].startswith(
                'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:'),
            'strict external solver bridge requires the pinned Cherchi 7bd6c269 source identity')
    variable = tuple(float(value) for value in independently_integrated)
    return {
        'raw_moments': tuple(float(value) for value in fixed_target + independently_integrated),
        'variable_raw_moments': variable,
        'exact_geometry_identity': (text_label(exact_geometry_identity, 'exact_geometry_identity')
                                    if exact_geometry_identity is not None else
                                    (digest(snapshot.identity) if len(vertices) else
                                     'EXPLICIT_EMPTY_MATERIAL')) ,
        **differential,
        'derivative_defined': derivative_defined,
        'variable_geometry_validation': variable_geometry_validation,
        'derivative_diagnostic': str(evaluation.get('derivative_diagnostic', '')),
        **labels,
    }


def rotational_objective_definition(mode, center_of_mass_weight=0.0,
                                    spinning_top_objective='INERTIA_AND_LOWERING'):
    """Versioned meaning of rotational objective values in saved run records.

    This identifies a problem; it does not approve its geometry or derivatives.
    c3 is measured in the input-normalized target frame, not world height.
    """
    weight = scalar(center_of_mass_weight, 'center_of_mass_weight')
    require(mode in ('SPINNING_TOP', 'YO_YO'), 'rotational objective requires top or yo-yo')
    require(weight >= 0. and (mode == 'SPINNING_TOP' or weight == 0.),
            'center_of_mass_weight requires a spinning top')
    require(spinning_top_objective in ('INERTIA_AND_LOWERING', 'LOWER_CENTER_OF_MASS'),
            'unknown spinning_top_objective')
    pure = spinning_top_objective == 'LOWER_CENTER_OF_MASS'
    require(not pure or (mode == 'SPINNING_TOP' and weight == 0.),
            'pure-height endpoint requires spinning top and zero additional weight')
    return {'schema': 'faster-rotational-objective-v1', 'mode': mode,
            'spinning_top_objective': spinning_top_objective,
            'center_of_mass_weight': weight,
            'expression': 'c3' if pure else 'f_yo+weight*c3' if mode == 'SPINNING_TOP' else 'f_yo',
            'f_yo': '(Ic11^2+Ic22^2+2*Ic12^2)/Ic33^2',
            'c3': 'l3/m', 'coordinate_units': 'INPUT_NORMALIZED_TARGET_FRAME',
            'constraints': ['l1', 'l2', 'Ic13', 'Ic23'] if mode == 'SPINNING_TOP'
                           else ['l1', 'l2', 'l3', 'Ic13', 'Ic23']}


def _validate_terminal_region_profile(profile):
    require(profile in ('INDEXED_BOUNDARY_REBUILD_V1', 'OWNED_INTERIOR_DIFFERENCE_V1'),
            'unknown terminal region profile')
    if profile == 'OWNED_INTERIOR_DIFFERENCE_V1':
        module = _native_module()
        require(hasattr(module, 'make_exact_interior_difference_region') and
                hasattr(module.ExactMaterialRegion, 'construction_operands') and
                hasattr(module.ExactMaterialRegion, 'boundary_coordinate_policy'),
                'owned terminal region profile requires a replay-capable native build')


def solve_cherchi_geometry_and_seal_faster(
        problem, *, geometry_evaluator, complement_evaluator, initial_quadrics,
        evaluator_settings, solver_settings, center_of_mass_weight=0.0,
        spinning_top_objective='INERTIA_AND_LOWERING', objective_settings=None,
        event_sink=None, stage_observer=None, support_basis=None,
        require_same_call_terminal=False,solver_boundary_profile='EXACT_NEF_REGION_V1',
        validation_profile='INDEPENDENT_EACH_EVALUATION_V1', material_validator=None,
        terminal_checkpoint_sink=None, _prepared_problem=None,
        input_validation_profile='RECHECK_INPUT_V1', terminal_region_profile='INDEXED_BOUNDARY_REBUILD_V1'):
    """Run Eq. 14 on a caller-owned Cherchi map, then independently seal it.

    The callback receives target-frame normalized coefficients and a locked
    context containing the exact outer/design meshes, densities, frame, and
    fixed-material moments. It must return the variable-material mesh and its
    matching full Jacobian from the same Cherchi Boolean call. With explicit
    ``support_basis``, it instead returns a native tangent basis and genuine
    directional Jacobian; explicit initial quadrics are then required. No
    ambient derivative is fabricated. Meshes are in
    the target frame; no moments supplied by the callback are trusted without
    integration against that mesh. ``complement_evaluator`` must construct the
    exact complementary side of the same terminal quadric mesh with the
    pinned Boolean backend. The complementary candidate is not accepted merely
    because it came from the opposite side: ``SealedSpinArtifact`` must still
    prove that retained and removed materials form an exact nonoverlapping
    partition of the locked design domain.
    An empty initial sequence is permitted to select the native deterministic
    initializer controlled by the locked request seed; nonempty sequences are
    explicit solver starts and remain bound into the result identity.
    """
    require(callable(geometry_evaluator), 'geometry_evaluator must be callable')
    require(callable(complement_evaluator), 'complement_evaluator must be callable')
    require(event_sink is None or callable(event_sink), 'event_sink must be callable or None')
    require(stage_observer is None or callable(stage_observer),
            'stage_observer must be callable or None')
    require(terminal_checkpoint_sink is None or callable(terminal_checkpoint_sink),
            'terminal_checkpoint_sink must be callable or None')
    require(input_validation_profile in ('RECHECK_INPUT_V1', 'CONSTRUCTOR_PROOF_REUSE_V1'),
            'unknown input validation profile')
    _validate_terminal_region_profile(terminal_region_profile)
    require(solver_boundary_profile in ('EXACT_NEF_REGION_V1','EXACT_RATIONAL_BOUNDARY_V1','EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1'),
            'unknown solver boundary profile')
    require(validation_profile in ('INDEPENDENT_EACH_EVALUATION_V1', 'LAYERED_SAME_CALL_V1'),
            'unknown material validation profile')
    require(validation_profile != 'LAYERED_SAME_CALL_V1' or callable(material_validator),
            'layered validation requires a same-call material validator')
    if validation_profile == 'LAYERED_SAME_CALL_V1':
        _require_same_call_producer(material_validator)
    if solver_boundary_profile in ('EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1'):
        _require_same_call_producer(material_validator)
    support = _frozen_support_basis(support_basis)
    initial = tuple(_finite_tuple(row, 10, 'initial quadric') for row in initial_quadrics)
    require(support is None or bool(initial),
            'supported material solve requires explicit initial quadrics')

    def run_stage(name, operation):
        return operation() if stage_observer is None else stage_observer(name, operation)

    if _prepared_problem is None:
        problem, request, identity = _prepare(problem, evaluator_settings, solver_settings,
                                             allow_reduced_stratum=support is not None,
            reuse_input_validation=input_validation_profile=='CONSTRUCTOR_PROOF_REUSE_V1')
    else:
        require(type(_prepared_problem) is _PreparedFasterProblem and support is None,
                'invalid prepared Faster problem')
        require(problem is _prepared_problem.problem and
                evaluator_settings == _prepared_problem.evaluator_settings and
                solver_settings == _prepared_problem.solver_settings and
                input_validation_profile == _prepared_problem.input_validation_profile and
                _native_module() is _prepared_problem.native_module,
                'prepared Faster problem dependency mismatch')
        request = _prepared_problem.request
        identity = json.loads(_prepared_problem.identity_json)
    if input_validation_profile != 'RECHECK_INPUT_V1':
        identity['input_validation_profile'] = input_validation_profile
    if terminal_region_profile != 'INDEXED_BOUNDARY_REBUILD_V1':
        identity['terminal_region_profile'] = terminal_region_profile
    if solver_boundary_profile!='EXACT_NEF_REGION_V1':
        identity['solver_boundary_profile']=solver_boundary_profile
    weight = scalar(center_of_mass_weight, 'center_of_mass_weight')
    require(weight >= 0.0, 'center_of_mass_weight must be nonnegative')
    request = replace(request, center_of_mass_weight=weight,
                      spinning_top_objective=spinning_top_objective)
    if objective_settings is not None:
        from .spin2014_objective_contract import Spin2014ObjectiveContract
        contract = Spin2014ObjectiveContract(objective_settings)
        contract.require_mode(request.mode)
        require(request.mode in ('SPINNING_TOP', 'YO_YO') and support is None,
                'original rotational objective requires full-space geometry')
        identity['objective_definition'] = contract.definition
        identity['objective_contract_identity'] = contract.identity
    elif request.mode in ('SPINNING_TOP', 'YO_YO'):
        identity['objective_definition'] = rotational_objective_definition(
            request.mode, weight, spinning_top_objective)
    if support is not None:
        identity['directional_mapping'] = {
            'schema': 'faster-material-supported-directions-v1',
            'support_basis': support.tolist(),
            'support_coordinates': 'DUAL_MONOMIAL',
            'derivative_coordinates': 'NATIVE_GEOMETRY_TANGENT',
            'direction_count': support.shape[1]-1,
            'initialization': 'EXPLICIT_CALLER_STARTS_NO_PROJECTION',
        }
    axes = np.asarray(request.target_axes, dtype=np.float64)
    origin = np.asarray(request.target_origin, dtype=np.float64)
    to_target = lambda vertices: (np.asarray(vertices, dtype=np.float64) - origin) @ axes
    context = {
        'outer_vertices': to_target(request.outer_vertices),
        'outer_faces': np.asarray(request.outer_faces, dtype=np.int32),
        'design_vertices': to_target(request.design_vertices),
        'design_faces': np.asarray(request.design_faces, dtype=np.int32),
        'outer_vertices_model': np.asarray(request.outer_vertices, dtype=np.float64),
        'design_vertices_model': np.asarray(request.design_vertices, dtype=np.float64),
        'material_density': request.material_density,
        'target_origin': (0.0, 0.0, 0.0),
        'target_axes': np.eye(3),
        'fixed_shell_moments_target': tuple(_transform_moments(
            request.fixed_shell_moments, axes.T, -axes.T @ origin)),
        'fixed_regions': tuple({
            'region_id': region.region_id,
            'vertices': to_target((region.mesh.arrays()[0] -
                                   np.asarray(problem.normalization['origin'])) /
                                  problem.normalization['L']),
            'faces': region.mesh.arrays()[1],
            'density': region.density / problem.normalization['R'],
        } for region in problem.fixed_regions if not region.mesh.explicit_empty),
        'normalization': dict(problem.normalization),
        'request_digest': digest(asdict(request)),
    }
    prepared_fixed = None
    if validation_profile == 'LAYERED_SAME_CALL_V1':
        prepared_fixed = _prepare_external_fixed_material(problem, request,
                                                         stage_observer=stage_observer)
        identity.update(validation_profile=validation_profile,
                        fixed_material_preparation_identity=prepared_fixed.identity)
    if solver_boundary_profile in ('EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1'):
        fixed_proof = material_validator.__self__.fixed_boundary_proof
        require(fixed_proof is not None,'fixed subset profile requires a prepared native fixed boundary proof')
        identity['fixed_boundary_proof_identity']=fixed_proof.reference_identity
    cache = {}

    def evaluate(coefficients):
        key = tuple(float(value) for value in coefficients)
        if key not in cache:
            raw = geometry_evaluator(key, context)
            mapped = run_stage(
                'external_geometry_moment_map',
                lambda: _external_geometry_moment_map(
                    problem, request, raw, expected_coefficients=key,
                    stage_observer=stage_observer, support_basis=support,
                    solver_boundary_profile=solver_boundary_profile,
                    prepared_fixed_material=prepared_fixed,
                    validation_profile=validation_profile, material_validator=material_validator),
            )
            cache[key] = (raw, mapped)
        return cache[key][1]

    run = run_stage(
        'native_newton_solver',
        lambda: solve_spin_it_faster_with_external_moment_evaluator(
            request, evaluate, initial_quadrics=initial, event_sink=event_sink,
            **({} if objective_settings is None else {'objective_settings': objective_settings}),
            **({} if support is None else {'support_basis': support})),
    )
    identity.update(
        operation='SOLVE_CHERCHI_EXTERNAL_GEOMETRY',
        native_request_digest=digest(asdict(request)),
        initial_quadrics=[list(row) for row in initial],
        native_run=_json_value(run),
        external_backend='pinned Cherchi callback; no fallback',
    )
    if terminal_checkpoint_sink is not None:
        terminal_checkpoint_sink('NATIVE_RETURNED', _json_value(run))
    if run.get('status') != 'CONVERGED' or run.get('quadric_coefficients') is None:
        return _result(str(run.get('status', 'UNKNOWN')), identity)

    coefficients = tuple(float(value) for value in run['quadric_coefficients'])
    retained_raw, retained_map = cache.get(coefficients, (None, None))
    if retained_raw is None:
        require(not require_same_call_terminal,
                'terminal material missing from accepted same-call evaluation cache')
        retained_raw = geometry_evaluator(coefficients, context)
        retained_map = _external_geometry_moment_map(
            problem, request, retained_raw, expected_coefficients=coefficients,
            support_basis=support)
    try:
        carved_raw = run_stage(
            'terminal_complement_boolean',
            lambda: complement_evaluator(retained_raw, context),
        )
        if terminal_checkpoint_sink is not None:
            terminal_checkpoint_sink('COMPLEMENT_RETURNED', carved_raw)
        source_identity = retained_raw.get('quadric_source_identity')
        if source_identity is not None:
            require(carved_raw.get('quadric_source_identity') == source_identity,
                    'terminal complement did not use the retained material quadric source')
            retained_source_vertices = retained_raw.get('quadric_source_vertices')
            retained_source_faces = retained_raw.get('quadric_source_faces')
            carved_source_vertices = carved_raw.get('quadric_source_vertices')
            carved_source_faces = carved_raw.get('quadric_source_faces')
            require(retained_source_vertices is not None and
                    retained_source_faces is not None and
                    carved_source_vertices is not None and
                    carved_source_faces is not None,
                    'terminal complement omitted the source mesh needed for exact identity validation')
            retained_source_vertices = np.ascontiguousarray(
                retained_source_vertices, dtype='<f8')
            carved_source_vertices = np.ascontiguousarray(
                carved_source_vertices, dtype='<f8')
            retained_source_faces = np.ascontiguousarray(
                retained_source_faces, dtype='<i4')
            carved_source_faces = np.ascontiguousarray(
                carved_source_faces, dtype='<i4')
            require(retained_source_vertices.ndim == 2 and
                    retained_source_vertices.shape[1:] == (3,) and
                    retained_source_faces.ndim == 2 and
                    retained_source_faces.shape[1:] == (3,) and
                    np.isfinite(retained_source_vertices).all() and
                    np.isfinite(carved_source_vertices).all() and
                    retained_source_vertices.shape == carved_source_vertices.shape and
                    retained_source_faces.shape == carved_source_faces.shape and
                    retained_source_vertices.tobytes(order='C') ==
                    carved_source_vertices.tobytes(order='C') and
                    retained_source_faces.tobytes(order='C') ==
                    carved_source_faces.tobytes(order='C'),
                    'terminal complement source mesh differs from the retained quadric source')
            inverse_operation = {'intersection': 'subtraction',
                                 'subtraction': 'intersection'}.get(
                                     retained_raw.get('boolean_operation'))
            require(inverse_operation is not None and
                    carved_raw.get('boolean_operation') == inverse_operation,
                    'terminal complement did not evaluate the opposite side of the same quadric')
        carved_map = run_stage('terminal_complement_audit', lambda:
            _external_geometry_moment_map(
                problem, request, carved_raw, require_derivative=False,
                support_basis=support, stage_observer=stage_observer))
    except Exception as exc:
        identity['terminal_complement'] = {
            'status': 'FAILED',
            'diagnostic': f'{type(exc).__name__}: {exc}',
            'backend': 'pinned Cherchi callback; no fallback',
        }
        return _result('FINAL_COMPLEMENT_AUDIT_FAILED', identity,
                       seal_error=f'{type(exc).__name__}: {exc}')
    require(retained_map['geometry_hash'] == run.get('trace', [{}])[-1].get(
        'material_geometry_hash', retained_map['geometry_hash']),
        'Cherchi terminal geometry differs from the Newton terminal event')
    require(len(retained_raw['variable_faces']) > 0,
            'converged Cherchi result has empty retained variable material')
    def target_to_model(mesh):
        vertices = np.asarray(mesh['variable_vertices'], dtype=np.float64)
        if 'variable_vertices_model' in mesh:
            model = np.asarray(mesh['variable_vertices_model'], dtype=np.float64)
        else:
            model = vertices @ axes.T + origin
        rational_vertices = mesh.get('variable_rational_vertices')
        if rational_vertices is not None:
            rational_vertices = _target_rationals_to_normalized_model(
                rational_vertices, request)
            model = np.asarray([[float(fraction_from_rational_text(value)) for value in row]
                                for row in rational_vertices], dtype=np.float64)
        return model, np.asarray(mesh['variable_faces'], dtype=np.int32), rational_vertices

    retained = target_to_model(retained_raw)
    carved = target_to_model(carved_raw)
    total_target = np.asarray(retained_map['raw_moments'], dtype=np.float64)
    if objective_settings is not None:
        # Independent original equations, with units applied by the same
        # declared configuration; neither native AD nor old Faster formulas.
        f, constraints = _independent_spin2014_rotational_physics(total_target, objective_settings)
        gradient = np.zeros(10)
        jacobian = np.zeros((len(constraints), 10))
        for i in range(10):
            perturbed = total_target.astype(complex)
            perturbed[i] += 1e-24j
            fi, ci = _independent_spin2014_rotational_physics(perturbed, objective_settings)
            gradient[i] = np.imag(fi)/1e-24
            jacobian[:, i] = np.imag(ci)/1e-24
        dual = np.asarray(coefficients).copy()
        dual[7:] *= 2
        dual /= np.linalg.norm(dual)
        mu = scalar(run['multiplier_mu'], 'original terminal mu', positive=True)
        residual = float(np.linalg.norm(np.r_[gradient-jacobian.T@np.asarray(run['multipliers'])+mu*dual, constraints]))
        require(np.isfinite(residual) and residual <= request.settings.equation_tolerance,
                'independent original objective terminal residual failed')
        require(run.get('objective_profile') == 'SPIN2014_OBJECTIVES_QUADRIC_V1',
                'native original objective profile mismatch')
        identity['original_objective_terminal_check'] = dict(objective=float(f),
            raw_equation_residual=residual, objective_contract_identity=contract.identity)
    try:
        artifact = run_stage('terminal_material_seal', lambda:
            _seal(problem, request, retained, carved,
                  (request.outer_vertices, request.outer_faces), total_target,
                  stage_observer=stage_observer,
                  reuse_problem_validation=input_validation_profile=='CONSTRUCTOR_PROOF_REUSE_V1',
                  terminal_region_profile=terminal_region_profile))
    except (ValueError, RuntimeError) as exc:
        return _result('FINAL_GEOMETRY_AUDIT_FAILED', identity,
                       seal_error=str(exc))
    identity['external_terminal_binding'] = {
        'retained_geometry': retained_map['geometry_hash'],
        'retained_topology': retained_map['topology_hash'],
        'complement_geometry': carved_map['geometry_hash'],
        'complement_topology': carved_map['topology_hash'],
        'fixed_moments_target': [float(x) for x in context['fixed_shell_moments_target']],
        'variable_moments_target': [float(x) for x in retained_map['variable_raw_moments']],
        'total_moments_target': [float(x) for x in total_target],
        'complement_construction': 'pinned-Cherchi-same-quadric-opposite-side-with-exact-partition-proof',
    }
    identity['terminal_complement'] = {
        'status': 'PASSED',
        'geometry_hash': carved_map['geometry_hash'],
        'topology_hash': carved_map['topology_hash'],
        'evaluation_map_identity': carved_map['evaluation_map_identity'],
        'backend': carved_map['boolean_backend_identity'],
        'quadric_source_identity': carved_raw.get('quadric_source_identity'),
        'boolean_operation': carved_raw.get('boolean_operation'),
    }
    return _result('CONVERGED', identity, artifact)


def _seal(problem, request, retained, carved, outer, total_moments, *, stage_observer=None,
          reuse_problem_validation=False, terminal_region_profile='INDEXED_BOUNDARY_REBUILD_V1'):
    _validate_terminal_region_profile(terminal_region_profile)
    def stage(name, operation):
        return operation() if stage_observer is None else stage_observer(name, operation)
    # Compare actual native outer output against the exact request arrays. The
    # separately sealed raw outer is the locked original, never a roundtrip mesh.
    for actual, expected, dtype in zip(outer, (request.outer_vertices, request.outer_faces), ('<f8', '<i4')):
        require(np.asarray(actual, dtype=dtype).shape == np.asarray(expected, dtype=dtype).shape and
                np.asarray(actual, dtype=dtype).tobytes() == np.asarray(expected, dtype=dtype).tobytes(),
                'native output outer bytes differ from normalized locked request')
    norm = problem.normalization
    def model_mesh(mesh):
        v, f = mesh[:2]
        if not len(v) and not len(f):
            return MeshSnapshot(None)
        rational_vertices = mesh[2] if len(mesh) > 2 else None
        if rational_vertices is None:
            return MeshSnapshot((np.asarray(v)*norm['L'] + np.array(norm['origin']),
                                 np.asarray(f, dtype=np.int32)))
        exact_world = []
        for row in rational_vertices:
            exact_row = []
            for axis, value in enumerate(row):
                world = fraction_from_rational_text(value) * Fraction.from_float(float(norm['L']))
                world += Fraction.from_float(float(norm['origin'][axis]))
                exact_row.append(fraction_to_rational_text(world))
            exact_world.append(tuple(exact_row))
        world_vertices = np.asarray([[float(fraction_from_rational_text(value)) for value in row]
                                     for row in exact_world], dtype=np.float64)
        return MeshSnapshot((world_vertices, np.asarray(f, dtype=np.int32)),
                            rational_vertices=exact_world)
    kept = stage('seal_retained_exact_snapshot', lambda: model_mesh(retained))
    removed = stage('seal_carved_exact_snapshot', lambda: model_mesh(carved))
    if removed.explicit_empty:
        material = problem.original
    elif terminal_region_profile == 'OWNED_INTERIOR_DIFFERENCE_V1':
        material = stage('seal_final_owned_difference', lambda:
            MeshSnapshot.from_exact_region(_native_module().make_exact_interior_difference_region(
                problem.original.exact_region(), removed.exact_region())))
    else:
        ov, of = problem.original.arrays()
        cv, cf = removed.arrays()
        material_vertices = np.vstack((ov, cv))
        material_faces = np.vstack((of, cf[:, ::-1]+len(ov)))
        exact_removed = removed.rational_vertices
        if exact_removed is None:
            exact_material = None
        else:
            exact_outer = problem.original.rational_vertices
            if exact_outer is None:
                exact_outer = tuple(tuple(
                    fraction_to_rational_text(Fraction.from_float(float(value))) for value in row)
                    for row in ov)
            exact_material = tuple(exact_outer) + tuple(exact_removed)
        material = stage('seal_final_exact_snapshot', lambda:
            MeshSnapshot((material_vertices, material_faces), rational_vertices=exact_material))
    claim = _transform_moments(total_moments, request.target_axes, request.target_origin)
    return stage('seal_independent_artifact_check', lambda:
        SealedSpinArtifact.create(problem, outer=problem.original,
            design_results=[(problem.design_regions[0].region_id, kept, removed)],
            final_material=material, newly_carved=removed,
            solver_claim={'frame': FRAME, 'order': list(ORDER), 'moments': claim},
            _reuse_problem_validation=reuse_problem_validation))


def evaluate_and_seal_faster_candidate(problem, coefficients, *, evaluator_settings, solver_settings):
    """Evaluate p and -p exactly once each; seal material without optimizing."""
    problem, request, identity = _prepare(
        problem, evaluator_settings, solver_settings, allow_reduced_stratum=True)
    coefficients = _finite_tuple(coefficients, 10, 'candidate coefficients')
    retained = evaluate_spin_it_faster_quadric(request, coefficients, evaluator_settings=evaluator_settings)
    if evaluator_settings.derivative_method != 'DISCRETE_REDUCED_STRATUM':
        require(retained.derivative_defined,
                'unsupported discrete derivative: ' + retained.derivative_diagnostic)
    carved = evaluate_spin_it_faster_quadric(request, tuple(-x for x in coefficients),
                                            evaluator_settings=evaluator_settings)
    if evaluator_settings.derivative_method != 'DISCRETE_REDUCED_STRATUM':
        require(carved.derivative_defined,
                'unsupported discrete complement derivative: ' + carved.derivative_diagnostic)
    require(retained.request_hash == carved.request_hash, 'native evaluation request identity mismatch')
    artifact = _seal(problem, request, (retained.material_vertices, retained.material_faces),
                     (carved.material_vertices, carved.material_faces),
                     (retained.outer_vertices, retained.outer_faces), retained.total_raw_moments)
    identity.update(operation='CANDIDATE', coefficients=list(coefficients), native_evaluation_count=2,
                    native_evaluations=[asdict(retained), asdict(carved)])
    return _result('NOT_OPTIMIZED_CANDIDATE', identity, artifact)


def solve_and_seal_faster(problem, *, initial_quadrics, evaluator_settings, solver_settings):
    """Preserve an actual native run; only converged exported geometry may seal."""
    problem, request, identity = _prepare(problem, evaluator_settings, solver_settings)
    initial = tuple(_finite_tuple(row, 10, 'initial quadric') for row in initial_quadrics)
    require(bool(initial), 'strict bridge requires an explicit nonempty initial sequence')
    run = solve_spin_it_faster(request, initial_quadrics=initial, evaluator_settings=evaluator_settings)
    identity.update(operation='SOLVE', initial_quadrics=initial, native_run=run.diagnostics,
                    derivative_metadata=run.derivative_metadata)
    if run.status != 'CONVERGED':
        return _result(run.status, identity)
    try:
        artifact = _seal(problem, request, (run.material_vertices, run.material_faces),
            (run.cavity_vertices, run.cavity_faces), (run.outer_vertices, run.outer_faces),
            run.final_raw_moments)
    except (ValueError, RuntimeError) as exc:
        return _result(run.status, identity, seal_error=str(exc))
    return _result(run.status, identity, artifact)


def solve_fullrank_and_seal_faster(
        problem, *, initial_quadrics, evaluator_settings, solver_settings,
        center_of_mass_weight=0.0):
    """Run the discrete full-rank solver and bind its result to a material seal.

    This is an explicit alternative to ``solve_and_seal_faster``: it calls the
    canonical full-rank Newton entry point, preserves its complete run record,
    re-evaluates both quadric orientations with that same discrete geometry
    map, then seals the exact material partition. The center-of-mass weight is
    an explicit request override because the locked geometric target does not
    itself encode this solver objective parameter.
    """
    problem, request, identity = _prepare(problem, evaluator_settings, solver_settings)
    center_of_mass_weight = scalar(center_of_mass_weight, 'center_of_mass_weight')
    require(center_of_mass_weight >= 0.,
            'center_of_mass_weight must be nonnegative')
    request = replace(request, center_of_mass_weight=center_of_mass_weight)
    identity['request_overrides'] = {
        'center_of_mass_weight': center_of_mass_weight,
    }
    identity['native_request_digest'] = digest(asdict(request))
    initial = tuple(_finite_tuple(row, 10, 'initial quadric') for row in initial_quadrics)
    require(bool(initial), 'strict full-rank bridge requires an explicit nonempty initial sequence')
    run = solve_spin_it_faster_fullrank_quadric(
        request, initial_quadrics=initial, evaluator_settings=evaluator_settings)
    native_run = asdict(run)
    native_run['solver_kind'] = 'FULL_RANK_DISCRETE'
    identity.update(operation='SOLVE_FULL_RANK', initial_quadrics=initial,
                    native_run=native_run)
    if run.status != 'CONVERGED':
        return _result(run.status, identity)
    require(run.final_geometry_verified and run.final_derivative_defined and
            run.quadric_coefficients is not None and run.final_geometry_hash,
            'full-rank solver reported convergence without its verified final state')
    retained = evaluate_spin_it_faster_fullrank_quadric(
        request, run.quadric_coefficients, evaluator_settings=evaluator_settings)
    carved = evaluate_spin_it_faster_fullrank_quadric(
        request, tuple(-value for value in run.quadric_coefficients),
        evaluator_settings=evaluator_settings)
    require(retained.derivative_defined and carved.derivative_defined,
            'full-rank final material/complement derivative is unavailable')
    require(retained.geometry_hash == run.final_geometry_hash,
            're-evaluated full-rank geometry differs from the solver final state')
    require(run.final_total_raw_moments is not None and
            np.allclose(retained.total_raw_moments, run.final_total_raw_moments,
                        rtol=1e-10, atol=1e-10),
            're-evaluated full-rank moments differ from the solver final state')
    try:
        artifact = _seal(
            problem, request,
            (retained.material_vertices, retained.material_faces),
            (carved.material_vertices, carved.material_faces),
            (request.outer_vertices, request.outer_faces),
            retained.total_raw_moments)
    except (ValueError, RuntimeError) as exc:
        return _result(run.status, identity, seal_error=str(exc))
    return _result(run.status, identity, artifact)


def _normalized_material_moments(mesh, density, normalization):
    """Integrate CAD reference material without replacing exact authority.

The expected retained material is not an OBJ operand. Its rational boundary
may have a degenerate binary64 display. Normalize the rational coordinates
before integration instead of constructing an invalid display MeshSnapshot.
"""
    length = scalar(normalization['L'], 'normalization length', positive=True)
    reference_density = scalar(normalization['R'], 'reference density', positive=True)
    density = scalar(density, 'material density', positive=True)
    origin = np.asarray(_finite_tuple(normalization['origin'], 3, 'normalization origin'))
    vertices, faces = mesh.arrays()
    exact = mesh.rational_vertices
    if exact is None:
        return _raw_moments(MeshSnapshot(((vertices-origin)/length, faces)),
                            density/reference_density)
    from .mass_moments import integrate_exact_rational_raw_moments
    exact_length = Fraction(length)
    exact_origin = [Fraction(float(x)) for x in origin]
    rows = [[(fraction_from_rational_text(x)-exact_origin[i])/exact_length
             for i,x in enumerate(row)] for row in exact]
    raw = integrate_exact_rational_raw_moments(rows, faces.tolist())
    factor = Fraction(density)/Fraction(reference_density)
    values = [float(x*factor) for x in raw]
    require(all(math.isfinite(x) for x in values) and values[0]>0.,
            'normalized exact CAD reference moments are invalid')
    return values


def export_spin_it_faster_quadric_material_step(
        result, *, occt_executable, brep_path, step_path):
    """Round-trip sealed Faster quadric material through exact OCCT CAD.

    This is a geometry handoff/method-validation path, not a paper-case claim.
    OCCT uses an exact half-space for planar candidates and converged
    prescribed-waterline plane results, and analytical ellipsoid, elliptic-
    paraboloid, elliptic-cylinder, or hyperboloid Booleans for supported
    quadrics, with closed faceted solids for fixed regions. The
    reloaded B-rep and STEP moments are independently compared with each
    authenticated material region in the seal. Fixed inputs must currently be
    single-shell OBJ solids; unsupported quadric classes fail explicitly.
    """
    require(isinstance(result, FasterMaterialResult),
            'expected a result returned by the strict Faster material bridge')
    identity = result.identity
    require(digest(material_run_scientific_identity(identity)) == result.scientific_digest,
            'Faster run identity digest mismatch')
    require(identity.get('artifact_scientific_digest') ==
            (None if result.artifact is None else result.artifact.scientific_digest),
            'Faster result is not bound to its sealed material artifact')
    require(result.artifact is not None and result.status in {
        'NOT_OPTIMIZED_CANDIDATE', 'CONVERGED'},
        'CAD export requires a sealed Faster candidate or converged solver result')
    artifact = result.artifact
    artifact.revalidate()
    problem = artifact.problem
    norm = problem.normalization
    origin = np.asarray(norm['origin'], dtype=np.float64)
    length = scalar(norm['L'], 'normalization length', positive=True)
    density_scale = scalar(norm['R'], 'normalization density', positive=True)

    operation = identity.get('operation')
    cad_target_origin = None
    cad_units_per_normalized_unit = None
    if result.status == 'NOT_OPTIMIZED_CANDIDATE':
        require(operation == 'CANDIDATE' and
                identity.get('native_evaluation_count') == 2,
                'candidate CAD export requires its bound p/-p native evaluations')
        coefficients = _finite_tuple(identity.get('coefficients', ()), 10,
                                      'candidate quadric coefficients')
        evaluations = identity.get('native_evaluations')
        require(isinstance(evaluations, list) and len(evaluations) == 2 and
                isinstance(evaluations[0].get('quadric_solid_is_positive_region'), bool),
                'candidate run lacks its native quadric-region sign witness')
        quadric_type = _cad_quadric_type(evaluations[0], coefficients)
        require(quadric_type in {'PLANE', 'ELLIPSOID', 'PARABOLOID',
                                 'HYPERBOLIC_PARABOLOID', 'ELLIPTIC_CYLINDER',
                                 'ONE_SHEET_HYPERBOLOID',
                                 'TWO_SHEET_HYPERBOLOID'} and
                _cad_quadric_type(evaluations[1], coefficients) == quadric_type,
                'OCCT CAD bridge does not support this candidate quadric type')
        positive_region = evaluations[0]['quadric_solid_is_positive_region']
    elif operation == 'SEAL_BUOYANT_PLANE_A':
        require(result.status == 'CONVERGED' and
                identity.get('prescribed_waterline_stable') is True and
                identity.get('equilibrium_constraints_verified') is True and
                identity.get('eq7_numerical_conditions_verified') is True,
                'buoyant CAD export requires the sealed stable solver result')
        payload = identity.get('request')
        require(isinstance(payload, dict),
                'buoyant CAD evidence omits its locked request snapshot')
        request_data = dict(payload)
        settings = request_data.get('settings')
        require(isinstance(settings, dict),
                'buoyant CAD request snapshot omits solver settings')
        request_data['settings'] = SpinSolverConfig(**settings)
        request = SpinItFasterRequest(**request_data)
        require(request.mode == 'BUOYANT_STABILITY' and
                request.request_digest == identity.get('request_digest'),
                'buoyant CAD request snapshot does not match its sealed digest')
        require(identity.get('problem_digest') == problem.input_digest,
                'buoyant CAD problem differs from the sealed A problem')
        require(len(problem.design_regions) == 1 and
                problem.design_regions[0].density == request.material_density,
                'buoyant CAD design region differs from the locked request')
        for snapshot, vertices, faces, label in (
                (problem.original, request.outer_vertices, request.outer_faces,
                 'buoyant request outer'),
                (problem.design_regions[0].mesh, request.design_vertices,
                 request.design_faces, 'buoyant request design')):
            actual_v, actual_f = snapshot.arrays()
            expected_v = np.asarray(vertices, dtype='<f8')
            expected_f = np.asarray(faces, dtype='<i4')
            require(actual_v.shape == expected_v.shape and
                    actual_f.shape == expected_f.shape and
                    actual_v.tobytes() == expected_v.tobytes() and
                    actual_f.tobytes() == expected_f.tobytes(),
                    f'{label} differs from the sealed CAD problem')
        axes = np.asarray(request.target_axes, dtype=np.float64)
        target_origin = np.asarray(request.target_origin, dtype=np.float64)
        require(np.array_equal(np.asarray(identity.get('target_axes'), dtype=np.float64), axes) and
                np.array_equal(np.asarray(identity.get('target_origin'), dtype=np.float64),
                               target_origin),
                'buoyant CAD target frame differs from the sealed request')
        plane_coefficients = _finite_tuple(
            identity.get('plane_coefficients_target', ()), 4,
            'sealed buoyant plane coefficients')
        split = evaluate_spin_it_faster_plane_halfspaces(
            request, coefficients=plane_coefficients)

        def verify_plane_region(target_mesh, sealed_region, label):
            target_vertices = np.asarray(target_mesh['vertices'], dtype=np.float64)
            target_faces = np.asarray(target_mesh['faces'], dtype=np.int32)
            model_vertices = target_vertices @ axes.T + target_origin
            sealed_vertices, sealed_faces = sealed_region.arrays()
            require(model_vertices.shape == sealed_vertices.shape and
                    target_faces.shape == sealed_faces.shape and
                    model_vertices.tobytes() == sealed_vertices.tobytes() and
                    target_faces.tobytes() == sealed_faces.tobytes(),
                    f'recomputed buoyant {label} differs from sealed material geometry')

        require(len(artifact.design_results) == 1,
                'buoyant plane artifact must contain one design result')
        _, retained, carved = artifact.design_results[0]
        verify_plane_region(split['positive'], retained, 'retained region')
        verify_plane_region(split['negative'], carved, 'carved region')
        raw_coefficients = np.asarray((plane_coefficients[0],
                                       plane_coefficients[1],
                                       plane_coefficients[2],
                                       plane_coefficients[3],
                                       0., 0., 0., 0., 0., 0.), dtype=np.float64)
        coefficient_norm = float(np.linalg.norm(raw_coefficients))
        require(math.isfinite(coefficient_norm) and coefficient_norm > 0.,
                'sealed buoyant plane has a zero/nonfinite coefficient vector')
        coefficients = tuple(float(value) for value in
                             raw_coefficients / coefficient_norm)
        quadric_type = 'PLANE'
        positive_region = True
        # Buoyant plane coefficients are in the request's physical model
        # coordinates, while this CAD bridge serializes the sealed geometry
        # in the artifact's normalized model frame.
        cad_target_origin = target_origin - origin
        cad_units_per_normalized_unit = 1. / length
    elif operation == 'SOLVE_CHERCHI_EXTERNAL_GEOMETRY':
        # Do not recover this sign/type by re-evaluating a CGAL material. The
        # analytic carrier is classified from the actual bound coefficients;
        # final OCCT moments still must match the independently sealed Cherchi
        # material. This does not assert analytic/discrete solids are equal.
        from .cherchi_cad import cherchi_cad_parameters
        parameters = cherchi_cad_parameters(result)
        coefficients = parameters['coefficients']
        quadric_type = parameters['quadric_type']
        positive_region = parameters['quadric_solid_is_positive_region']
    else:
        native_run = identity.get('native_run', {})
        require(operation in {'SOLVE', 'SOLVE_FULL_RANK'} and
                native_run.get('status') == 'CONVERGED',
                'optimized CAD export requires a converged bound native run')
        coefficients = _finite_tuple(
            native_run.get('quadric_coefficients', ()), 10,
            'converged quadric coefficients')
        # Re-evaluate the exact locked request to recover the bounded-ellipsoid
        # sign used by native Boolean selection. This does not alter the solve.
        evaluator_settings = SpinEvaluatorConfig(**identity['evaluator_settings'])
        solver_settings = SpinSolverConfig(**identity['solver_settings'])
        locked_problem, request, rebuilt_identity = _prepare(
            problem, evaluator_settings, solver_settings)
        request_overrides = identity.get('request_overrides', {})
        if operation == 'SOLVE_FULL_RANK':
            require(request_overrides.keys() == {'center_of_mass_weight'} and
                    native_run.get('solver_kind') == 'FULL_RANK_DISCRETE',
                    'full-rank CAD evidence lacks its explicit solver request identity')
            weight = scalar(request_overrides['center_of_mass_weight'],
                            'center_of_mass_weight')
            require(weight >= 0., 'center_of_mass_weight must be nonnegative')
            request = replace(request, center_of_mass_weight=weight)
            rebuilt_identity['native_request_digest'] = digest(asdict(request))
        else:
            require(not request_overrides,
                    'legacy solver CAD evidence contains unsupported request overrides')
        require(rebuilt_identity['native_request_digest'] == identity.get('native_request_digest'),
                'CAD request does not match the Faster run request identity')
        sign_witness = evaluate_spin_it_faster_quadric(
            request, coefficients, evaluator_settings=evaluator_settings)
        require(sign_witness.request_hash and sign_witness.base_request_hash,
                'native CAD sign evaluation lacks request binding')
        quadric_type = _cad_quadric_type(asdict(sign_witness), coefficients)
        require(quadric_type in {'ELLIPSOID', 'PARABOLOID',
                                 'HYPERBOLIC_PARABOLOID', 'ELLIPTIC_CYLINDER',
                                 'ONE_SHEET_HYPERBOLOID',
                                 'TWO_SHEET_HYPERBOLOID'},
                'OCCT CAD bridge does not support this solved quadric type: ' +
                repr((sign_witness.interface_type, sign_witness.interface_subtype)))
        positive_region = bool(sign_witness.quadric_solid_is_positive_region)
        problem = locked_problem

    axes = np.asarray(identity.get('target_axes'), dtype=np.float64)
    target_origin = np.asarray(identity.get('target_origin'), dtype=np.float64)
    require(axes.shape == (3, 3) and target_origin.shape == (3,) and
            np.isfinite(axes).all() and np.isfinite(target_origin).all(),
            'Faster run is missing its locked target frame')
    if cad_target_origin is None:
        cad_target_origin = target_origin
        cad_units_per_normalized_unit = 1.0
    cad_coefficients = quadric_coefficients_to_cad_frame(
        coefficients, target_axes=axes, target_origin=cad_target_origin,
        original_units_per_normalized_unit=cad_units_per_normalized_unit)

    retained_by_id = {region_id: retained for region_id, retained, _ in artifact.design_results}
    records = []
    # The optimized/design entity must come first. Feed the locked source
    # domain to OCCT, then compare its Boolean result with the sealed retained
    # mesh; feeding the retained mesh would apply the quadric a second time.
    for region in problem.design_regions:
        retained = retained_by_id[region.region_id]
        if not retained.explicit_empty:
            records.append(('design', region.region_id, region.mesh, retained, region.density))
    for region in problem.fixed_regions:
        if not region.mesh.explicit_empty:
            records.append(('fixed', region.region_id, region.mesh, region.mesh, region.density))
    require(bool(records) and records[0][0] == 'design',
            'sealed Faster artifact has no nonempty editable material region')

    artifact_regions = artifact.facts['regions']
    expected_by_id = {}
    for region_record in artifact_regions:
        if region_record['explicit_empty']:
            continue
        expected_by_id[(region_record['id'], region_record['kind'])] = np.asarray(
            region_record['normalized']['moments'], dtype=np.float64)

    def normalized_mesh(snapshot):
        vertices, faces = snapshot.arrays()
        mapped = (vertices - origin) / length
        return mapped, faces

    def obj_text(vertices, faces):
        lines = [f'v {x:.17g} {y:.17g} {z:.17g}' for x, y, z in vertices]
        lines.extend('f ' + ' '.join(str(int(index) + 1) for index in face)
                     for face in faces)
        return '\n'.join(lines) + '\n'

    protocol_lines = ['FASTER_3D_QUADRIC_MATERIAL_V2',
                      '1' if positive_region else '0',
                      ' '.join(format(float(value), '.17g') for value in cad_coefficients),
                      quadric_type,
                      str(len(records))]
    entity_ids, roles, densities, expected_entity_moments = [], [], [], []
    with tempfile.TemporaryDirectory(prefix='spin-faster-3d-cad-') as temporary:
        temporary_path = Path(temporary)
        for index, (role, region_id, input_mesh, expected_mesh, physical_density) in enumerate(records):
            vertices, faces = normalized_mesh(input_mesh)
            require(len(vertices) > 0 and len(faces) > 0,
                    f'{role} material {region_id} has no closed mesh')
            density = float(physical_density / density_scale)
            require(math.isfinite(density) and density > 0.0,
                    f'{role} material {region_id} has invalid normalized density')
            record_kind = 'retained' if role == 'design' else 'fixed'
            expected = expected_by_id.get((region_id, record_kind))
            require(expected is not None and expected.shape == (10,),
                    f'sealed artifact omits material moments for {region_id}')
            expected_vertices, expected_faces = normalized_mesh(expected_mesh)
            recomputed = np.asarray(_normalized_material_moments(
                expected_mesh, physical_density, norm), dtype=np.float64)
            expected_mass = max(abs(float(expected[0])), 1e-14)
            expected_length = max(float(np.ptp(expected_vertices, axis=0).max()), 1e-14)
            independent_scales = expected_mass * np.asarray(
                (1., *([expected_length] * 3), *([expected_length**2] * 6)))
            independent_error = float(np.max(np.abs(recomputed - expected) /
                                             independent_scales))
            require(independent_error <= 1e-10,
                    'independent CAD input moments differ from sealed region '
                    f'{region_id} (normalized error {independent_error:.17g})')
            safe_id = digest({'role': role, 'region_id': region_id})
            obj_path = temporary_path / f'entity-{index}.obj'
            obj_path.write_text(obj_text(vertices, faces), encoding='ascii')
            protocol_lines.append(
                f'ENTITY {role} {safe_id} {density:.17g} {json.dumps(str(obj_path))}')
            entity_ids.append(safe_id)
            roles.append(role)
            densities.append(density)
            expected_entity_moments.append(expected)
        protocol_path = temporary_path / 'faster-material.txt'
        protocol_path.write_text('\n'.join(protocol_lines) + '\n', encoding='ascii')

        executable = Path(occt_executable).resolve()
        require(executable.is_file() and os.access(executable, os.X_OK),
                'occt_executable must be the pinned OCCT bridge')
        brep = Path(brep_path).absolute()
        step = Path(step_path).absolute()
        sidecar = step.with_suffix(step.suffix + '.identity.json')
        require(brep != step, 'B-rep and STEP paths must differ')
        if brep.exists() or step.exists() or sidecar.exists():
            raise FileExistsError('CAD paths must be new; refusing to overwrite prior evidence')
        require(brep.parent.is_dir() and step.parent.is_dir(),
                'CAD output parent directories must already exist')
        # Keep unverified CAD files in private sibling directories. The exact
        # Boolean and moment checks below must finish before any final path is
        # published or given a trusted sidecar.
        brep_stage_directory = tempfile.TemporaryDirectory(
            prefix='.spin-faster-brep-', dir=str(brep.parent))
        step_stage_directory = tempfile.TemporaryDirectory(
            prefix='.spin-faster-step-', dir=str(step.parent))
        staged_brep = Path(brep_stage_directory.name) / 'material.brep'
        staged_step = Path(step_stage_directory.name) / 'material.step'
        completed = subprocess.run(
            [str(executable), '--faster-3d-quadric-material', str(protocol_path),
             str(staged_brep), str(staged_step)], check=False, capture_output=True,
            text=True, timeout=60)
    def cleanup_staged_outputs():
        brep_stage_directory.cleanup()
        step_stage_directory.cleanup()

    if completed.returncode != 0:
        cleanup_staged_outputs()
        raise RuntimeError('OCCT Faster 3D material export failed: ' + completed.stderr.strip())
    try:
        cad_result = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        cleanup_staged_outputs()
        raise RuntimeError('OCCT Faster 3D export returned invalid evidence JSON') from error
    require(cad_result.get('status') == 'FASTER_3D_MATERIAL_STEP_BREP_VALIDATED',
            'OCCT did not validate its 3D B-rep/STEP round trip')
    require(cad_result.get('material_entity_ids') == entity_ids and
            cad_result.get('entity_roles') == roles and
            np.array_equal(np.asarray(cad_result.get('entity_densities'), dtype=np.float64),
                           np.asarray(densities, dtype=np.float64)),
            'OCCT CAD material identities, roles, or densities differ from the seal')
    observed_entities = np.asarray(cad_result.get('entity_mass_moments'), dtype=np.float64)
    require(observed_entities.shape == (len(records), 10) and
            np.isfinite(observed_entities).all(),
            'OCCT did not return finite per-region mass moments')
    all_errors = []
    for index, expected in enumerate(expected_entity_moments):
        mass = max(abs(float(expected[0])), 1e-14)
        vertices, _ = normalized_mesh(records[index][3])
        characteristic_length = max(float(np.ptp(vertices, axis=0).max()), 1e-14)
        scales = mass * np.asarray((1., *([characteristic_length] * 3),
                                    *([characteristic_length**2] * 6)))
        all_errors.extend((np.abs(observed_entities[index] - expected) / scales).tolist())
    maximum_error = max(all_errors, default=math.inf)
    if all_errors:
        worst_flat_index = int(np.argmax(all_errors))
        worst_region_index, worst_moment_index = divmod(worst_flat_index, 10)
        worst_region = records[worst_region_index][1]
        worst_delta = float(observed_entities[worst_region_index, worst_moment_index] -
                            expected_entity_moments[worst_region_index][worst_moment_index])
    else:
        worst_region, worst_moment_index, worst_delta = '<none>', -1, math.nan
    if not (math.isfinite(maximum_error) and maximum_error <= 1e-6):
        cleanup_staged_outputs()
        raise ValueError(
            'analytic OCCT material differs from the sealed triangulated material '
            f'(normalized moment error {maximum_error:.17g}; region={worst_region!r}, '
            f'moment_index={worst_moment_index}, raw_delta={worst_delta:.17g})')
    expected_total = np.asarray(artifact.facts['normalized']['moments'], dtype=np.float64)
    step_total = np.asarray(cad_result['raw_mass_moments']['step'], dtype=np.float64)
    outer_vertices, _ = normalized_mesh(problem.original)
    total_length = max(float(np.ptp(outer_vertices, axis=0).max()), 1e-14)
    total_scales = max(abs(float(expected_total[0])), 1e-14) * np.asarray(
        (1., *([total_length] * 3), *([total_length**2] * 6)))
    total_error = float(np.max(np.abs(step_total - expected_total) / total_scales))
    if not (math.isfinite(total_error) and total_error <= 1e-6):
        cleanup_staged_outputs()
        raise ValueError('reloaded STEP total mass moments differ from the sealed material')
    classification = ('METHOD_VALIDATION_ONLY' if result.status == 'NOT_OPTIMIZED_CANDIDATE'
                      else 'SOLVER_RESULT_CAD_VALIDATION')
    def file_sha256(path):
        hasher = hashlib.sha256()
        with Path(path).open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                hasher.update(block)
        return hasher.hexdigest()

    cad_input_identity = digest({
        'schema': 'faster-3d-quadric-cad-input-v2',
        'coefficients_input_normalized_target': coefficients,
        'coefficients_input_normalized_model': cad_coefficients,
        'quadric_type': quadric_type,
        'quadric_solid_is_positive_region': positive_region,
        'entities': [{'role': row[0], 'region_id': row[1],
                      'density_normalized': densities[index],
                      'input_mesh': row[2].identity,
                      'expected_material_mesh': row[3].identity}
                     for index, row in enumerate(records)],
        'artifact_scientific_digest': artifact.scientific_digest,
    })
    enriched = {**cad_result, 'classification': classification,
                'source_operation': operation,
                'source_scope': identity.get('scope'),
                'faster_result_status': result.status,
                'faster_result_scientific_digest': result.scientific_digest,
                'artifact_scientific_digest': artifact.scientific_digest,
                'artifact_evidence_digest': artifact.evidence_digest,
                'cad_input_digest': cad_input_identity,
                'occt_executable_sha256': file_sha256(executable),
                'brep_sha256': file_sha256(staged_brep),
                'step_sha256': file_sha256(staged_step),
                'material_region_ids': [record[1] for record in records],
                'maximum_normalized_material_moment_error': maximum_error,
                'normalized_total_moment_error': total_error,
                'region_moment_errors': all_errors}
    staged_sidecar = Path(step_stage_directory.name) / 'material.step.identity.json'
    staged_sidecar.write_text(canonical(enriched) + '\n', encoding='ascii')
    published = []
    try:
        # Hard links are atomic and fail rather than overwriting a path that
        # appeared after the initial no-overwrite check. Publish the sidecar
        # last so its presence means all three files passed validation.
        for staged, final in ((staged_brep, brep), (staged_step, step),
                              (staged_sidecar, sidecar)):
            os.link(staged, final)
            published.append(final)
    except OSError:
        for final in published:
            try:
                final.unlink()
            except FileNotFoundError:
                pass
        cleanup_staged_outputs()
        raise
    cleanup_staged_outputs()
    return enriched
