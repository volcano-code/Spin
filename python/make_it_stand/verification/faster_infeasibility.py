"""Exact continuous moment-capacity certificates for Faster box fixtures.

This checker deliberately supports only closed, triangulated axis-aligned boxes
and a fixed shell made from their nested outer/inner boundaries. Unsupported
geometry returns UNKNOWN. OBJ decimal coordinates and JSON numbers are treated
as exact rationals; no optimizer result is used as proof.
"""
from __future__ import annotations

from collections import defaultdict, deque
from decimal import Decimal
from fractions import Fraction as Q
import hashlib
import argparse
import json
from pathlib import Path
import sys
from typing import Any

from .feasibility import (ContinuousMomentCapacityCertificate, ContinuousMomentSymmetryWitness,
                          FeasibilityReport)

_SCOPE = 'CONTINUOUS_GEOMETRIC_MOMENT_CAPACITY'


def _q(value: Any) -> Q:
    if isinstance(value, bool):
        raise ValueError('boolean is not a numeric geometry value')
    if isinstance(value, Q):
        return value
    if isinstance(value, Decimal):
        return Q(value)
    if type(value) is int:
        return Q(value)
    if type(value) is float:
        return Q(*value.as_integer_ratio())
    raise ValueError('unsupported exact numeric value')


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_hash(value: Any) -> str:
    def encode(item):
        if isinstance(item, Q):
            return [item.numerator, item.denominator]
        if isinstance(item, dict):
            return {key: encode(part) for key, part in item.items()}
        if isinstance(item, (tuple, list)):
            return [encode(part) for part in item]
        return item
    encoded = json.dumps(encode(value), sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def _read_obj(path: Path) -> tuple[tuple[Q, Q, Q], tuple[tuple[int, int, int], ...]]:
    vertices: list[tuple[Q, Q, Q]] = []
    faces: list[tuple[int, int, int]] = []
    for line_number, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
        fields = line.split('#', 1)[0].split()
        if not fields:
            continue
        if fields[0] == 'v':
            if len(fields) < 4:
                raise ValueError(f'{path.name}:{line_number}: malformed vertex')
            vertices.append(tuple(Q(field) for field in fields[1:4]))
        elif fields[0] == 'f':
            if len(fields) != 4:
                raise ValueError(f'{path.name}:{line_number}: only triangle faces are supported')
            indices = []
            for field in fields[1:]:
                index = int(field.split('/', 1)[0])
                if index < 0:
                    index = len(vertices) + index
                else:
                    index -= 1
                if index < 0 or index >= len(vertices):
                    raise ValueError(f'{path.name}:{line_number}: face index is out of range')
                indices.append(index)
            faces.append(tuple(indices))
    if not vertices or not faces:
        raise ValueError(f'{path.name}: empty triangle mesh')
    _validate_oriented_closed_mesh(vertices, faces, path.name)
    return tuple(vertices), tuple(faces)


def _validate_oriented_closed_mesh(vertices, faces, label: str) -> None:
    edge_uses: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    face_adjacency = [set() for _ in faces]
    face_edges: dict[tuple[int, int], list[int]] = defaultdict(list)
    seen = set()
    for face_index, face in enumerate(faces):
        if len(set(face)) != 3 or any(index < 0 or index >= len(vertices) for index in face):
            raise ValueError(f'{label}: invalid face {face_index}')
        key = tuple(sorted(face))
        if key in seen:
            raise ValueError(f'{label}: duplicate face {face_index}')
        seen.add(key)
        a, b, c = (vertices[index] for index in face)
        cross = ((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]),
                 (b[2]-a[2])*(c[0]-a[0])-(b[0]-a[0])*(c[2]-a[2]),
                 (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))
        if not any(cross):
            raise ValueError(f'{label}: degenerate face {face_index}')
        for start, end in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge = tuple(sorted((start, end)))
            edge_uses[edge].append((start, end))
            face_edges[edge].append(face_index)
    for edge, uses in edge_uses.items():
        if len(uses) != 2 or uses[0] != (uses[1][1], uses[1][0]):
            raise ValueError(f'{label}: edge {edge} is not a consistently oriented 2-manifold edge')
        first, second = face_edges[edge]
        face_adjacency[first].add(second)
        face_adjacency[second].add(first)
    for vertex_index in range(len(vertices)):
        incident = {i for i, face in enumerate(faces) if vertex_index in face}
        if not incident:
            raise ValueError(f'{label}: unreferenced vertex {vertex_index}')
        pending = [next(iter(incident))]
        reached = set()
        while pending:
            current = pending.pop()
            if current in reached:
                continue
            reached.add(current)
            pending.extend(face_adjacency[current] & incident)
        if reached != incident:
            raise ValueError(f'{label}: vertex {vertex_index} has multiple surface fans')


def _face_components(faces):
    edge_to_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, face in enumerate(faces):
        for edge in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge_to_faces[tuple(sorted(edge))].append(index)
    adjacency = [set() for _ in faces]
    for sharing in edge_to_faces.values():
        if len(sharing) == 2:
            a, b = sharing
            adjacency[a].add(b)
            adjacency[b].add(a)
    components = []
    unseen = set(range(len(faces)))
    while unseen:
        root = min(unseen)
        pending = [root]
        component = set()
        while pending:
            item = pending.pop()
            if item not in unseen:
                continue
            unseen.remove(item)
            component.add(item)
            pending.extend(adjacency[item] & unseen)
        components.append(tuple(sorted(component)))
    return tuple(components)


def _moments(vertices, faces, face_indices=None):
    selected = range(len(faces)) if face_indices is None else face_indices
    volume = Q(0)
    first = [Q(0), Q(0), Q(0)]
    for face_index in selected:
        a, b, c = (vertices[index] for index in faces[face_index])
        cross_bc = (b[1]*c[2]-b[2]*c[1], b[2]*c[0]-b[0]*c[2], b[0]*c[1]-b[1]*c[0])
        signed_volume = sum((a[k]*cross_bc[k] for k in range(3)), Q(0)) / 6
        volume += signed_volume
        for k in range(3):
            first[k] += signed_volume * (a[k]+b[k]+c[k]) / 4
    return volume, tuple(first)


def _validate_box(vertices, faces, expected_sign: int, label: str):
    if len(vertices) != 8 or len(faces) != 12:
        raise ValueError(f'{label}: expected an 8-vertex, 12-triangle box boundary')
    lower = tuple(min(vertex[k] for vertex in vertices) for k in range(3))
    upper = tuple(max(vertex[k] for vertex in vertices) for k in range(3))
    if any(lower[k] >= upper[k] for k in range(3)):
        raise ValueError(f'{label}: degenerate box extent')
    corners = {(x, y, z) for x in (lower[0], upper[0])
                         for y in (lower[1], upper[1])
                         for z in (lower[2], upper[2])}
    if set(vertices) != corners:
        raise ValueError(f'{label}: vertices are not the eight box corners')
    side_counts = defaultdict(list)
    for face in faces:
        points = [vertices[index] for index in face]
        fixed_axes = [axis for axis in range(3)
                      if points[0][axis] == points[1][axis] == points[2][axis]
                      and points[0][axis] in (lower[axis], upper[axis])]
        if len(fixed_axes) != 1:
            raise ValueError(f'{label}: face is not on exactly one box plane')
        axis = fixed_axes[0]
        side_counts[(axis, points[0][axis])].append(points)
    if len(side_counts) != 6 or any(len(group) != 2 for group in side_counts.values()):
        raise ValueError(f'{label}: every box side must consist of two triangles')
    for (axis, _), triangles in side_counts.items():
        transverse = [index for index in range(3) if index != axis]
        full_area = (upper[transverse[0]]-lower[transverse[0]]) * (upper[transverse[1]]-lower[transverse[1]])
        for points in triangles:
            u = [point[transverse[0]] for point in points]
            v = [point[transverse[1]] for point in points]
            twice_area = abs((u[1]-u[0])*(v[2]-v[0])-(u[2]-u[0])*(v[1]-v[0]))
            if twice_area != full_area:
                raise ValueError(f'{label}: triangle does not cover half of its box side')
    volume, first = _moments(vertices, faces)
    expected_volume = expected_sign * (upper[0]-lower[0])*(upper[1]-lower[1])*(upper[2]-lower[2])
    center = tuple((lower[k]+upper[k])/2 for k in range(3))
    if volume != expected_volume or first != tuple(expected_volume*center[k] for k in range(3)):
        raise ValueError(f'{label}: orientation or exact box moments do not match')
    return lower, upper, volume, first


def _mesh_asset(manifest_path: Path, field: str, manifest: dict) -> Path:
    raw = manifest[field]
    path = Path(raw)
    path = path if path.is_absolute() else manifest_path.parent / path
    return path.resolve(strict=True)


def _region_asset(manifest_path: Path, region: dict) -> Path:
    path = Path(region['mesh'])
    return (path if path.is_absolute() else manifest_path.parent / path).resolve(strict=True)


def verify_axis_aligned_yoyo_moment_capacity(manifest_path: str | Path, *,
                                              expected_input_hash: str | None = None) -> FeasibilityReport:
    """Certify a negative first-moment upper bound for a nested box fixture.

    The bound permits every measurable retained subset of the editable box,
    which is a superset of the quadric-hollowing design space. Other modes,
    axes, meshes, or material layouts are reported UNKNOWN.
    """
    source = Path(manifest_path).resolve()
    empty = ContinuousMomentCapacityCertificate()
    try:
        manifest_bytes = source.read_bytes()
        manifest = json.loads(manifest_bytes, parse_float=Decimal, parse_int=int)
        if manifest.get('schema') != 'faster-cherchi-problem-manifest-v1':
            raise ValueError('unsupported problem manifest schema')
        if manifest.get('target', {}).get('mode') != 'YO_YO':
            raise ValueError('certificate only applies to YO_YO center-of-mass constraints')
        target = tuple(_q(value) for value in manifest['target']['origin'])
        axis_vector = tuple(_q(value) for value in manifest['target']['axis'])
        if len(target) != 3 or len(axis_vector) != 3:
            raise ValueError('target frame must be three-dimensional')
        nonzero = [(index, value) for index, value in enumerate(axis_vector) if value]
        if len(nonzero) != 1 or abs(nonzero[0][1]) != 1:
            raise ValueError('only a cardinal unit axis is supported')
        axis, sign_q = nonzero[0]
        sign = int(sign_q)
        if len(manifest.get('design_regions', ())) != 1 or len(manifest.get('fixed_regions', ())) != 1:
            raise ValueError('certificate requires one editable region and one fixed shell')
        design_region = manifest['design_regions'][0]
        fixed_region = manifest['fixed_regions'][0]
        editable_density = _q(design_region['density'])
        fixed_density = _q(fixed_region['density'])
        if editable_density <= 0 or fixed_density <= 0:
            raise ValueError('material densities must be strictly positive')
        outer_path = _mesh_asset(source, 'outer_mesh', manifest)
        design_path = _region_asset(source, design_region)
        fixed_path = _region_asset(source, fixed_region)
        outer_vertices, outer_faces = _read_obj(outer_path)
        design_vertices, design_faces = _read_obj(design_path)
        fixed_vertices, fixed_faces = _read_obj(fixed_path)
        outer_lo, outer_hi, outer_volume, _ = _validate_box(outer_vertices, outer_faces, 1, 'outer mesh')
        design_lo, design_hi, design_volume, _ = _validate_box(design_vertices, design_faces, 1, 'editable box')
        if not all(outer_lo[k] < design_lo[k] < design_hi[k] < outer_hi[k] for k in range(3)):
            raise ValueError('editable box is not strictly inside the outer box')
        components = _face_components(fixed_faces)
        if len(components) != 2 or any(len(component) != 12 for component in components):
            raise ValueError('fixed shell must contain exactly the outer and inner box components')
        component_boxes = []
        for component in components:
            used = sorted({index for face_index in component for index in fixed_faces[face_index]})
            local_vertices = tuple(fixed_vertices[index] for index in used)
            index_map = {old: new for new, old in enumerate(used)}
            local_faces = tuple(tuple(index_map[index] for index in fixed_faces[face_index]) for face_index in component)
            lower = tuple(min(vertex[k] for vertex in local_vertices) for k in range(3))
            sign_expected = 1 if lower == outer_lo else -1 if lower == design_lo else 0
            if sign_expected == 0:
                raise ValueError('fixed shell component does not match either nested box')
            component_boxes.append(_validate_box(local_vertices, local_faces, sign_expected, 'fixed shell component'))
        if {item[0] for item in component_boxes} != {outer_lo, design_lo}:
            raise ValueError('fixed shell components do not match outer and editable box bounds')
        shell_volume, shell_first = _moments(fixed_vertices, fixed_faces)
        if shell_volume != outer_volume-design_volume:
            raise ValueError('fixed shell volume is not outer-box minus editable-box volume')
        expected_shell_first = tuple(outer_first-design_first for outer_first, design_first in zip(
            _moments(outer_vertices, outer_faces)[1], _moments(design_vertices, design_faces)[1]))
        if shell_first != expected_shell_first:
            raise ValueError('fixed shell first moments are not outer-box minus editable-box moments')

        relative_fixed = fixed_density * (shell_first[axis]-target[axis]*shell_volume) * sign
        transverse_area = Q(1)
        for index in range(3):
            if index != axis:
                transverse_area *= design_hi[index]-design_lo[index]
        directed_lo, directed_hi = sorted((design_lo[axis]*sign, design_hi[axis]*sign))
        directed_target = target[axis]*sign
        # Integral of max(x-target, 0), not half the square of the retained
        # width. When the entire box lies beyond the target, the lower-limit
        # contribution is nonzero and a width-only formula underestimates the
        # available first moment (potentially creating a false certificate).
        positive_hi = max(Q(0), directed_hi-directed_target)
        positive_lo = max(Q(0), directed_lo-directed_target)
        capacity = editable_density*transverse_area*(positive_hi**2-positive_lo**2)/2
        upper_bound = relative_fixed+capacity

        assets = ((str(source), hashlib.sha256(manifest_bytes).hexdigest()),
                  (str(outer_path), _sha256(outer_path)),
                  (str(design_path), _sha256(design_path)),
                  (str(fixed_path), _sha256(fixed_path)))
        input_hash = _canonical_hash({'schema': 'continuous-box-capacity-input-v1',
                                      'assets': assets, 'target': target, 'axis': axis_vector,
                                      'mode': manifest['target']['mode'],
                                      'fixed_density': fixed_density,
                                      'editable_density': editable_density})
        if expected_input_hash is not None and expected_input_hash != input_hash:
            return FeasibilityReport(reason='INPUT_IDENTITY_MISMATCH', input_hash=input_hash,
                                     system_hash=input_hash, checker_sha256=_sha256(Path(__file__)),
                                     scope=_SCOPE)
        certificate = ContinuousMomentCapacityCertificate(
            constrained_axis=axis, direction_sign=sign, target_coordinate=target[axis],
            fixed_density=fixed_density, fixed_volume=shell_volume,
            fixed_first_moment=shell_first[axis], fixed_relative_first_moment=relative_fixed,
            editable_density=editable_density, box_lower=design_lo, box_upper=design_hi,
            maximum_editable_positive_moment=capacity, total_moment_upper_bound=upper_bound,
            input_sha256=input_hash, asset_sha256=tuple((str(path), digest) for path, digest in assets),
            assumptions=('YO_YO requires zero first moment along the selected axis',
                         'all editable material is a measurable subset of the certified box',
                         'editable density is constant and strictly positive',
                         'the exact fixed mesh equals the outer box minus the editable box'))
        reason = 'EXACT_GEOMETRIC_MOMENT_CAPACITY_CERTIFICATE' if certificate.valid else 'NO_STRICT_NEGATIVE_CAPACITY_BOUND'
        return FeasibilityReport(reason=reason, input_hash=input_hash, system_hash=input_hash,
                                 continuous_certificate=certificate, checker_sha256=_sha256(Path(__file__)),
                                 scope=_SCOPE)
    except (OSError, UnicodeError, json.JSONDecodeError, KeyError, TypeError, ValueError, ZeroDivisionError) as error:
        return FeasibilityReport(reason='INVALID_OR_UNSUPPORTED_INPUT: ' + str(error),
                                 checker_sha256=_sha256(Path(__file__)), scope=_SCOPE)


def verify_centered_axis_ellipsoid_witness(manifest_path: str | Path, *,
                                            expected_input_hash: str | None = None) -> FeasibilityReport:
    """Prove one nonempty analytic retained-material witness for a centered box fixture.

    The ellipsoid is a mathematical existence witness only; its parameters are
    never passed to the optimizer as an initial solution.
    """
    source = Path(manifest_path).resolve()
    scope = 'CONTINUOUS_GEOMETRIC_SYMMETRY_WITNESS'
    try:
        manifest_bytes = source.read_bytes()
        manifest = json.loads(manifest_bytes, parse_float=Decimal, parse_int=int)
        if manifest.get('schema') != 'faster-cherchi-problem-manifest-v1' or manifest.get('target', {}).get('mode') != 'YO_YO':
            raise ValueError('witness applies only to a Faster YO_YO manifest')
        target = tuple(_q(value) for value in manifest['target']['origin'])
        axis_vector = tuple(_q(value) for value in manifest['target']['axis'])
        if len(target) != 3 or len(axis_vector) != 3:
            raise ValueError('target frame must be three-dimensional')
        nonzero = [(index, value) for index, value in enumerate(axis_vector) if value]
        if len(nonzero) != 1 or abs(nonzero[0][1]) != 1:
            raise ValueError('witness requires a cardinal unit spindle axis')
        axis, sign_q = nonzero[0]
        sign = int(sign_q)
        if len(manifest.get('design_regions', ())) != 1 or len(manifest.get('fixed_regions', ())) != 1:
            raise ValueError('witness requires one editable box and one fixed shell')
        design_region, fixed_region = manifest['design_regions'][0], manifest['fixed_regions'][0]
        editable_density, fixed_density = _q(design_region['density']), _q(fixed_region['density'])
        if editable_density <= 0 or fixed_density <= 0:
            raise ValueError('material densities must be positive')
        outer_path = _mesh_asset(source, 'outer_mesh', manifest)
        design_path = _region_asset(source, design_region)
        fixed_path = _region_asset(source, fixed_region)
        outer_vertices, outer_faces = _read_obj(outer_path)
        design_vertices, design_faces = _read_obj(design_path)
        fixed_vertices, fixed_faces = _read_obj(fixed_path)
        outer_lo, outer_hi, outer_volume, _ = _validate_box(outer_vertices, outer_faces, 1, 'outer mesh')
        design_lo, design_hi, design_volume, _ = _validate_box(design_vertices, design_faces, 1, 'editable box')
        if not all(outer_lo[k] < design_lo[k] < design_hi[k] < outer_hi[k] for k in range(3)):
            raise ValueError('editable box is not strictly inside the outer box')
        target_center = tuple((design_lo[k]+design_hi[k])/2 for k in range(3))
        if target != target_center or tuple((outer_lo[k]+outer_hi[k])/2 for k in range(3)) != target:
            raise ValueError('fixture boxes and spindle origin are not centrally aligned')
        if len({outer_hi[k]-outer_lo[k] for k in range(3)}) != 1 or len({design_hi[k]-design_lo[k] for k in range(3)}) != 1:
            raise ValueError('fixed shell must be formed by concentric cubes for isotropic fixed inertia')
        components = _face_components(fixed_faces)
        if len(components) != 2 or any(len(component) != 12 for component in components):
            raise ValueError('fixed shell must have outer and inner cube components')
        component_bounds = []
        for component in components:
            used = sorted({index for face_index in component for index in fixed_faces[face_index]})
            local_vertices = tuple(fixed_vertices[index] for index in used)
            index_map = {old: new for new, old in enumerate(used)}
            local_faces = tuple(tuple(index_map[index] for index in fixed_faces[face_index]) for face_index in component)
            lower = tuple(min(vertex[k] for vertex in local_vertices) for k in range(3))
            expected_sign = 1 if lower == outer_lo else -1 if lower == design_lo else 0
            if expected_sign == 0:
                raise ValueError('fixed shell component is not one of the nested cubes')
            box_lo, box_hi, _, _ = _validate_box(local_vertices, local_faces, expected_sign, 'fixed shell component')
            component_bounds.append((box_lo, box_hi))
        if set(component_bounds) != {(outer_lo, outer_hi), (design_lo, design_hi)}:
            raise ValueError('fixed shell does not equal the nested cube boundaries')
        shell_volume, shell_first = _moments(fixed_vertices, fixed_faces)
        if shell_volume != outer_volume-design_volume or any(shell_first[k] != shell_volume*target[k] for k in range(3)):
            raise ValueError('fixed shell is not centered with the validated nested-box moments')
        fixed_second_diagonal = tuple(fixed_density * (
            outer_volume*(outer_hi[k]-outer_lo[k])**2/12
            - design_volume*(design_hi[k]-design_lo[k])**2/12) for k in range(3))
        fixed_second = tuple(tuple(fixed_second_diagonal[i] if i == j else Q(0) for j in range(3))
                             for i in range(3))
        if not fixed_second_diagonal[0] == fixed_second_diagonal[1] == fixed_second_diagonal[2] > 0:
            raise ValueError('fixed shell inertia is not isotropic and positive')
        # Frozen independent witness: semiaxes (1/2, 3/4, 3/4), centered at
        # the target and wholly inside the editable unit cube.
        semiaxes = (Q(1, 2), Q(3, 4), Q(3, 4))
        if any(semiaxes[k] > min(target[k]-design_lo[k], design_hi[k]-target[k]) for k in range(3)):
            raise ValueError('registered witness ellipsoid does not fit inside the design box')
        axis_inertia_factor = sum((semiaxes[k]**2 for k in range(3) if k != axis), Q(0))
        other_factors = tuple(sum((semiaxes[k]**2 for k in range(3) if k != other), Q(0))
                              for other in range(3) if other != axis)
        inertia_gap = min(axis_inertia_factor-value for value in other_factors)/5
        assets = ((str(source), hashlib.sha256(manifest_bytes).hexdigest()),
                  (str(outer_path), _sha256(outer_path)),
                  (str(design_path), _sha256(design_path)),
                  (str(fixed_path), _sha256(fixed_path)))
        input_hash = _canonical_hash({'schema': 'centered-ellipsoid-witness-input-v1',
                                      'assets': assets, 'target': target, 'axis': axis_vector,
                                      'fixed_density': fixed_density, 'editable_density': editable_density,
                                      'semiaxes': semiaxes})
        if expected_input_hash is not None and expected_input_hash != input_hash:
            return FeasibilityReport(reason='INPUT_IDENTITY_MISMATCH', input_hash=input_hash,
                                     system_hash=input_hash, checker_sha256=_sha256(Path(__file__)), scope=scope)
        witness = ContinuousMomentSymmetryWitness(
            constrained_axis=axis, target=target, design_lower=design_lo, design_upper=design_hi,
            ellipsoid_semiaxes=semiaxes, ellipsoid_density=editable_density,
            fixed_density_weighted_second_moment=fixed_second, inertia_gap_per_mass=inertia_gap,
            input_sha256=input_hash, asset_sha256=tuple((str(path), digest) for path, digest in assets))
        reason = 'EXACT_ANALYTIC_SYMMETRY_WITNESS' if witness.valid else 'INVALID_ANALYTIC_SYMMETRY_WITNESS'
        return FeasibilityReport(reason=reason, input_hash=input_hash, system_hash=input_hash,
                                 continuous_witness=witness, checker_sha256=_sha256(Path(__file__)), scope=scope)
    except (OSError, UnicodeError, json.JSONDecodeError, KeyError, TypeError, ValueError, ZeroDivisionError) as error:
        return FeasibilityReport(reason='INVALID_OR_UNSUPPORTED_INPUT: ' + str(error),
                                 checker_sha256=_sha256(Path(__file__)), scope=scope)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--centered-witness', action='store_true',
                        help='check the registered centered-ellipsoid existence witness')
    parser.add_argument('--expected-input-hash')
    parser.add_argument('--output', type=Path, help='write the report to a new or existing parent directory')
    args = parser.parse_args(argv)
    checker = (verify_centered_axis_ellipsoid_witness if args.centered_witness
               else verify_axis_aligned_yoyo_moment_capacity)
    report = checker(args.manifest, expected_input_hash=args.expected_input_hash)
    text = json.dumps(report.to_dict(), sort_keys=True, indent=2, allow_nan=False) + '\n'
    if args.output is None:
        sys.stdout.write(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding='utf-8')
        sys.stdout.write(f'{report.status}: {args.output}\n')
    return 0 if report.status in ('INFEASIBLE_CERTIFIED', 'FEASIBLE_WITNESS') else 2


if __name__ == '__main__':
    raise SystemExit(main())
