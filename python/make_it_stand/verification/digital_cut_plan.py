"""Exact-on-grid verification for explicitly supplied digital cut plans.

Spin-It (2014) used manually selected physical cuts.  This module deliberately
does not search for or synthesize those cuts.  It verifies what a supplied cut
set does on one *discrete* material model: face-connected void reachability,
pre-existing cavity exposure, remaining material components, and fixed-material
exclusions.  It makes no continuous-mesh, CAD, wall-thickness, or paper-case
claim.
"""
from __future__ import annotations

from collections import deque
from collections.abc import Iterable, Sequence
import math
from typing import Any


Cell = tuple[int, int, int]
_NEIGHBORS: tuple[Cell, ...] = (
    (-1, 0, 0), (1, 0, 0),
    (0, -1, 0), (0, 1, 0),
    (0, 0, -1), (0, 0, 1),
)


def _dimensions(value: Sequence[int]) -> Cell:
    if len(value) != 3 or any(type(axis) is not int or axis <= 0 for axis in value):
        raise ValueError("dimensions must contain exactly three positive integers")
    return int(value[0]), int(value[1]), int(value[2])


def _finite_vector3(value: Sequence[float], *, name: str) -> tuple[float, float, float]:
    if len(value) != 3:
        raise ValueError(name + " must contain exactly three finite values")
    try:
        vector = float(value[0]), float(value[1]), float(value[2])
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(name + " must contain exactly three finite values") from error
    if not all(math.isfinite(component) for component in vector):
        raise ValueError(name + " must contain exactly three finite values")
    return vector


def _cut_plane_geometry(
    grid_origin: Sequence[float],
    voxel_spacing: Sequence[float],
    plane_point: Sequence[float],
    plane_normal: Sequence[float],
    kerf_width: float,
) -> tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
    float,
]:
    origin = _finite_vector3(grid_origin, name="grid_origin")
    spacing = _finite_vector3(voxel_spacing, name="voxel_spacing")
    point = _finite_vector3(plane_point, name="plane_point")
    normal = _finite_vector3(plane_normal, name="plane_normal")
    if any(axis <= 0.0 for axis in spacing):
        raise ValueError("voxel_spacing must be strictly positive")
    normal_length = math.sqrt(sum(axis * axis for axis in normal))
    if not math.isfinite(normal_length) or normal_length == 0.0:
        raise ValueError("plane_normal must be finite and nonzero")
    normalized = tuple(axis / normal_length for axis in normal)
    try:
        kerf = float(kerf_width)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("kerf_width must be finite and strictly positive") from error
    if not math.isfinite(kerf) or kerf <= 0.0:
        raise ValueError("kerf_width must be finite and strictly positive")
    return origin, spacing, point, normalized, kerf


def rasterize_discrete_cut_plane(
    dimensions: Sequence[int],
    *,
    grid_origin: Sequence[float],
    voxel_spacing: Sequence[float],
    plane_point: Sequence[float],
    plane_normal: Sequence[float],
    kerf_width: float,
) -> list[Cell]:
    """Return cells whose boxes intersect a supplied infinite kerf slab.

    ``grid_origin`` is the lower corner of cell (0, 0, 0). The supplied normal
    is normalized internally. A cell is selected iff its closed axis-aligned
    box intersects the slab of ``kerf_width`` centered on the given plane;
    this is an exact plane/box interval test, not center sampling. The plane is
    infinite across the finite grid, so this function does not model a finite
    saw blade, toolpath, or automatically selected paper cut.
    """
    grid = _dimensions(dimensions)
    origin, spacing, point, normal, kerf = _cut_plane_geometry(
        grid_origin, voxel_spacing, plane_point, plane_normal, kerf_width
    )
    half_kerf = 0.5 * kerf
    support_radius = 0.5 * sum(
        abs(normal[axis]) * spacing[axis] for axis in range(3)
    )
    selected: list[Cell] = []
    for x in range(grid[0]):
        for y in range(grid[1]):
            for z in range(grid[2]):
                cell = (x, y, z)
                center = tuple(
                    origin[axis] + (cell[axis] + 0.5) * spacing[axis]
                    for axis in range(3)
                )
                signed_distance = sum(
                    normal[axis] * (center[axis] - point[axis])
                    for axis in range(3)
                )
                tolerance_scale = max(
                    1.0, abs(signed_distance), half_kerf, support_radius,
                    *(abs(coordinate) for coordinate in center),
                    *(abs(coordinate) for coordinate in point),
                )
                tolerance = 16.0 * math.ulp(tolerance_scale)
                if abs(signed_distance) <= half_kerf + support_radius + tolerance:
                    selected.append(cell)
    return selected


def _cell(value: object, dimensions: Cell, *, name: str) -> Cell:
    if not isinstance(value, tuple | list) or len(value) != 3 or any(type(axis) is not int for axis in value):
        raise ValueError(name + " must contain integer xyz cells")
    cell = int(value[0]), int(value[1]), int(value[2])
    if any(axis < 0 or axis >= limit for axis, limit in zip(cell, dimensions)):
        raise ValueError(name + " contains a cell outside the grid")
    return cell


def _cells(values: Iterable[object], dimensions: Cell, *, name: str) -> set[Cell]:
    try:
        return {_cell(value, dimensions, name=name) for value in values}
    except TypeError as error:
        raise ValueError(name + " must be an iterable of cells") from error


def _all_cells(dimensions: Cell) -> set[Cell]:
    return {
        (x, y, z)
        for x in range(dimensions[0])
        for y in range(dimensions[1])
        for z in range(dimensions[2])
    }


def _neighbors(cell: Cell, dimensions: Cell) -> tuple[Cell, ...]:
    result: list[Cell] = []
    for offset in _NEIGHBORS:
        other = tuple(cell[axis] + offset[axis] for axis in range(3))
        if all(0 <= coordinate < dimensions[axis] for axis, coordinate in enumerate(other)):
            result.append(other)  # type: ignore[arg-type]
    return tuple(result)


def _components(cells: set[Cell], dimensions: Cell) -> list[set[Cell]]:
    remaining = set(cells)
    components: list[set[Cell]] = []
    while remaining:
        first = min(remaining)
        component: set[Cell] = {first}
        queue: deque[Cell] = deque((first,))
        remaining.remove(first)
        while queue:
            cell = queue.popleft()
            for other in _neighbors(cell, dimensions):
                if other in remaining:
                    remaining.remove(other)
                    component.add(other)
                    queue.append(other)
        components.append(component)
    return sorted(components, key=lambda component: min(component))


def _boundary_cells(cells: set[Cell], dimensions: Cell) -> set[Cell]:
    return {
        cell for cell in cells
        if any(cell[axis] in (0, dimensions[axis] - 1) for axis in range(3))
    }


def _reachable_from_exterior(void_cells: set[Cell], dimensions: Cell) -> set[Cell]:
    reachable = _boundary_cells(void_cells, dimensions)
    queue: deque[Cell] = deque(sorted(reachable))
    while queue:
        cell = queue.popleft()
        for other in _neighbors(cell, dimensions):
            if other in void_cells and other not in reachable:
                reachable.add(other)
                queue.append(other)
    return reachable


def _listed(cells: Iterable[Cell]) -> list[list[int]]:
    return [list(cell) for cell in sorted(cells)]


def verify_discrete_cut_plan(
    dimensions: Sequence[int],
    material_cells: Iterable[object],
    *,
    cut_cells: Iterable[object],
    fixed_material_cells: Iterable[object] = (),
) -> dict[str, Any]:
    """Audit an explicit voxel cut plan without adding any cut geometry.

    ``cut_cells`` is the actual swept/removal set, so it may include air cells
    on its route to an existing cavity.  Only its intersection with material is
    removed.  Connectivity is strictly 6-neighbor face connectivity: diagonal
    contact never creates a false opening or a false material connection.
    """
    grid = _dimensions(dimensions)
    all_cells = _all_cells(grid)
    material = _cells(material_cells, grid, name="material_cells")
    cuts = _cells(cut_cells, grid, name="cut_cells")
    fixed = _cells(fixed_material_cells, grid, name="fixed_material_cells")
    if not fixed <= material:
        raise ValueError("fixed material must be a subset of material_cells")
    prohibited = fixed & cuts
    if prohibited:
        raise ValueError("cut plan removes fixed material cells: " + repr(_listed(prohibited)))

    void_before = all_cells - material
    exterior_before = _reachable_from_exterior(void_before, grid)
    cavities_before = _components(void_before - exterior_before, grid)

    removed = material & cuts
    material_after = material - removed
    void_after = all_cells - material_after
    exterior_after = _reachable_from_exterior(void_after, grid)
    cavity_records: list[dict[str, Any]] = []
    for index, cavity in enumerate(cavities_before):
        exposed = bool(cavity & exterior_after)
        cavity_records.append({
            "cavity_index": index,
            "cells": _listed(cavity),
            "cell_count": len(cavity),
            "exposed": exposed,
        })
    material_components = _components(material_after, grid)
    post_cut_void_components = _components(void_after - exterior_after, grid)
    exposed = sum(record["exposed"] for record in cavity_records)
    return {
        "schema": "discrete-digital-cut-plan-v1",
        "status": "VERIFIED_DISCRETE_CUT_PLAN",
        "scope": "DISCRETE_VOXEL_CONNECTIVITY_ONLY_NOT_CONTINUOUS_CAD_OR_PAPER_CASE",
        "connectivity": "FACE_6_NEIGHBOR",
        "dimensions": list(grid),
        "input_material_cell_count": len(material),
        "fixed_material_cell_count": len(fixed),
        "planned_cut_cell_count": len(cuts),
        "removed_material_cell_count": len(removed),
        "removed_material_cells": _listed(removed),
        "preexisting_cavity_count": len(cavity_records),
        "exposed_preexisting_cavity_count": exposed,
        "unexposed_preexisting_cavity_count": len(cavity_records) - exposed,
        "cavities": cavity_records,
        "cut_material_component_count": len(material_components),
        "post_cut_enclosed_void_component_count": len(post_cut_void_components),
        "post_cut_exterior_void_cell_count": len(exterior_after),
        "continuous_geometry_status": "UNKNOWN",
        "cut_plane_or_toolpath_status": "EXPLICIT_INPUT_REQUIRED",
        "paper_case_status": "NOT_GRANTED",
    }


def verify_discrete_cut_plane(
    dimensions: Sequence[int],
    material_cells: Iterable[object],
    *,
    grid_origin: Sequence[float],
    voxel_spacing: Sequence[float],
    plane_point: Sequence[float],
    plane_normal: Sequence[float],
    kerf_width: float,
    fixed_material_cells: Iterable[object] = (),
) -> dict[str, Any]:
    """Rasterize and audit one explicitly supplied plane cut.

    This applies a digital cut to a discrete voxel material occupancy and then
    runs the same face-connectivity/fixed-material checks as an explicit cell
    plan. It makes no continuous surface, CAD, or author-case claim.
    """
    grid = _dimensions(dimensions)
    origin, spacing, point, normal, kerf = _cut_plane_geometry(
        grid_origin, voxel_spacing, plane_point, plane_normal, kerf_width
    )
    cut_cells = rasterize_discrete_cut_plane(
        grid,
        grid_origin=origin,
        voxel_spacing=spacing,
        plane_point=point,
        plane_normal=normal,
        kerf_width=kerf,
    )
    material = _cells(material_cells, grid, name="material_cells")
    base = verify_discrete_cut_plan(
        grid,
        material,
        cut_cells=cut_cells,
        fixed_material_cells=fixed_material_cells,
    )
    remaining = material - _cells(cut_cells, grid, name="cut_cells")
    result = dict(base)
    result.update({
        "schema": "discrete-digital-cut-plane-v1",
        "status": "VERIFIED_DISCRETE_CUT_PLANE",
        "scope": "DISCRETE_VOXEL_CUT_ONLY_NOT_CONTINUOUS_CAD_OR_PAPER_CASE",
        "cut_plane": {
            "point": list(point),
            "normalized_normal": list(normal),
            "kerf_width": kerf,
            "grid_origin": list(origin),
            "voxel_spacing": list(spacing),
            "rasterization": "EXACT_CLOSED_CELL_BOX_INTERSECTION_WITH_INFINITE_KERF_SLAB",
        },
        "remaining_material_cell_count": len(remaining),
        "remaining_material_cells": _listed(remaining),
        "continuous_geometry_status": "UNKNOWN",
        "cut_plane_or_toolpath_status": "FINITE_TOOL_GEOMETRY_NOT_MODELED",
        "paper_case_status": "NOT_GRANTED",
    })
    return result
