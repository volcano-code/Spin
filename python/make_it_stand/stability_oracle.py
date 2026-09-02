"""Deterministic fixed-support quasi-static gravity/table tilt oracle."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class QuasiStaticTiltResult:
    protocol: str
    contact_policy: str
    tilt_degrees: float
    azimuth_degrees: float
    inside_support: bool
    signed_margin: float
    nominal_signed_margin: float
    center_height: float
    critical_tilt_degrees: float
    projected_center: tuple[float, float, float]
    tilted_gravity: tuple[float, float, float]


def _unit(value: Sequence[float], name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (3,) or not np.isfinite(result).all():
        raise ValueError(f"{name} must be a finite 3-vector")
    length = float(np.linalg.norm(result))
    if length <= 1.0e-15:
        raise ValueError(f"{name} must be nonzero")
    return result / length


def _basis(up: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    reference = np.asarray((1.0, 0.0, 0.0))
    if abs(float(up @ reference)) > 0.9:
        reference = np.asarray((0.0, 0.0, 1.0))
    first = np.cross(up, reference)
    first /= np.linalg.norm(first)
    second = np.cross(up, first)
    return first, second


def _convex_hull(points: np.ndarray) -> np.ndarray:
    unique = sorted(set((float(point[0]), float(point[1])) for point in points))
    if len(unique) < 3:
        raise ValueError("support polygon must span a nonzero 2D area")

    def cross(origin: tuple[float, float], a: tuple[float, float], b: tuple[float, float]) -> float:
        return (a[0] - origin[0]) * (b[1] - origin[1]) - (
            a[1] - origin[1]
        ) * (b[0] - origin[0])

    lower: list[tuple[float, float]] = []
    for point in unique:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0.0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[float, float]] = []
    for point in reversed(unique):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0.0:
            upper.pop()
        upper.append(point)
    hull = np.asarray(lower[:-1] + upper[:-1], dtype=np.float64)
    if hull.shape[0] < 3:
        raise ValueError("support polygon convex hull is degenerate")
    return hull


def _signed_margin(hull: np.ndarray, point: np.ndarray) -> float:
    margins = []
    for index in range(hull.shape[0]):
        first = hull[index]
        second = hull[(index + 1) % hull.shape[0]]
        edge = second - first
        length = float(np.linalg.norm(edge))
        if length <= 1.0e-15:
            continue
        offset = point - first
        margins.append(
            float(edge[0] * offset[1] - edge[1] * offset[0]) / length
        )
    if not margins:
        raise ValueError("support polygon has no valid edges")
    return min(margins)


def evaluate_quasi_static_tilt(
    support_polygon: Sequence[Sequence[float]],
    center_of_mass: Sequence[float],
    gravity: Sequence[float],
    *,
    tilt_degrees: float,
    azimuth_degrees: float,
) -> QuasiStaticTiltResult:
    """Intersect the tilted gravity ray with one fixed projected support hull."""

    polygon = np.asarray(support_polygon, dtype=np.float64)
    center = np.asarray(center_of_mass, dtype=np.float64)
    if (
        polygon.ndim != 2
        or polygon.shape[1] != 3
        or polygon.shape[0] < 3
        or not np.isfinite(polygon).all()
    ):
        raise ValueError("support_polygon must contain at least three finite 3D points")
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("center_of_mass must be a finite 3-vector")
    angle = float(tilt_degrees)
    azimuth = float(azimuth_degrees)
    if not math.isfinite(angle) or angle < 0.0 or angle >= 90.0:
        raise ValueError("tilt_degrees must be finite in [0, 90)")
    if not math.isfinite(azimuth):
        raise ValueError("azimuth_degrees must be finite")

    down = _unit(gravity, "gravity")
    up = -down
    first, second = _basis(up)
    plane_origin = polygon.mean(axis=0)
    coordinates = np.column_stack((polygon @ first, polygon @ second))
    hull = _convex_hull(coordinates)
    center_height = float((center - plane_origin) @ up)
    if center_height <= 0.0:
        raise ValueError("center_of_mass must lie above the support plane")

    radians = math.radians(angle)
    azimuth_radians = math.radians(azimuth)
    horizontal = (
        math.cos(azimuth_radians) * first
        + math.sin(azimuth_radians) * second
    )
    tilted_down = math.cos(radians) * down + math.sin(radians) * horizontal
    denominator = float(tilted_down @ up)
    if denominator >= -1.0e-15:
        raise ValueError("tilted gravity does not intersect the support plane")
    distance = -center_height / denominator
    projected = center + distance * tilted_down
    projected_2d = np.asarray((projected @ first, projected @ second))
    center_2d = np.asarray((center @ first, center @ second))
    margin = _signed_margin(hull, projected_2d)
    nominal_margin = _signed_margin(hull, center_2d)
    critical = math.degrees(math.atan2(max(0.0, nominal_margin), center_height))
    tolerance = 128.0 * np.finfo(np.float64).eps * max(
        1.0,
        float(np.linalg.norm(hull, ord=np.inf)),
    )
    return QuasiStaticTiltResult(
        protocol="quasi_static_table_tilt",
        contact_policy="fixed_projected_support_polygon",
        tilt_degrees=angle,
        azimuth_degrees=azimuth % 360.0,
        inside_support=margin >= -tolerance,
        signed_margin=margin,
        nominal_signed_margin=nominal_margin,
        center_height=center_height,
        critical_tilt_degrees=critical,
        projected_center=tuple(float(item) for item in projected),
        tilted_gravity=tuple(float(item) for item in tilted_down),
    )


def evaluate_tilt_azimuths(
    support_polygon: Sequence[Sequence[float]],
    center_of_mass: Sequence[float],
    gravity: Sequence[float],
    *,
    tilt_degrees: float,
    azimuths_degrees: Iterable[float],
) -> tuple[QuasiStaticTiltResult, ...]:
    return tuple(
        evaluate_quasi_static_tilt(
            support_polygon,
            center_of_mass,
            gravity,
            tilt_degrees=tilt_degrees,
            azimuth_degrees=float(azimuth),
        )
        for azimuth in azimuths_degrees
    )


__all__ = [
    "QuasiStaticTiltResult",
    "evaluate_quasi_static_tilt",
    "evaluate_tilt_azimuths",
]
