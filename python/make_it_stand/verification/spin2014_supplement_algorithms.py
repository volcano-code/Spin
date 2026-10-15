"""Independent transcription of Spin-It (2014) Supplement Algorithms 1 and 2.

The reference evaluates the published per-triangle surface polynomials directly
and deliberately does *not* import the production surface-moment implementation.
It uses the project raw-moment order ``(M,LX,LY,LZ,QXX,QYY,QZZ,QYZ,QXZ,QXY)``;
the paper's order is documented below and the conversion is explicit.

This is a local numerical reference, not a claim that its arithmetic sequence
or floating-point results are author-bitwise identical.
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np


ORDER = ("M", "LX", "LY", "LZ", "QXX", "QYY", "QZZ", "QYZ", "QXZ", "QXY")
PAPER_ORDER = ("M", "LX", "LY", "LZ", "QXY", "QYZ", "QXZ", "QXX", "QYY", "QZZ")
ALGORITHM = "spin-it-2014-supplement-algorithms-1-2-surface-polynomial"

# Paper rows [1,x,y,z,xy,yz,xz,x^2,y^2,z^2] to the native raw-moment order.
_PAPER_TO_NATIVE = np.asarray((0, 1, 2, 3, 9, 7, 8, 4, 5, 6), dtype=int)
_FACTORS = np.asarray((1 / 6, 1 / 24, 1 / 24, 1 / 24,
                       1 / 120, 1 / 120, 1 / 120, 1 / 60, 1 / 60, 1 / 60))


def _coordinates(vertices: Iterable[Iterable[float]]) -> np.ndarray:
    result = np.asarray(tuple(vertices), dtype=float)
    if result.ndim != 2 or result.shape[1:] != (3,) or not len(result):
        raise ValueError("vertices must be a nonempty Nx3 numeric array")
    if not np.all(np.isfinite(result)):
        raise ValueError("vertices must be finite")
    return result


def _triangles(faces: Iterable[Iterable[int]], vertex_count: int) -> np.ndarray:
    result = np.asarray(tuple(faces), dtype=int)
    if result.ndim != 2 or result.shape[1:] != (3,) or not len(result):
        raise ValueError("faces must be a nonempty Mx3 integer array")
    if np.any(result < 0) or np.any(result >= vertex_count):
        raise ValueError("face index is outside vertices")
    if any(len(set(face.tolist())) != 3 for face in result):
        raise ValueError("face repeats a vertex")
    return result


def _skew(vector: np.ndarray) -> np.ndarray:
    x, y, z = vector
    return np.asarray(((0.0, -z, y), (z, 0.0, -x), (-y, x, 0.0)))


def _cycle(vector: np.ndarray) -> np.ndarray:
    """Paper's bar operator: (x,y,z) -> (y,z,x)."""
    return np.asarray((vector[1], vector[2], vector[0]))


def triangle_trace(a: Iterable[float], b: Iterable[float], c: Iterable[float]) -> dict[str, object]:
    """Evaluate both published algorithms for one oriented triangle.

    The returned arrays before scaling use the paper's row order.  The final
    arrays use :data:`ORDER`, making any convention conversion inspectable.
    """

    a, b, c = (np.asarray(value, dtype=float) for value in (a, b, c))
    if any(value.shape != (3,) or not np.all(np.isfinite(value)) for value in (a, b, c)):
        raise ValueError("triangle coordinates must be three finite values")
    u, v = b - a, c - a
    normal = np.cross(u, v)
    if not np.all(np.isfinite(normal)) or np.all(normal == 0.0):
        raise ValueError("triangle must be nondegenerate with finite normal")

    h1 = a + b + c
    h2 = a * a + b * (a + b)
    h3 = h2 + c * h1
    h4 = a * a * a + b * h2 + c * h3
    h5, h6, h7 = h3 + a * (h1 + a), h3 + b * (h1 + b), h3 + c * (h1 + c)
    h8 = _cycle(a) * h5 + _cycle(b) * h6 + _cycle(c) * h7
    paper_moments = np.concatenate(((normal[0] * h1[0],), normal * h3, normal * h8, normal * h4))

    dn = (_skew(v) - _skew(u), -_skew(v), _skew(u))
    d3 = (2 * a + b + c, a + 2 * b + c, a + b + 2 * c)
    d5a, d6b, d7c = 6 * a + 2 * b + 2 * c, 2 * a + 6 * b + 2 * c, 2 * a + 2 * b + 6 * c
    d5b, d5c, d6c = 2 * a + 2 * b + c, 2 * a + b + 2 * c, a + 2 * b + 2 * c
    d8 = (
        _cycle(a) * d5a + _cycle(b) * d5b + _cycle(c) * d5c,
        _cycle(a) * d5b + _cycle(b) * d6b + _cycle(c) * d6c,
        _cycle(a) * d5c + _cycle(b) * d6c + _cycle(c) * d7c,
    )
    d4 = (h5, h6, h7)
    paper_jacobian = np.zeros((10, 9), dtype=float)
    for vertex in range(3):
        block = np.zeros((10, 3), dtype=float)
        block[0] = h1[0] * dn[vertex][0]
        block[0, 0] += normal[0]
        for axis in range(3):
            block[1 + axis] = h3[axis] * dn[vertex][axis]
            block[1 + axis, axis] += normal[axis] * d3[vertex][axis]
            block[4 + axis] = h8[axis] * dn[vertex][axis]
            block[4 + axis, axis] += normal[axis] * d8[vertex][axis]
            # [n h]_* in Algorithm 2 has its nonzero entry in the next
            # cyclic coordinate; spelling it out avoids an implicit convention.
            block[4 + axis, (axis + 1) % 3] += normal[axis] * d4[vertex][axis]
            block[7 + axis] = h4[axis] * dn[vertex][axis]
            block[7 + axis, axis] += normal[axis] * d4[vertex][axis]
        paper_jacobian[:, 3 * vertex:3 * vertex + 3] = block

    moments = np.zeros(10, dtype=float)
    jacobian = np.zeros((10, 9), dtype=float)
    for paper_row, native_row in enumerate(_PAPER_TO_NATIVE):
        moments[native_row] = _FACTORS[paper_row] * paper_moments[paper_row]
        jacobian[native_row] = _FACTORS[paper_row] * paper_jacobian[paper_row]
    if not np.all(np.isfinite(moments)) or not np.all(np.isfinite(jacobian)):
        raise OverflowError("published polynomial evaluation overflowed")
    return {
        "algorithm": ALGORITHM,
        "paper_order": PAPER_ORDER,
        "order": ORDER,
        "normal": normal,
        "h": (h1, h2, h3, h4, h5, h6, h7, h8),
        "paper_unscaled_moments": paper_moments,
        "paper_unscaled_jacobian": paper_jacobian,
        "moments": moments,
        "jacobian": jacobian,
    }


def integrate_surface_moments(vertices: Iterable[Iterable[float]], faces: Iterable[Iterable[int]], density: float = 1.0) -> np.ndarray:
    """Algorithm 1, summed over an oriented triangle surface and density-scaled."""

    points = _coordinates(vertices)
    triangles = _triangles(faces, len(points))
    if isinstance(density, bool) or not isinstance(density, (int, float)) or not math.isfinite(density):
        raise ValueError("density must be finite")
    result = np.zeros(10, dtype=float)
    for face in triangles:
        result += triangle_trace(*(points[index] for index in face))["moments"]
    result *= float(density)
    if not np.all(np.isfinite(result)):
        raise OverflowError("published moment integration overflowed")
    return result


def differentiate_surface_moments(vertices: Iterable[Iterable[float]], faces: Iterable[Iterable[int]], density: float = 1.0) -> np.ndarray:
    """Algorithm 2, summed into the 10-by-(3N) vertex derivative matrix."""

    points = _coordinates(vertices)
    triangles = _triangles(faces, len(points))
    if isinstance(density, bool) or not isinstance(density, (int, float)) or not math.isfinite(density):
        raise ValueError("density must be finite")
    result = np.zeros((10, 3 * len(points)), dtype=float)
    for face in triangles:
        local = triangle_trace(*(points[index] for index in face))["jacobian"]
        for local_vertex, global_vertex in enumerate(face):
            result[:, 3 * global_vertex:3 * global_vertex + 3] += local[:, 3 * local_vertex:3 * local_vertex + 3]
    result *= float(density)
    if not np.all(np.isfinite(result)):
        raise OverflowError("published derivative integration overflowed")
    return result
