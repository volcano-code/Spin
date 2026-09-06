"""Local consistency of the published surface kernel against signed tetrahedra.

All gates below are IMPLEMENTATION_DEFAULT engineering choices, not paper
tolerances. This does not establish geometry legality, source/build attestation,
optimization feasibility, a release gate, or reproduction of a paper case.

The verification frame is selected once from INPUT: origin is the joint bbox
midpoint, L its largest span, and R the largest region density. Coordinates
q=(x-origin)/L and densities rho/R remain in that frame for all perturbations.
The independent integrator's internal conditioning may change, but it returns
moments in this SAME q frame; we do not use its internal frame as a new unit.

For monomial degree d, reported moments have units R*L**(3+d), and derivatives
with respect to physical vertex coordinates have units R*L**(2+d). The latter
are derivatives of centered monomials, keeping origin/L/R constant. To recover
world raw moments: M=M_c, L_i=L_ci+o_i*M, and
Q_ij=Q_cij+o_i*L_cj+o_j*L_ci+o_i*o_j*M. Apply the same linear map to Jacobian
rows AFTER the stated derivative unit conversion. Physical-coordinate scale,
density and translation covariance is tested separately from this normalized
numerical check. Native matrices use all symmetric terms in ORDER, 3N columns
interleaved xyz, with each region's columns concatenated in input order.
"""

from __future__ import annotations

from collections import Counter
import math
from pathlib import Path

import numpy as np

from ..spin2014 import _native_module
from . import mass_moments
from .mass_moments import DensityRegion, TriangleMesh
from .run_evidence import _digest, _record


ORDER = ("M", "LX", "LY", "LZ", "QXX", "QYY", "QZZ", "QYZ", "QXZ", "QXY")
ALGORITHM = "spin-it-2014-supplement-algorithms-1-2-surface-polynomial"
STEPS = (2e-5, 1e-5, 5e-6)
RELATIVE_TOLERANCE = 1e-6
NEAR_ZERO_CUTOFF = 1e-5
ABSOLUTE_TOLERANCE = 1e-10
MOMENT_TOLERANCE = 1e-10
ROUNDING_MULTIPLIER = 64
DEGREES = (0, 1, 1, 1, 2, 2, 2, 2, 2, 2)


def _independent_vector(mesh: TriangleMesh, density: float) -> np.ndarray:
    """Density once, using ONLY the independent Python geometric integrator."""
    m = mass_moments.integrate_geometric_moments(mesh)
    q = m.second_moment
    values = density * np.array((m.volume, *m.first_moment,
                                q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1]))
    if values.shape != (10,) or not np.all(np.isfinite(values)):
        raise ValueError("independent moments must be ten finite values")
    return values


def assess_error_trend(coarse_error, fine_error, coarse_floor, fine_floor) -> str:
    """Classify adjacent halving steps; this never substitutes for accuracy."""
    values = (coarse_error, fine_error, coarse_floor, fine_floor)
    if not all(math.isfinite(v) and v >= 0 for v in values):
        raise ValueError("trend inputs must be finite nonnegative measurements")
    coarse_noise, fine_noise = coarse_error <= coarse_floor, fine_error <= fine_floor
    if coarse_noise and fine_noise:
        return "NOISE_LIMITED"
    if coarse_noise or fine_noise:
        return "FLOOR_TRANSITION"
    ratio = coarse_error / fine_error
    return "SECOND_ORDER" if 2 <= ratio <= 6 else "INCONSISTENT"


def _finite_array(value, shape, label):
    array = np.asarray(value)
    if array.dtype.kind not in "fiu" or array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(label + " has invalid shape, dtype or nonfinite entries")
    return np.asarray(array, dtype=float)


def verify_surface_derivatives(regions) -> dict:
    """Measure all entries at all three fixed steps using the actual backend.

    Each region must have positive net oriented volume, as required by the
    independent integrator. Inward cavity components belong within that region's
    mesh. The native primitive also permits standalone signed inward components.
    Malformed data/backend failure raises explicitly; numerical disagreement
    returns status='inconsistent'. There is no caller-provided pass flag.
    """
    regions = tuple(regions)
    if not regions or any(not isinstance(r, DensityRegion) for r in regions):
        raise ValueError("expected nonempty DensityRegion inputs")
    if any(isinstance(r.density, bool) or not math.isfinite(r.density) or r.density <= 0
           for r in regions):
        raise ValueError("region density must be finite and positive")
    vertices = []
    for region in regions:
        if not isinstance(region.mesh, TriangleMesh) or not region.mesh.vertices:
            raise ValueError("region requires nonempty TriangleMesh")
        vertices.append(_finite_array(region.mesh.vertices, (len(region.mesh.vertices), 3), "vertices"))
    combined = np.concatenate(vertices)
    minima, maxima = combined.min(axis=0), combined.max(axis=0)
    spans = maxima-minima
    length = float(np.max(spans))
    origin = minima + spans*.5
    density_scale = float(max(r.density for r in regions))
    if not np.all(np.isfinite(origin)) or not math.isfinite(length) or length <= 0:
        raise ValueError("input normalization is unrepresentable or degenerate")
    moment_units = [density_scale * length**(3+d) for d in DEGREES]
    derivative_units = [density_scale * length**(2+d) for d in DEGREES]
    if not all(math.isfinite(v) and v > 0 for v in (*moment_units, *derivative_units)):
        raise ValueError("normalization units are unrepresentable")
    normalized = []
    for r, v in zip(regions, vertices):
        density = r.density/density_scale
        if density <= 0:
            raise ValueError("normalized density underflow")
        normalized.append(DensityRegion(TriangleMesh((v-origin)/length, r.mesh.faces), density))

    native = _native_module()
    function = getattr(native, "differentiate_surface_raw_moments", None)
    if not callable(function):
        raise RuntimeError("published 2014 surface differential backend is unavailable")
    # Reuse the existing observed-file identity helper, including its hash-race
    # checks. These records identify observed bytes, never attest their build.
    identity = {
        "native_module": _record(Path(native.__file__).resolve()),
        "checker": _record(Path(__file__).resolve()),
        "independent_integrator": _record(Path(mass_moments.__file__).resolve()),
        "input_sha256": _digest([{"vertices": v.tolist(), "faces": [list(f) for f in r.mesh.faces],
                                   "density": float(r.density)} for r, v in zip(regions, vertices)]),
        "build_source_attestation": "NOT_ESTABLISHED",
    }
    native_moments, independent_moments, jacobians, region_columns = [], [], [], []
    column = 0
    for index, region in enumerate(normalized):
        count = 3*len(region.mesh.vertices)
        independent_moments.append(_independent_vector(region.mesh, region.density))
        value = function(region.mesh.vertices, region.mesh.faces, region.density)
        if not isinstance(value, dict) or tuple(value.get("order", ())) != ORDER or \
                value.get("algorithm") != ALGORITHM or value.get("coordinates") != "interleaved-vertex-xyz":
            raise ValueError("native algorithm, ordering or coordinate identity mismatch")
        native_moments.append(_finite_array(value.get("moments"), (10,), "native moments"))
        jacobians.append(_finite_array(value.get("jacobian"), (10, count), "native Jacobian"))
        region_columns.append({"region": index, "first_column": column, "column_count": count,
                               "normalized_density": region.density})
        column += count
    analytic = np.concatenate(jacobians, axis=1)
    # Check regional errors as well as their aggregate, so cancellation between
    # density regions cannot conceal an incorrect native region.
    regional_moment_errors = np.abs(np.array(native_moments)-np.array(independent_moments))
    total_native, total_independent = np.sum(native_moments, axis=0), np.sum(independent_moments, axis=0)
    moment_error = max(float(np.max(regional_moment_errors)), float(np.max(np.abs(total_native-total_independent))))
    steps = []
    for h in STEPS:
        finite = np.zeros_like(analytic)
        floors = np.zeros_like(analytic)
        plus_values, minus_values = np.zeros_like(analytic), np.zeros_like(analytic)
        column = 0
        for region in normalized:
            for vertex in range(len(region.mesh.vertices)):
                for axis in range(3):
                    plus = np.array(region.mesh.vertices, dtype=float)
                    minus = plus.copy()
                    plus[vertex, axis] += h
                    minus[vertex, axis] -= h
                    p = _independent_vector(TriangleMesh(plus, region.mesh.faces), region.density)
                    m = _independent_vector(TriangleMesh(minus, region.mesh.faces), region.density)
                    plus_values[:, column], minus_values[:, column] = p, m
                    finite[:, column] = (p-m)/(2*h)
                    floors[:, column] = ROUNDING_MULTIPLIER*np.finfo(float).eps*np.maximum(1, np.maximum(abs(p), abs(m)))/h
                    column += 1
        if not np.all(np.isfinite(finite)):
            raise ValueError("finite differences are nonfinite")
        errors = abs(analytic-finite)
        magnitude = np.maximum(abs(analytic), abs(finite))
        near = magnitude <= NEAR_ZERO_CUTOFF
        relative = np.divide(errors, magnitude, out=np.zeros_like(errors), where=~near)
        consistent = np.where(near, errors <= ABSOLUTE_TOLERANCE, relative <= RELATIVE_TOLERANCE)
        steps.append({
            "h": h, "finite_difference": finite.tolist(),
            "regional_plus_moments": plus_values.tolist(), "regional_minus_moments": minus_values.tolist(),
            "absolute_errors": errors.tolist(), "smooth_relative_errors": relative.tolist(),
            "near_zero_mask": near.tolist(), "roundoff_floors": floors.tolist(),
            "entry_consistent": consistent.tolist(), "all_entries_consistent": bool(np.all(consistent)),
            "max_absolute_error": float(errors.max()),
            "max_smooth_relative_error": float(relative.max()),
            "max_near_zero_absolute_error": float(np.max(errors[near], initial=0)),
        })
    trends = []
    for coarse, fine in zip(steps, steps[1:]):
        findings = [[assess_error_trend(coarse["absolute_errors"][r][c], fine["absolute_errors"][r][c],
                                      coarse["roundoff_floors"][r][c], fine["roundoff_floors"][r][c])
                     for c in range(column)] for r in range(10)]
        counts = dict(Counter(item for row in findings for item in row))
        trends.append({"coarse_h": coarse["h"], "fine_h": fine["h"],
                       "findings": findings, "counts": counts})
    consistent = moment_error <= MOMENT_TOLERANCE and all(s["all_entries_consistent"] for s in steps) and \
        all("INCONSISTENT" not in t["counts"] for t in trends)
    return {
        "scope": "local-surface-tetrahedral-differential-consistency", "status": "consistent" if consistent else "inconsistent",
        "algorithm": ALGORITHM, "order": list(ORDER), "identity": identity,
        "normalization": {"origin": origin.tolist(), "length_scale": length, "density_scale": density_scale,
                          "moment_units": moment_units, "physical_centered_derivative_units": derivative_units,
                          "degrees": list(DEGREES), "frame_selection": "fixed-input-bbox"},
        "controls": {"classification": "IMPLEMENTATION_DEFAULT", "steps": list(STEPS),
                     "smooth_relative_tolerance": RELATIVE_TOLERANCE, "near_zero_cutoff": NEAR_ZERO_CUTOFF,
                     "near_zero_absolute_tolerance": ABSOLUTE_TOLERANCE, "moment_absolute_tolerance": MOMENT_TOLERANCE,
                     "roundoff_floor": "64*float64_epsilon*max(1,abs(m_plus),abs(m_minus))/h",
                     "second_order_ratio_interval": [2, 6], "accuracy_required_in_noise_regime": True},
        "region_columns": region_columns, "analytic_jacobian": analytic.tolist(),
        "native_moments": total_native.tolist(), "independent_moments": total_independent.tolist(),
        "regional_moment_absolute_errors": regional_moment_errors.tolist(), "moment_max_absolute_error": moment_error,
        "steps": steps, "trends": trends,
    }
