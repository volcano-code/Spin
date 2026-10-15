"""Independent scalar transcription of Spin-It (2014) paper Givens equations.

This module intentionally contains only the published objective/constraint
values.  It imports no native extension, optimization code, or automatic
differentiation implementation.  Derivatives remain a separately measured
native interface property.
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np


PAPER_TEXT_SHA256 = "5cf5ab047189575df48b01796de73c72ceadc602a3380b74f33d2dda81459161"
ORDER = ("M", "LX", "LY", "LZ", "QXX", "QYY", "QZZ", "QYZ", "QXZ", "QXY")
OBJECTIVES = frozenset(("YOYO", "TOP_INERTIA_ONLY", "TOP_LOWERING_ONLY",
                        "TOP_MASS_AND_LOWERING", "STATIC_LOWERING"))


def _finite_scalar(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(name + " must be finite")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(name + " must be finite") from error
    if not math.isfinite(result):
        raise ValueError(name + " must be finite")
    return result


def stability_threshold_squared(
    mass: float,
    gravity: float,
    center_height: float,
    axial_inertia: float,
    middle_inertia: float,
) -> float:
    """Paper Eq. (1): the strict lower bound on stable-spin ``omega**2``.

    The formula is defined only for a positive mass and distinct stable-axis
    ordering ``Ic > Ib``.  This scalar condition does not claim contact
    stability, numerical convergence, or an author experiment outcome.
    """

    m = _finite_scalar(mass, "mass")
    g = _finite_scalar(gravity, "gravity")
    ell = _finite_scalar(center_height, "center_height")
    ic = _finite_scalar(axial_inertia, "axial_inertia")
    ib = _finite_scalar(middle_inertia, "middle_inertia")
    if m <= 0 or g < 0 or ell < 0 or ic <= ib:
        raise ValueError("Eq. (1) requires M>0, g>=0, ell>=0 and Ic>Ib")
    return m * g * ell / (ic - ib)


def evaluate_paper_givens(
    objective: str,
    raw_moments: Iterable[float],
    inertia_weight: float = 1.0,
    center_weight: float = 1.0,
    phi: float = 0.0,
) -> dict[str, object]:
    """Evaluate paper Eq. (2)/(3), Eqs. (5)--(8), or the static objective.

    Raw moment order is ``M,Lx,Ly,Lz,Qxx,Qyy,Qzz,Qyz,Qxz,Qxy``.  For the
    top objectives, the parallel-axis subtraction is ``Lz**2/M``.  The
    returned constraint vector has the exact order used in the paper's
    displayed set constraints, followed by the diagonalization condition.
    """

    if objective not in OBJECTIVES:
        raise ValueError("unsupported Spin-It 2014 objective")
    raw = np.asarray(tuple(raw_moments), dtype=float)
    if raw.shape != (10,) or not np.all(np.isfinite(raw)):
        raise ValueError("raw_moments must be ten finite values")
    gamma_i = _finite_scalar(inertia_weight, "inertia_weight")
    gamma_c = _finite_scalar(center_weight, "center_weight")
    angle = _finite_scalar(phi, "phi")
    if gamma_i < 0 or gamma_c < 0:
        raise ValueError("weights must be nonnegative")
    if raw[0] <= 0:
        raise ValueError("mass must be positive")
    if objective == "STATIC_LOWERING":
        if angle != 0.0:
            raise ValueError("static objective requires phi=0")
        return {
            "objective": gamma_c * (raw[3] / raw[0]) ** 2,
            "constraints": np.asarray((raw[1], raw[2])),
            "formula_scope": "STATIC_LOWERING",
        }

    m, lx, ly, lz, qxx, qyy, qzz, qyz, qxz, qxy = raw
    ic = qxx + qyy
    if ic <= 0:
        raise ValueError("axial inertia must be positive")
    cosine, sine = math.cos(angle), math.sin(angle)
    delta = 0.0 if objective == "YOYO" else lz * lz / m
    ia = cosine * cosine * (qyy + qzz) + sine * sine * (qxx + qzz) + 2.0 * cosine * sine * qxy - delta
    ib = sine * sine * (qyy + qzz) + cosine * cosine * (qxx + qzz) - 2.0 * cosine * sine * qxy - delta
    value = gamma_i * ((ia / ic) ** 2 + (ib / ic) ** 2)
    if objective == "TOP_LOWERING_ONLY":
        value += gamma_c * (lz / m) ** 2
    elif objective == "TOP_MASS_AND_LOWERING":
        value += gamma_c * lz * lz

    constraints = [lx, ly]
    if objective == "YOYO":
        constraints.append(lz)
    constraints.extend((qxz, qyz,
                        cosine * sine * (qxx - qyy) + (cosine * cosine - sine * sine) * qxy))
    return {
        "objective": value,
        "constraints": np.asarray(constraints),
        "principal_moments": np.asarray((ia, ib, ic)),
        "formula_scope": "YOYO_EQ2_EQ5_EQ6" if objective == "YOYO" else "TOP_EQ3_EQ7_EQ8",
    }
