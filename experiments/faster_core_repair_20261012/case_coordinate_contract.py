"""Equivalent units only: input moments must already be in the same target frame.

The quadric helper only changes length units in the same target frame; it
does not translate an origin, reconstruct material or certify a production map.
"""
import math
import numpy as np


def positive(value, name):
    result=float(value)
    if not math.isfinite(result) or result<=0.:
        raise ValueError('INVALID_COORDINATE_'+name)
    return result


def normalized_faster_weight(model_weight, length):
    length=positive(length,'LENGTH')
    weight=float(model_weight)
    normalized=weight*length
    if not math.isfinite(weight) or weight<0. or not math.isfinite(normalized):
        raise ValueError('INVALID_COORDINATE_WEIGHT')
    return normalized


def normalize_target_moments(model_moments, length, density):
    length=positive(length,'LENGTH')
    density=positive(density,'DENSITY')
    values=np.array(model_moments,dtype=float,copy=True)
    if values.shape!=(10,) or not np.isfinite(values).all() or values[0]<=0.:
        raise ValueError('INVALID_TARGET_MOMENTS')
    try:
        scales=np.array([density*length**power for power in (3,4,4,4,5,5,5,5,5,5)])
    except OverflowError as error:
        raise ValueError('UNREPRESENTABLE_COORDINATE_SCALE') from error
    if not np.isfinite(scales).all() or np.any(scales<=0.):
        raise ValueError('UNREPRESENTABLE_COORDINATE_SCALE')
    with np.errstate(over='ignore',invalid='ignore'):
        normalized=values/scales
    if not np.isfinite(normalized).all() or normalized[0]<=0.:
        raise ValueError('UNREPRESENTABLE_NORMALIZED_MOMENTS')
    return normalized


def normalize_target_quadric(coefficients, length):
    """Return unit coefficients and d(unit coefficients)/d(input coefficients).

    For x_model=L*x_normalized, p=D*q/||D*q||. Positive normalization
    preserves the selected material side. Canonical cross terms are matrix
    entries (twice yz/xz/xy in the polynomial), unchanged by isotropic scaling.
    This Jacobian describes the coordinate conversion, not material transport.
    """
    length=positive(length,'LENGTH')
    q=np.array(coefficients,dtype=float,copy=True)
    if q.shape!=(10,) or not np.isfinite(q).all():
        raise ValueError('INVALID_QUADRIC_COEFFICIENTS')
    with np.errstate(over='ignore',invalid='ignore',under='ignore'):
        scale=np.array([1.,length,length,length,*([length*length]*6)])
        transformed=scale*q
        norm=float(np.linalg.norm(transformed))
    if not np.isfinite(scale).all() or np.any(scale<=0.) or not math.isfinite(norm) or norm<=0.:
        raise ValueError('UNREPRESENTABLE_QUADRIC_SCALE')
    p=transformed/norm
    jac=(np.eye(10)-np.outer(p,p))*scale[None,:]/norm
    if not np.isfinite(jac).all():
        raise ValueError('UNREPRESENTABLE_QUADRIC_JACOBIAN')
    return p,jac
