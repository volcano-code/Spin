"""Preregistered comparison and explicit mass-moment-to-SI conversion.

Input moments already include material density. This adapter applies only the
declared model mass and length units; it cannot apply density a second time.
"""
import itertools
import math
import numpy as np


def scenarios(protocol):
    expected = dict(materials=['original', 'optimized'], collision_proxies=['default', 'robust'],
        azimuth_degrees=list(range(0, 360, 45)), time_steps_seconds=[.0005, .00025])
    if any(protocol.get(key) != value for key, value in expected.items()):
        raise ValueError('CONTACT_MATRIX_DIFFERS_FROM_PREREGISTERED_PROTOCOL')
    for key, value in dict(duration_seconds=8., nominal_spin_hz=35., initial_tilt_degrees=2.,
        lateral_friction=.8, rolling_friction=0., spinning_friction=0., linear_damping=0.,
        angular_damping=0., restitution=0., gravity_m_s2=9.81).items():
        if type(protocol.get(key)) not in (int, float) or protocol[key] != value:
            raise ValueError('CONTACT_PARAMETER_CHANGED: '+key)
    if protocol.get('persistent_drive_torque') is not False or protocol.get('world_fixed_pivot') is not False:
        raise ValueError('UNEXPECTED_DRIVE_OR_FIXED_PIVOT')
    cases = [dict(material=m, proxy=p, azimuth_degrees=a, dt_seconds=dt,
                  case_id=f'{m}-{p}-dir{a:03d}-dt{dt:g}')
             for m, p, a, dt in itertools.product(expected['materials'], expected['collision_proxies'],
                 expected['azimuth_degrees'], expected['time_steps_seconds'])]
    if len(cases) != 64 or protocol.get('unique_scenarios') != 64:
        raise ValueError('CONTACT_SCENARIO_COUNT_MISMATCH')
    return cases


def physical_properties(raw_mass_moments, *, mass_unit_kg, length_unit_m):
    raw = np.asarray(raw_mass_moments, dtype=float)
    if (raw.shape != (10,) or not np.isfinite(raw).all() or raw[0] <= 0 or
            not math.isfinite(mass_unit_kg) or mass_unit_kg <= 0 or
            not math.isfinite(length_unit_m) or length_unit_m <= 0):
        raise ValueError('INVALID_MATERIAL_MOMENTS_OR_UNITS')
    mass = raw[0]*mass_unit_kg
    center = raw[1:4]/raw[0]*length_unit_m
    q = np.array([[raw[4], raw[9], raw[8]], [raw[9], raw[5], raw[7]],
                  [raw[8], raw[7], raw[6]]]) * mass_unit_kg*length_unit_m**2
    inertia = np.trace(q)*np.eye(3)-q-mass*((center@center)*np.eye(3)-np.outer(center, center))
    eigenvalues, axes = np.linalg.eigh(inertia)
    if not np.isfinite(eigenvalues).all() or eigenvalues[0] <= 0:
        raise ValueError('NONPOSITIVE_CENTRAL_INERTIA')
    if np.linalg.det(axes) < 0:
        axes[:, 0] *= -1
    # Also enforce the principal-moment triangle inequalities of physical bodies.
    if eigenvalues[2] > eigenvalues[0]+eigenvalues[1] + 1e-12*eigenvalues[2]:
        raise ValueError('NONPHYSICAL_INERTIA_TRIANGLE_INEQUALITY')
    return mass, center, inertia, eigenvalues, axes
