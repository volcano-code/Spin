"""Exact local criteria for a DECLARED ideal upright principal-axis top.

References: Lewis et al. (1992), Sections 2.3 and 4.1; Spin-It (2014),
Equation 1. This is an independent validation oracle, not an optimization
objective and not a contact simulator. Inertia is shifted from CoM to the
fixed pivot before using the heavy-top criteria. No geometry alignment is
inferred or silently projected: the caller declares an exactly axial CoM
and a principal frame with spin along axis 3. A mesh needs its own checks.
"""
from fractions import Fraction as Q


def _rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, float, str, Q)):
        raise ValueError('finite real scalar or rational string required')
    try:
        return Q(value)
    except (ValueError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError('finite real scalar or rational string required') from exc


def principal_axis_critical_report(mass, height, central_principal_moments,
                                   omega, *, gravity='9.81'):
    """Return conditional algebraic facts, never full material/physics admission.

    Scalar units must be consistent (e.g. kg, m, s, radians). Rational strings
    express exact decimals; floats are interpreted as their exact binary64
    values, not rounded onto presumed equalities. Equality/repeated roots are
    not accepted as strict spectral stability. Below a sufficient energy bound
    is not, by itself, an instability proof.
    """
    m, h, w, g = map(_rational, (mass, height, omega, gravity))
    if min(m,h,g) <= 0:
        raise ValueError('positive mass, axial height and gravity required')
    if not isinstance(central_principal_moments, (tuple,list)) or len(central_principal_moments) != 3:
        raise ValueError('three central principal moments required')
    ca, cb, cc = map(_rational, central_principal_moments)
    if min(ca,cb,cc) <= 0 or 2*max(ca,cb,cc) > ca+cb+cc:
        raise ValueError('nonphysical central inertia')
    a, b, c = ca+m*h*h, cb+m*h*h, cc
    torque = m*g*h
    speed2 = w*w
    # Characteristic polynomial of the transverse four-state Euler/vertical
    # linearization, independently evaluated without constructing that matrix.
    pa = ((2*a*b-c*(a+b)+c*c)*speed2-(a+b)*torque)/(a*b)
    pb = ((c-a)*speed2-torque)*((c-b)*speed2-torque)/(a*b)
    discriminant = pa*pa-4*pb
    if pb < 0 or discriminant < 0 or pa < 0:
        linear = 'EXPONENTIALLY_UNSTABLE'
    elif pb > 0 and pa > 0 and discriminant > 0:
        linear = 'SPECTRALLY_STABLE_SIMPLE'
    else:
        linear = 'MARGINAL_REQUIRES_FURTHER_ANALYSIS'

    threshold = None
    if a == b:
        condition = 'NOT_APPLICABLE_REPEATED_TRANSVERSE_MOMENTS'
    elif c <= max(a,b):
        condition = 'NOT_APPLICABLE_AXIS_NOT_MAXIMAL'
    else:
        threshold = torque/(c-max(a,b))
        condition = 'SATISFIED_FOR_DECLARED_IDEAL_MODEL' if speed2 > threshold else 'NOT_SATISFIED'
    symmetric = 4*a*torque/(c*c) if a == b else None
    return dict(schema='fixed-pivot-principal-axis-criteria-v1',
        scope='CONDITIONAL_IDEAL_UPRIGHT_MODEL_NOT_MESH_ADMISSION',
        assumptions=['EXACT_AXIAL_COM_ABOVE_FIXED_POINT', 'PRINCIPAL_FRAME_SPIN_AXIS_3',
                     'NO_CONTACT_SLIP_FRICTION_OR_DISSIPATION'],
        input_origin='CENTER_OF_MASS', analyzed_origin='FIXED_PIVOT',
        mass=str(m), height=str(h), gravity=str(g), angular_speed=str(w),
        central_principal_moments=list(map(str,(ca,cb,cc))),
        pivot_principal_moments=list(map(str,(a,b,c))),
        asymmetric_sufficient_omega_squared=None if threshold is None else str(threshold),
        asymmetric_condition=condition,
        axisymmetric_linear_omega_squared=None if symmetric is None else str(symmetric),
        polynomial={'a':str(pa),'b':str(pb),'discriminant':str(discriminant)},
        linear_state=linear,
        material_admission='NOT_GRANTED', geometry_alignment='NOT_CHECKED',
        global_contact_stability='NOT_GRANTED', paper_reproduction='NOT_GRANTED')
