"""Independent rational-to-binary64 plane check, not a source mesh audit.

Compute an exact rational origin + k*side/2**level, then round once.
No native grid code or libm fma is called. Enumerated numerical roots are
test probes, not original-paper inputs or material-validity certificates.
"""
from fractions import Fraction
import json
import math


def planes(origin, side, level):
    o, s = Fraction(origin), Fraction(side)
    return [float(o + k * s / (1 << level)) for k in range((1 << level) + 1)]


def check():
    binade = math.ldexp(1.0, 48)
    examples = []
    for origin, expected in ((math.nextafter(binade, 0.0), 7),
                             (-binade - .5, 0), (binade, 0)):
        p = planes(origin, 1.0, 4)
        collapsed = sum(a == b for a, b in zip(p, p[1:]))
        assert collapsed == expected
        assert p[0] < p[1]  # The old first-cell check cannot detect case0.
        coarse = planes(origin, 1.0, 3)
        assert all(a < b for a, b in zip(coarse, coarse[1:]))
        examples.append({"origin_hex": origin.hex(), "collapsed": collapsed})

    accepted = rejected = checked = 0
    # Broad sign/binade/scale coverage, bounded small enumerations. Passing
    # this sweep is not the mathematical proof for all supported depths.
    for exponent in (-100, -20, 0, 48, 100):
        b = math.ldexp(1.0, exponent)
        ulp = math.ulp(b)
        for sign in (-1.0, 1.0):
            for origin in (sign*b, math.nextafter(sign*b, 0.0),
                           math.nextafter(sign*b, sign*math.inf), sign*b-8*ulp):
                for side in (8*ulp, 16*ulp, 32*ulp, 64*ulp, b/2):
                    for level in range(7):
                        p = planes(origin, side, level)
                        magnitude = max(abs(p[0]), abs(p[-1]))
                        spacing = math.nextafter(magnitude, math.inf)-magnitude
                        step = math.ldexp(side, -level)
                        if step > spacing:
                            assert all(a < z for a, z in zip(p, p[1:])), (
                                origin.hex(), side.hex(), level)
                            accepted += 1
                        else:
                            rejected += 1
                        checked += 1
    return {"scope": "SYNTHETIC_NUMERICAL_GUARD_CHECK",
            "exact_rational_rounding": True, "examples": examples,
            "roots_checked": checked, "strict_bound_admitted": accepted,
            "strict_bound_rejected": rejected,
            "admitted_collapsed_intervals": 0}


if __name__ == "__main__":
    print(json.dumps(check(), sort_keys=True, indent=2, allow_nan=False))
