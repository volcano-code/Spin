"""Exact sufficient inequality; callers must separately bind geometric proofs."""
from fractions import Fraction


def transfer(squared_distance, outer_linf_bound, design_linf_bound, required):
    distance2, outer, design, wall = map(Fraction,
        (squared_distance, outer_linf_bound, design_linf_bound, required))
    if min(distance2, outer, design) < 0 or wall <= 0:
        raise ValueError('INVALID_DISTANCE_TRANSFER_INPUT')
    # Bidirectional Hausdorff perturbation changes set distance by at most
    # the sum of Euclidean bounds. In R^3, ||x||2 <= sqrt(3)||x||inf <= 2||x||inf.
    sufficient2 = (wall + 2 * (outer + design)) ** 2
    return dict(at_least_required=distance2 >= sufficient2,
                required_exact=str(wall), sufficient_squared_distance_exact=str(sufficient2),
                available_squared_distance_exact=str(distance2),
                method='EXACT_HAUSDORFF_LINF_TO_L2_FACTOR_2_SUFFICIENT_V1')
