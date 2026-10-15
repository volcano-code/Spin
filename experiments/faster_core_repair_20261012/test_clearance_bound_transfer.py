from fractions import Fraction
import pytest


def test_exact_conservative_distance_transfer_and_insufficient_margin():
    from clearance_bound_transfer import transfer
    # sqrt(1) - sqrt(3)*(.1+.1) > .5; use 2 rather than sqrt(3).
    assert transfer('1', '1/10', '1/10', '1/2')['at_least_required']
    assert not transfer('1', '1/10', '1/10', '4/5')['at_least_required']
    assert transfer('1/4', '0', '0', '1/2')['at_least_required']
    with pytest.raises(ValueError):
        transfer('1', '-1/10', '0', '1/2')
