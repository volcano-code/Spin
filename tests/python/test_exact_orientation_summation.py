"""Exact orientation sums must survive cancellation and growing denominators."""
from fractions import Fraction
import numpy as np
import pytest

from make_it_stand.verification import cherchi_motion_differential as motion


def test_balanced_exact_sum_preserves_large_cancellation_and_empty_sum():
    accumulator = motion._ExactFractionSum()
    assert accumulator.total() == 0
    values = [Fraction(2**400), Fraction(1,3), -Fraction(2**400), Fraction(-1,7)]
    for value in values:
        accumulator.add(value)
    assert accumulator.total() == Fraction(4,21)


def test_balanced_exact_sum_preserves_all_nonpower_of_two_prefixes():
    accumulator = motion._ExactFractionSum()
    expected = Fraction(0)
    for i in range(1,101):
        value = Fraction((-1)**i,2*i+1)
        expected += value
        accumulator.add(value)
        assert accumulator.total() == expected


@pytest.mark.parametrize('flip', [False,True])
def test_orientation_uses_the_existing_balanced_reduction(monkeypatch,flip):
    points = [(Fraction(0),)*3,(Fraction(1),Fraction(0),Fraction(0)),
              (Fraction(0),Fraction(1),Fraction(0)),(Fraction(0),Fraction(0),Fraction(1))]
    faces = np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]],dtype=np.int64)
    if flip:
        faces[1] = faces[1][[0,2,1]]
    calls=[]
    original=motion._ExactFractionSum.add
    def observe(self,value):
        calls.append(value)
        return original(self,value)
    monkeypatch.setattr(motion._ExactFractionSum,'add',observe)
    _,audit=motion.normalize_closed_component_orientation(np.asarray(points,dtype=float),faces,
                                                         rational_vertices=points)
    assert audit['components'][0]['exact_signed_volume']=='1/6'
    assert len(calls)==(8 if flip else 4)
