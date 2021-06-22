import pytest
from pyxmath.elliptic_curve import *
from pyxmath.number_theory import *
from pyxmath.number_theory import Field as F


def test_add():
    ec = EC([1, 1, 0, 1])
    p = 101
    f = F(p)
    add1 = ec.add(PT(f(47), f(12)), PT(f(47), f(12)))
    add2 = ec.add(PT(f(47), f(12)), add1)
    assert add1 == PT(f(6), f(83))
    assert add2 == PT(f(23), f(77))
    assert ec.mul(PT(f(47), f(12)), 2) == add1
    assert ec.mul(PT(f(47), f(12)), 3) == add2
    # Points exist
    assert modular_sqrt(ec.find_y_square(47), p) == 12
    assert p - 12 == 89
    assert modular_sqrt(ec.find_y_square(0), p) == 1
    assert modular_sqrt(ec.find_y_square(87), p) == 24
    # Points does not exist
    assert modular_sqrt(ec.find_y_square(4), p) == 0
    assert modular_sqrt(ec.find_y_square(7), p) == 0
    # find_points
    assert len(find_points(ec, p)) == 105