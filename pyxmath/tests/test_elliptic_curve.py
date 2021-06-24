import pytest
from pyxmath.elliptic_curve import *
from pyxmath.number_theory import *
from pyxmath.number_theory import Field as F
from pyxmath.finite_mono_polynomial import *


#@pytest.mark.skip(reason=".")
def test_field():
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


def test_poly():
    ec = EC([1, 0, 0, 1])
    p = 7691
    poly = FiniteMonoPolynomial([1, 0, 1], p)

    Q = PT(poly([6145, 633]), poly([109, 7372]))
    assert ec.add(Q, Q) == PT(poly([2096, 4448]), poly([6039, 7386]))
    assert ec.add(Q, PT(poly([2096, 4448]), poly([6039, 7386]))) == PT(poly([619, 6652]), poly([1174, 2367]))
    assert ec.mul(Q, 135) == PT(poly([1403, 5806]), poly([2370, 6091]))
