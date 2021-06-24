import pytest
from pyxmath.number_theory import *
from pyxmath.finite_mono_polynomial import *


# @pytest.mark.skip(reason=".")
def test_add_sub():
    poly = FiniteMonoPolynomial([4, 1, 0, 1], 11)
    with pytest.raises(ValueError):
        poly.add([1, 10, 10, 10], [2, 2, 2])
    assert poly.add([10, 10, 10], [2, 2, 2]) == [1, 1, 1]

    with pytest.raises(ValueError):
        poly([10, 10, 10]) + poly([2, 2, 2, 3])
    assert poly([10, 10, 10]) + poly([2, 2, 2]) == [1, 1, 1]
    assert poly([10, 10, 10]) + [2, 2, 2] == [1, 1, 1]
    assert [2, 2, 2] + poly([10, 10, 10]) == [1, 1, 1]

    with pytest.raises(ValueError):
        poly.sub([1, 2, 2, 2], [10, 10, 10])
    assert poly.sub([2, 2, 2], [10, 10, 10]) == [3, 3, 3]

    assert poly([2, 2, 2]) - poly([10, 10, 10]) == [3, 3, 3]
    assert poly([2, 2, 2]) - [10, 10, 10] == [3, 3, 3]
    assert [2, 2, 2] - poly([10, 10, 10]) == [3, 3, 3]


def test_mul():
    poly = FiniteMonoPolynomial([4, 1, 0, 1], 11)
    assert poly.mul([1, 1, 1], [2, 2, 2]) == [8, 3, 4]
    '''
    1u^3+0u^3+1u^1+4u^0 => u^3 = 10u^+7

                         2 2 2
                        x1 1 1
                        ------
         2u^4+4u^3+6u^2+4u^1+2
      2u^4=2u*u^3=20u*(10u^+7)
      4u^3=4u^3  =  4*(10u^+7)
                        10 7
                       x 2 4
                       ------
               20u^2+54u^1+28
    6u^2+4u^1+2 + 20u^2+54u^1+28 = 26u^2+58u^1+30 = 4u^2+3u^1+8 (mod 11)
    => [8, 3, 4]
    '''
    assert poly([1, 1, 1]) * poly([2, 2, 2]) == [8, 3, 4]
    assert poly([1, 1, 1]) * [2, 2, 2] == [8, 3, 4]
    assert [1, 1, 1] * poly([2, 2, 2]) == [8, 3, 4]


def test_inv_div():
    # x^3+2x+1 역순
    poly = FiniteMonoPolynomial([1, 2, 0, 1], 3)
    assert poly.inv([1, 0, 1]) == [2, 1, 2]
    assert poly.mul([1, 0, 1], [2, 1, 2]) == [1, 0, 0]
    assert poly.div([1, 0, 1], [1, 0, 1]) == [1, 0, 0]

    assert poly([1, 0, 1]) / poly([1, 0, 1]) == [1, 0, 0]
    assert poly([1, 0, 1]) / [1, 0, 1] == [1, 0, 0]
    assert [1, 0, 1] / poly([1, 0, 1]) == [1, 0, 0]

    poly = FiniteMonoPolynomial([1, 1, 0, 1], 2)
    assert poly.inv([1, 1, 1]) == [0, 0, 1]
    assert poly.mul([1, 1, 1], [0, 0, 1]) == [1, 0, 0]
    assert poly.div([1, 1, 1], [1, 1, 1]) == [1, 0, 0]

    poly = FiniteMonoPolynomial([1, 0, 0, 1, 1], 2)
    assert poly.inv([1, 1, 1, 0]) == [0, 1, 1, 1]
    assert poly.mul([1, 1, 1, 0], [0, 1, 1, 1]) == [1, 0, 0, 0]
    assert poly.div([1, 1, 1, 0], [1, 1, 1, 0]) == [1, 0, 0, 0]

    # 5이상인곳에서 inv test
    poly = FiniteMonoPolynomial([1, 2, 0, 1], 5)
    assert poly.inv([2, 4, 3]) == [0, 2, 2]
    assert poly.mul([2, 4, 3], [0, 2, 2]) == [1, 0, 0]
    assert poly.inv([4, 4, 3]) == [2, 0, 3]
    assert poly.mul([4, 4, 3], [2, 0, 3]) == [1, 0, 0]
    # https://math.stackexchange.com/questions/124300/finding-inverse-of-polynomial-in-a-field


def test_etc():
    # pow
    poly = FiniteMonoPolynomial([1, 2, 0, 1], 67)
    m2 = poly.mul([1, 0, 1], [1, 0, 1])
    m3 = poly.mul([1, 0, 1], m2)
    assert poly.pow([1, 0, 1], 3) == m3

    # Frobenius endomorphism
    poly = FiniteMonoPolynomial([1, 0, 1], 67)
    assert poly.pow([16, 2], 67 ** 2) == [16, 2]
    assert frob_end_pi(poly([16, 2]), 67, 2) == [16, 2]
    assert frob_end_pi(poly([39, 30]), 67, 2) == [39, 30]

    poly = FiniteMonoPolynomial([2, 0, 0, 1], 67)
    assert frob_end_pi(poly([8, 4, 15]), 67, 3) == [8, 4, 15]

    poly = FiniteMonoPolynomial([1, 2, 0, 1], 3)
    assert poly([1, 0, 1]).neg() == [2, 0, 2]
    assert -poly([1, 0, 1]) == [2, 0, 2]
