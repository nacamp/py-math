import pytest
from pyxmath.number_theory import *
from pyxmath.finite_mono_polynomial import *


# @pytest.mark.skip(reason=".")
def test_field():
    assert xgcd(3, 9)[0] == 3
    if sys.version_info.major >= 3 and sys.version_info.minor >= 8:
        assert inv(3, 7) == pow(3, -1, 7)
    f = Field(5)
    assert f.add(3, 4) == 2
    assert Field(5, 3) + Field(5, 4) == 2
    assert Field(5, 3) + 4 == 2
    assert 4 + Field(5, 3) == 2

    assert f.sub(4, 3) == 1
    assert Field(5, 4) - Field(5, 3) == 1
    assert Field(5, 4) - 3 == 1
    assert 4 - Field(5, 3) == 1

    assert f.mul(3, 4) == 2
    assert Field(5, 3) * Field(5, 4) == 2
    assert Field(5, 3) * 4 == 2
    assert 3 * Field(5, 4) == 2

    assert f.inv(2) == 3
    assert f.div(3, 2) == 4
    assert Field(5, 3) / Field(5, 2) == 4
    assert 3 / Field(5, 2) == 4
    assert Field(5, 3) / 2 == 4

    assert f.pow(3, 2) == 4
    assert Field(5, 3) ** 2 == 4
    assert Field(5, 3) ** Field(5, 2) == 4

    assert f.neg(1) == 4
    assert -Field(5, 1) == 4
    assert -Field(5, 1) == 4

    assert -f(1) == 4


def test_frob_end_pi():
    q = 67
    assert frob_end_pi(15, q, 1) == 15
    assert frob_end_pi(50, q, 1) == 50


def test_yield_operation():
    p = 67
    # NotImplemented
    with pytest.raises(TypeError):
        Field(p, 3) + 4.4

    poly = FiniteMonoPolynomial([1, 2, 0, 1], p)
    assert Field(p, 3) + poly([1, 1, 1]) == [4, 1, 1]


def test_rmod():
    p = 71
    assert rmod(0.5, p) == 36
    assert rmod(72, p) == 1
    assert rmod(72.5, p) == 37


def test_modular_sqrt_for_power():
    # r2 = a (mod p)  4**2 % 7 == 9 % 7
    # assert modular_sqrt(9, 7) == 31
    assert modular_sqrt(9, 7) == 4
    assert modular_sqrt_for_square_of_prime(30, 7**2) == 31
    assert modular_sqrt_for_square_of_prime(13290059, 127 ** 2) == 14109
    # https://math.stackexchange.com/questions/95443/how-to-compute-modular-square-roots-when-modulus-is-non-prime

def test_chinese_remainder():
    m = [2, 3, 7, 23]
    a = [1, 0, 1, 20]
    assert chinese_remainder(m, a) == 687
