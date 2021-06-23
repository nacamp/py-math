import pytest
from pyxmath.number_theory import *


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

    # __call__
    assert -f(1) == 4
    assert -f(1) == 4


def test_frob_end_pi():
    q = 67
    assert frob_end_pi(15, q, 1) == 15
    assert frob_end_pi(50, q, 1) == 50
