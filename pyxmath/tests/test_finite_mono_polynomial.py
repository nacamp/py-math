import pytest
from pyxmath.finite_mono_polynomial import *


def test_add_sub():
    field = FiniteMonoPolynomial([4, 1, 0, 1], 11)
    with pytest.raises(FiniteMonoPolynomialError):
        field.add([1, 10, 10, 10], [2, 2, 2])
    assert field.add([10, 10, 10], [2, 2, 2]) == [1, 1, 1]

    with pytest.raises(FiniteMonoPolynomialError):
        field.sub([1, 2, 2, 2], [10, 10, 10])
    assert field.sub([2, 2, 2], [10, 10, 10]) == [3, 3, 3]
