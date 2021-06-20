import pytest
from pyxmath.finite_mono_polynomial import *


def test_add_sub():
    poly = FiniteMonoPolynomial([4, 1, 0, 1], 11)
    with pytest.raises(FiniteMonoPolynomialError):
        poly.add([1, 10, 10, 10], [2, 2, 2])
    assert poly.add([10, 10, 10], [2, 2, 2]) == [1, 1, 1]

    with pytest.raises(FiniteMonoPolynomialError):
        poly.sub([1, 2, 2, 2], [10, 10, 10])
    assert poly.sub([2, 2, 2], [10, 10, 10]) == [3, 3, 3]

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