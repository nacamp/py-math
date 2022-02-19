import pytest
from pyxmath.elliptic_curve import *
from pyxmath.number_theory import *
from pyxmath.number_theory import Field as F
from pyxmath.finite_mono_polynomial import *


# @pytest.mark.skip(reason=".")
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


# @pytest.mark.skip(reason=".")
def test_at_projective_jacobian_coordinates():
    ec = EC([1, 1, 0, 1])
    ec_p = EC([1, 1, 0, 1], 'projective')
    ec_j = EC([1, 1, 0, 1], 'jacobian')
    p = 101
    f = F(p)
    add1 = ec._add_at_projective(PT(f(47), f(12)).to_projective(), PT(f(47), f(12)).to_projective())
    add2 = ec._add_at_projective(PT(f(47), f(12)).to_projective(), add1)
    assert add1.to_affine_from_projective() == PT(f(6), f(83))
    assert ec_p.add(PT(f(47), f(12)).to_projective(), PT(f(47), f(12)).to_projective()) == add1
    assert add2.to_affine_from_projective() == PT(f(23), f(77))
    assert ec_p.mul(PT(f(47), f(12)).to_projective(), 2).to_affine_from_projective() == add1.to_affine_from_projective()
    assert ec_p.mul(PT(f(47), f(12)).to_projective(), 3).to_affine_from_projective() == add2.to_affine_from_projective()
    # not same
    # assert ec_p.mul(PT(f(47), f(12)).to_projective(), 3) == add2

    add1 = ec._add_at_jacobian(PT(f(47), f(12)).to_projective(), PT(f(47), f(12)).to_projective())
    add2 = ec._add_at_jacobian(PT(f(47), f(12)).to_projective(), add1)
    assert add1.to_affine_from_jacobian() == PT(f(6), f(83))
    assert ec_j.add(PT(f(47), f(12)).to_projective(), PT(f(47), f(12)).to_projective()) == add1
    assert add2.to_affine_from_jacobian() == PT(f(23), f(77))
    assert ec_j.mul(PT(f(47), f(12)).to_projective(), 2).to_affine_from_jacobian() == add1.to_affine_from_jacobian()
    assert ec_j.mul(PT(f(47), f(12)).to_projective(), 3).to_affine_from_jacobian() == add2.to_affine_from_jacobian()
    # not same
    # assert ec_j.mul(PT(f(47), f(12)).to_projective(), 3) == add2

    ec = EC([1, 0, 0, 1])
    ec_j = EC([1, 0, 0, 1], 'jacobian')
    p = 7691
    poly = FiniteMonoPolynomial([1, 0, 1], p)

    Q = PT(poly([6145, 633]), poly([109, 7372]))
    assert ec.add(Q, Q) == PT(poly([2096, 4448]), poly([6039, 7386]))
    assert ec.add(Q, Q) == ec_j.add(Q.to_projective(), Q.to_projective()).to_affine_from_jacobian()
    # assert ec._add_at_projective(Q.to_projective(),Q.to_projective()).to_affine_from_projective() == PT(poly([2096, 4448]), poly([6039, 7386]))
    # #

    assert ec_j.add(Q.to_projective(),
                    PT(poly([2096, 4448]), poly([6039, 7386])).to_projective()).to_affine_from_jacobian() == PT(
        poly([619, 6652]), poly([1174, 2367]))
    assert ec_j.mul(Q.to_projective(), 135).to_affine_from_jacobian() == PT(poly([1403, 5806]), poly([2370, 6091]))


def test_poly():
    ec = EC([1, 0, 0, 1])
    p = 7691
    poly = FiniteMonoPolynomial([1, 0, 1], p)

    Q = PT(poly([6145, 633]), poly([109, 7372]))
    assert ec.add(Q, Q) == PT(poly([2096, 4448]), poly([6039, 7386]))
    assert ec.add(Q, PT(poly([2096, 4448]), poly([6039, 7386]))) == PT(poly([619, 6652]), poly([1174, 2367]))
    assert ec.mul(Q, 135) == PT(poly([1403, 5806]), poly([2370, 6091]))


# def test_torsion():
#     q = 11
#     ec = EC([4, 0, 0, 1])
#     points = find_points(ec, q)
#     assert r_torsion(ec, points, 3) == [[0, 9], [0, 2], None]
#
#     poly = FiniteMonoPolynomial([1, 0, 1], q)
#     points = find_points2(ec, poly.elements())
#     print(points)
#     # assert len(points)==2
#     assert  poly([0,1])**2  ==   [10, 0]
#     # assert ec.mul(PT(poly([7,9]),poly([0,1])), 2)  ==  [[8, 0], [0, 10]]
#     # assert ec.mul(PT(poly([8,0]),poly([0,1])), 2)  ==  [[8, 0], [0, 10]]
#     # assert ec.mul(PT(poly([7,9]),poly([0,1])), 2)  ==  [[8, 0], [0, 10]]
#     #[8, 0], [10, 0]
#     # assert r_torsion(ec, poly.elements(), 3) == [[0, 9], [0, 2], None]

def test_miller_weil():
    ec = EC([34, 30, 0, 1])
    p = 631
    f = F(p)
    P = PT(f(36), f(60))
    Q = PT(f(121), f(387))
    S = PT(f(0), f(36))
    Q_S = ec.add(Q, S)
    assert Q_S == PT(f(176), f(486))

    # miller
    m = 5
    m_bs = [int(x) for x in bin(m)[2:]][::-1]
    assert ec.miller(P, Q_S, m_bs) == 103
    assert ec.miller(P, S, m_bs) == 219
    assert ec.miller(Q, ec.sub(P, S), m_bs) == 284
    assert ec.miller(Q, -S, m_bs) == 204

    # weil
    # e5(P,Q)
    assert ec.weil_pairing(P, Q, S, 5) == 242


def test_miller_tate():
    ec = EC([-3, 0, 0, 1])
    p = 5
    poly = FiniteMonoPolynomial([2, 0, 1], p)
    # p = 631
    f = F(p)
    P = PT(f(3), f(2))
    Q = PT(poly([1, 1]), poly([2, 4]))
    # S = PT(poly([0, 2]), poly([2, 1]))
    # Q_S = ec.add(Q, S)
    # assert Q_S == PT(poly([1, 3]), poly([2, 0]))
    R = PT(poly([0, 2]), poly([2, 1]))
    Q_R = ec.add(Q, R)
    assert Q_R == PT(poly([1, 3]), poly([2, 0]))
    m = 3
    ms = [int(x) for x in bin(m)[2:]][::-1]
    assert ec.miller(P, Q_R, ms) == [1, 1]
    assert ec.miller(P, R, ms) == [4, 0]
    assert ec.miller(P, Q_R, ms) / ec.miller(P, R, ms) == [4, 4]
    assert ec.tate_pairing(P, Q, R, m) == [4, 4]
    # DQ = ([2]Q) − (Q)
    assert ec.tate_pairing(P, ec.mul(Q, 2), R, m) == [4, 2]
    assert ec.tate_pairing(ec.mul(P, 2), Q, R, m) == [2, 3]


# p≡3 mod 4 Supersingular elliptic curves
def test_tate_m_exceed_p():
    ec = EC([0, 1, 0, 1])
    p = 71
    points = find_points(ec, p)
    poly = FiniteMonoPolynomial([1, 0, 1], p)
    f = F(p)
    P = PT(f(11), f(8))
    Q = PT(poly([60, 0]), poly([0, 8]))
    P_Members = torsion(ec, P, p + 1)
    m = len(P_Members)
    # m cannot exceed p
    if m > p:
        m = p
    # print(len(P_Members), P_Members)
    Q_Members = torsion(ec, Q, p ** 2)
    # print(len(Q_Members), Q_Members)

    print(ec.tate_pairing(P, ec.mul(Q, 2), Q, m))
    print(ec.tate_pairing(P, Q, Q, m) * 2)

    # e(P+Y,Q) = e(P,Q)*e(Y*Q) = [30,33]
    Y = PT(f(2), f(9))
    assert ec.tate_pairing(ec.add(P, Y), Q, Q, m) == ec.tate_pairing(P, Q, Q, m) * ec.tate_pairing(Y, Q, Q, m) \
           == [30, 33]

    # e(P,Q)**2 = e([2]P,Q) = e(P,[2]Q)
    # print(ec.tate_pairing(P, Q, Q, 71) ** 2)
    assert ec.tate_pairing(P, Q, Q, m) ** 2 == ec.tate_pairing(ec.mul(P, 2), Q, Q, m) \
           == ec.tate_pairing(P, ec.mul(Q, 2), Q, m) == [18, 23]

    # e([P,Q)**(2*3) = e([3]P,[2]Q)
    # print(ec.tate_pairing(P, Q, Q, m) ** 6)
    assert ec.tate_pairing(P, Q, Q, m) ** 6 == ec.tate_pairing(ec.mul(P, 3), ec.mul(Q, 2), Q, m) \
           == [57, 36]
    # print((ec.tate_pairing(P, Q, Q, m) ** 2) * 3)
    # print(ec.tate_pairing(P, Q, Q, m) ** 3)
    # print(ec.tate_pairing(P, Q, Q, m) ** 3 + ec.tate_pairing(P, Q, Q, m) ** 3)
    assert p % 4 == 3


def test_etc():
    # is_supersingular
    ec = EC([0, 1, 0, 1])
    # not prime
    p = 75
    assert p % 4 == 3
    assert not is_supersingular(ec, p)
    # prime
    p = 71
    assert p % 4 == 3
    assert is_supersingular(ec, p)
    p = 5
    assert p % 3 == 2
    assert not is_supersingular(ec, p)

    ec = EC([1, 0, 0, 1])
    p = 5
    assert p % 3 == 2
    assert is_supersingular(ec, p)

    # 7*7 field
    # https://graui.de/code/elliptic2/
    ec = EC([0, 1, 0, 1])
    p = 7 * 7
    points = find_points(ec, p)
    assert len(points) == 50
