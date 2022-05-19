import math

from ec import (G1FromBytes, G1Generator, G1Infinity, G2FromBytes, G2Generator,
                G2Infinity, JacobianPoint, default_ec, default_ec_twist,
                sign_Fq2, twist, untwist, y_for_x)
from fields import Fq, Fq2, Fq6, Fq12, setDebug

G1Element = JacobianPoint
G2Element = JacobianPoint


def test_paraemters():
    # E(Fq) = q + 1 -t
    t = default_ec.x + 1
    assert default_ec.q + 1 - t == default_ec.n * default_ec.h

    def q_x(x):
        return ((x - 1) ** 2) / 3 * r_x(x) + x

    def r_x(x):
        return x ** 4 - x ** 2 + 1

    def t_x(x):
        return x + 1

    # https://eprint.iacr.org/2006/372.pdf
    # https://hackmd.io/@yelhousni/bls12_subgroup_check
    # https://eprint.iacr.org/2006/372.pdf p28
    # r = p + 1 - (x+1)
    assert q_x(default_ec.x) / default_ec.q == 1
    assert r_x(default_ec.x) / default_ec.n == 1
    assert r_x(default_ec.x) * default_ec.h / (q_x(default_ec.x) + 1 - t_x(default_ec.x)) == 1
    print('rho : ', math.log10(default_ec.q) / math.log10(default_ec.n))
    print(r_x(default_ec.x))
    print(default_ec.h)
    print(0x396c8c005555e1568c00aaab0000aaa)
    print(
        0x5d543a95414e7f1091d50792876a202cd91de4547085abaa68a205b2e5a7ddfa628f1cb4d9e82ef21537e293a6691ae1616ec6e786f0c70cf1c38e31c7238e5)
    E = 31
    s = 7
    a = s + (s % 4 * 11)
    b = ((s * 3) % 4) * 11
    print(a, b)

    print(7 * E % 44, 8 * E % 44)
    print(a * E % 44, b * E % 44)


def test_fq():
    setDebug(True)
    Fq6(11, Fq2(7, 1, 1), Fq2(7, 2, 2), Fq2(7, 3, 3)) * Fq6(17, Fq2(7, 4, 4), Fq2(7, 5, 5), Fq2(7, 6, 6))


def test_twist():
    setDebug(False)
    q = default_ec.q
    g = G1Generator()
    g2 = G2Generator()
    s = g2 + g2
    assert g2.is_on_curve()
    assert untwist(twist(s.to_affine())) == s.to_affine()
    assert untwist(5 * twist(s.to_affine())) == (5 * s).to_affine()

    '''
    Fq12 : c0 + c1w, w^2=r=v => v + 0w , w^3=w^2w=vw => 0 + vw
    wsq = Fq12(ec.q, f.root, Fq6.zero(ec.q))
    wcu = Fq12(ec.q, Fq6.zero(ec.q), f.root)

    #twist
    new_x = point.x * wsq
    new_y = point.y * wcu
    # untwist
    AffinePoint(point.x / wsq, point.y / wcu, False, ec)

    https://eprint.iacr.org/2022/352.pdf
    https://hackmd.io/@yelhousni/bls12_subgroup_check
    D-type twist: y^2=x^3+b/ξ E'->E(Fq2->Fq12) : (x,y) -> (ξ^1/3x,ξ^1/2y), (w2x, w3x)
    M-type twist: y^2=x3^+bξ  E'->E(Fq2->Fq12) : (x,y) -> (ξ^2/3x/ξ,ξ^1/2y/ξ), (x/w2, x/w3)

    bls12-381 is D-type twist ? M-type twist ?
    '''


def test_final_exponentiation():
    # https://mathworld.wolfram.com/CyclotomicPolynomial.html
    # pairing for beginners 113 page
    """
    k = 12, d = k/2
    (q^k -1)/r = [(q^d-1)(q^d+1)/phi_k(q)]*[phi_k(q)/r]
    phi_12 = q^4-q^2+1
    (q^6+1)/phi_12(q) = q^2+1
    element^(q^6-1)(q^2+1)]*[q^4-q^2+1/r]
    =element^q^6 * element^-1 * element^q^2 * element^1 * element^(q^4-q^2+1/r)
    """
    pass

# def test_qi_power():
#     # Frobenius endomorphism 별도로 아래 작업한 이유는 속도인데 방법은?
#     one = Fq(default_ec.q, 1)
#     two = one + one
#     a = Fq2(default_ec.q, two, two)
#     b = Fq6(default_ec.q, a, a, a)
#     c = Fq12(default_ec.q, b, b)
#     for base in (a, b, c):
#         for expo in range(1, base.extension):
#             assert base.qi_power(expo) == pow(base, pow(default_ec.q, expo))
