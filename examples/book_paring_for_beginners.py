# https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjhiuXAzbzrAhWzL6YKHQohB1oQFjAAegQIBBAB&url=http%3A%2F%2Fwww.craigcostello.com.au%2Fpairings%2FPairingsForBeginners.pdf&usg=AOvVaw1H5dLtelG00vWsvWRGxBNZ
# ParingsForBeginners.pdf
import inspect
import math
from pyxmath import *
from pyxmath.elliptic_curve import *
from pyxmath.finite_mono_polynomial import *


def find_r(order):
    rs = []
    for i in range(1, order + 1):
        if order % i == 0:
            rs.append(i)
    return rs


def solve_poly(coef, x):
    return sum(c * x ** i for i, c in enumerate(coef))


# x^2  ~ q (mod n)
def qr(x_2, mod):
    if x_2 >= mod:
        x_2 = x_2 % mod
    results = []
    for x in range(mod):
        if x ** 2 % mod == x_2:
            results.append(x)
    return results


def g_p_q(p, q, s):
    y_p = p.y
    x_p = p.x
    x_q = q.x
    print(f'y +{-s}x + {-y_p + s * x_p}')  # y -y_p -s(x-x_p)
    print('/')
    print(f'x+ {x_p + x_q - (s ** 2)}')


def p22():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([1, 1, 0, 1])
    p = 101
    points = find_points(ec, p)
    order = len(points)
    print('points : ', points)
    print(f'order = {order} = {get_prime_factors(order)}')
    f = Field(101)
    print('-------------------')
    # P is a generator
    P = PT(f(47), f(12))
    print(f'{order // 1}-torsion')
    print(torsion(ec, P, order // 1))
    print('-------------------')
    print('r : ')
    rs = find_r(order)
    print(find_r(order))
    print('-------------------')
    print(f'{order // 5}-torsion')
    order5_p = ec.mul(P, 5)
    print(torsion(ec, order5_p, order // 5))
    print('-------------------')
    print(f'{order // 35}-torsion')
    order35_p = ec.mul(P, 35)
    print(torsion(ec, order35_p, order // 35))


def p23():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([100, 905, 0, 1])
    q = 1021
    points = find_points(ec, q)
    order = len(points)
    print('points : ', points)
    print('order : ', order)
    f = Field(q)
    P = PT(f(1006), f(416))
    Q = PT(f(612), f(827))
    print('P2, Q2')
    P_2 = ec.mul(P, 966 / 2)
    assert ec.mul(Q, 483) == ec.mul(P_2, 1)
    print(P_2, ec.mul(Q, 483), ec.mul(P_2, 1))
    for i in range(1, 2 + 1):
        print(i, ec.mul(P_2, i))

    print('P3, Q3')
    P_3 = ec.mul(P, 966 / 3)
    assert ec.mul(Q, 322) == ec.mul(P_3, 3)
    print(P_3, ec.mul(Q, 322), ec.mul(P_3, 3))
    for i in range(1, 3 + 1):
        print(i, ec.mul(P_3, i))
        if ec.mul(Q, 322) == ec.mul(P_3, i):
            break

    print('P7, Q7')
    P_7 = ec.mul(P, 966 / 7)
    assert ec.mul(Q, 138) == ec.mul(P_7, 1)
    print(P_7, ec.mul(Q, 138), ec.mul(P_7, 1))
    for i in range(1, 7 + 1):
        print(i, ec.mul(P_7, i))

    print('P23, Q23')
    P_23 = ec.mul(P, 966 / 23)
    print(P_23, ec.mul(Q, 42), ec.mul(P_23, 1))
    for i in range(1, 23 + 1):
        print(i, ec.mul(P_23, i))
        if ec.mul(Q, 42) == ec.mul(P_23, i):
            break

    k = chinese_remainder([2, 3, 7, 23], [1, 0, 1, 20])
    assert Q == ec.mul(P, k)


def p24():
    print(inspect.stack()[0][3], '>>>>>')
    q = 115792089210356248762697446949407573530086143415290314195533631308867097853951
    r = 115792089210356248762697446949407573529996955224135760342422259061068512044369
    b = 41058363725152142129326129780047268409114441015993725554835256314039467401291
    xG = 48439561293906451759052585252797914202762949526041747995844080717082404635286
    yG = 36134250956749795798585127919587881956611106672985015071877198253568414405109
    xH = 53987601597021778433910548064987973235945515666715026302948657055639179420355
    yH = 53690949263410447908824456005055253553237881490194075871737490561466076234637
    ec = EC([b, -3, 0, 1])
    f = Field(q)
    G = PT(f(xG), f(yG))
    print(ec.mul(G, r - 1))
    for i in range(1, 10):  # for i in range(1, r+1):
        pt = ec.mul(G, i)
        if pt.x == xH and pt.y == yH:
            print(f'k={i}')
            break
        else:
            print(pt)


def p26():
    print(inspect.stack()[0][3], '>>>>>')
    q = 67
    # q^1
    assert frob_end_pi(15, q, 1) == 15
    assert frob_end_pi(50, q, 1) == 50
    # q^2
    poly = FiniteMonoPolynomial([1, 0, 1], q)
    assert poly.pow([16, 2], 67 ** 2) == [16, 2]
    assert frob_end_pi(poly([16, 2]), q, 2) == [16, 2]
    assert frob_end_pi(poly([39, 30]), q, 2) == [39, 30]
    # q^3
    poly = FiniteMonoPolynomial([2, 0, 0, 1], q)
    assert frob_end_pi(poly([8, 4, 15]), q, 3) == [8, 4, 15]


def p27():
    print(inspect.stack()[0][3], '>>>>>')
    P = (-1, 2)

    # (x,y) -> (7x, y)
    q = 19
    ec = EC([5, 0, 0, 1])
    e = find_points(ec, q)
    for p in e:
        if p is not None:
            p.x = p.x * 7 % q
            if not ec.is_on_curve(p):
                print('not exist')

    # endomorphism
    f = Field(q)
    ksi_3 = 7
    ksi = PT(ksi_3 * P[0], f(2))
    print(ksi, ec.is_on_curve(ksi))

    # automorphism
    ksi_2 = 11
    ksi_inv = PT(ksi_2 * ksi.x, f(2))
    ksi_inv.x %= q
    assert ksi_inv.x == P[0] % q
    print((ksi_inv.x - q, ksi_inv.y), ec.is_on_curve(ksi_inv))

    # (x,y) -> ((8u+11)x, y)
    # endomorphism
    q = 23
    poly = FiniteMonoPolynomial([1, 0, 1], q)
    ksi_3 = poly([11, 8])
    ksi = PT(ksi_3 * P[0], poly([2, 0]))
    print(ksi, ec.is_on_curve(ksi))

    # automorphism
    ksi_2 = poly([11, 15])
    ksi_inv = PT(ksi_2 * ksi.x, poly([2, 0]))
    assert ksi_inv.x == P[0] % q
    ksi_inv.x.val_coefs[0] -= q
    print(ksi_inv, ec.is_on_curve(ksi))
    # https://www.mathphysicsbook.com/mathematics/mathematical-structures/defining-mathematical-structures-and-mappings/


def p28():
    print(inspect.stack()[0][3], '>>>>>')
    q = 67
    ec = EC([3, 4, 0, 1])
    points = find_points(ec, q)
    order = len(points)
    t = order - (q + 1)

    f = Field(q)
    P1 = PT(f(15), f(50))
    for t in range(int(-2 * math.sqrt(q)), int(2 * math.sqrt(q))):
        if ec.mul(P1, (1 - (-t) + q)) is None:
            print(f't={t} is found')
            break
    print(ec.mul(P1, (1 - (-t) + q)))
    assert frob_end_pi(P1.x, q, 1) == P1.x
    assert frob_end_pi(P1.y, q, 1) == P1.y
    assert frob_end_pi(P1.x, q, 2) == P1.x
    assert frob_end_pi(P1.y, q, 2) == P1.y

    poly = FiniteMonoPolynomial([1, 0, 1], q)
    P2 = PT(poly([16, 2]), poly([39, 30]))
    print(ec.add(
        ec.mul(PT(frob_end_pi(P2.x, q, 1), frob_end_pi(P2.y, q, 1)), t),
        ec.mul(P2, 68)))

    poly = FiniteMonoPolynomial([2, 0, 0, 1], q)
    P3 = PT(poly([8, 4, 15]), poly([21, 30, 44]))
    print(ec.add(ec.add(
        ec.mul(PT(frob_end_pi(P3.x, q, 2), frob_end_pi(P3.y, q, 2)), 1),
        ec.mul(PT(frob_end_pi(P3.x, q, 1), frob_end_pi(P3.y, q, 1)), t)),
        ec.mul(P3, 67)
    ))


def p29():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([1, 1, 0, 1])
    q = 101
    points = find_points(ec, q)
    order = len(points)
    # print('points : ', points)
    print(f'order = {order} = {get_prime_factors(order)}')
    print(f'2-tortion: {r_torsion(ec, points, 2)}')

    print('ψ2(x) = 4x3 + 4x + 4')
    for x in range(101):
        y = solve_poly([4, 4, 0, 4], x) % q
        if y == 0:
            print(x, end='')
    print('\nψ2(x) = 3x4 +6x2 +12x+100')
    print('ψ3(x) = (x+73)(x+84)(x2 +45x+36)')
    print(f'{(101 - 73) % q}, {(101 - 84) % q}')
    for x in range(101):
        y = solve_poly([100, 12, 6, 0, 3], x) % q
        if y == 0:
            print(x, end=' ')
    print(f'\n2-tortion: {r_torsion(ec, points, 3)}')


def p37():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([20, 20, 0, 1])
    q = 103
    f = Field(q)
    P = PT(f(26), f(20))
    Q = PT(f(63), f(78))
    R = PT(f(59), f(95))
    T = PT(f(77), f(84))

    print('f = 6y+71x2+91x+91 / x2+70x2+11')
    print(f'x2 + {-(P.x + Q.x)}x+{P.x * Q.x}')
    print(f'x2 + {-(R.x + T.x)}x+{R.x * T.x}')

    print('Gp,q>>')
    s = ec.slope(P, Q)
    g_p_q(P, Q, s)

    print('Gr,t>>')
    s = ec.slope(R, T)
    g_p_q(R, T, s)

    print('Gu,u>>')
    U = ec.add(R, T)
    print(U)
    print('Gr,t>>')
    s = ec.slope(U, U)
    g_p_q(U, U, s)

    print('y -y_p -s(x-x_p)>>')
    s = ec.slope(P, Q)
    y_p = P.y
    x_p = P.x
    print(f'y +{-s}x + {-y_p + s * x_p}')  # y -y_p -s(x-x_p)

    print('x-u.x>>')
    print(f'x +{-U.x} ')


def p43():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([-2, -1, 0, 1])
    q = 163
    f = Field(q)
    P = PT(f(43), f(154))
    Q = PT(f(46), f(38))
    y_p = P.y
    x_p = P.x
    # x_q = Q.x
    s = ec.slope(P, Q)
    print(f'y +{-s}x + {-y_p + s * x_p}')  # y -y_p -s(x-x_p)
    y_p = P.y
    x_p = P.x
    # x_q = Q.x
    s = ec.slope(P, P)
    print(f'y +{-s}x + {-y_p + s * x_p}')  # y -y_p -s(x-x_p)


def p48():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([1, 0, 0, 1])
    q = 7691
    points = find_points(ec, q)
    order = len(points)
    # print('points : ', points)
    print(f'order = {order} = {get_prime_factors(order)}')
    # print(points)

    f = Field(q)
    P = PT(f(2693), f(4312))
    print(len(torsion(ec, P, q)))

    poly = FiniteMonoPolynomial([1, 0, 1], q)
    Q = PT(poly([6145, 633]), poly([109, 7372]))
    print(len(torsion(ec, Q, q)))


def p51():
    print(inspect.stack()[0][3], '>>>>>')
    q = 11
    ec = EC([4, 0, 0, 1])
    points = find_points(ec, q)
    order = len(points)
    print('points : ', points)
    print(f'order = {order} = {get_prime_factors(order)}')
    print(f'3-tortion: {r_torsion(ec, points, 3)}')

    poly = FiniteMonoPolynomial([1, 0, 1], q)
    # P = PT(poly([8, 0]), poly([0, 1]))
    print(f'3-tortion: {torsion(ec, PT(poly([8, 0]), poly([0, 1])), q ** 2)}')
    print(f'3-tortion: {torsion(ec, PT(poly([7, 2]), poly([0, 10])), q ** 2)}')
    print(f'3-tortion: {torsion(ec, PT(poly([7, 9]), poly([0, 1])), q ** 2)}')


def p51():
    print(inspect.stack()[0][3], '>>>>>')
    q = 11
    ec = EC([4, 0, 0, 1])
    points = find_points(ec, q)
    order = len(points)
    print('points : ', points)
    print(f'order = {order} = {get_prime_factors(order)}')
    print(f'3-tortion: {r_torsion(ec, points, 3)}')

    poly = FiniteMonoPolynomial([1, 0, 1], q)
    print(f'3-tortion: {torsion(ec, PT(poly([8, 0]), poly([0, 1])), q ** 2)}')
    print(f'3-tortion: {torsion(ec, PT(poly([7, 2]), poly([0, 10])), q ** 2)}')
    print(f'3-tortion: {torsion(ec, PT(poly([7, 9]), poly([0, 1])), q ** 2)}')


def p56():
    print(inspect.stack()[0][3], '>>>>>')
    q = 59
    ec = EC([1, 0, 0, 1])
    points = find_points(ec, q)
    order = len(points)
    print('points : ', points)
    print(f'order = {order} = {get_prime_factors(order)}')
    print(f'5-tortion: {r_torsion(ec, points, 5)}')

    f = Field(q)
    assert frob_end_pi(f(18), q, 1) == f(18)

    poly = FiniteMonoPolynomial([1, 0, 1], q)
    assert frob_end_pi(poly([0, 37]), q, 1) != [0, 37]
    assert frob_end_pi(poly([0, 37]), q, 2) == poly([0, 37])

    # _24i+29, distortion map
    assert poly([29, 24]) * f(18) == [50, 19]
    assert poly([29, 24]) ** 2 * f(18) == [50, 40]
    assert poly([29, 24]) ** 3 * f(18) == [18]

    assert poly([29, 24]) * poly([41, 38]) == [41, 21]
    assert poly([29, 24]) ** 2 * poly([41, 38]) == [36, 0]
    assert poly([29, 24]) ** 3 * poly([41, 38]) == [41, 38]


def p57():
    print(inspect.stack()[0][3], '>>>>>')
    q = 59
    ec = EC([1, 0, 0, 1])
    assert is_supersingular(ec, q) == True
    points = find_points(ec, q)
    order = len(points)
    print('points : ', points)
    print(f'order = {order} = {get_prime_factors(order)}')
    print('r = ', end='')
    for i in get_prime_factors(order):
        print(i[0], end=',')
    print('')
    print(f'5-tortion: {r_torsion(ec, points, 5)}')

    # distortion
    f = Field(q)
    poly = FiniteMonoPolynomial([1, 0, 1], q)
    pt1 = [25, 30]
    pt2 = [poly([34, 0]), poly([0, 30])]
    assert pt1[0] * (-1) % q == pt2[0]
    assert pt1[1] * poly([0, 1]) == pt2[1]

    assert pt2[0] * (-1) == poly(pt1[0])
    assert poly([0, 29]) * poly([0, 1]) == poly([30, 0])
    assert poly([0, 30]) * poly([0, 1]) == poly([29, 0])
    assert poly([0, 28]) * poly([0, 1]) == poly([31, 0])
    assert poly([0, 31]) * poly([0, 1]) == poly([28, 0])


def p61():
    print(inspect.stack()[0][3], '>>>>>')
    q = 11
    ec = EC([4, 0, 0, 1])
    # assert is_supersingular(ec, q) == True
    points = find_points(ec, q)
    order = len(points)
    print('E points : ', points)
    print(f'order = {order} = {get_prime_factors(order)}')
    print('r = ', end='')
    for i in get_prime_factors(order):
        print(i[0], end=',')
    print('')

    q = 11
    ec1 = EC([-4, 0, 0, 1])  # y2=x3+4w6, w2+1=0 w2= -1, w6=(u2)**3=-1
    # TODO EC([4*poly([-1,0]), 0, 0, 1]) # y2=x3+4poly([-1,0])
    points = find_points(ec1, q)
    # print('E\', points : ', points)
    print(f'E\',points : {r_torsion(ec1, points, 3)}')

    f = Field(q)
    poly = FiniteMonoPolynomial([1, 0, 1], q)
    print('\nfind E/q2')
    pt = PT(f(3), f(1))
    for i in range(1, 1000):
        npt = ec1.mul(pt, i)
        if npt is None:
            break
        print(npt, end=' => ')
        print(npt.x * poly([-1, 0]), npt.y * -poly([0, 1]))

    'distortion'
    poly = FiniteMonoPolynomial([1, 0, 1], q)
    pt1 = [3, 10]
    pt2 = [poly([8, 0]), poly([0, 1])]  # twist

    # x=w=[0,1], x2=w2=[0,1]*[0,1]=[-1,0]
    # E' -> E,  (x/w2, y/w3)
    assert pt2[0] * (-1) == poly(pt1[0])  # x/w2
    assert pt2[1] / -poly([0, 1]) == poly([pt1[1], 0])  # y/w3
    # E -> E' (w2x, w3y)
    assert pt1[0] * (-1) % q == pt2[0]  # w2x
    assert pt1[1] * -poly([0, 1]) == pt2[1]  # w3y


def p62():
    print(inspect.stack()[0][3], '>>>>>')
    q = 103
    ec = EC([72, 0, 0, 1])
    # assert is_supersingular(ec, q) == True
    points = find_points(ec, q)
    order = len(points)
    print('points : ', points)
    print(f'order = {order} = {get_prime_factors(order)}')
    print('r = ', end='')
    for i in get_prime_factors(order):
        print(i[0], end=',')
    print('')

    q = 103
    ec1 = EC([-72 * 2, 0, 0, 1])  # y2=x3+72w6, w6+2=0 w6= -2
    # TODO EC([72*2*poly([-1,0,0,0,0,0]), 0, 0, 1 ]) y2=x3+72*2poly([-1,0,0,0,0,0])
    points = find_points(ec1, q)
    # print('E\',points : ', points)
    print(f'E\',points : {r_torsion(ec, points, 7)}')

    f = Field(q)
    poly = FiniteMonoPolynomial([2, 0, 0, 0, 0, 0, 1], q)
    ec2 = EC([72 * 2 * poly([-1, 0, 0, 0, 0, 0]), 0, 0, 1])
    assert ec1.is_on_curve(PT(f(33), f(19)))
    assert ec2.is_on_curve(PT(f(33), f(19)))
    assert ec2.is_on_curve(PT(poly([0, 0, 101, 0, 0, 0]), poly([0, 0, 0, 8, 0, 0])))

    print('\nfind E/q6')
    pt = PT(f(33), f(19))
    for i in range(1, 1000):
        npt = ec1.mul(pt, i)
        if npt is None:
            break
        print(npt, end=' => ')
        print(npt.x / poly([0, 0, 1, 0, 0, 0]), npt.y / poly([0, 0, 0, 1, 0, 0]))

    print('\ndistortion')
    pt = PT(poly([0, 0, 0, 0, 35, 0]), poly([0, 0, 0, 42, 0, 0]))
    assert ec.is_on_curve(pt)
    for i in range(1, 1000):
        npt = ec.mul(pt, i)
        if npt is None:
            break
        print(npt)
        # x2=w2=[0,0,1,0,0,0], x3=w3=[0,0,0,1,0,0],
        # E -> E' (w2x, w3y)
        x_1 = npt.x * poly([0, 0, 1, 0, 0, 0])
        y_1 = npt.y * poly([0, 0, 0, 1, 0, 0])
        print('E\' ', x_1, y_1)
        # E' -> E,  (x/w2, y/w3)
        print('E  ', x_1.val_coefs[0] / poly([0, 0, 1, 0, 0, 0]), y_1.val_coefs[0] / poly([0, 0, 0, 1, 0, 0]))
    # not use b in calculating mul,add in y2=x3+ax+b


def p69():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([0, -1, 0, 1])
    q = 23
    f = Field(q)
    poly = FiniteMonoPolynomial([1, 0, 1], q)

    P = PT(f(2), f(11))
    Q = PT(poly([21, 0]), poly([0, 12]))
    S = PT(poly([18, 10]), poly([13, 13]))
    assert ec.weil_pairing(P, Q, S, 3) == [11, 15]
    assert ec.weil_pairing(ec.mul(P, 2), Q, S, 3) == [11, 8]
    assert ec.weil_pairing(P, ec.mul(Q, 2), S, 3) == [11, 8]
    # R = PT(poly([0, 17]), poly([21, 2]))
    # print(ec.weil_pairing2(P, Q, S, R, 3))


def p73():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([-3, 0, 0, 1])
    q = 5
    poly = FiniteMonoPolynomial([2, 0, 1], q)
    f = Field(q)
    P = PT(f(3), f(2))
    Q = PT(poly([1, 1]), poly([2, 4]))
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

    # (q^k-1)/r
    e = (q ** 2 - 1) // 3
    # tr -> Tr
    assert ec.tate_pairing(P, Q, R, m) ** 2 ** e == poly([2, 4])
    assert ec.tate_pairing(P, ec.mul(Q, 2), R, m) ** e == poly([2, 4])
    assert ec.tate_pairing(ec.mul(P, 2), Q, R, m) ** e == poly([2, 4])


def p78():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([15, 21, 0, 1])
    q = 47
    poly = FiniteMonoPolynomial([5, 0, -4, 0, 1], q)
    f = Field(q)
    P = PT(f(45), f(23))
    Q = PT(poly([29, 0, 31, 0]), poly([0, 11, 0, 35]))
    # Q=R
    R = PT(poly([29, 0, 31, 0]), poly([0, 11, 0, 35]))

    m = 17
    ms = [int(x) for x in bin(m)[2:]][::-1]
    assert ec.tate_pairing(P, Q, R, m) == [22, 10, 6, 17]
    # 시간많이걸림
    # Tr
    assert ec.tate_pairing(P, Q, R, m) ** 287040 == [39, 45, 43, 33]


def p82():
    print(inspect.stack()[0][3], '>>>>>')
    # ec = EC([14, 0, 0, 1])
    q = 4219433269001672285392043949141038139415112441532591511251381287775317505016692408034796798044263154903329667
    r = 2236970611075786789503882736578627885610300038964062133851391137376569980702677867
    print(q.bit_length())
    # print(len(bin(q))-2)
    print(r.bit_length())
    k = 9
    assert (q ** k - 1) % r == 0
    print('rho : ', math.log10(q) / math.log10(r))


def p86():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([314159, 0, 0, 1])
    q = 425708693169757088196017853607835113595127103859429924930531263283244403251872949802982860038531930965867890444658222153407204383584492024637762799391807569669124814253270947366226515064812665901907204494611177526596015257984009814596057160388672298355821309046798841446111721495601835913381835880170934319890420895521320439930666405003725309562669243847766834546592867695533445054256132471093279787853214492986394176521193456205570309658462204234557728373615304193316916440130004424612327
    r = 57896044618658097711785492504343953926634992332820282019728792003956564820063
    # assert is_supersingular(ec, q) # overflow
    print(q.bit_length())
    print(r.bit_length())
    print(q.bit_length() - r.bit_length())
    print('rho : ', math.log10(q) / math.log10(r))


def p88():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([314159, 0, 0, 1])
    q = 425708693169757088196017853607835113595127103859429924930531263283244403251872949802982860038531930965867890444658222153407204383584492024637762799391807569669124814253270947366226515064812665901907204494611177526596015257984009814596057160388672298355821309046798841446111721495601835913381835880170934319890420895521320439930666405003725309562669243847766834546592867695533445054256132471093279787853214492986394176521193456205570309658462204234557728373615304193316916440130004424612327
    r = 57896044618658097711785492504343953926634992332820282019728792003956564820063
    # assert is_supersingular(ec, q) # overflow
    print(q.bit_length())
    print(r.bit_length())
    print(q.bit_length() - r.bit_length())
    print('rho : ', math.log10(q) / math.log10(r))
    # k = 9
    # assert (q ** k - 1) % r == 0
    # rho = math.log10(q) / math.log10(r)
    # print('rho : ', rho)


p28()

# TODO p32
'''
https://en.wikipedia.org/wiki/Schoof%27s_algorithm
https://github.com/pdinges/python-schoof
'''
