# https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjhiuXAzbzrAhWzL6YKHQohB1oQFjAAegQIBBAB&url=http%3A%2F%2Fwww.craigcostello.com.au%2Fpairings%2FPairingsForBeginners.pdf&usg=AOvVaw1H5dLtelG00vWsvWRGxBNZ
# ParingsForBeginners.pdf
import inspect
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
    p = 1021
    points = find_points(ec, p)
    order = len(points)
    print('points : ', points)
    print('order : ', order)


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
    # E_3
    # q distortion map
    q = 19
    f = Field(q)
    assert f(7) == 7
    assert f(7) * f(7) == 11
    assert f(7) * f(7) * f(7) == 1

    # q^2 distortion map
    q = 23
    poly = FiniteMonoPolynomial([1, 0, 1], q)
    assert poly([11, 8]) * poly([11, 8]) == [11, 15]
    assert poly([11, 8]) * poly([11, 8]) * poly([11, 8]) == [1, 0]


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
    ec = EC([-4, 0, 0, 1])  # y2=x3+4w6, w2+1=0 w2= -1, w6=(u2)**3=-1
    # TODO EC([4*poly([-1,0]), 0, 0, 1]) # y2=x3+4poly([-1,0])
    points = find_points(ec, q)
    print('E\', points : ', points)

    # distortion
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
    print('E\',points : ', points)
    print(len(points))
    f = Field(q)
    pt = PT(f(33), f(19))
    for i in range(1, 1000):
        npt = ec1.mul(pt, i)
        if npt is None:
            break
        print(npt)

    poly = FiniteMonoPolynomial([2, 0,0,0,0,0, 1], q)
    ec2 = EC([72*2*poly([-1,0,0,0,0,0]), 0, 0, 1 ])
    assert ec2.is_on_curve(PT(f(33), f(19)))
    assert ec2.is_on_curve(PT(poly([0, 0, 101, 0, 0, 0]), poly([0, 0, 0, 8, 0, 0])))

    # distortion
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


p62()

# TODO p32
'''
https://en.wikipedia.org/wiki/Schoof%27s_algorithm
https://github.com/pdinges/python-schoof
'''
