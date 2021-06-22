# https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjhiuXAzbzrAhWzL6YKHQohB1oQFjAAegQIBBAB&url=http%3A%2F%2Fwww.craigcostello.com.au%2Fpairings%2FPairingsForBeginners.pdf&usg=AOvVaw1H5dLtelG00vWsvWRGxBNZ
# ParingsForBeginners.pdf
import inspect
from pyxmath.elliptic_curve import *


def find_r(order):
    rs = []
    for i in range(1, order + 1):
        if order % i == 0:
            rs.append(i)
    return rs


def torsion(ec, p, r):
    print(f'{r}-torsion')
    for i in range(1, r + 1):
        r = ec.mul(p, i)
        if r is None:
            print(None)
            print(f'\r{p} order :', i)
        else:
            print(r, end=" ")


def p22():
    print(inspect.stack()[0][3], '>>>>>')
    ec = EC([1, 1, 0, 1])
    p = 101
    points = find_points(ec, p)
    order = len(points)
    print('points : ', points)
    print('order : ', order)
    f = Field(101)
    print('-------------------')
    # P is a generator
    P = PT(f(47), f(12))
    torsion(ec, P, order // 1)
    print('-------------------')
    print('r : ')
    rs = find_r(order)
    print(find_r(order))
    print('-------------------')
    order5_p = ec.mul(P, 5)
    torsion(ec, order5_p, order // 5)
    print('-------------------')
    order35_p = ec.mul(P, 35)
    torsion(ec, order35_p, order // 35)


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
    print(ec.mul(G, r-1))
    for i in range(1, 10): #for i in range(1, r+1):
        pt = ec.mul(G, i)
        if pt.x == xH and pt.y == yH:
            print(f'k={i}')
            break
        else:
            print(pt)

p24()
