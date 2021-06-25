from pyxmath import *
from pyxmath.number_theory import *


class EC:
    def __init__(self, coefs):
        self.coefs = coefs

    def neg(self, p1):
        return -p1

    def add(self, p1, p2):
        if p1 is None or p2 is None:
            return p1 if p2 is None else p2
        if p1 is None or p2 is None:
            return p1 if p2 is None else p2

        if p2.x == p1.x and p2.y == p1.y:
            return self.double(p1)
        # TODO: 이상하다 p2.y이로 비교해야 되는데
        elif p2.x == p1.x:
            return None
        else:
            try:
                l = (p2.y - p1.y) / (p2.x - p1.x)
            except (ValueError, ZeroDivisionError):
                return None
        x = l ** 2 - p1.x - p2.x
        y = -l * x + l * p1.x - p1.y
        return PT(x, y)

    def sub(self, p1, p2):
        p2_copy = p2.copy()
        p2_copy.y = -p2_copy.y
        return self.add(p1, p2_copy)

    def mul(self, p1, n):
        if n == 0:
            return None
        elif n == 1:
            return p1
        elif not n % 2:
            return self.mul(self.double(p1), n // 2)
        else:
            return self.add(self.mul(self.double(p1), int(n // 2)), p1)

    def double(self, p1):
        a = self.coefs[1]
        y2 = 2 * p1.y
        if y2 == 0:
            return None
        l = (3 * p1.x ** 2 + a) / y2
        x = l ** 2 - 2 * p1.x
        y = -l * x + l * p1.x - p1.y
        return PT(x, y)

    def find_y_square(self, p_x):
        y = 0
        for deg, c in enumerate(self.coefs):
            y = y + c * p_x ** deg
        return y

    def slope(self, p1, p2):
        a = self.coefs[1]
        if p1 is None or p2 is None:
            return p1 if p2 is None else p2

        if p1.x == p2.x:
            if p1.y == -p2.y:
                return None
            return (3 * p1.x ** 2 + a) / (2 * p2.y)
        else:
            try:
                return (p2.y - p1.y) / (p2.x - p1.x)
            except (ValueError, ZeroDivisionError):
                return None

def find_points(ec, p):
    ps = []
    ps.append(None)
    for i in range(p):
        r = modular_sqrt(ec.find_y_square(i), p)
        if r != 0:
            ps.append([i, r])
            ps.append([i, p - r])
    return ps
