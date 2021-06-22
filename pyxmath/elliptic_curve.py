from pyxmath.number_theory import *


class PT():
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        else:
            return self.x == other[0] and self.y == other[1]


class EC:
    def __init__(self, coefs):
        self.coefs = coefs

    def add(self, p1, p2):
        if p1 is None or p2 is None:
            return p1 if p2 is None else p2
        if p1 is None or p2 is None:
            return p1 if p2 is None else p2

        if p2.x == p1.x and p2.y == p1.y:
            return self.double(p1)
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


def find_points(ec, p):
    ps = []
    ps.append([0, 0])
    for i in range(p):
        r = modular_sqrt(ec.find_y_square(i), p)
        if r is not 0:
            ps.append([i, r])
            ps.append([i, p - r])
    return ps
