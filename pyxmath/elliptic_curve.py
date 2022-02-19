import math
from pyxmath import *
from pyxmath.number_theory import *


class EC:
    def __init__(self, coefs, coordinates='affine'):
        self.coefs = coefs
        self.coordinates = coordinates

    def neg(self, p1):
        return -p1

    def add(self, p1, p2):
        if self.coordinates == 'jacobian':
            return self._add_at_jacobian(p1, p2)
        elif self.coordinates == 'projective':
            return self._add_at_projective(p1, p2)

        if p1 is None or p2 is None:
            return p1 if p2 is None else p2
        if p1 is None or p2 is None:
            return p1 if p2 is None else p2

        if p2.x == p1.x and p2.y == p1.y:
            return self.double(p1)
        elif p2.x == p1.x and p2.y == (-p1.y):
            return None
        else:
            try:
                l = (p2.y - p1.y) / (p2.x - p1.x)
            except (ValueError, ZeroDivisionError):
                return None
        x = l ** 2 - p1.x - p2.x
        y = -l * x + l * p1.x - p1.y
        return PT(x, y)

    def _add_at_projective(self, p1, p2):
        if p1 is None or p2 is None:
            return p1 if p2 is None else p2
        if p1 is None or p2 is None:
            return p1 if p2 is None else p2

        if p2.x == p1.x and p2.y == p1.y and p2.z == p1.z:
            return self._double_at_projective(p1)
        elif p2.x == p1.x and p2.y == (-p1.y):
            return None
        else:
            u = p2.y * p1.z - p1.y * p2.z
            v = p2.x * p1.z - p1.x * p2.z
            w = u ** 2 * p1.z * p2.z - v ** 3 - 2 * v ** 2 * p1.x * p2.z
            x = v * w
            y = u * (v ** 2 * p1.x * p2.z - w) - v ** 3 * (p1.y * p2.z)
            z = v ** 3 * p1.z * p2.z
            return PT(x, y, z)

    def _add_at_jacobian(self, p1, p2):
        if p1 is None or p2 is None:
            return p1 if p2 is None else p2
        if p1 is None or p2 is None:
            return p1 if p2 is None else p2

        if p2.x == p1.x and p2.y == p1.y and p2.z == p1.z:
            return self._double_at_jacobian(p1)
        elif p2.x == p1.x and p2.y == (-p1.y):
            return None
        else:
            r = p1.x * p2.z ** 2
            s = p2.x * p1.z ** 2
            t = p1.y * p2.z ** 3
            u = p2.y * p1.z ** 3
            v = s - r
            w = u - t
            x = -v ** 3 - 2 * r * v ** 2 + w ** 2
            y = -t * v ** 3 + (r * v ** 2 - x) * w
            z = v * p1.z * p2.z
            return PT(x, y, z)

    def sub(self, p1, p2):
        p2_copy = p2.copy()
        p2_copy.y = -p2_copy.y
        return self.add(p1, p2_copy)

    def mul(self, p1, n):
        if p1 is None:
            return None
        elif n == 0:
            return None
        elif n == 1:
            return p1
        elif not n % 2:
            return self.mul(self.double(p1), n // 2)
        else:
            return self.add(self.mul(self.double(p1), int(n // 2)), p1)

    def double(self, p1):
        if self.coordinates == 'jacobian':
            return self._double_at_jacobian(p1)
        elif self.coordinates == 'projective':
            return self._double_at_projective(p1)
        a = self.coefs[1]
        y2 = 2 * p1.y
        if y2 == 0:
            return None
        l = (3 * p1.x ** 2 + a) / y2
        x = l ** 2 - 2 * p1.x
        y = -l * x + l * p1.x - p1.y
        return PT(x, y)

    def _double_at_projective(self, p1):
        a = self.coefs[1]
        t = a * p1.z ** 2 + 3 * p1.x ** 2
        u = p1.y * p1.z
        v = u * p1.x * p1.y
        w = t ** 2 - 8 * v
        x = 2 * u * w
        y = t * (4 * v - w) - 8 * p1.y ** 2 * u ** 2
        z = 8 * u ** 3
        return PT(x, y, z)

    def _double_at_jacobian(self, p1):
        a = self.coefs[1]
        v = 4 * p1.x * p1.y ** 2
        w = 3 * p1.x ** 2 + a * p1.z ** 4
        x = -2 * v + w ** 2
        y = -8 * p1.y ** 4 + (v - x) * w
        z = 2 * p1.y * p1.z
        return PT(x, y, z)

    def find_y_square(self, p_x):
        y = 0
        for deg, c in enumerate(self.coefs):
            y = y + c * p_x ** deg
        return y

    def find_y_square2(self, x, y, ef):
        points = []
        for sqrt_y in ef:
            c = sqrt_y ** 2
            # if c % x.q == y:
            if c == y:
                points.append((x, sqrt_y))
        return points

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

    def miller_g(self, p, q, g):
        s = self.slope(p, q)
        if s is None:
            return g.x - p.x
        else:
            return (g.y - p.y - s * (g.x - p.x)) / (g.x + p.x + q.x - s ** 2)

    def miller(self, p, g, ms):
        t = p
        f = 1
        for i in range(len(ms) - 2, -1, -1):
            f = f ** 2 * self.miller_g(t, t, g)
            t = self.mul(t, 2)  # t = 2*t
            if ms[i] == 1:
                f = f * self.miller_g(t, p, g)
                t = self.add(t, p)  # t = t + p
        return f

    def weil_pairing(self, p, q, s, m):
        ms = [int(x) for x in bin(m)[2:]][::-1]
        # fp(q+s)/fp(s)  / fq(p-s)/fq(-s)
        return (self.miller(p, self.add(q, s), ms) / self.miller(p, s, ms)) / (
                self.miller(q, self.sub(p, s), ms) / self.miller(q, -s, ms))

    def weil_pairing2(self, p, q, s, r, m):
        ms = [int(x) for x in bin(m)[2:]][::-1]
        # fp(q+s)/fp(s)  / fq(p+r)/fq(r)
        return (self.miller(p, self.add(q, s), ms) / self.miller(p, s, ms)) / (
                self.miller(q, self.add(p, r), ms) / self.miller(q, r, ms))

    def tate_pairing(self, p, q, s, m):
        ms = [int(x) for x in bin(m)[2:]][::-1]
        # fp(q+s)/fp(s)
        return self.miller(p, self.add(q, s), ms) / self.miller(p, s, ms)

    def is_on_curve(self, p1):
        if p1 is None:
            return True
        return p1.y ** 2 - p1.x ** 3 - self.coefs[0] == 0


# https://en.wikipedia.org/wiki/Supersingular_elliptic_curve
def is_supersingular(ec, p):
    f = get_prime_factors(p)
    if len(f) > 1 or f[0][1] > 1:
        return False
    if p % 3 == 2:
        if ec.coefs[1] == 0:
            return True
    if p % 4 == 3:
        if ec.coefs[0] == 0:
            return True
    return False
    # points = find_points(ec, p)
    # if len(points) == p + 1:
    #     return True
    # else:
    #     return False


def find_points(ec, q):
    ps = []
    ps.append(None)
    f = Field(q)
    for i in range(q):
        iq = ec.find_y_square(i) % q
        if iq == 0:
            ps.append(PT(f(i), f(0)))
            factors = get_prime_factors(q)
            if factors[0][1] > 1:
                for j in range(1, q):
                    if j * j % q == 0:
                        ps.append(PT(f(0), f(j)))
        else:
            try:
                r = modular_sqrt_for_square_of_prime(iq, q)
                if r != 0:
                    ps.append(PT(f(i), f(r)))
                    ps.append(PT(f(i), f(q - r)))
            except (ValueError,):
                pass
    return ps


# TODO QR, find_points2
'''
https://math.stackexchange.com/questions/3204686/polynomial-evaluates-to-quadratic-residue-in-p-cases
https://crypto.stanford.edu/pbc/notes/numbertheory/qr.html
https://pypi.org/project/labmath/
https://web.northeastern.edu/dummit/teaching_sp20_3527/3527_lecture_32_applications_of_quadratic_reciprocity.pdf
'''


# https://en.wikipedia.org/wiki/Quadratic_residue
# https://rkm0959.tistory.com/20
# https://eli.thegreenplace.net/2009/03/07/computing-modular-square-roots-in-python
# https://gist.github.com/nakov/60d62bdf4067ea72b7832ce9f71ae079
# x^2  ~ q (mod n)
# def qr(x_2, mod):
#     if x_2 >= mod:
#         x_2 = x_2 % mod
#     results = []
#     for x in range(mod):
#         if x ** 2 % mod == x_2:
#             results.append(x)
#     return results
# def find_points2(ec, polys):
#     ps = []
#     ps.append(None)
#     for i in polys:
#         iq = ec.find_y_square(i)
#         print(qr(iq, i.q))
#         ps.append(PT(i, iq))
#         # # 0**2 = 0
#         # if iq == 0:
#         #     ps.append(PT(f(i), f(0)))
#         # else:
#         #     r = modular_sqrt(iq, q)
#         #     if r != 0:
#         #         ps.append(PT(f(i), f(r)))
#         #         ps.append(PT(f(i), f(q - r)))
#     return ps


def torsion(ec, p, r):
    rs = []
    for i in range(1, r + 1):
        r = ec.mul(p, i)
        if r is None:
            rs.append(None)
            break
        else:
            rs.append(r)
    return rs


def r_torsion(ec, ps, r):
    rs = []
    for p in ps:
        for j in range(1, r + 1):
            if j == 1 and p is None:
                break
            rp = ec.mul(p, j)
            if j == r:
                if rp is None:
                    rs.append(rp)
                    return rs
                rs = []
                break
            else:
                rs.append(rp)
