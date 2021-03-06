import sys
from fractions import Fraction
from pyxmath import *
from functools import reduce


def xgcd(a, b):
    """
    a * x + b * y = gcd
    """
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = b, a

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    return old_r, old_s, old_t


def inv(n, p):
    # if sys.version_info.major >= 3 and sys.version_info.minor >= 8:
    #     return pow(n, -1, p)
    if n == 0:
        raise ZeroDivisionError('division by zero')
    gcd, x, y = xgcd(n, p)
    if gcd != 1:
        raise ValueError(
            f'{n} is not invertible for the given modulus {p}')
    else:
        return x % p


def rmod(a, p):
    if round(a) == a:
        if a == 0:
            return 0
        return a % p
    n = Fraction(str(a))
    return n.numerator * inv(n.denominator, p) % p


class Field():
    def __init__(self, q, n=None):
        self.q = q
        self.n = n

    def __repr__(self):
        return repr(self.n)

    def __call__(self, a):
        return self.__class__(self.q, a % self.q)

    def add(self, a, b):
        return (a + b) % self.q

    def __add__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.q, (self.n + other.n) % self.q)
        else:
            if type(other) != int:
                return NotImplemented
            return self.__class__(self.q, (self.n + other) % self.q)

    __radd__ = __add__

    def sub(self, a, b):
        return (a - b) % self.q

    def __sub__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.q, (self.n - other.n) % self.q)
        else:
            if type(other) != int:
                return NotImplemented
            return self.__class__(self.q, (self.n - other) % self.q)

    def __rsub__(self, other):
        if type(other) != int:
            return NotImplemented
        return self.__class__(self.q, (other - self.n) % self.q)

    def mul(self, a, b):
        return (a * b) % self.q

    def __mul__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.q, (self.n * other.n) % self.q)
        else:
            if type(other) != int:
                return NotImplemented
            return self.__class__(self.q, (self.n * other) % self.q)

    __rmul__ = __mul__

    def div(self, a, b):
        return (a * inv(b, self.q)) % self.q

    def __truediv__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.q, (self.n * inv(other.n, self.q)) % self.q)
        else:
            if type(other) != int:
                return NotImplemented
            return self.__class__(self.q, (self.n * inv(other, self.q)) % self.q)

    def __rtruediv__(self, other):
        if type(other) != int:
            return NotImplemented
        return self.__class__(self.q, (other * inv(self.n, self.q)) % self.q)

    def inv(self, a):
        return inv(a, self.q)

    def pow(self, a, exp):
        return (a ** exp) % self.q

    def __pow__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.q, (self.n ** other.n) % self.q)
        else:
            return self.__class__(self.q, (self.n ** other) % self.q)

    def __mod__(self, other):
        return self.__class__(self.q, self.n % other)

    def neg(self, a):
        return (a * (-1)) % self.q

    def __neg__(self):
        return self.__class__(self.q, (self.n * (-1)) % self.q)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.n == other.n
        else:
            return self.n == other


# Frobenius endomorphism ??
# https://en.wikipedia.org/wiki/Frobenius_endomorphism
# https://ko.wikipedia.org/wiki/%ED%94%84%EB%A1%9C%EB%B2%A0%EB%8B%88%EC%9A%B0%EC%8A%A4_%EC%82%AC%EC%83%81
def frob_end_pi(r, q, exp=1):
    return r ** (q ** exp) % q


def modular_sqrt_for_square_of_prime(a, p):
    p_2 = p
    f = get_prime_factors(p_2)
    if len(f) > 1:
        raise ValueError(f'{p} is not prime or square of prime')
    if f[0][1] == 1:
        return modular_sqrt(a, p)
    elif f[0][1] == 2:
        p_1 = f[0][0]
        mod_p_1 = a % p_1
        mod_p_2 = a % p_2
        b_1 = modular_sqrt(mod_p_1, p_1)
        if p_1 - b_1 < b_1:
            b_1 = p_1 - b_1
        b_2 = mod_p_2 % p_2

        lhs = (2 * b_1 * p_1) / p_1
        rhs = (b_2 - b_1 ** 2) / p_1
        try:
            m = inv(lhs, p_1) * rhs % p_1
            return int(b_1 + m * p_1)
        except Exception as e:
            raise ValueError(
                f'{a} is not invertible for the given modulus {p}')
    else:
        raise ValueError(f'{p} is not prime or square of prime')


# https://gist.github.com/nakov/60d62bdf4067ea72b7832ce9f71ae079
# Tonelli???Shanks
# todo study later
def modular_sqrt(a, p):
    def legendre_symbol(a, p):
        """ Compute the Legendre symbol a|p using
            Euler's criterion. p is a prime, a is
            relatively prime to p (if p divides
            a, then a|p = 0)

            Returns 1 if a has a square root modulo
            p, -1 otherwise.
        """
        # https://brilliant.org/wiki/legendre-symbol/
        # https://ko.wikipedia.org/wiki/????????????_??????
        ls = pow(a, (p - 1) // 2, p)
        return -1 if ls == p - 1 else ls

    """ Find a quadratic residue (mod p) of 'a'. p
        must be an odd prime.

        Solve the congruence of the form:
            x^2 = a (mod p)
        And returns x. Note that p - x is also a root.

        0 is returned is no square root exists for
        these a and p.

        The Tonelli-Shanks algorithm is used (except
        for some simple cases in which the solution
        is known from an identity). This algorithm
        runs in polynomial time (unless the
        generalized Riemann hypothesis is false).
    """
    # Simple cases
    #
    if legendre_symbol(a, p) != 1:
        return 0
    elif a == 0:
        return 0
    elif p == 2:
        return p
    elif p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    # Q2^s => s2^e
    # Partition p-1 to s * 2^e for an odd s (i.e.
    # reduce all the powers of 2 from p-1)
    #
    s = p - 1
    e = 0
    while s % 2 == 0:
        s //= 2
        e += 1

    # Find some 'n' with a legendre symbol n|p = -1.
    # Shouldn't take long.
    #
    n = 2  # z
    while legendre_symbol(n, p) != -1:
        n += 1

    # Here be dragons!
    # Read the paper "Square roots from 1; 24, 51,
    # 10 to Dan Shanks" by Ezra Brown for more
    # information
    #

    # x is a guess of the square root that gets better
    # with each iteration.
    # b is the "fudge factor" - by how much we're off
    # with the guess. The invariant x^2 = ab (mod p)
    # is maintained throughout the loop.
    # g is used for successive powers of n to update
    # both a and b
    # r is the exponent - decreases with each update
    #
    x = pow(a, (s + 1) // 2, p)  # R
    b = pow(a, s, p)  # t
    g = pow(n, s, p)  # c
    r = e  # S

    while True:
        t = b
        m = 0
        for m in range(r):
            if t == 1:
                break
            t = pow(t, 2, p)

        if m == 0:
            return x

        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m


# https://medium.com/analytics-vidhya/chinese-remainder-theorem-using-python-25f051e391fc
def chinese_remainder(m, a):
    sum = 0
    prod = reduce(lambda acc, b: acc * b, m)
    for n_i, a_i in zip(m, a):
        p = prod // n_i
        sum += a_i * inv(p, n_i) * p
    return sum % prod
