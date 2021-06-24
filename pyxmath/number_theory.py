import sys


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


class Field():
    def __init__(self, p, n=None):
        self.p = p
        self.n = n

    def __repr__(self):
        return repr(self.n)

    def __call__(self, a):
        return self.__class__(self.p, a % self.p)

    def add(self, a, b):
        return (a + b) % self.p

    def __add__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.p, (self.n + other.n) % self.p)
        else:
            return self.__class__(self.p, (self.n + other) % self.p)

    __radd__ = __add__

    def sub(self, a, b):
        return (a - b) % self.p

    def __sub__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.p, (self.n - other.n) % self.p)
        else:
            return self.__class__(self.p, (self.n - other) % self.p)

    def __rsub__(self, other):
        return self.__class__(self.p, (other - self.n) % self.p)

    def mul(self, a, b):
        return (a * b) % self.p

    def __mul__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.p, (self.n * other.n) % self.p)
        else:
            return self.__class__(self.p, (self.n * other) % self.p)

    __rmul__ = __mul__

    def div(self, a, b):
        return (a * inv(b, self.p)) % self.p

    def __truediv__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.p, (self.n * inv(other.n, self.p)) % self.p)
        else:
            return self.__class__(self.p, (self.n * inv(other, self.p)) % self.p)

    def __rtruediv__(self, other):
        return self.__class__(self.p, (other * inv(self.n, self.p)) % self.p)

    def inv(self, a):
        return inv(a, self.p)

    def pow(self, a, exp):
        return (a ** exp) % self.p

    def __pow__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.p, (self.n ** other.n) % self.p)
        else:
            return self.__class__(self.p, (self.n ** other) % self.p)

    def neg(self, a):
        return (a * (-1)) % self.p

    def __neg__(self):
        return self.__class__(self.p, (self.n * (-1)) % self.p)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.n == other.n
        else:
            return self.n == other


# Frobenius endomorphism π
# https://en.wikipedia.org/wiki/Frobenius_endomorphism
# https://ko.wikipedia.org/wiki/%ED%94%84%EB%A1%9C%EB%B2%A0%EB%8B%88%EC%9A%B0%EC%8A%A4_%EC%82%AC%EC%83%81
def frob_end_pi(r, q, exp=1):
    return r ** (q ** exp) % q


# https://gist.github.com/nakov/60d62bdf4067ea72b7832ce9f71ae079
# Tonelli–Shanks
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
    n = 2
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
    x = pow(a, (s + 1) // 2, p)
    b = pow(a, s, p)
    g = pow(n, s, p)
    r = e

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
