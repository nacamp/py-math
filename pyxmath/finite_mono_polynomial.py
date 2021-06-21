import itertools
from pyxmath import number_theory as nth

# class FiniteMonoPolynomialError(Exception):
#     def __init__(self, message):
#         self.message = message


class FiniteMonoPolynomial():
    @staticmethod
    def neg(a):
        return [x * (-1) for x in a]

    def __init__(self, coefs, mod):
        '''

        Parameters
        ----------
        :param coef: list x^3+2x+1=0 : [1,2,0,1]
        :param mod: int
        '''

        self.coefs = coefs
        self.deg = len(coefs) - 1
        self.mod = mod
        self.fields = []
        # ex: x^3+2x+1=0 (mod3) : x^3=0x^2-2x-1 : [2,1,0] :  left_coef=1, right_coefs=[2,1,0]
        self.right_coefs = coefs[:]
        self.right_coefs.pop()
        self.right_coefs = [x * (-1) % self.mod for x in self.right_coefs]

    def add(self, r_a, r_b):
        '''add

        :param r_a: right_coefs
        :param r_b: right_coefs
        :return:
        '''
        self._validate(r_a, r_b)
        return [sum(x) % self.mod for x in itertools.zip_longest(r_a, r_b, fillvalue=0)]

    def sub(self, r_a, r_b):
        '''sub

        :param r_a: right_coefs
        :param r_b: right_coefs
        :return:
        '''
        self._validate(r_a, r_b)
        aa = list(r_a[:])
        bb = list(r_b[:])
        for i in range(max(len(aa), len(bb))):
            if len(aa) == i:
                aa.append(0)
            if len(bb) == i:
                bb.append(0)
            aa[i] = (aa[i] - bb[i]) % self.mod
        return aa

    def mul(self, r_a, r_b):
        '''mul

        :param r_a: right_coefs
        :param r_b: right_coefs
        :return:
        '''
        self._validate(r_a, r_b)
        r_coefs_len = len(self.right_coefs)
        over_coefs = self._mul_coef(r_a, r_b)
        coefs = over_coefs[0:r_coefs_len]
        l = len(over_coefs)
        while l > r_coefs_len:
            over_coefs = self._mul_coef(self.right_coefs, over_coefs[r_coefs_len:])
            l = len(over_coefs)
            coefs = [sum(x) for x in itertools.zip_longest(over_coefs[0:r_coefs_len], coefs, fillvalue=0)]
        return ([x % self.mod for x in coefs])

    def div(self, r_a, r_b):
        self._validate(r_a, r_b)
        return self.mul(r_a, self.inv(r_b))

    def poly_round_div(self, numerator, denominator):
        # n//d = q + r
        quotient = []
        n = numerator[:]
        while len(n) >= len(denominator):
            diff_deg = len(n) - len(denominator)
            d = [0] * diff_deg + denominator
            c = (n[-1] * nth.inv(d[-1], self.mod)) % self.mod
            quotient.insert(0, c)
            d = [x * (-c) for x in d]
            n = [sum(x) for x in itertools.zip_longest(n, d, fillvalue=0)]
            n.pop()
        return quotient

    def xgcd(self, a):
        """
        a * x + b * y = gcd
        """
        s, old_s = [0], [1]
        t, old_t = [1], [0]
        r, old_r = self.coefs, a

        while sum(r) != 0:
            quotient = self.poly_round_div(old_r, r)
            old_r, r = r, [sum(x) % self.mod for x in
                           itertools.zip_longest(old_r, [x * (-1) for x in self._mul_coef(quotient, r)], fillvalue=0)]
            old_s, s = s, [sum(x) % self.mod for x in
                           itertools.zip_longest(old_s, [x * (-1) for x in self._mul_coef(quotient, s)], fillvalue=0)]
            old_t, t = t, [sum(x) % self.mod for x in
                           itertools.zip_longest(old_t, [x * (-1) for x in self._mul_coef(quotient, t)], fillvalue=0)]
            while len(r) and r[-1] == 0:
                r.pop()
        return old_r, old_s, old_t
        # # old_r[0]이 1이 아닌경우는 old_r[0]으로 나눠야 한다.
        # old_s_inv = mul_inverse_mod(old_r[0], self.mod)
        # # result = [ x % 3 for x in self.poly_mul2(old_s, [old_s_inv])]
        # result = [x % self.mod for x in self._mul_coef(old_s, [old_s_inv])]
        # return result + ([0] * (len(self.irr_coef) - len(result)))

    def inv(self, a):
        gcd, s, t = self.xgcd(a)
        # gcd[0]이 1이 아닌경우는 old_r[0]으로 나눠야 한다.
        s_inv = nth.inv(gcd[0], self.mod)
        result = [x % self.mod for x in self._mul_coef(s, [s_inv])]
        return result + ([0] * (len(self.right_coefs) - len(result)))


    def _mul_coef(self, r_a, r_b):
        aa = list(r_a[:])
        bb = list(r_b[:])
        while len(aa) > 0 and aa[-1] == 0:
            aa.pop()
        while len(bb) > 0 and bb[-1] == 0:
            bb.pop()
        r = [0] * (len(aa) + len(bb) - 1)
        for i, c in enumerate(bb):
            for j in range(len(aa)):
                r[i + j] += c * aa[j]
        return r

    def _validate(self, r_a, r_b):
        if isinstance(r_a, int):
            r_a = [r_a]
        if isinstance(r_b, int):
            r_b = [r_b]
        if len(r_a) > self.deg or len(r_b) > self.deg:
            raise ValueError(f'exceed the maximum degree {self.deg}')