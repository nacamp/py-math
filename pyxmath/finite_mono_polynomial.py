import itertools
from pyxmath import number_theory as nth


class FiniteMonoPolynomial():

    def __call__(self, val_coefs):
        if type(val_coefs) == int or type(val_coefs) == float:
            val_coefs = [val_coefs]
        return self.__class__(self.coefs, self.q, val_coefs)

    def __init__(self, coefs, q, val_coefs=None):
        '''

        Parameters
        ----------
        :param coef: list x^3+2x+1=0 : [1,2,0,1]
        :param p: int
        '''

        self.coefs = coefs
        self.deg = len(coefs) - 1
        self.q = q
        self.fields = []
        # ex: x^3+2x+1=0 (mod3) : x^3=0x^2-2x-1 : [2,1,0] :  left_coef=1, right_coefs=[2,1,0]
        self.right_coefs = coefs[:]
        self.right_coefs.pop()
        self.right_coefs = [x * (-1) % self.q for x in self.right_coefs]
        self.val_coefs = val_coefs

    def __repr__(self):
        return repr(self.val_coefs)

    def add(self, r_a, r_b):
        '''add

        :param r_a: right_coefs
        :param r_b: right_coefs
        :return:
        '''
        self._validate(r_a, r_b)
        return [sum(x) % self.q for x in itertools.zip_longest(r_a, r_b, fillvalue=0)]

    def __add__(self, other):
        if isinstance(other, self.__class__):
            self._validate(self.val_coefs, other.val_coefs)
            return self.__class__(self.coefs, self.q, self.add(self.val_coefs, other.val_coefs))
        else:
            return self.__class__(self.coefs, self.q, self.add(self.val_coefs, self.change_to_coefs(other)))

    __radd__ = __add__

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
            aa[i] = (aa[i] - bb[i]) % self.q
        return aa

    def __sub__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.coefs, self.q, self.sub(self.val_coefs, other.val_coefs))
        else:
            return self.__class__(self.coefs, self.q, self.sub(self.val_coefs, self.change_to_coefs(other)))

    def __rsub__(self, other):
        return self.__class__(self.coefs, self.q, self.sub(self.change_to_coefs(other), self.val_coefs))

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
        return ([x % self.q for x in coefs])

    def __mul__(self, other):
        if isinstance(other, self.__class__):
            self._validate(self.val_coefs, other.val_coefs)
            return self.__class__(self.coefs, self.q, self.mul(self.val_coefs, other.val_coefs))
        else:
            return self.__class__(self.coefs, self.q, self.mul(self.val_coefs, self.change_to_coefs(other)))

    __rmul__ = __mul__

    def div(self, r_a, r_b):
        self._validate(r_a, r_b)
        return self.mul(r_a, self.inv(r_b))

    def __truediv__(self, other):
        if isinstance(other, self.__class__):
            self._validate(self.val_coefs, other.val_coefs)
            return self.__class__(self.coefs, self.q, self.div(self.val_coefs, other.val_coefs))
        else:
            return self.__class__(self.coefs, self.q, self.div(self.val_coefs, self.change_to_coefs(other)))

    def __rtruediv__(self, other):
        return self.__class__(self.coefs, self.q, self.div(self.change_to_coefs(other), self.val_coefs))

    def poly_round_div(self, numerator, denominator):
        # n//d = q + r
        quotient = []
        n = numerator[:]
        while len(n) >= len(denominator):
            diff_deg = len(n) - len(denominator)
            d = [0] * diff_deg + denominator
            c = (n[-1] * nth.inv(d[-1], self.q)) % self.q
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
            old_r, r = r, [sum(x) % self.q for x in
                           itertools.zip_longest(old_r, [x * (-1) for x in self._mul_coef(quotient, r)], fillvalue=0)]
            old_s, s = s, [sum(x) % self.q for x in
                           itertools.zip_longest(old_s, [x * (-1) for x in self._mul_coef(quotient, s)], fillvalue=0)]
            old_t, t = t, [sum(x) % self.q for x in
                           itertools.zip_longest(old_t, [x * (-1) for x in self._mul_coef(quotient, t)], fillvalue=0)]
            while len(r) and r[-1] == 0:
                r.pop()
        return old_r, old_s, old_t
        # # old_r[0]??? 1??? ??????????????? old_r[0]?????? ????????? ??????.
        # old_s_inv = mul_inverse_mod(old_r[0], self.mod)
        # # result = [ x % 3 for x in self.qoly_mul2(old_s, [old_s_inv])]
        # result = [x % self.mod for x in self._mul_coef(old_s, [old_s_inv])]
        # return result + ([0] * (len(self.irr_coef) - len(result)))

    def inv(self, a):
        gcd, s, t = self.xgcd(a)
        # gcd[0]??? 1??? ??????????????? old_r[0]?????? ????????? ??????.
        s_inv = nth.inv(gcd[0], self.q)
        result = [x % self.q for x in self._mul_coef(s, [s_inv])]
        return result + ([0] * (len(self.right_coefs) - len(result)))

    def pow(self, a, n):
        if n == 0:
            return [1] + [0] * (len(a) - 1)
        x = a[:]
        for i in range(n - 1):
            x = self.mul(x, a)
        return (x)

    def __pow__(self, other, mod=None):
        return self.__class__(self.coefs, self.q, self.pow(self.val_coefs, other))
        # if mod is not None:
        #     return self.__class__(self.coefs, self.q, self.pow(self.val_coefs, other))
        # else:
        #     return self.__class__(self.coefs, self.q, self.pow(self.val_coefs, other))

    def mod(self, r_a):
        return ([x % self.q for x in r_a])

    # https://math.stackexchange.com/questions/844231/polynomial-multiplication-modulo-polynomial
    # TODO: when use mod, coefs or poly
    def __mod__(self, other):
        return self.__class__(self.coefs, self.q, [x % other for x in self.val_coefs])

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            aa = list(self.val_coefs[:])
            bb = list(other.val_coefs[:])
            for i in range(max(len(self.val_coefs), len(other.val_coefs))):
                if len(aa) == i:
                    aa.append(0)
                if len(bb) == i:
                    bb.append(0)
            return aa == bb
            # self.val_coefs == other.val_coefs
        else:
            return self.val_coefs == self.change_to_coefs(other)

    def neg(self):
        return [x * (-1) % self.q for x in self.val_coefs]

    def __neg__(self):
        return self.__class__(self.coefs, self.q, self.neg())

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

    def change_to_coefs(self, other):
        if type(other) == int:
            return [other] + [0] * (self.deg - 1)
        if isinstance(other, nth.Field):
            return [other.n] + [0] * (self.deg - 1)
        if type(other) == list:
            return other
        raise ValueError('int, Field, list only supported')

    def div_q_r(self, i):
        q = i // self.q
        r = i % self.q
        return q, r

    def elements(self):
        coefs = []
        coef_0 = [0] * self.deg
        for i in range(self.q ** self.deg):
            coef = coef_0[:]
            q = i
            for j in range(self.q):
                q, r = self.div_q_r(q)
                coef[j] = r
                if q == 0:
                    break
            coefs.append(self(coef))
        return coefs

    def make_power_representation(self):
        right_coefs_len = len(self.right_coefs)
        element = [[1] + [0] * (right_coefs_len - 1)]
        for i in range(self.q ** right_coefs_len - 2):
            f = element[i][:]
            f.insert(0, 0)
            ms = f.pop(right_coefs_len)
            if ms != 0:
                # ????????? deg??? ?????? 0??? ????????? ?????? deg ????????? ?????? ??? ?????? ?????????.
                # [ms * irr_coef[0], ms * irr_coef[1], ...] + f
                new_field = []
                zip_object = zip([x * ms for x in self.right_coefs], f)
                for a, b in zip_object:
                    new_field.append(a + b)
                element.append([x % self.q for x in new_field])
            else:
                element.append([x % self.q for x in f])
        return element;