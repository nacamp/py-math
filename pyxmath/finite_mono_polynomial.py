import itertools


class FiniteMonoPolynomialError(Exception):
    def __init__(self, message):
        self.message = message


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
            raise FiniteMonoPolynomialError(f'exceed the maximum degree {self.deg}')
