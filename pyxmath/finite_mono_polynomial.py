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
        # ex: x^3+2x+1=0 (mod3) :  x^3=0x^2-2x-1 : [2,1,0]
        self.finite_coefs = coefs[:]
        self.finite_coefs.pop()
        self.finite_coefs = [x * (-1) % self.mod for x in self.finite_coefs]

    def add(self, a, b):
        if len(a) > self.deg or len(b) > self.deg:
            raise FiniteMonoPolynomialError(f'exceed the maximum degree {self.deg}')
        if isinstance(a, int):
            a = [a]
        if isinstance(b, int):
            b = [b]
        return [sum(x) % self.mod for x in itertools.zip_longest(a, b, fillvalue=0)]

    def sub(self, a, b):
        if len(a) > self.deg or len(b) > self.deg:
            raise FiniteMonoPolynomialError(f'exceed the maximum degree {self.deg}')
        if isinstance(a, int):
            a = [a]
        if isinstance(b, int):
            b = [b]
        aa = list(a[:])
        bb = list(b[:])
        for i in range(max(len(aa), len(bb))):
            if len(aa) == i:
                aa.append(0)
            if len(bb) == i:
                bb.append(0)
            aa[i] = (aa[i] - bb[i]) % self.mod
        return aa