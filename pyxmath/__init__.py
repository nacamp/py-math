import math


class PT():
    def __init__(self, x, y, z=None):
        self.x = x
        self.y = y
        self.z = z

    def to_projective(self):
        return self.__class__(self.x, self.y, 1)

    def to_projective(self):
        return self.__class__(self.x, self.y, 1)

    def to_affine_from_projective(self):
        return self.__class__(self.x / self.z, self.y / self.z)

    def to_affine_from_jacobian(self):
        return self.__class__(self.x / self.z ** 2, self.y / self.z ** 3)

    def __repr__(self):
        return repr([self.x, self.y, self.z])

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        elif other is None:
            return False
        else:
            return self.x == other[0] and self.y == other[1]

    def __neg__(self):
        return self.__class__(self.x, -self.y)

    def copy(self):
        return self.__class__(self.x, self.y)


# https://gist.github.com/dzhou/2632362
# https://ratsgo.github.io/data%20structure&algorithm/2017/10/07/prime/
def prime_sieve(sieveSize):
    # creating Sieve (0~n까지의 slot)
    sieve = [True] * (sieveSize + 1)
    # 0과 1은 소수가 아니므로 제외
    sieve[0] = False
    sieve[1] = False
    # 2부터 (루트 n) + 1까지의 숫자를 탐색
    for i in range(2, int(math.sqrt(sieveSize)) + 1):
        # i가 소수가 아니면 pass
        if sieve[i] == False:
            continue
        # i가 소수라면 i*i~n까지 숫자 가운데 i의 배수를
        # 소수에서 제외
        for pointer in range(i ** 2, sieveSize + 1, i):
            sieve[pointer] = False
    primes = []
    # sieve 리스트에서 True인 것이 소수이므로
    # True인 값의 인덱스를 결과로 저장
    for i in range(sieveSize + 1):
        if sieve[i] == True:
            primes.append(i)
    return primes


# 소인수분해, prime factorization
def get_prime_factors(n):
    # n 범위 내의 소수를 구한다
    primelist = prime_sieve(n)
    # 이 소수들 중 n으로 나누어 떨어지는
    # 소수를 구하고, 몇 번 나눌 수 있는지 계산
    # 예 : n = 8, factors = [(2, 3)]
    # 예 : n = 100, fcount = [(2: 2), (5: 2)]
    factors = []
    for p in primelist:
        count = 0
        while n % p == 0:
            n /= p
            count += 1
        if count > 0:
            factors.append((p, count))
    return factors


__all__ = ['PT', 'prime_sieve', 'get_prime_factors']
