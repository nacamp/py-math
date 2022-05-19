'''
https://hackmd.io/@benjaminion/bls12-381
https://loup-vaillant.fr/tutorials/cofactor
https://github.com/zcash/librustzcash/blob/6e0364cd42a2b3d2b958a54771ef51a8db79dd29/pairing/src/bls12_381/README.md#generators
https://medium.com/atomrigslab/pairing-based-cryptography%EC%99%80-bls-signature%EC%9D%98-%EC%9D%B4%ED%95%B4-part-1-f4ec67d69940

parameter
https://datatracker.ietf.org/doc/html/draft-yonezawa-pairing-friendly-curves-01#section-4.2

implementations
https://github.com/Chia-Network/bls-signatures/tree/main/python-impl
https://github.com/algorand/bls_sigs_ref/

용어들
https://datatracker.ietf.org/doc/html/rfc8017#page-11
I2OSP - Integer-to-Octet-String primitive
OS2IP - Octet-String-to-Integer primitive

https://en.wikipedia.org/wiki/HKDF
HKDF extracts a pseudorandom key (PRK) using an HMAC hash function (e.g. HMAC-SHA256) on an optional salt (acting as a key) and any potentially weak input key material (IKM) (acting as data).
It then generates similarly cryptographically strong output key material (OKM) of any desired length by repeatedly generating PRK-keyed hash-blocks and then appending them into the output key material, finally truncating to the desired length.

https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-14#section-6.6.2
Simplified Shallue-van de Woestijne-Ulas method : "simplified SWU" map
6.6.1.  Shallue-van de Woestijne method

https://www.johndcook.com/blog/2019/04/20/isogeny-based-cryptography/
https://www.johndcook.com/blog/2019/04/21/what-is-an-isogeny/
isogeny-based encryption
https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-14#section-8.8
BLS12-381 G2

https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-14#page-77
3-isogeny map for BLS12-381 G2

DST : a domain separation tag
'''
from pyxmath import *
from pyxmath.number_theory import *
from pyxmath.elliptic_curve import *
from pyxmath.finite_mono_polynomial import *

x = -0xD201000000010000
q = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
r2 = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
h = 0x396C8C005555E1568C00AAAB0000AAAB
h2 = 0x5d543a95414e7f1091d50792876a202cd91de4547085abaa68a205b2e5a7ddfa628f1cb4d9e82ef21537e293a6691ae1616ec6e786f0c70cf1c38e31c7238e5
f = Field(q)
G1 = PT(f(3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507),
        f(1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569))
poly = FiniteMonoPolynomial([1, 0, 1], q)
G2 = PT(poly([352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160, 3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758]),
        poly([1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905, 927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582]))

def cofactor1():
    def _h1(x, q):
        # return (x - 1) ** 2 / 3
        return (x - 1) ** 2 * inv(3, q) % q

    def _h2(x, q):
        # return (x ** 8 - 4 * 4 ** 7 + 5 * x ** 6 - 4 * x ** 4 + 6 * x ** 3 - 4 * x ** 2 - x ** 1 + 13)/9
        return (x ** 8 - 4 * x ** 7 + 5 * x ** 6 - 4 * x ** 4 + 6 * x ** 3 - 4 * x ** 2 - 4 * x ** 1 + 13) * inv(9, q) % q

    assert _h1(x, q) == h


    ec = EC([4, 0, 0, 1])
    #f = Field(q)
    P1 = PT(f(0x04), f(0x0a989badd40d6212b33cffc3f3763e9bc760f988c9926b26da9dd85e928483446346b8ed00e1de5d5ea93e354abe706c))
    assert G1 == ec.mul(P1, h)
    P2 = PT(poly([0x02, 0x00]),
            poly([0x013a59858b6809fca4d9a3b6539246a70051a3c88899964a42bc9a69cf9acdd9dd387cfa9086b894185b9a46a402be73,
                  0x02d27e0ec3356299a346a09ad7dc4ef68a483c3aed53f9139d2f929a3eecebf72082e5e58c6da24ee32e03040c406d4f]))
    assert G2 == ec.mul(P2, h2)

def gen_pk():
    # KeyGen
    # 1. PRK = HKDF-Extract("BLS-SIG-KEYGEN-SALT-", IKM || I2OSP(0, 1))
    # 2. OKM = HKDF-Expand(PRK, keyInfo || I2OSP(L, 2), L)
    # 3. SK = OS2IP(OKM) mod r
    # 4. return SK
    '''
    pk = OKM*G1Generator()
    '''
    pass

def gen_sign():
    '''
    sk.value * g2_map(message, dst)
    Hp2->hash_to_field
        return hash_to_field(msg, count, dst, q, 2, 64, expand_message_xmd, hashlib.sha256)
    https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-14#page-19
    hash_to_field
    expand_message_xmd

    opt_swu2_map(hash_to_field)
    http://ed25519.cr.yp.to/ed25519-20110926.pdf
    https://eprint.iacr.org/2019/403.pdf
    https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-14#page-17
    '''
    pass


# https://loup-vaillant.fr/tutorials/cofactor
def cofactor2():
    # 여인수,cofactor
    # for i in range(1, 45):
    #     print(i, 25*i % 44)
    q = 44

    def _d(a, b):
        r = 0
        for i in range(1, b + 1):
            r += a
        return r % q

    # def x(a, b):
    #     r = 0
    #     for i in range(1, b + 1):
    #         r += a
    #     return r % q
    #
    # def e(a,b):
    #     return (x(4,a) + x(11,b))%44
    #
    # print(x(25, 4))
    # print(e(2, 2))

    '''
    X25519 is designed in such a way that:
    -The chosen generator of the curve has prime order.
    -The private key is a multiple of the cofactor.


    B   = 4                  -- Generator of the prime-order group
    sa  = 20                 -- Alice's private key
    sb  = 28                 -- Bob's   private key
    SA  = B.sa  = 4.20  = 36 -- Alice's public key
    SB  = B.sb  = 4.28  = 24 -- Bob's   public key
    ssa = SB.sa = 24.20 = 40 -- Shared secret (computed by Alice)
    ssb = SA.sb = 36.28 = 40 -- Shared secret (computed by Bob)


    LO  = 11                 -- random low order point
    HA  = SA+LO = 36+11 =  3 -- Alice's "hidden" key
    ssb = HA.sb =  3.28 = 40
    '''

    B = 4  # Generator of the prime-order group
    sa = 20  # Alice's private key
    sb = 28  # Bob's   private key
    SA = _d(B, sa)  # Alice's public key
    SB = _d(B, sb)  # Bob's   public key
    ssa = _d(SB, sa)  # Shared secret (computed by Alice)
    assert ssa == 40  # 4.28.20
    ssb = _d(SA, sb)  # Shared secret (computed by Bob)
    assert ssb == 40  # 4.20.28

    LO = 11  # random low order point
    HA = SA + LO  # Alice's "hidden" key  # 36 + 11
    assert ssb == _d(HA, sb)  # (36+11).(4x5)%44, 11은 소멸

    O = 11  # low order point (order 4)
    B = 4  # base point      (prime order)
    D = _d(B, LO) == 15  # "dirty" base point
    # SA = B.sa = 4.20  = 36 # Alice's public key
    # HA = D.sa = 15.20 = 36 # ???

    d = 33  # random multiple of the prime order
    da = _d(sa, d)  # Alice's "dirty" secret key
    HA = _d(D, da)  # Alice's hidden key
    assert SA != HA
