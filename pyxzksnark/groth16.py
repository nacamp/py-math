# https://codeocean.com/capsule/8850121/tree/v3
from fractions import Fraction
import math
from random import *
from pyxmath.number_theory import *
from pyxmath.finite_mono_polynomial import *
from pyxmath.elliptic_curve import *
from pyxmath.number_theory import Field as F
from pyxzksnark.qap import *

##### start... ######

# matlab에 있는 초기 데이타
R = [1, 2, 4, 8, 4, -12]
A = [
    [0, 1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, -4, -4, 1, 1, 0]]
B = [
    [0, 0, 1, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0]]

C = [
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1]
]

NumWires = len(A[0])
NumGates = len(A)
Roots_of_Z = list(i for i in range(1, NumGates + 1))
s = set()
p = 71
q = p ** 2 - 1
Ap, Bp, Cp, Zp = r1cs_to_qap(A, B, C)


# scale and do modulus
def _find_denominator(Xp, s):
    for a in Xp:
        for aa in a:
            n = Fraction(str(aa))
            s.add(n.denominator)


_find_denominator(Ap, s)
_find_denominator(Bp, s)
_find_denominator(Cp, s)

lcm_val = math.lcm(*s)


def _scale_and_mod(Xp, lcm_val):
    for a in Xp:
        for i, aa in enumerate(a):
            a[i] = int(lcm_val * aa) % q


_scale_and_mod(Ap, lcm_val)
_scale_and_mod(Bp, lcm_val)
_scale_and_mod(Cp, lcm_val ** 2)
for i, aa in enumerate(Zp):
    Zp[i] = rmod(aa, q)
for i, aa in enumerate(R):
    R[i] = rmod(aa, q)

# P(x)=H(x)*Z(x)
Apoly, Bpoly, Cpoly, Ppoly = create_solution_polynomials2(R, Ap, Bp, Cp, q)
Px = Ppoly
Hx = Hpoly = create_divisor_polynomial2(Ppoly, Zp, q)

##### CSR  start
# find generator
ec = EC([0, 1, 0, 1])


def _max_gen(ec, p):
    points = find_points(ec, p)
    torsions = []
    tmp_len = 0
    max_order_id = 0
    for i, p in enumerate(points):
        ps = torsion(ec, p, len(points) + 1)
        torsions.append(ps)
        if tmp_len < len(ps):
            max_order_id = i
            tmp_len = len(ps)
    return torsions[max_order_id][0], tmp_len


g, g_order = _max_gen(ec, p)
# f = F(p)
# g = PT(f(11), f(8))

poly = FiniteMonoPolynomial([1, 0, 1], p)
# TODO find torsion in poly
h = PT(poly([60, 0]), poly([0, 8]))

# tau
# # 실패 metlab도 동일
# alpha = 4685
# beta = 2741
# gamma = 1289
# delta = 1513
# x_val = 3729

alpha = 1 + randrange(q - 1)
beta = 1 + randrange(q - 1)
gamma = 1 + randrange(q - 1)
while math.gcd(gamma, q) != 1:
    gamma = 1 + randrange(q - 1)
delta = randrange(q)
while math.gcd(delta, q) != 1:
    delta = 1 + randrange(q - 1)
x_val = 1
while x_val in Roots_of_Z:
    x_val = 1 + randrange(q - 1)
# 당분간 고정값으로 설정
alpha = 716
beta = 2126
gamma = 181
delta = 4279
x_val = 4707
tau = [alpha, beta, gamma, delta, x_val]
print('tau', tau)

Ax_val = list(eval_poly(v, x_val) % q for v in Ap)
Bx_val = list(eval_poly(v, x_val) % q for v in Bp)
Cx_val = list(eval_poly(v, x_val) % q for v in Cp)
Zx_val = eval_poly(Zp, x_val) % q
Hx_val = eval_poly(Hpoly, x_val) % q

sigma1_1 = []
sigma1_1.append(ec.mul(g, alpha))
sigma1_1.append(ec.mul(g, beta))
sigma1_1.append(ec.mul(g, delta))

sigma1_2 = []
for i in range(NumGates):
    sigma1_2.append(ec.mul(g, x_val ** i % q))

sigma1_3 = []
# TODO: Ind_pri, Ind_pub find using cod
Ind_pri = [2, 3, 4, 5]
Ind_pub = [1, 6]
PUBLIC_VAL = VAL = [0] * 6
for i in range(len(Ind_pub)):
    VAL[Ind_pub[i] - 1] = (beta * Ax_val[Ind_pub[i] - 1] + alpha * Bx_val[Ind_pub[i] - 1] + Cx_val[
        Ind_pub[i] - 1]) * inv(gamma, q) % q;
    sigma1_3.append(ec.mul(g, VAL[Ind_pub[i] - 1]))

sigma1_4 = []
for i in range(len(Ind_pri)):
    val = (beta * Ax_val[Ind_pri[i] - 1] + alpha * Bx_val[Ind_pri[i] - 1] + Cx_val[Ind_pri[i] - 1]) * inv(delta, q) % q;
    sigma1_4.append(ec.mul(g, val))

sigma1_5 = []
for i in range(NumGates - 1):
    val = x_val ** i * inv(delta, q) * Zx_val % q;
    sigma1_5.append(ec.mul(g, val))

sigma2_1 = []
sigma2_1.append(ec.mul(h, beta))
sigma2_1.append(ec.mul(h, gamma))
sigma2_1.append(ec.mul(h, delta))

sigma2_2 = []
for i in range(NumGates):
    val = x_val ** i % q;
    sigma2_2.append(ec.mul(h, val))

Length_CRS = len(sigma1_1) + len(sigma1_2) + len(sigma1_3) + len(sigma1_4) + len(sigma1_5) + len(sigma2_1) + len(
    sigma2_2);
print('Length_CRS:', Length_CRS)

# CRS is valid
lval = (sum(a * b for a, b in zip(R, Ax_val)) * sum(a * b for a, b in zip(R, Bx_val)) - sum(
    a * b for a, b in zip(R, Cx_val))) % q
rval = Zx_val * Hx_val % q
assert lval == rval

# prove
r = 3421
s = 3819
r = 1 + randrange(q - 1)
s = 1 + randrange(q - 1)

Proof_A = sigma1_1[0]
for i in range(NumWires):
    temp = None
    for j in range(NumGates):
        temp = ec.add(temp, ec.mul(sigma1_2[j], Ap[i][j]))
    Proof_A = ec.add(Proof_A, ec.mul(temp, R[i]))
Proof_A = ec.add(Proof_A, ec.mul(sigma1_1[2], r))

Proof_B = sigma2_1[0]
for i in range(NumWires):
    temp = None
    for j in range(NumGates):
        temp = ec.add(temp, ec.mul(sigma2_2[j], Bp[i][j]))
    Proof_B = ec.add(Proof_B, ec.mul(temp, R[i]))
Proof_B = ec.add(Proof_B, ec.mul(sigma2_1[2], s));

temp_Proof_B = sigma1_1[1];
for i in range(NumWires):
    temp = None
    for j in range(NumGates):
        temp = ec.add(temp, ec.mul(sigma1_2[j], Bp[i][j]))
    temp_Proof_B = ec.add(temp_Proof_B, ec.mul(temp, R[i]))
temp_Proof_B = ec.add(temp_Proof_B, ec.mul(sigma1_1[2], s))
Proof_C = ec.add(ec.mul(Proof_A, s), ec.mul(temp_Proof_B, r))
Proof_C = ec.sub(Proof_C, ec.mul(sigma1_1[2], s * r))
for i in range(len(Ind_pri)):
    Proof_C = ec.add(Proof_C, ec.mul(sigma1_4[i], R[Ind_pri[i] - 1]))
for i in range(NumGates - 1):  # i=0:NumGates-2
    Proof_C = ec.add(Proof_C, ec.mul(sigma1_5[i], Hpoly[i]));
proof = [Proof_A, Proof_B, Proof_C]
print('proof : ', proof)

# Check the completeness of proof
A = (alpha + sum(a * b for a, b in zip(R, Ax_val)) + r * delta) % q
B = (beta + sum(a * b for a, b in zip(R, Bx_val)) + s * delta) % q
C = (inv(delta, q) * (sum(R[v - 1] * (beta * Ax_val[v - 1] + alpha * Bx_val[v - 1] + Cx_val[v - 1]) for v in
                          Ind_pri) + Hx_val * Zx_val) + A * s + B * r + (-r * s * delta) % q) % q
lhs = A * B % q
rhs = ((alpha * beta) % q + sum(gamma * R[v - 1] * PUBLIC_VAL[v - 1] for v in Ind_pub)) % q
rhs = (rhs + C * delta) % q
assert lhs == rhs

assert Proof_A == ec.mul(g, A)
assert Proof_B == ec.mul(h, B)
assert Proof_C == ec.mul(g, C)

# verify
# m cannot exceed p
if g_order > p:
    m = p
LHS = ec.tate_pairing(Proof_A, Proof_B, Proof_B, m)

RHS = ec.tate_pairing(sigma1_1[0], sigma2_1[0], sigma2_1[0], m)
temp = None;
for i in range(len(Ind_pub)):
    temp = ec.add(temp, ec.mul(sigma1_3[i], R[Ind_pub[i] - 1]))
RHS = RHS * ec.tate_pairing(temp, sigma2_1[1], sigma2_1[1], m);
RHS = RHS * ec.tate_pairing(proof[2], sigma2_1[2], sigma2_1[2], m)
assert LHS == RHS
