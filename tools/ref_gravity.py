#!/usr/bin/env python3
"""Independent Python (sympy / numpy / scipy) check of sections 11-14 of 7+1.nb:
gravity on the 7+1 spacetime (one time coordinate x0, seven space coordinates x1 ... x7), the Einstein and
Lovelock tensors, the generalized frame matrices, the canonical spin connection with a vanishing total
covariant derivative of the frame matrices, the covariant Dirac equation of dirac16complex, and the
three cosmological epochs (primordial near x0 = 0, modern, far future) with the pair creation of dirac16complex.

Conventions (as in 7+1.nb):  one time-like dimension x0 and seven space-like dimensions x1 ... x7 (declared in xAct, as in
the reference notebook Einstein-Lovelock-7+1-Nash.nb, with the signature list {7, 1, 0}: seven space-like, one time-like,
none degenerate);  metric components with the signs of eta17 = diag(+1, -1, ..., -1), which is simply the squares of the
octonion basis (e_0^2 = +1, e_i^2 = -1), so g_00 > 0;  Christoffel symbols chr[l][m][n] = {l; m n};
Riemann tensor (Landau-Lifshitz)  R^r_{s m n} = d_m chr^r_{s n} - d_n chr^r_{s m} + chr^r_{m l} chr^l_{s n} - chr^r_{n l} chr^l_{s m};
Ricci R_{s n} = R^r_{s r n};  Einstein G^m_n = R^m_n - delta^m_n R/2;  frame indices A, B = 0..7 with eta17.

Run from the repository root:   python3 tools/ref_gravity.py
Every line of output is  name : value ; the script exits with status 1 if any check is not True.
"""
import itertools, math, os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from ref71 import G17, eta17

R = {}                                   # name -> result (True/False or a number)
G = {A: G17[A] for A in range(8)}        # Gamma17[A], lower frame index
Gup = {A: eta17[A, A]*G[A] for A in range(8)}
I16, Z16 = sp.eye(16), sp.zeros(16)

from lovelock71 import geometry, ricci, rud, perm_sign, pair_sets, ordered_seqs, lovelockE, lovelockL, D, X

# ======================================================================= 0. eta17: the squares of the octonion basis
# The spacetime has one time-like dimension (x0) and seven space-like dimensions (x1..x7); xAct declares it with the
# signature {7, 1, 0} (section 12.2 of 7+1.nb). The matrix eta17 = diag(+1, -1, ..., -1) is simply the squares of the
# octonion basis: e_0^2 = +1, e_i^2 = -1, i.e. eta_AB = Re(e_A e_B), and Re(x x) = x0^2 - x1^2 - ... - x7^2.
from ref71 import omul
from octalg import basis
R['eta17 is simply the squares of the octonion basis: eta_AB = Re(e_A e_B), e_0^2 = +1, e_i^2 = -1 (i = 1..7)'] = (
    sp.Matrix(8, 8, lambda A, B: omul(basis(A), basis(B))[0]) == eta17)
_xq = sp.symbols('q0:8', real=True)
R['the quadratic form of the octonions has the signs of the squares of the octonion basis: Re(x x) = x0^2 - x1^2 - ... - x7^2 = x.eta17.x'] = (
    sp.expand(omul(list(_xq), list(_xq))[0] - sum(eta17[A, A]*_xq[A]**2 for A in range(8))) == 0)

# ======================================================================= 1. the cosmological metric
t = X[0]
a = sp.Function('a', real=True, positive=True)(t)
b = sp.Function('b', real=True, positive=True)(t)
gC = sp.diag(1, -a**2, -a**2, -a**2, -b**2, -b**2, -b**2, -b**2)
giC, chrC, RiemC = geometry(gC)
RicC, RsC, GC = ricci(RiemC, giC)
RudC = rud(RiemC, giC)
Ha, Hb = sp.symbols('H_a H_b', real=True)
dHa, dHb = sp.symbols('dH_a dH_b', real=True)               # dH/dx0
def hubble(e):
    """rewrite in terms of Ha = a'/a, Hb = b'/b and their derivatives"""
    e = e.subs({sp.Derivative(a, (t, 2)): (dHa + Ha**2)*a, sp.Derivative(b, (t, 2)): (dHb + Hb**2)*b})
    e = e.subs({sp.Derivative(a, t): Ha*a, sp.Derivative(b, t): Hb*b})
    return sp.expand(sp.simplify(e))
K = {(0, 1): -(dHa + Ha**2), (0, 4): -(dHb + Hb**2), (1, 2): -Ha**2, (4, 5): -Hb**2, (1, 4): -Ha*Hb}
R["cosmological metric: sectional curvatures R^{AB}_{AB} = -(Ha'+Ha^2), -(Hb'+Hb^2), -Ha^2, -Hb^2, -Ha Hb"] = all(
    sp.expand(hubble(RudC[(i, j, i, j)]) - K[(i, j)]) == 0 for (i, j) in K)
R['cosmological metric: Riemann R^{cd}_{ab} non-zero only for {c,d} = {a,b}'] = all(set(k[:2]) == set(k[2:]) for k in RudC)
G00 = hubble(GC[0, 0]); G11 = hubble(GC[1, 1]); G44 = hubble(GC[4, 4])
R['G^0_0 = 3Ha^2 + 12 Ha Hb + 6 Hb^2'] = sp.expand(G00 - (3*Ha**2 + 12*Ha*Hb + 6*Hb**2)) == 0
R['G^1_1 = 2Ha\' + 3Ha^2 + 4Hb\' + 10Hb^2 + 8HaHb'] = sp.expand(G11 - (2*dHa + 3*Ha**2 + 4*dHb + 10*Hb**2 + 8*Ha*Hb)) == 0
R['G^4_4 = 3Ha\' + 6Ha^2 + 3Hb\' + 6Hb^2 + 9HaHb'] = sp.expand(G44 - (3*dHa + 6*Ha**2 + 3*dHb + 6*Hb**2 + 9*Ha*Hb)) == 0
R['scalar curvature R = -2(3Ha\' + 6Ha^2 + 4Hb\' + 10Hb^2 + 12 HaHb)'] = sp.expand(hubble(RsC) + 2*(3*dHa + 6*Ha**2 + 4*dHb + 10*Hb**2 + 12*Ha*Hb)) == 0

# ======================================================================= 2. Lovelock tensors of the cosmological metric
EC = {k: sp.Matrix(D, D, lambda m, n: sp.simplify(lovelockE(RudC, k, m, n))) for k in range(0, 5)}
LC = {k: sp.simplify(lovelockL(RudC, k)) for k in range(0, 5)}
R['E^(1) = Einstein tensor G (cosmological metric)'] = sp.simplify(EC[1] - GC) == sp.zeros(D)
R['L_1 = R'] = sp.simplify(LC[1] - RsC) == 0
R['E^(0) = -delta/2'] = EC[0] == -sp.eye(D)/2
R['trace identity E^(k)mu_mu = -(8-2k)/2 L_k, k = 0..4'] = all(sp.simplify(EC[k].trace() + sp.Rational(D - 2*k, 2)*LC[k]) == 0 for k in range(5))
R['E^(4) = 0 identically in 8 dimensions, but L_4 (Euler density) is not 0'] = EC[4] == sp.zeros(D) and LC[4] != 0
R['E^(k) diagonal: E^(k)mu_nu = 0 for mu != nu, and E^(k)1_1 = E^(k)2_2 = E^(k)3_3, E^(k)4_4 = ... = E^(k)7_7'] = all(
    EC[k][m, n] == 0 for k in range(5) for m in range(D) for n in range(D) if m != n) and all(
    sp.simplify(EC[k][1, 1] - EC[k][i, i]) == 0 for k in range(5) for i in (2, 3)) and all(
    sp.simplify(EC[k][4, 4] - EC[k][i, i]) == 0 for k in range(5) for i in (5, 6, 7))
E2_00 = -12*Hb*(4*Ha**3 + 18*Ha**2*Hb + 12*Ha*Hb**2 + Hb**3)
E3_00 = 360*Ha**2*Hb**3*(4*Ha + 3*Hb)
R['E^(2)0_0 = -12 Hb (4Ha^3 + 18Ha^2Hb + 12HaHb^2 + Hb^3)'] = sp.expand(hubble(EC[2][0, 0]) - E2_00) == 0
R['E^(3)0_0 = 360 Ha^2 Hb^3 (4Ha + 3Hb)'] = sp.expand(hubble(EC[3][0, 0]) - E3_00) == 0
Rl_ = {}
for (p_, s_, m_, n_), v_ in RiemC.items():
    Rl_[(p_, s_, m_, n_)] = gC[p_, p_]*v_                                         # R_{rsmn}, diagonal metric
RicU_ = giC*RicC*giC
L2GB = sp.simplify(RsC**2 - 4*sum(RicC[i, j]*RicU_[i, j] for i in range(D) for j in range(D))
                   + sum(v_**2*giC[k_[0], k_[0]]*giC[k_[1], k_[1]]*giC[k_[2], k_[2]]*giC[k_[3], k_[3]] for k_, v_ in Rl_.items()))
HL = sp.zeros(D)
for m_ in range(D):
    for n_ in range(D):
        v_ = RsC*RicC[m_, n_] - 2*sum(RicC[m_, x_]*giC[x_, x_]*RicC[x_, n_] for x_ in range(D))
        v_ += -2*sum(RicU_[x_, y_]*Rl_.get((m_, x_, n_, y_), 0) for x_ in range(D) for y_ in range(D))
        v_ += sum(Rl_.get((m_, x_, y_, z_), 0)*giC[x_, x_]*giC[y_, y_]*giC[z_, z_]*Rl_.get((n_, x_, y_, z_), 0) for x_ in range(D) for y_ in range(D) for z_ in range(D))
        HL[m_, n_] = 2*v_ - gC[m_, n_]*L2GB/2
R['E^(2) = Lanczos tensor 2(R R_mn - 2R_ma R^a_n - 2R^ab R_manb + R_m^abc R_nabc) - g_mn L_2/2, and L_2 = R^2 - 4 Ric.Ric + Riem.Riem (cosmological metric)'] = (
    (giC*HL - EC[2]).applyfunc(sp.simplify) == sp.zeros(D) and sp.simplify(L2GB - LC[2]) == 0)
R['closed forms'] = {(k, c): hubble(EC[k][c, c]) for k in (1, 2, 3) for c in (0, 1, 4)}
# conservation (contracted Bianchi) for each k: d0 E^0_0 + 3Ha (E^0_0 - E^1_1) + 4Hb (E^0_0 - E^4_4) = 0
def cons(k):
    e0, e1, e4 = EC[k][0, 0], EC[k][1, 1], EC[k][4, 4]
    return sp.simplify(sp.diff(e0, t) + 3*sp.diff(a, t)/a*(e0 - e1) + 4*sp.diff(b, t)/b*(e0 - e4))
R['divergence-free E^(k), k = 1, 2, 3 (cosmological metric, arbitrary a, b)'] = all(cons(k) == 0 for k in (1, 2, 3))
# variational origin: minisuperspace with lapse N
N = sp.Function('N', real=True, positive=True)(t)
gN = sp.diag(N**2, -a**2, -a**2, -a**2, -b**2, -b**2, -b**2, -b**2)
giN, chrN, RiemN = geometry(gN)
RudN = rud(RiemN, giN)
def EL(Lag, qf):
    return sp.diff(Lag, qf) - sp.diff(sp.diff(Lag, sp.diff(qf, t)), t) + sp.diff(sp.diff(Lag, sp.diff(qf, t, 2)), t, 2)
mini = True
for k in range(1, 5):
    Lag = sp.simplify(N*a**3*b**4*lovelockL(RudN, k))
    eN = sp.simplify(EL(Lag, N).subs(N, 1).doit()); ea = sp.simplify(EL(Lag, a).subs(N, 1).doit()); eb = sp.simplify(EL(Lag, b).subs(N, 1).doit())
    mini = mini and sp.simplify(eN + 2*a**3*b**4*EC[k][0, 0]) == 0 and sp.simplify(ea + 6*a**2*b**4*EC[k][1, 1]) == 0 \
        and sp.simplify(eb + 8*a**3*b**3*EC[k][4, 4]) == 0
R['variation of sqrt|g| L_k (lapse N, a, b): EL_N = -2 a^3 b^4 E^(k)0_0, EL_a = -6 a^2 b^4 E^(k)1_1, EL_b = -8 a^3 b^3 E^(k)4_4, k = 1..4'] = mini

# ======================================================================= 3. an inhomogeneous metric: exact divergence at a point
class Dual:
    """number + first derivatives along x0 and x4 (exact rationals)"""
    __slots__ = ('v', 'd')
    def __init__(self, v, d): self.v = v; self.d = tuple(d)
    def __add__(s, o):
        if isinstance(o, Dual): return Dual(s.v + o.v, [x + y for x, y in zip(s.d, o.d)])
        return Dual(s.v + o, s.d)
    __radd__ = __add__
    def __mul__(s, o):
        if isinstance(o, Dual): return Dual(s.v*o.v, [s.v*y + x*o.v for x, y in zip(s.d, o.d)])
        return Dual(s.v*o, [x*o for x in s.d])
    __rmul__ = __mul__
    def __eq__(s, o):
        if not isinstance(o, Dual) and o == 0: return s.v == 0 and all(x == 0 for x in s.d)
        return NotImplemented
    def __ne__(s, o): return not s.__eq__(o)
x0, x4 = X[0], X[4]
Nf = 1 + x0*x4/5 + x0**2/11
Af = 1 + x0/3 + x4**2/7
Bf = 2 - x0*x4/4 + x0**2/9 + x4/13
gI = sp.diag(Nf**2, -Af**2, -Af**2, -Af**2, -Bf**2, -Bf**2, -Bf**2, -Bf**2)
giI, chrI, RiemI = geometry(gI, sp.cancel)
RudI = rud(RiemI, giI, sp.cancel)
PT = {x0: sp.Rational(1, 2), x4: sp.Rational(1, 3)}
val = lambda e: Fr(str(sp.nsimplify(sp.sympify(e).subs(PT))))
RudD = {key: Dual(val(v), [val(sp.diff(v, x0)), val(sp.diff(v, x4))]) for key, v in RudI.items()}
chrP = [[[val(chrI[l][m][n]) for n in range(D)] for m in range(D)] for l in range(D)]
gP = [[val(gI[i, j]) for j in range(D)] for i in range(D)]
RicI = [[sum(RiemI.get((r, s, r, n), 0) for r in range(D)) for n in range(D)] for s in range(D)]
RsI = sum(giI[s, s]*RicI[s][s] for s in range(D))
GI = [[val(sum(giI[m, s]*RicI[s][n] for s in range(D)) - (RsI/2 if m == n else 0)) for n in range(D)] for m in range(D)]
divs, syms, e1g = {}, {}, None
for k in (1, 2, 3):
    Ek = [[lovelockE(RudD, k, m, n) for n in range(D)] for m in range(D)]
    Ek = [[(e if isinstance(e, Dual) else Dual(Fr(e), [0, 0])) for e in row] for row in Ek]
    div = []
    for n in range(D):
        s_ = Ek[0][n].d[0] + Ek[4][n].d[1]
        s_ += sum(chrP[m][m][l]*Ek[l][n].v for m in range(D) for l in range(D))
        s_ -= sum(chrP[l][m][n]*Ek[m][l].v for m in range(D) for l in range(D))
        div.append(s_)
    divs[k] = div
    syms[k] = all(sum(gP[m][l]*Ek[l][n].v for l in range(D)) == sum(gP[n][l]*Ek[l][m].v for l in range(D)) for m in range(D) for n in range(D))
    if k == 1: e1g = all(Ek[m][n].v == GI[m][n] for m in range(D) for n in range(D))
    if k == 3: R['inhomogeneous metric: E^(3)0_4 at the point (exact)'] = str(Ek[0][4].v)
R['inhomogeneous metric (N, A, B of x0 and x4): E^(1) = G at the point'] = e1g
R['inhomogeneous metric: nabla_mu E^(k)mu_nu = 0 exactly at x0 = 1/2, x4 = 1/3, k = 1, 2, 3'] = all(all(v == 0 for v in divs[k]) for k in (1, 2, 3))
R['inhomogeneous metric: E^(k)_mu nu symmetric, k = 1, 2, 3'] = all(syms.values())

# ======================================================================= 4. a plane gravitational wave (pp-wave)
u = X[0] - X[1]
Fpp = sp.Function('F', real=True)(u, X[2], X[3])
gW = sp.Matrix(sp.diag(1, -1, -1, -1, -1, -1, -1, -1))
gW[0, 0] += Fpp; gW[1, 1] += Fpp; gW[0, 1] -= Fpp; gW[1, 0] -= Fpp
giW, chrW, RiemW = geometry(gW)
RudW = rud(RiemW, giW)
EW = {k: sp.Matrix(D, D, lambda m, n: sp.simplify(lovelockE(RudW, k, m, n))) for k in (1, 2, 3)}
lapF = sp.diff(Fpp, X[2], 2) + sp.diff(Fpp, X[3], 2)
lW = sp.Matrix([1, -1, 0, 0, 0, 0, 0, 0])                # l_mu
EWlow = (gW*EW[1]).applyfunc(sp.simplify)
R['pp-wave: E^(1)_mu nu = +(1/2) (d2^2 + d3^2) F l_mu l_nu, l = dx0 - dx1'] = sp.simplify(EWlow - lapF/2*lW*lW.T) == sp.zeros(D)
R['pp-wave: E^(2) = E^(3) = 0 for every profile F'] = EW[2] == sp.zeros(D) and EW[3] == sp.zeros(D)

# ======================================================================= 5. Kasner solutions (primordial epoch)
p, q = sp.symbols('p q', real=True)
gK = sp.diag(1, -t**(2*p), -t**(2*p), -t**(2*p), -t**(2*q), -t**(2*q), -t**(2*q), -t**(2*q))
giK, chrK, RiemK = geometry(gK)
RudK = rud(RiemK, giK)
EK = {k: [sp.factor(sp.simplify(lovelockE(RudK, k, c, c)*t**(2*k))) for c in (0, 1, 4)] for k in (1, 2, 3)}
R['Kasner E^(k) components x0^(2k) (00, 11, 44)'] = EK
pE1, qE1 = (1 + 2*sp.sqrt(2))/7, (2 - 3*sp.sqrt(2))/14
pE2, qE2 = (1 - 2*sp.sqrt(2))/7, (2 + 3*sp.sqrt(2))/14
R['Einstein Kasner: p = (1+-2 sqrt2)/7, q = (2-+3 sqrt2)/14 solve E^(1) = 0; 3p+4q = 1, 3p^2+4q^2 = 1'] = all(
    sp.simplify(e.subs({p: pp, q: qq})) == 0 for (pp, qq) in ((pE1, qE1), (pE2, qE2)) for e in EK[1]) and all(
    sp.simplify(3*pp + 4*qq - 1) == 0 and sp.simplify(3*pp**2 + 4*qq**2 - 1) == 0 for (pp, qq) in ((pE1, qE1), (pE2, qE2)))
cub = sp.expand((4*p**3 + 18*p**2*q + 12*p*q**2 + q**3).subs(q, (3 - 3*p)/4)*64)
R['Gauss-Bonnet Kasner cubic (64 x e4 at q = (3-3p)/4)'] = cub
gbroots = sorted([r for r in sp.Poly(cub, p).nroots(n=30) if abs(sp.im(r)) < 1e-25], key=lambda z: float(sp.re(z)))
gbroots = [sp.re(r) for r in gbroots]
R['Gauss-Bonnet Kasner exponents p (three real roots)'] = tuple(gbroots)
R['Gauss-Bonnet Kasner: the three (p, (3-3p)/4) solve E^(2) = 0, sum of the 7 exponents = 3'] = len(gbroots) == 3 and all(
    abs(sp.N(e.subs({p: r, q: (3 - 3*r)/4}), 25)) < 1e-18 for r in gbroots for e in EK[2])
R['cubic Lovelock Kasner: (p, q) = (-15/7, 20/7) solves E^(3) = 0, 3p + 4q = 5; degenerate families q = 0 and p = 0'] = all(
    sp.simplify(e.subs({p: sp.Rational(-15, 7), q: sp.Rational(20, 7)})) == 0 for e in EK[3]) and all(
    sp.simplify(e.subs(q, 0)) == 0 for e in EK[3]) and all(sp.simplify(e.subs(p, 0)) == 0 for e in EK[3])
R['Einstein Kasner is not a cubic-Lovelock vacuum (E^(3) != 0)'] = any(sp.simplify(e.subs({p: pE1, q: qE1})) != 0 for e in EK[3])

# ======================================================================= 6. de Sitter (far future)
Hs = sp.symbols('H', positive=True)
dS8 = {(k, c): sp.simplify(hubble(EC[k][c, c]).subs({Ha: Hs, Hb: Hs, dHa: 0, dHb: 0})) for k in (1, 2, 3) for c in (0, 1, 4)}
R['8D de Sitter (a = b = e^{H x0}): E^(1) = 21 H^2, E^(2) = -420 H^4, E^(3) = 2520 H^6 (times delta)'] = all(
    dS8[(1, c)] == 21*Hs**2 and dS8[(2, c)] == -420*Hs**4 and dS8[(3, c)] == 2520*Hs**6 for c in (0, 1, 4))
dS4 = {(k, c): sp.simplify(hubble(EC[k][c, c]).subs({Ha: Hs, Hb: 0, dHa: 0, dHb: 0})) for k in (1, 2, 3) for c in (0, 1, 4)}
R['4D de Sitter on a static T^4 (a = e^{H x0}, b = b0): E^(1) = diag(3,3,6) H^2, E^(2) = diag(0,0,-12) H^4, E^(3) = 0'] = (
    (dS4[(1, 0)], dS4[(1, 1)], dS4[(1, 4)]) == (3*Hs**2, 3*Hs**2, 6*Hs**2) and (dS4[(2, 0)], dS4[(2, 1)], dS4[(2, 4)]) == (0, 0, -12*Hs**4)
    and all(dS4[(3, c)] == 0 for c in (0, 1, 4)))
mod = {(k, c): sp.simplify(hubble(EC[k][c, c]).subs({Hb: 0, dHb: 0})) for k in (1, 2, 3) for c in (0, 1, 4)}
R['static T^4 (Hb = 0): E^(1) = (3Ha^2, 2Ha\'+3Ha^2, 3Ha\'+6Ha^2), E^(2) = (0, 0, -12 Ha^2 (Ha\'+Ha^2)), E^(3) = 0'] = (
    sp.expand(mod[(1, 0)] - 3*Ha**2) == 0 and sp.expand(mod[(1, 1)] - 2*dHa - 3*Ha**2) == 0 and sp.expand(mod[(1, 4)] - 3*dHa - 6*Ha**2) == 0
    and mod[(2, 0)] == 0 and mod[(2, 1)] == 0 and sp.expand(mod[(2, 4)] + 12*Ha**2*(dHa + Ha**2)) == 0 and all(mod[(3, c)] == 0 for c in (0, 1, 4)))

# ======================================================================= 6b. exponential vacuum solutions (constant Hubble rates)
# Vacuum equations in the form of Einstein-Lovelock-7+1-Nash.nb: 0 = -Lambda + H^-2 w1 E1 + H^-4 w2 E2 + H^-6 w3 E3 (mixed components),
# for a = e^{H x0}, b = e^{q H x0}: three distinct equations (00, 11, 44) with one linear relation (divergence identity).
from lovelock71 import cosmological_closed_forms
CFx = cosmological_closed_forms()
HaS, HbS, dHaS, dHbS = sp.symbols('H_a H_b dH_a dH_b', real=True)
Hc = sp.Symbol('H', positive=True); qS = sp.Symbol('q', real=True)
w1S, w2S, w3S, LamS = sp.symbols('w1 w2 w3 Lambda')
tExp = X[0]
gExp = sp.diag(1, *([-sp.exp(2*Hc*tExp)]*3), *([-sp.exp(2*qS*Hc*tExp)]*4))
giE, chrE, RiemE = geometry(gExp)
RudE = rud(RiemE, giE)
EExp = {k: [sp.simplify(lovelockE(RudE, k, c, c)) for c in (0, 1, 4)] for k in (1, 2, 3)}
R['exponential expansion a = e^{H x0}, b = e^{q H x0}: E^(k) (00, 11, 44) = closed forms with Ha = H, Hb = q H, Ha\' = Hb\' = 0 (k = 1, 2, 3)'] = all(
    sp.simplify(EExp[k][i] - CFx[(k, c)].subs({HaS: Hc, HbS: qS*Hc, dHaS: 0, dHbS: 0})) == 0 for k in (1, 2, 3) for i, c in enumerate((0, 1, 4)))
eqsX = [sp.expand(-LamS + w1S*EExp[1][i]/Hc**2 + w2S*EExp[2][i]/Hc**4 + w3S*EExp[3][i]/Hc**6) for i in range(3)]
R['exponential expansion: 3 Ha (eq00 - eq11) + 4 Hb (eq00 - eq44) = 0 identically (the three vacuum equations have rank 2)'] = sp.expand(
    3*Hc*(eqsX[0] - eqsX[1]) + 4*qS*Hc*(eqsX[0] - eqsX[2])) == 0
solX = sp.solve(eqsX, [w2S, w3S], dict=True)
R['exponential vacuum solutions: one solution (w2, w3) for given w1, Lambda, q, satisfying all three equations'] = (
    len(solX) == 1 and all(sp.simplify(e.subs(solX[0])) == 0 for e in eqsX))
solID = {k_: sp.simplify(v.subs({qS: -1})) for k_, v in solX[0].items()}
R['eternal inflation a = e^{H x0} with deflation b = e^{-H x0}: w2 = Lambda/4 - w1/2, w3 = w1/24 - Lambda/36 (Lambda = 0: w2/w1 = -1/2, w3/w1 = 1/24)'] = (
    sp.simplify(solID[w2S] - (LamS/4 - w1S/2)) == 0 and sp.simplify(solID[w3S] - (w1S/24 - LamS/36)) == 0)
R['NUM exponential vacuum solution (w2, w3) for general q'] = (sp.factor(solX[0][w2S]), sp.factor(solX[0][w3S]))
R['NUM eq00, eq11, eq44 for q = -1'] = tuple(sp.expand(e.subs({qS: -1})) for e in eqsX)

# ======================================================================= 7. frame matrices, spin connection, total covariant derivative
def frame_checks(g, e, name, simp=sp.simplify):
    gi, chr_, Riem = geometry(g, simp)
    Ric, Rs, Gm = ricci(Riem, gi, simp)
    E = e.inv().applyfunc(simp)                          # E[mu, A] = e_A^mu
    res = {}
    res['eta e e = g'] = (e.T*eta17*e - g).applyfunc(simp) == sp.zeros(D)
    Gam = {m: sum((e[A, m]*G[A] for A in range(D)), Z16) for m in range(D)}
    Gamu = {m: sum((E[m, A]*Gup[A] for A in range(D)), Z16) for m in range(D)}
    res['{Gamma_mu, Gamma_nu} = 2 g_mu nu'] = all((Gam[m]*Gam[n] + Gam[n]*Gam[m] - 2*g[m, n]*I16).applyfunc(simp) == Z16
                                                 for m in range(D) for n in range(m, D))
    res['Gamma_mu real'] = all(not x.has(sp.I) for m in range(D) for x in Gam[m])
    om = {}
    for m in range(D):
        om[m] = sp.Matrix(D, D, lambda A, B: simp(sum(e[A, n]*(sp.diff(E[n, B], X[m]) + sum(chr_[n][m][l]*E[l, B] for l in range(D)))
                                                    for n in range(D))))
    omlow = {m: (eta17*om[m]).applyfunc(simp) for m in range(D)}
    res['omega_mu AB = -omega_mu BA'] = all((omlow[m] + omlow[m].T).applyfunc(simp) == sp.zeros(D) for m in range(D))
    res['vielbein postulate d_mu e^A_nu - chr^l_mu nu e^A_l + omega_mu^A_B e^B_nu = 0'] = all(
        simp(sp.diff(e[A, n], X[m]) - sum(chr_[l][m][n]*e[A, l] for l in range(D)) + sum(om[m][A, B]*e[B, n] for B in range(D))) == 0
        for m in range(D) for n in range(D) for A in range(D))
    Om = {m: sum((omlow[m][A, B]*Gup[A]*Gup[B]/4 for A in range(D) for B in range(D) if omlow[m][A, B] != 0), Z16).applyfunc(simp) for m in range(D)}
    res['Omega_mu real, Gamma17[0] Omega_mu^T Gamma17[0] = -Omega_mu'] = all(not x.has(sp.I) for m in range(D) for x in Om[m]) and all(
        (G[0]*Om[m].T*G[0] + Om[m]).applyfunc(simp) == Z16 for m in range(D))
    res['total covariant derivative of the frame matrices: d_mu Gamma_nu - chr^l_mu nu Gamma_l + [Omega_mu, Gamma_nu] = 0'] = all(
        (sp.diff(Gam[n], X[m]) - sum((chr_[l][m][n]*Gam[l] for l in range(D)), Z16) + Om[m]*Gam[n] - Gam[n]*Om[m]).applyfunc(simp) == Z16
        for m in range(D) for n in range(D))
    res['constant frame matrices: omega_mu^A_B Gamma^B + [Omega_mu, Gamma^A] = 0'] = all(
        (sum((om[m][A, B]*Gup[B] for B in range(D)), Z16) + Om[m]*Gup[A] - Gup[A]*Om[m]).applyfunc(simp) == Z16
        for m in range(D) for A in range(D))
    F = {(m, n): (sp.diff(Om[n], X[m]) - sp.diff(Om[m], X[n]) + Om[m]*Om[n] - Om[n]*Om[m]).applyfunc(simp) for m in range(D) for n in range(D)}
    ok = True
    for m in range(D):
        for n in range(m + 1, D):
            M = Z16
            for A in range(D):
                for B in range(D):
                    v = simp(sum(eta17[A, A]*e[A, r]*Riem.get((r, s, m, n), 0)*E[s, B] for r in range(D) for s in range(D)))
                    if v != 0: M = M + v*Gup[A]*Gup[B]/4
            if (F[(m, n)] - M).applyfunc(simp) != Z16: ok = False
    res['spin curvature F_mu nu = (1/4) R_AB mu nu Gamma^A Gamma^B'] = ok
    res['Bianchi identity d_l F_mn + [Omega_l, F_mn] + cyclic = 0'] = all(
        (sum((sp.diff(F[(y, z)], X[x]) + Om[x]*F[(y, z)] - F[(y, z)]*Om[x] for (x, y, z) in ((l, m, n), (m, n, l), (n, l, m))), Z16)).applyfunc(simp) == Z16
        for (l, m, n) in itertools.combinations(range(D), 3))
    lic = (sum((Gamu[m]*Gamu[n]*F[(m, n)] for m in range(D) for n in range(D)), Z16)/2).applyfunc(simp)
    res['(1/2) Gamma^mu Gamma^nu F_mu nu = -(R/4) 1_16'] = (lic + Rs*I16/4).applyfunc(simp) == Z16
    for k_, v_ in res.items(): R[name + ': ' + k_] = v_
    return gi, chr_, Riem, Rs, e, E, Gam, Gamu, om, Om

eC = sp.diag(1, a, a, a, b, b, b, b)
giC, chrC, RiemC, RsC, eC, EfC, GamC, GamuC, omC, OmC = frame_checks(gC, eC, 'cosmological metric')
R['cosmological metric: Gamma^mu Omega_mu = (3a\'/(2a) + 2b\'/b) Gamma17[0]'] = (sum((GamuC[m]*OmC[m] for m in range(D)), Z16)
    - (3*sp.diff(a, t)/(2*a) + 2*sp.diff(b, t)/b)*G[0]).applyfunc(sp.simplify) == Z16
R['cosmological metric: Omega_i = (a\'/2) Gamma17[i] Gamma17[0] (i = 1..3), Omega_k = (b\'/2) Gamma17[k] Gamma17[0] (k = 4..7), Omega_0 = 0'] = (
    OmC[0] == Z16 and all((OmC[i] - sp.diff(a, t)/2*G[i]*G[0]).applyfunc(sp.simplify) == Z16 for i in (1, 2, 3))
    and all((OmC[k] - sp.diff(b, t)/2*G[k]*G[0]).applyfunc(sp.simplify) == Z16 for k in (4, 5, 6, 7)))
Fs = sp.Function('F', real=True)(u, X[2], X[3])
eW = sp.Matrix(sp.eye(8))
eW[0, 0] = 1 + Fs/2; eW[0, 1] = -Fs/2; eW[1, 0] = Fs/2; eW[1, 1] = 1 - Fs/2
frame_checks(gW, eW, 'pp-wave')
Nn = sp.Function('N', real=True, positive=True)(X[0], X[4])
An = sp.Function('A', real=True, positive=True)(X[0], X[4])
Bn = sp.Function('B', real=True, positive=True)(X[0], X[4])
gIn = sp.diag(Nn**2, -An**2, -An**2, -An**2, -Bn**2, -Bn**2, -Bn**2, -Bn**2)
eIn = sp.diag(Nn, An, An, An, Bn, Bn, Bn, Bn)
frame_checks(gIn, eIn, 'inhomogeneous metric N, A, B (x0, x4)')

# ======================================================================= 7b. Cartan's first structure equation; a static metric with a radial coordinate
def cartan_ok(e, om):
    """de^A + omega^A_B ^ e^B = 0 in components: d_m e^A_n - d_n e^A_m + omega_m^A_B e^B_n - omega_n^A_B e^B_m = 0"""
    return all(sp.simplify(sp.diff(e[A, n], X[m]) - sp.diff(e[A, m], X[n]) + sum(om[m][A, B]*e[B, n] - om[n][A, B]*e[B, m] for B in range(D))) == 0
               for A in range(D) for m in range(D) for n in range(m + 1, D))
R["Cartan's first structure equation de^A + omega^A_B ^ e^B = 0 (torsion-free) for the cosmological frame"] = cartan_ok(eC, omC)
rr = X[1]
fR = sp.Function('f', positive=True)(rr); hR = sp.Function('h', positive=True)(rr)
gR = sp.diag(fR, -hR, *([-rr**2]*6))
eR = sp.diag(sp.sqrt(fR), sp.sqrt(hR), *([rr]*6))
giR_, chrR_, RiemR_, RsR_, eR_, ER_, GamR_, GamuR_, omR_, OmR_ = frame_checks(gR, eR, 'static radial metric f(r) dx0^2 - h(r) dr^2 - r^2 (dx2^2 + ... + dx7^2)')
R["Cartan's first structure equation for the static radial metric"] = cartan_ok(eR, omR_)
Cst = {(A, B, Cc): sp.simplify(sum(ER_[m, B]*ER_[n, Cc]*(sp.diff(eR[A, n], X[m]) - sp.diff(eR[A, m], X[n])) for m in range(D) for n in range(D)))
       for A in range(D) for B in range(D) for Cc in range(D)}
Cstl = {k_: eta17[k_[0], k_[0]]*v_ for k_, v_ in Cst.items()}
R['structure coefficients: de^A = (1/2) C^A_BC e^B ^ e^C and omega_C,AB (form index C) = -(1/2)(C_ACB + C_BAC - C_CBA) (static radial metric)'] = all(
    sp.simplify(sum(ER_[m, Cc]*sum(eta17[A, K]*omR_[m][K, B] for K in range(D)) for m in range(D))
                + sp.Rational(1, 2)*(Cstl[(A, Cc, B)] + Cstl[(B, A, Cc)] - Cstl[(Cc, B, A)])) == 0 for A in range(D) for B in range(D) for Cc in range(D))
R['static radial metric: Omega_0 = (f\'/(4 sqrt(f h))) Gamma17[1] Gamma17[0], Omega_1 = 0, Omega_k = (1/(2 sqrt h)) Gamma17[1] Gamma17[k] (k = 2..7)'] = (
    (OmR_[0] - sp.diff(fR, rr)/(4*sp.sqrt(fR*hR))*G[1]*G[0]).applyfunc(sp.simplify) == Z16 and OmR_[1] == Z16
    and all((OmR_[k] - 1/(2*sp.sqrt(hR))*G[1]*G[k]).applyfunc(sp.simplify) == Z16 for k in range(2, 8)))
R['static radial metric: Gamma^mu Omega_mu = (1/sqrt h)(f\'/(4f) + 3/r) Gamma17up[1], so (i Gamma^mu D_mu - m) psi = [i Gamma17up[0] d0/sqrt f + i Gamma17up[1] (d_r + f\'/(4f) + 3/r)/sqrt h + sum_k i Gamma17up[k] d_k/r - m] psi'] = (
    (sum((GamuR_[m]*OmR_[m] for m in range(D)), Z16) - (sp.diff(fR, rr)/(4*fR) + 3/rr)/sp.sqrt(hR)*Gup[1]).applyfunc(sp.simplify) == Z16)
chiR = sp.Matrix([sp.Function('r%d' % i)(X[0], X[1], X[5]) for i in range(16)])
DslR = lambda f_: sum((GamuR_[m]*(f_.diff(X[m]) + OmR_[m]*f_) for m in range(D)), sp.zeros(16, 1))
wR = fR**sp.Rational(-1, 4)*rr**-3
R['static radial metric: psi = f^(-1/4) r^(-3) chi removes the spin connection: Gamma^mu D_mu psi = f^(-1/4) r^(-3) Gamma^mu d_mu chi'] = (
    DslR(wR*chiR) - wR*sum((GamuR_[m]*chiR.diff(X[m]) for m in range(D)), sp.zeros(16, 1))).applyfunc(sp.simplify) == sp.zeros(16, 1)
# minimal coupling of dirac16complex to an abelian gauge field (cosmological metric)
qc = sp.Symbol('q', real=True)
YY = (X[0], X[1], X[5])
Af = [sp.Function('A%d' % m_)(*YY) for m_ in range(D)]
alf = sp.Function('alpha')(*YY)
psiG = sp.Matrix([sp.Function('u%d' % i)(*YY) for i in range(16)])
DslA = lambda f_, Avec: sum((GamuC[m]*(f_.diff(X[m]) + OmC[m]*f_ - sp.I*qc*Avec[m]*f_) for m in range(D)), sp.zeros(16, 1))
Ash = [Af[m] + sp.diff(alf, X[m]) for m in range(D)]
R['U(1) gauge covariance: Gamma^mu (D_mu - i q A\'_mu)(e^{i q alpha} psi) = e^{i q alpha} Gamma^mu (D_mu - i q A_mu) psi with A\' = A + d alpha'] = (
    (DslA(sp.exp(sp.I*qc*alf)*psiG, Ash) - sp.exp(sp.I*qc*alf)*DslA(psiG, Af)).applyfunc(sp.simplify) == sp.zeros(16, 1))
psbG = sp.Matrix([[sp.Function('w%d' % i)(*YY) for i in range(16)]])
DpsiA = {m: psiG.diff(X[m]) + OmC[m]*psiG - sp.I*qc*Af[m]*psiG for m in range(D)}
DpsbA = {m: psbG.diff(X[m]) - psbG*OmC[m] + sp.I*qc*Af[m]*psbG for m in range(D)}
lhsA = sum((sp.diff(a**3*b**4*(psbG*GamuC[m]*psiG)[0, 0], X[m]) for m in range(D)), 0)
rhsA = a**3*b**4*sum(((DpsbA[m]*GamuC[m]*psiG)[0, 0] + (psbG*GamuC[m]*DpsiA[m])[0, 0] for m in range(D)), 0)
R['current identity with the gauge field: d_mu(sqrt|g| psibar Gamma^mu psi) = sqrt|g| [(D_mu psibar) Gamma^mu psi + psibar Gamma^mu D_mu psi], D including -i q A_mu'] = sp.simplify(sp.expand(lhsA - rhsA)) == 0

# ======================================================================= 8. the covariant Dirac equation of dirac16complex
Y = (X[0], X[1], X[5])
psi = sp.Matrix([sp.Function('f%d' % i)(*Y) for i in range(16)])
psb = sp.Matrix([[sp.Function('h%d' % i)(*Y) for i in range(16)]])
Dpsi = {m: psi.diff(X[m]) + OmC[m]*psi for m in range(D)}
Dpsb = {m: psb.diff(X[m]) - psb*OmC[m] for m in range(D)}
sqg = a**3*b**4
lhs = sum((sp.diff(sqg*(psb*GamuC[m]*psi)[0, 0], X[m]) for m in range(D)), 0)
rhs = sqg*sum(((Dpsb[m]*GamuC[m]*psi)[0, 0] + (psb*GamuC[m]*Dpsi[m])[0, 0] for m in range(D)), 0)
R['current identity d_mu(sqrt|g| psibar Gamma^mu psi) = sqrt|g| [(D_mu psibar) Gamma^mu psi + psibar Gamma^mu D_mu psi]'] = sp.simplify(sp.expand(lhs - rhs)) == 0
Dsl = lambda f: sum((GamuC[m]*(f.diff(X[m]) + OmC[m]*f) for m in range(D)), sp.zeros(16, 1))
boxs = sp.zeros(16, 1)
for m in range(D):
    if giC[m, m] != 0:
        boxs += giC[m, m]*(Dpsi[m].diff(X[m]) + OmC[m]*Dpsi[m] - sum((chrC[l][m][m]*Dpsi[l] for l in range(D)), sp.zeros(16, 1)))
R['Lichnerowicz: (Gamma^mu D_mu)^2 psi = (g^mu nu nabla_mu D_nu - R/4) psi'] = (Dsl(Dsl(psi)) - boxs + RsC*psi/4).applyfunc(
    lambda x: sp.simplify(sp.expand(x))) == sp.zeros(16, 1)
chi = sp.Matrix([sp.Function('c%d' % i)(*Y) for i in range(16)])
R['rescaling psi = (a^3 b^4)^(-1/2) chi removes the spin connection: Gamma^mu D_mu psi = (a^3 b^4)^(-1/2) Gamma^mu d_mu chi'] = (
    Dsl(chi/sp.sqrt(sqg)) - sum((GamuC[m]*chi.diff(X[m]) for m in range(D)), sp.zeros(16, 1))/sp.sqrt(sqg)).applyfunc(sp.simplify) == sp.zeros(16, 1)
# Hamiltonian form of the mode equation: chi = e^{i(p.x + q.y)} c(x0)
Pm = sp.symbols('p1:4', real=True); Qm = sp.symbols('q4:8', real=True); mD = sp.symbols('m', positive=True)
alpha = {i: G[0]*Gup[i] for i in range(1, 8)}
hK = sum((alpha[i]*Pm[i - 1]/a for i in (1, 2, 3)), Z16) + sum((alpha[k]*Qm[k - 4]/b for k in (4, 5, 6, 7)), Z16) + G[0]*mD
cvec = sp.Matrix([sp.Function('c%d' % i)(t) for i in range(16)])
phase = sp.exp(sp.I*(sum(Pm[i]*X[i + 1] for i in range(3)) + sum(Qm[k]*X[k + 4] for k in range(4))))
dir_ = (sp.I*Dsl(phase*cvec/sp.sqrt(sqg)) - mD*phase*cvec/sp.sqrt(sqg))
R['mode equation: (i Gamma^mu D_mu - m) psi = 0  <=>  i dc/dx0 = h(x0) c, h = alpha.(p/a, q/b) + beta m'] = (
    (G[0]*dir_*sp.sqrt(sqg)/phase - (sp.I*cvec.diff(t) - hK*cvec))).applyfunc(sp.simplify) == sp.zeros(16, 1)
# conformal invariance (8D conformally flat metric C(x0)^2 eta in conformal time x0)
Cf = sp.Function('C', real=True, positive=True)(t)
gCf = sp.diag(Cf**2, *([-Cf**2]*7))
giCf, chrCf, RiemCf = geometry(gCf)
ECf = sp.diag(*([1/Cf]*8))
omCf = {m: sp.Matrix(D, D, lambda A, B: sp.simplify(sum(Cf*(1 if A == n else 0)*(sp.diff(ECf[n, B], X[m]) + sum(chrCf[n][m][l]*ECf[l, B] for l in range(D)))
                                                       for n in range(D)))) for m in range(D)}
OmCf = {m: sum((((eta17*omCf[m])[A, B])*Gup[A]*Gup[B]/4 for A in range(D) for B in range(D)), Z16).applyfunc(sp.simplify) for m in range(D)}
fflat = sp.Matrix([sp.Function('k%d' % i)(*Y) for i in range(16)])
w = Cf**sp.Rational(-7, 2)*fflat
lhs = sum(((Gup[m]/Cf)*(w.diff(X[m]) + OmCf[m]*w) for m in range(D)), sp.zeros(16, 1))
rhs = Cf**sp.Rational(-9, 2)*sum((Gup[A]*fflat.diff(X[A]) for A in range(D)), sp.zeros(16, 1))
R['conformal weight 7/2: Gamma^mu D_mu (C^(-7/2) f) = C^(-9/2) Gamma^A d_A f for g = C(x0)^2 eta'] = (lhs - rhs).applyfunc(sp.simplify) == sp.zeros(16, 1)

# ======================================================================= 9. modes and pair creation in the three epochs
npf = lambda M: np.array(M.tolist(), dtype=float)
al = {i: npf(alpha[i]) for i in alpha}; be = npf(G[0])
na = np.array([1.0, 2.0, 2.0])/3.0; nb = np.array([1.0, -1.0, 1.0, 1.0])/2.0
An_ = sum(na[i]*al[i + 1] for i in range(3)); Bn_ = sum(nb[k]*al[k + 4] for k in range(4))
Zc = An_ @ Bn_ @ be
R['2x2 reduction: A = alpha.n_a, B = alpha.n_b, beta: squares 1, mutually anticommuting; Z = A B beta: Z^2 = -1, Z^T = -Z, commutes with A, B, beta'] = bool(
    np.allclose(An_ @ An_, np.eye(16)) and np.allclose(Bn_ @ Bn_, np.eye(16)) and np.allclose(An_ @ Bn_ + Bn_ @ An_, 0)
    and np.allclose(An_ @ be + be @ An_, 0) and np.allclose(Bn_ @ be + be @ Bn_, 0) and np.allclose(Zc @ Zc, -np.eye(16))
    and np.allclose(Zc.T, -Zc) and np.allclose(Zc @ An_, An_ @ Zc) and np.allclose(Zc @ Bn_, Bn_ @ Zc) and np.allclose(Zc @ be, be @ Zc))
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
def evolve(hfun, t0, t1, c0):
    sol = solve_ivp(lambda tt_, y: -1j*(hfun(tt_) @ y), (t0, t1), c0.astype(complex), method='DOP853', rtol=1e-11, atol=1e-13)
    return sol.y[:, -1]
def bogoliubov(hfun, t0, t1, n=2):
    """|beta|^2 averaged over the n/2 positive-energy in-states, and |alpha|^2"""
    w0, v0 = np.linalg.eigh(hfun(t0)); w1, v1 = np.linalg.eigh(hfun(t1))
    h = n//2; b2 = a2 = 0.0
    for j in range(h, n):
        c = evolve(hfun, t0, t1, v0[:, j])
        b2 += np.sum(abs(v1[:, :h].conj().T @ c)**2); a2 += np.sum(abs(v1[:, h:].conj().T @ c)**2)
    return b2/h, a2/h
pK = (1 + 2*math.sqrt(2))/7; qK = (2 - 3*math.sqrt(2))/14
kas_b = {}
for (ka, kb) in ((1.0, 1.0), (2.0, 1.0), (1.0, 3.0)):
    lam = (kb/ka)**(1/(qK - pK)); Ksc = ka*lam**(1 - pK)
    hf = lambda tt_, ka=ka, kb=kb: ka*tt_**(-pK)*s1 + kb*tt_**(-qK)*s2
    b2, a2 = bogoliubov(hf, 1e-15*min(1, lam), 300.0*max(1, lam))
    kas_b[(ka, kb)] = (b2, a2, Ksc)
R['Kasner (Einstein branch), massless: |beta|^2 for (|p|, |q|) = (1,1), (2,1), (1,3)'] = tuple(kas_b[k][0] for k in ((1.0, 1.0), (2.0, 1.0), (1.0, 3.0)))
R['Kasner scaling variable K = |p| lambda^(1-p), lambda = (|q|/|p|)^(1/(q-p))'] = tuple(kas_b[k][2] for k in ((1.0, 1.0), (2.0, 1.0), (1.0, 3.0)))
R['Kasner: |alpha|^2 + |beta|^2 = 1 (to 1e-8)'] = all(abs(v[0] + v[1] - 1) < 1e-8 for v in kas_b.values())
ex = (1 - pK)/(qK - pK); kb2 = 2.5; ka2 = (1.0/kb2**ex)**(1/(1 - ex)); lam2 = (kb2/ka2)**(1/(qK - pK))
b2s, _ = bogoliubov(lambda tt_: ka2*tt_**(-pK)*s1 + kb2*tt_**(-qK)*s2, 1e-15*min(1, lam2), 300.0*max(1, lam2))
R['Kasner scaling: (|p|, |q|) = (1.430..., 2.5) has K = 1 and the same |beta|^2 as (1, 1) (to 1e-6)'] = abs(b2s - kas_b[(1.0, 1.0)][0]) < 1e-6
# 16x16 versus 2x2 (Kasner, massive, generic momentum direction)
mK = 0.3; pv = np.array([0.6, 0.0, 0.8]); qv = np.array([0.5, 0.5, -0.5, 0.5])*1.5
h16 = lambda tt_: sum(pv[i]*tt_**(-pK)*al[i + 1] for i in range(3)) + sum(qv[k]*tt_**(-qK)*al[k + 4] for k in range(4)) + mK*be
h2p = lambda tt_: np.linalg.norm(pv)*tt_**(-pK)*s1 + np.linalg.norm(qv)*tt_**(-qK)*s2 + mK*s3
h2m = lambda tt_: np.linalg.norm(pv)*tt_**(-pK)*s1 + np.linalg.norm(qv)*tt_**(-qK)*s2 - mK*s3
b16, a16 = bogoliubov(h16, 0.05, 40.0, 16); b2p, a2p = bogoliubov(h2p, 0.05, 40.0, 2); b2m, a2m = bogoliubov(h2m, 0.05, 40.0, 2)
R['Kasner, m = 0.3, |p| = 1, |q| = 1.5, x0: 0.05 -> 40: |beta|^2 of the 16x16 evolution (average), of the sector Z = +i (+m sigma3) and of the sector Z = -i (-m sigma3)'] = (b16, b2p, b2m)
R['16x16 evolution = 4 copies of the 2x2 evolution with +m sigma3 and 4 copies with -m sigma3: average agrees (to 1e-8), unitarity; the two sectors differ'] = (
    abs(b16 - (b2p + b2m)/2) < 1e-8 and abs(a16 + b16 - 1) < 1e-8 and abs(b2p - b2m) > 1e-4)
R['isotropic expansion, massless: no pair creation (|beta|^2 < 1e-20)'] = bogoliubov(lambda tt_: tt_**(-0.4)*(0.6*s1 + 0.8*s2), 1e-6, 1e2)[0] < 1e-20
dS = {}
for mH in (0.25, 0.5, 1.0):
    hf = lambda tt_, mH=mH: math.exp(-tt_)*s1 + mH*s3
    dS[mH] = bogoliubov(hf, -math.log(2000.0), math.log(1e6/mH**2))[0]
R['far future (de Sitter, H = 1): |beta|^2 for m/H = 0.25, 0.5, 1'] = tuple(dS[k] for k in (0.25, 0.5, 1.0))
R['de Sitter: |beta|^2 = 1/(exp(2 pi m/H) + 1) (Fermi-Dirac at the Gibbons-Hawking temperature H/2pi), to 1e-7'] = all(
    abs(dS[k] - 1/(math.exp(2*math.pi*k) + 1)) < 1e-7 for k in dS)
lz = {}
for k2 in (0.1, 0.5, 1.0):
    lz[k2] = bogoliubov(lambda e_, k=math.sqrt(k2): k*s1 + e_*s3, -150.0, 150.0)[0]
R['integrator benchmark: Landau-Zener |beta|^2 = exp(-pi kappa^2) (to 1e-6)'] = all(abs(lz[k] - math.exp(-math.pi*k)) < 1e-6 for k in lz)

if __name__ == '__main__':
    bad = [k for k, v in R.items() if isinstance(v, (bool, np.bool_)) and not v]
    nbool = sum(isinstance(v, (bool, np.bool_)) for v in R.values())
    for k, v in R.items():
        print(f'{k} : {v}')
    print(f'boolean checks: {nbool}   failed: {len(bad)}')
    print('VERDICT:', 'ALL CHECKS PASSED' if not bad else 'FAILED: ' + '; '.join(bad))
    sys.exit(1 if bad else 0)
