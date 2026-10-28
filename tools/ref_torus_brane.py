#!/usr/bin/env python3
"""Independent Python checks of section 25 of 7+1.nb: the auxiliary torus and the world brane (definitions and worked examples).

    python3 tools/ref_torus_brane.py

Computes, without reading the notebook:
  the auxiliary torus T^4 = R^4/(Lc Z)^4 -- its size b^8 = det(-g_kl) and proper volume V_perp = b^4 Lc^4, the invariance under
  y = lambda y', the 384 symmetries that fix y = 0 and their orbits on the lattice shells, orthonormal Fourier modes and Parseval's
  identity, the Laplacian, the Kaluza-Klein masses of a scalar (Laplace-Beltrami operator) and of dirac16complex, Jacobi's r_4(k) and
  Weyl's law, the numbers for M8 = 10 TeV, the degrees of freedom of the graviton on the torus;
  Newton's law on the torus -- the eight-dimensional potential, the four-dimensional Fourier transform of (r^2 + y^2)^(-5/2), the Poisson
  duality "images = Kaluza-Klein modes", the Epstein value Z(5/2), the near and far fields, and graviton + radion = 5/3;
  Kaluza-Klein momentum in the bulk and at the brane, the delta function of the torus;
  the world brane W = {x4 = ... = x7 = 0} -- induced metrics, the volume element, Christoffel symbols at the brane (extrinsic curvature),
  the explicit Gamma matrices and the intertwiner U (Schur average; Schur's lemma by the null space), the action of the brane Lorentz
  group and of the transverse rotations, the Spin(4) structure of the transverse index, the chirality content, an explicit spinor, the
  Kaluza-Klein operator; the brane tension and smeared sources, bulk radiation, the normalcy temperature, the tension as vacuum energy,
  the Green's function of the four transverse directions.
Every line 'key : value' is printed; boolean checks are summarized at the end (exit code 0 if all are True)."""
import itertools, math, os, sys
from collections import OrderedDict
import numpy as np
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
R = OrderedDict()

# ---------------------------------------------------------------------------------------------------- constants (as in sections 17-24)
MPl = 2.435e18; M8 = 1.0e4; Vperp = MPl**2/M8**6; lp = Vperp**0.25
hbarc_m = 1.973269804e-16                     # GeV m
r_proton_m = 0.8414e-15                       # CODATA 2018 proton rms charge radius
a0_m = 5.29177210903e-11                      # CODATA 2018 Bohr radius
mH = 125.20; mr = 10.0
hbar = 6.582119569e-25; Mpc_km = 3.0856775814913673e19; GeV_per_K = 8.617333262e-14; Tcmb = 2.7255*GeV_per_K
H0 = 67.4/Mpc_km*hbar; Om = 0.315
Og = (math.pi**2/15)*Tcmb**4/(3*H0**2*MPl**2); Or = Og*(1 + 7/8*(4/11)**(4/3)*3.044); OL = 1 - Om - Or
rhoL = 3*H0**2*MPl**2*OL

# ==================================================================================================== A. the torus: definition
x0 = sp.Symbol('x0', real=True)
X = sp.symbols('x0:8', real=True)
a = sp.Function('a')(X[0]); b = sp.Function('b')(X[0])
Lc, lam, bb = sp.symbols('L_c lambda b', positive=True)
gCos = sp.diag(1, -a**2, -a**2, -a**2, -b**2, -b**2, -b**2, -b**2)
gk = gCos[4:, 4:]
vol = sp.integrate(sp.sqrt(sp.factor((-gk).det())).subs(b, bb), (X[4], 0, Lc), (X[5], 0, Lc), (X[6], 0, Lc), (X[7], 0, Lc))
R['torus: det(-g_kl) = b^8 for the cosmological metric; the proper volume of [0, Lc)^4 is b^4 Lc^4'] = (
    sp.simplify((-gk).det() - b**8) == 0 and sp.simplify(vol - bb**4*Lc**4) == 0)
dy = sp.Matrix(sp.symbols('dy4:8', real=True))
R['only b Lc is physical: y = lambda y\' maps (Lc, b) to (Lc/lambda, lambda b); b^2 dy.dy, b Lc and b^4 Lc^4 unchanged'] = (
    sp.simplify(bb**2*(lam*dy).dot(lam*dy) - (lam*bb)**2*dy.dot(dy)) == 0 and sp.simplify((Lc/lam)*(lam*bb) - Lc*bb) == 0
    and sp.simplify((Lc/lam)**4*(lam*bb)**4 - Lc**4*bb**4) == 0)
# the symmetries fixing y = 0: signed permutation matrices
SP = []
for p in itertools.permutations(range(4)):
    for s in itertools.product((1, -1), repeat=4):
        M = np.zeros((4, 4), dtype=int)
        for i in range(4): M[i, p[i]] = s[i]
        SP.append(M)
SPset = {tuple(M.flatten()) for M in SP}
shell = lambda k: [v for v in itertools.product(range(-3, 4), repeat=4) if sum(t*t for t in v) == k]
def orbits(k):
    cl = OrderedDict()
    for v in shell(k):
        cl.setdefault(tuple(sorted(abs(t) for t in v)), []).append(v)
    return list(cl.values())
def orbit_ok(k):
    for c in orbits(k):
        o = {tuple(M @ np.array(c[0])) for M in SP}
        if o != set(c): return False
    return True
R['symmetries fixing y = 0: 384 = 2^4 4! distinct signed permutation matrices, orthogonal integer matrices; shells 1, 2, 3, 4, 9 have orbits {8}, {24}, {32}, {8, 16}, {8, 96}'] = (
    len(SP) == 384 and len(SPset) == 384 and all((M @ M.T == np.eye(4, dtype=int)).all() for M in SP)
    and all(orbit_ok(k) for k in (1, 2, 3, 4, 9))
    and [sorted(len(c) for c in orbits(k)) for k in (1, 2, 3, 4, 9)] == [[8], [24], [32], [8, 16], [8, 96]])
R['NUM orbit sizes of the shells k = 1..12'] = [sorted(len(c) for c in orbits(k)) for k in range(1, 13)]

# ==================================================================================================== B. fields on the torus
Y = sp.Matrix(X[4:])
def q(n): return sp.Matrix(n)*2*sp.pi/Lc
V = bb**4*Lc**4
def phi(n): return sp.exp(sp.I*(q(n).T*Y)[0])/sp.sqrt(V)
def phic(n): return sp.exp(-sp.I*(q(n).T*Y)[0])/sp.sqrt(V)
def tint(f):
    return sp.simplify(sp.integrate(bb**4*f, (X[4], 0, Lc), (X[5], 0, Lc), (X[6], 0, Lc), (X[7], 0, Lc)))
on = [tint(phic(n)*phi(m)) for n, m in (((1, 0, -2, 3), (1, 0, -2, 3)), ((1, 0, -2, 3), (1, 1, -2, 3)), ((0, 0, 0, 0), (0, 2, 0, 0)))]
R['Fourier modes phi_n = exp(i q_n.y)/sqrt(V_perp), q_n = 2 pi n/Lc: orthonormal for b^4 d^4y (three exact pairs)'] = on == [1, 0, 0]
f = sp.prod([sp.cos(sp.pi*Y[k]/Lc)**2 for k in range(4)])
w1 = {0: sp.Rational(1, 2), 1: sp.Rational(1, 4), -1: sp.Rational(1, 4)}
def coef(n): return sp.sqrt(V)*sp.prod([w1[t] for t in n])
coef_ok = all(sp.simplify(tint(phic(n)*f) - coef(n)) == 0 for n in ((0, 0, 0, 0), (1, 0, -1, 0), (1, 1, 1, -1)))
pars = sp.simplify(sum(coef(n)**2 for n in itertools.product((-1, 0, 1), repeat=4)))
R['Parseval: f = prod cos^2(pi y_k/Lc) has 81 coefficients sqrt(V_perp) prod w(n_k) (3 checked by integration), Sum |c_n|^2 = Int b^4 f^2 = (3/8)^4 V_perp'] = (
    coef_ok and sp.simplify(pars - tint(f**2)) == 0 and sp.simplify(pars - sp.Rational(3, 8)**4*V) == 0)
n0 = (2, -1, 0, 3)
lap = -sum(sp.diff(phi(n0), Y[k], 2) for k in range(4))/bb**2
R['Laplacian of the torus: -(1/b^2) Sum d^2/dy_k^2 phi_n = (|q_n|^2/b^2) phi_n'] = sp.simplify(lap - (q(n0).dot(q(n0))/bb**2)*phi(n0)) == 0
E, m = sp.symbols('E m', positive=True); k1, k2, k3 = sp.symbols('k1:4', real=True)
gS = sp.diag(1, -1, -1, -1, -bb**2, -bb**2, -bb**2, -bb**2); gSi = gS.inv(); sg = sp.sqrt(-gS.det())
nq = (1, -2, 0, 2)
wave = sp.exp(-sp.I*E*X[0] + sp.I*(k1*X[1] + k2*X[2] + k3*X[3]) + sp.I*(q(nq).T*Y)[0])
box = sum(sp.diff(sg*gSi[A, A]*sp.diff(wave, X[A]), X[A]) for A in range(8))/sg
R['scalar Kaluza-Klein masses: (Box + m^2) exp(-iEx0 + ik.x + iq.y) = (-E^2 + |k|^2 + |q|^2/b0^2 + m^2) exp(...) on diag(1,-1,-1,-1,-b0^2,...)'] = (
    sp.simplify((box + m**2*wave)/wave - (-E**2 + k1**2 + k2**2 + k3**2 + q(nq).dot(q(nq))/bb**2 + m**2)) == 0)
# graviton degrees of freedom: symmetric traceless tensors of the little group
sym = lambda d: d*(d + 1)//2
massless = (sym(2) - 1, 2*4, sym(4))          # helicity 2; four vectors (2 states each); ten scalars (one trace combination)
massive = (sym(3) - 1, 3*3, (sym(3) - 1) + 1)  # spin 2; three spin-1; five + one spin-0
R['graviton on the torus: D(D-3)/2 = 20 (D = 8) = 2 + 8 + 10 (zero mode: graviton, 4 graviphotons, 10 moduli) = 5 + 9 + 6 (massive level)'] = (
    8*5//2 == 20 and sym(6) - 1 == 20 and sum(massless) == 20 and massless == (2, 8, 10) and sum(massive) == 20 and massive == (5, 9, 6))

# ==================================================================================================== C. counting the levels
kmax = int(math.floor((mr*lp/(4*math.pi))**2))           # |n|^2 <= kmax  <=>  mu_n <= m_r/2 (as in section 22.2)
K = kmax                                                   # the shells k = 0..kmax
sig = np.zeros(K + 1, dtype=np.int64)
for d_ in range(1, K + 1): sig[d_::d_] += d_
r4 = 8*sig.copy(); r4[4::4] -= 32*sig[1:K//4 + 1]; r4[0] = 1
g = np.arange(-10, 11); grid = (g[:, None, None, None]**2 + g[None, :, None, None]**2 + g[None, None, :, None]**2 + g[None, None, None, :]**2).ravel()
cnt = np.bincount(grid, minlength=401)
R['Jacobi: #{n in Z^4 : |n|^2 = k} = r_4(k) = 8 sigma(k) - 32 sigma(k/4), direct count for k = 1..100'] = bool((cnt[1:101] == r4[1:101]).all())
R['NUM r_4(k), k = 0..16'] = [int(x) for x in r4[:17]]
Rs = (10, 30, 100, 300, 1000)
NR = [int(r4[:Rv*Rv + 1].sum()) for Rv in Rs]
dev = [NR[i]/(math.pi**2*Rv**4/2) - 1 for i, Rv in enumerate(Rs)]
R['NUM N(R) for R = 10, 30, 100, 300, 1000'] = NR
R['NUM N(R)/(pi^2 R^4/2) - 1 for R = 10, 30, 100, 300, 1000'] = dev
R['Weyl\'s law: N(10) = 49689, N(100) = 493490641; |N(R)/(pi^2 R^4/2) - 1| < 1e-2, 1e-3, 3e-5, 5e-6, 1e-6'] = (
    NR[0] == 49689 and NR[2] == 493490641 and all(abs(dv) < t for dv, t in zip(dev, (1e-2, 1e-3, 3e-5, 5e-6, 1e-6))))
R['NUM V_perp (GeV^-4), l_perp (GeV^-1), l_perp (m), l_perp (fm), 2 pi/l_perp (MeV)'] = (Vperp, lp, lp*hbarc_m, lp*hbarc_m*1e15, 2*math.pi/lp*1e3)
MD = (2*math.pi)**(2/3)*M8
R['NUM l_perp/proton radius, Bohr radius/l_perp, R = l_perp/(2 pi) (GeV^-1, m), M_D = (2 pi)^(2/3) M8 (GeV)'] = (
    lp*hbarc_m/r_proton_m, a0_m/(lp*hbarc_m), lp/(2*math.pi), lp/(2*math.pi)*hbarc_m, MD)
Nmu = lambda mu: Vperp*mu**4/(32*math.pi**2)
lat5 = int(r4[:kmax + 1].sum())
R['NUM N(mu) = V_perp mu^4/(32 pi^2) for mu = 5 GeV, m_h, M8; lattice count below 5 GeV'] = (Nmu(5.0), Nmu(mH), Nmu(M8), lat5)
R['numbers today: V_perp = 5.929e12 GeV^-4, l_perp = 1560.4 GeV^-1 = 3.079e-13 m, 2 pi/l_perp = 4.027 MeV; N(5 GeV) = 1.17e13 (lattice to 1e-5); N(M8) = M_Pl^2/(32 pi^2 M8^2)'] = (
    abs(Vperp/5.929e12 - 1) < 1e-3 and abs(lp/1560.4 - 1) < 1e-4 and abs(lp*hbarc_m/3.079e-13 - 1) < 1e-3 and abs(2*math.pi/lp*1e3/4.027 - 1) < 1e-3
    and abs(lat5/Nmu(5.0) - 1) < 1e-5 and abs(Nmu(5.0)/1.17e13 - 1) < 0.01 and abs(Nmu(M8)/(MPl**2/(32*math.pi**2*M8**2)) - 1) < 1e-12)
lev = [(k, int(r4[k]), int(4*r4[k]), 2*math.pi*math.sqrt(k)/lp*1e3) for k in range(1, 9)]
R['NUM first levels: k, r_4(k), Dirac fermions 4 r_4(k), mass 2 pi sqrt(k)/l_perp (MeV), k = 1..8'] = lev

# ---------------------------------------------------------------------------------------------------- C'. from eight to four dimensions
def scalar_curvature_diag(gd, xs):
    """scalar curvature of a diagonal metric diag(gd) in the coordinates xs, with the conventions of 7+1.nb (Landau-Lifshitz):
    R^r_smn = d_m {r sn} - d_n {r sm} + {r ml}{l sn} - {r nl}{l sm}, R_sn = R^r_srn, R = g^sn R_sn"""
    n_ = len(xs)
    gi = [1/gd[i] for i in range(n_)]
    def chr_(r, s_, m_):
        return sp.simplify(gi[r]*(sp.diff(gd[r], xs[m_])*(1 if r == s_ else 0) + sp.diff(gd[r], xs[s_])*(1 if r == m_ else 0)
                                  - sp.diff(gd[s_], xs[r])*(1 if s_ == m_ else 0))/2)
    Ch = [[[chr_(r, s_, m_) for m_ in range(n_)] for s_ in range(n_)] for r in range(n_)]
    def ric(s_, n2):
        return sum(sp.diff(Ch[r][s_][n2], xs[r]) - sp.diff(Ch[r][s_][r], xs[n2])
                   + sum(Ch[r][r][l]*Ch[l][s_][n2] - Ch[r][n2][l]*Ch[l][s_][r] for l in range(n_)) for r in range(n_))
    return sp.simplify(sum(gi[i]*ric(i, i) for i in range(n_)))
Ha_ = sp.Function('a')(X[0])
R4 = scalar_curvature_diag([1, -Ha_**2, -Ha_**2, -Ha_**2], X[:4])
b0s = sp.Symbol('b_0', positive=True); apos = sp.Symbol('a_pos', positive=True)
R8 = scalar_curvature_diag([1, -Ha_**2, -Ha_**2, -Ha_**2, -b0s**2, -b0s**2, -b0s**2, -b0s**2], X)
Hs, dHs = sp.symbols('H dH')
toH_ = lambda e: sp.simplify(e.subs(sp.diff(Ha_, X[0], 2), (dHs + Hs**2)*Ha_).subs(sp.diff(Ha_, X[0]), Hs*Ha_))
R['from eight to four dimensions: for diag(1, -a^2, -a^2, -a^2, -b0^2, ..., -b0^2) (static torus) R = R4 = -6(Ha\' + 2 Ha^2), the scalar curvature of the FRW metric diag(1, -a^2, -a^2, -a^2); sqrt|g| = a^3 b0^4'] = (
    sp.simplify(toH_(R8) - toH_(R4)) == 0 and sp.simplify(toH_(R4) + 6*(dHs + 2*Hs**2)) == 0
    and sp.simplify(sp.sqrt(-sp.diag(1, -apos**2, -apos**2, -apos**2, -b0s**2, -b0s**2, -b0s**2, -b0s**2).det()) - apos**3*b0s**4) == 0)
# the quaternion subalgebra of the octonions (the brane directions): its products in the octonion table of 7+1.nb
from ref71 import CO, mul, basis
qprod = {(i, j): mul(CO, basis(i), basis(j)) for i in range(4) for j in range(4)}
closed = all(all(qprod[(i, j)][k] == 0 for k in range(4, 8)) for i in range(4) for j in range(4))
re_q = [[qprod[(i, j)][0] for j in range(4)] for i in range(4)]
R['quaternion subalgebra of the brane directions: e1 e2 = -e3, e2 e3 = -e1, e3 e1 = -e2 (e_i e_j = -delta_ij - eps_ijk e_k), span{e0, ..., e3} closed under the octonion product, Re(e_A e_B) = diag(1, -1, -1, -1)'] = (
    closed and re_q == [[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]
    and list(qprod[(1, 2)]) == [0, 0, 0, -1, 0, 0, 0, 0] and list(qprod[(2, 3)]) == [0, -1, 0, 0, 0, 0, 0, 0] and list(qprod[(3, 1)]) == [0, 0, -1, 0, 0, 0, 0, 0])
R['NUM quaternion products e1 e2, e2 e3, e3 e1: coefficients of e3, e1, e2'] = (int(qprod[(1, 2)][3]), int(qprod[(2, 3)][1]), int(qprod[(3, 1)][2]))

# ==================================================================================================== D. Newton's law on the torus
rr, kap, Mm, GN, Vs, MP = sp.symbols('r kappa M G_N V_perp M_Pl', positive=True)
Om6 = 2*sp.pi**sp.Rational(7, 2)/sp.gamma(sp.Rational(7, 2)); Om2 = 2*sp.pi**sp.Rational(3, 2)/sp.gamma(sp.Rational(3, 2))
R['Newton in 8D: Omega_6 = 16 pi^3/15, (1/r^6) d/dr (r^6 d/dr r^-5) = 0, Phi = -(5/6) kappa M/(5 Omega_6 r^5) = -5 kappa M/(32 pi^3 r^5); D = 4: -G M/r'] = (
    sp.simplify(Om6 - 16*sp.pi**3/15) == 0 and sp.simplify(sp.diff(rr**6*sp.diff(rr**-5, rr), rr)) == 0
    and sp.simplify(kap*sp.Rational(5, 6)/(5*Om6) - 5*kap/(32*sp.pi**3)) == 0
    and sp.simplify(8*sp.pi*GN*sp.Rational(1, 2)/(1*Om2) - GN) == 0)
# four-dimensional Fourier transform of (r^2 + y^2)^(-5/2): (4 pi^2/p) Int_0^oo y^2 J_1(p y) f(y) dy = (4 pi^2/(3 r)) e^(-p r)
mp.mp.dps = 30
def ft(p, r):
    return 4*mp.pi**2/p*mp.quadosc(lambda y: y**2*mp.besselj(1, p*y)*(r*r + y*y)**-2.5, [0, mp.inf], omega=p)
fts = [(p_, r_, ft(mp.mpf(p_), mp.mpf(r_)), 4*mp.pi**2/(3*mp.mpf(r_))*mp.exp(-mp.mpf(p_)*mp.mpf(r_))) for p_, r_ in ((1.3, 0.7), (2*math.pi, 0.5), (0.4, 2.0))]
R['4D Fourier transform: Int d^4y exp(-ip.y) (r^2 + |y|^2)^(-5/2) = (4 pi^2/(3 r)) exp(-p r) (mpmath, three (p, r), to 1e-20)'] = all(abs(x/y_ - 1) < 1e-20 for *_, x, y_ in fts)
kk = np.arange(K + 1, dtype=np.float64)
def S_img(r):
    u0 = r*r + K + 0.5
    return float(np.sum(r4*(r*r + kk)**-2.5)) + math.pi**2*(2/math.sqrt(u0) - 2*r*r/(3*u0**1.5))
def S_kk(r):
    m_ = 20000
    return 4*math.pi**2/(3*r)*float(np.sum(r4[:m_]*np.exp(-2*math.pi*r*np.sqrt(kk[:m_]))))
pois = []
for r_ in (0.25, 0.5, 1.0, 2.0):
    si, sk = S_img(r_), S_kk(r_)
    pois.append((si, sk, si/sk - 1))
R['NUM images S(r), Kaluza-Klein (4 pi^2/(3 r)) Sum exp(-2 pi r |m|), relative difference, at r = 1/4, 1/2, 1, 2'] = pois
R['Poisson duality: Sum_n (r^2 + |n|^2)^(-5/2) = (4 pi^2/(3 r)) Sum_m exp(-2 pi r |m|) to 1e-8 at r = 1/4, 1/2, 1, 2'] = all(abs(t[2]) < 1e-8 for t in pois)
Z52 = float(7*mp.zeta(2.5)*mp.zeta(1.5))
Z52s = float(np.sum(r4[1:]*kk[1:]**-2.5)) + 2*math.pi**2/math.sqrt(K + 0.5)
R['NUM Z(5/2) = 7 zeta(5/2) zeta(3/2), shell sum with tail'] = (Z52, Z52s)
near = [(r_**5*S_img(r_) - 1)/(r_**5*Z52) for r_ in (0.05, 0.1)]
far = [(3*r_*S_kk(r_)/(4*math.pi**2) - 1)/(8*math.exp(-2*math.pi*r_)) for r_ in (1.0, 2.0, 3.0)]
R['NUM near field (r^5 S - 1)/(r^5 Z(5/2)) at r = 1/20, 1/10; far field (3 r S/(4 pi^2) - 1)/(8 exp(-2 pi r)) at r = 1, 2, 3'] = near + far
R['Z(5/2) = 7 zeta(5/2) zeta(3/2) = Sum r_4(k) k^-5/2 (to 1e-8); near field ratio -> 1 (0.98 < . < 1 at r = 1/10, 1/20); far field ratio -> 1 (to 2e-2 at r = 2, 3)'] = (
    abs(Z52s/Z52 - 1) < 1e-8 and all(0.98 < x < 1 for x in near) and near[0] > near[1] and all(abs(x - 1) < 2e-2 for x in far[1:]))
G4 = kap/(8*sp.pi*Vs)
R['graviton + radion: (5 kappa M/(32 pi^3))(4 pi^2/(3 V_perp r)) = (5/3) G4 M/r; radion (1/(sqrt3 M_Pl))^2/(4 pi) = (2/3)/(8 pi M_Pl^2)'] = (
    sp.simplify((5*kap/(32*sp.pi**3))*(4*sp.pi**2/(3*Vs)) - sp.Rational(5, 3)*G4) == 0
    and sp.simplify((1/(sp.sqrt(3)*MP))**2/(4*sp.pi) - sp.Rational(2, 3)/(8*sp.pi*MP**2)) == 0)

# ==================================================================================================== E. Kaluza-Klein momentum, bulk and brane
nA, mA = (1, 0, 0, 2), (0, 1, -1, 0)
lA = tuple(i + j for i, j in zip(nA, mA)); lB = (1, 1, -1, 1)
v1 = tint(phi(nA)*phi(mA)*phic(lA)); v2 = tint(phi(nA)*phi(mA)*phic(lB))
at0 = sp.simplify((phi(nA)*phi(mA)*phic(lB)).subs({Y[k]: 0 for k in range(4)}))
R['Kaluza-Klein momentum: bulk vertex Int b^4 phi_n phi_m phi_l^* = 1/sqrt(V_perp) for l = n + m, 0 otherwise; at y = 0 phi_n phi_m phi_l^* = V_perp^(-3/2)'] = (
    sp.simplify(v1 - 1/sp.sqrt(V)) == 0 and v2 == 0 and sp.simplify(at0 - V**sp.Rational(-3, 2)) == 0)
yv, nn = sp.symbols('y n', real=True)
dc = sp.integrate(sp.exp(-sp.I*2*sp.pi*nn*yv/Lc)*sp.DiracDelta(yv), (yv, -Lc/2, Lc/2))
rho4 = sp.Symbol('rho_4', positive=True)
sm = sp.integrate(bb**4*rho4*sp.DiracDelta(X[4])*sp.DiracDelta(X[5])*sp.DiracDelta(X[6])*sp.DiracDelta(X[7])/bb**4,
                  (X[4], -Lc/2, Lc/2), (X[5], -Lc/2, Lc/2), (X[6], -Lc/2, Lc/2), (X[7], -Lc/2, Lc/2))
R['delta function of the torus: Fourier coefficients Int exp(-i q_n y) delta(y) dy = 1; Int b^4 d^4y rho4 delta^4(y)/b^4 = rho4 (average rho4/V_perp)'] = (
    sp.simplify(dc - 1) == 0 and sp.simplify(sm - rho4) == 0)

# ==================================================================================================== F. the world brane: geometry
J = sp.Matrix(8, 4, lambda A, mu: 1 if A == mu else 0)        # dX^A/dxi^mu for X = (xi, 0, 0, 0, 0)
on_brane = {X[4]: 0, X[5]: 0, X[6]: 0, X[7]: 0}
hCos = sp.simplify(J.T*gCos.subs(on_brane)*J)
nP = 1 + X[0]*X[4]/5 + X[0]**2/11; aP = 1 + X[0]/3 + X[4]**2/7; bP = 2 - X[0]*X[4]/4 + X[0]**2/9 + X[4]/13
gPoly = sp.diag(nP**2, -aP**2, -aP**2, -aP**2, -bP**2, -bP**2, -bP**2, -bP**2)
hPoly = sp.simplify(J.T*gPoly.subs(on_brane)*J)
R['induced metric h = J^T g J: FRW diag(1, -a^2, -a^2, -a^2) for the cosmological metric; diag(n^2, -a^2, -a^2, -a^2) at x4 = 0 with n = 1 + x0^2/11, a = 1 + x0/3 for the inhomogeneous metric'] = (
    sp.simplify(hCos - sp.diag(1, -a**2, -a**2, -a**2)) == sp.zeros(4) and
    sp.simplify(hPoly - sp.diag((1 + X[0]**2/11)**2, -(1 + X[0]/3)**2, -(1 + X[0]/3)**2, -(1 + X[0]/3)**2)) == sp.zeros(4))
ap, bp = sp.symbols('a_p b_p', positive=True)
volC = sp.simplify((sp.sqrt(-gCos.det()) - sp.sqrt(-hCos.det())*sp.sqrt((-gCos[4:, 4:]).det())).subs({a: ap, b: bp}))
pt = {X[0]: sp.Rational(1, 2), X[4]: 0}
gP0 = gPoly.subs(pt); hP0 = hPoly.subs(pt)
volP = sp.sqrt(-gP0.det()) - sp.sqrt(-hP0.det())*sp.sqrt((-gP0[4:, 4:]).det())
Fy = sp.exp(X[4] + 2*X[5])*sp.cos(X[6] - X[7])
tf = sp.integrate(bp**4*(sp.DiracDelta(X[4])*sp.DiracDelta(X[5])*sp.DiracDelta(X[6])*sp.DiracDelta(X[7])/bp**4)*Fy,
                  (X[4], -1, 1), (X[5], -1, 1), (X[6], -1, 1), (X[7], -1, 1))
R['brane volume element: sqrt|g| = sqrt|h| b^4 at the brane (cosmological metric; inhomogeneous metric at x0 = 1/2, exact), so Int_W d^4x sqrt|h| L = Int d^8X sqrt|g| delta^4(y)/b^4 L; test Int b^4 (delta^4/b^4) F = F(0)'] = (
    volC == 0 and sp.simplify(volP) == 0 and sp.simplify(tf - 1) == 0)
def christoffel(gm, k, mu, nu, at):
    """{k mu nu} = (1/2) g^ks (d_mu g_s nu + d_nu g_s mu - d_s g_mu nu) for a DIAGONAL metric gm, evaluated at the point 'at'"""
    assert gm.is_diagonal()
    return sp.simplify(((sp.diff(gm[k, mu], X[nu]) + sp.diff(gm[k, nu], X[mu]) - sp.diff(gm[mu, nu], X[k]))/(2*gm[k, k])).subs(at))
chrC = [christoffel(gCos, k, mu, nu, on_brane) for k in range(4, 8) for mu in range(4) for nu in range(4)]
chrP = {(k, mu, nu): sp.nsimplify(christoffel(gPoly, k, mu, nu, pt)) for k in range(4, 8) for mu in range(4) for nu in range(mu, 4)}
nz = {key: v for key, v in chrP.items() if v != 0}
R['Christoffel symbols {k mu nu} at the brane: all 0 for the cosmological metric (totally geodesic); for the inhomogeneous metric at x0 = 1/2 only {4 00} = 1458/58619 = n d4n/b^2'] = (
    all(sp.simplify(c) == 0 for c in chrC) and list(nz.keys()) == [(4, 0, 0)] and nz[(4, 0, 0)] == sp.Rational(1458, 58619)
    and sp.simplify((nP*sp.diff(nP, X[4])/bP**2).subs(pt) - sp.Rational(1458, 58619)) == 0)
R['NUM Christoffel {4 00} at the brane (inhomogeneous metric, x0 = 1/2)'] = float(sp.Rational(1458, 58619))

# ==================================================================================================== G. dirac16complex at the brane
from ref71 import G17
I = sp.I
s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -I], [I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]]); Z2 = sp.zeros(2)
def blk(A_, B_, C_, D_): return sp.Matrix(sp.BlockMatrix([[A_, B_], [C_, D_]]))
gam = [sp.diag(1, 1, -1, -1)] + [blk(Z2, s, -s, Z2) for s in (s1, s2, s3)]
g5 = sp.simplify(I*gam[0]*gam[1]*gam[2]*gam[3])
rho = [I*gam[0], -gam[1], -gam[2], -gam[3]]
kp = sp.kronecker_product
I4, I16 = sp.eye(4), sp.eye(16)
Rb = [kp(gam[mu], I4) for mu in range(4)] + [kp(g5, rho[k]) for k in range(4)]
eta = sp.diag(1, -1, -1, -1, -1, -1, -1, -1)
G = [G17[A] for A in range(8)]
R['R_A = (gamma_mu (x) 1_4, gamma5 (x) rho_k) and Gamma17[A] both satisfy X_A X_B + X_B X_A = 2 eta_AB 1_16'] = all(
    Rb[A]*Rb[B] + Rb[B]*Rb[A] == 2*eta[A, B]*I16 and G[A]*G[B] + G[B]*G[A] == 2*eta[A, B]*I16 for A in range(8) for B in range(8))
def sperm(M):
    return all(sum(1 for i in range(16) if M[i, j] != 0) == 1 and sum(1 for j2 in range(16) if M[j, j2] != 0) == 1 for j in range(16)) and set(M) <= {0, 1, -1}
sp_codes = [[int(sum(G[A][i, j]*(i + 1) for i in range(16))) for j in range(16)] for A in range(8)]
R['every Gamma17[A] is a signed permutation matrix'] = all(sperm(G[A]) for A in range(8))
R['NUM Gamma17[A] as signed permutations: column j has the entry sign(c) in row |c|, c listed for A = 0..7'] = sp_codes
# the intertwiner U: Schur average over the 256 Clifford products (seed: the first matrix unit that gives a non-zero sum)
Gn = [np.array(M.tolist(), dtype=complex) for M in G]; Rn = [np.array(M.tolist(), dtype=complex) for M in Rb]
subs_ = [c for r_ in range(9) for c in itertools.combinations(range(8), r_)]
def cprod(ms, Kidx):
    P = np.eye(16, dtype=complex)
    for A in Kidx: P = P @ ms[A]
    return P
GP = [cprod(Gn, Kidx) for Kidx in subs_]; RPi = [np.linalg.inv(cprod(Rn, Kidx)) for Kidx in subs_]
for j in range(256):
    s = np.zeros((16, 16), dtype=complex); s[j//16, j % 16] = 1
    u = sum(GP[i] @ s @ RPi[i] for i in range(256))
    if np.abs(u).max() > 1e-9: break
seed = (j//16 + 1, j % 16 + 1)
u = u/np.sqrt(np.trace(u.conj().T @ u).real/16)
U8n = 2*math.sqrt(2)*u
U8 = sp.Matrix(16, 16, lambda i, k: int(round(U8n[i, k].real)) + I*int(round(U8n[i, k].imag)))
R['intertwiner: Schur average with the seed matrix unit (%d, %d); 2 sqrt(2) U has entries 0, +-1, +-i (rounding error < 1e-12) with 8 non-zero entries in each row and column; Gamma17[A].U8 = U8.R_A and U8^dag.U8 = 8 1_16 exactly' % seed] = (
    np.abs(U8n - np.array(U8.tolist(), dtype=complex)).max() < 1e-12 and set(U8) <= {0, 1, -1, I, -I}
    and all(sum(1 for k in range(16) if U8[i, k] != 0) == 8 and sum(1 for k in range(16) if U8[k, i] != 0) == 8 for i in range(16))
    and all(G[A]*U8 == U8*Rb[A] for A in range(8)) and sp.expand(U8.H*U8) == 8*I16)
# Schur's lemma: the intertwiners X (Gamma_A X = X R_A for all A) form a one-dimensional space
L = np.vstack([np.kron(np.eye(16), Gn[A]) - np.kron(Rn[A].T, np.eye(16)) for A in range(8)])
sv = np.linalg.svd(L, compute_uv=False)
R['Schur\'s lemma: the linear map X -> (Gamma_A X - X R_A)_A has a one-dimensional null space (one singular value < 1e-10, the next > 0.1)'] = (
    int(np.sum(sv < 1e-10)) == 1 and sorted(sv)[1] > 0.1)
R['NUM Re(2 sqrt(2) U)'] = [[int(sp.re(U8[i, k])) for k in range(16)] for i in range(16)]
R['NUM Im(2 sqrt(2) U)'] = [[int(sp.im(U8[i, k])) for k in range(16)] for i in range(16)]
U = U8/(2*sp.sqrt(2)); Ud = U.H
ok_mn = all(sp.simplify(G[mu]*G[nu] - U*kp(gam[mu]*gam[nu], I4)*Ud) == sp.zeros(16) for mu in range(4) for nu in range(4))
ok_kl = all(sp.simplify(G[4 + k]*G[4 + l] - U*kp(I4, rho[k]*rho[l])*Ud) == sp.zeros(16) for k in range(4) for l in range(4))
ok_mk = all(sp.simplify(G[mu]*G[4 + k] - U*kp(gam[mu]*g5, rho[k])*Ud) == sp.zeros(16) for mu in range(4) for k in range(4))
R['brane Lorentz group on the brane spinor index, transverse rotations on a: Gamma_mu Gamma_nu = U (gamma_mu gamma_nu (x) 1) U^dag, Gamma_k Gamma_l = U (1 (x) rho_k rho_l) U^dag, Gamma_mu Gamma_k = U (gamma_mu gamma5 (x) rho_k) U^dag'] = ok_mn and ok_kl and ok_mk
r4567 = rho[0]*rho[1]*rho[2]*rho[3]
P = {s_: (I4 + s_*r4567)/2 for s_ in (1, -1)}
even = [(0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2)]
R['Spin(4) on the transverse index: rho4567 = -gamma5 (4x4) commutes with rho_k rho_l, anticommutes with rho_k; (rho_k rho_l + s rho_m rho_n) P_s = 0 and (rho_k rho_l - s rho_m rho_n) P_s has rank 2 for even (k,l,m,n): 4 = (2,1) + (1,2)'] = (
    r4567 == -g5 and all(r4567*rho[k]*rho[l] == rho[k]*rho[l]*r4567 for k in range(4) for l in range(4))
    and all(r4567*rho[k] == -rho[k]*r4567 for k in range(4))
    and all(sp.expand((rho[p_[0]]*rho[p_[1]] + s_*rho[p_[2]]*rho[p_[3]])*P[s_]) == sp.zeros(4) for p_ in even for s_ in (1, -1))
    and all(((rho[p_[0]]*rho[p_[1]] - s_*rho[p_[2]]*rho[p_[3]])*P[s_]).rank() == 2 for p_ in even for s_ in (1, -1)))
chi = I*sp.prod(G)
PL17 = (I16 - chi)/2
PL4, PR4 = (I4 - g5)/2, (I4 + g5)/2
R['chirality: chi = i Gamma_0...Gamma_7 = U (gamma5 (x) rho4567) U^dag; P_L17 = U (P_L (x) P_+ + P_R (x) P_-) U^dag, each part of rank 4'] = (
    sp.simplify(chi - U*kp(g5, r4567)*Ud) == sp.zeros(16) and sp.simplify(PL17 - U*(kp(PL4, P[1]) + kp(PR4, P[-1]))*Ud) == sp.zeros(16)
    and kp(PL4, P[1]).rank() == 4 and kp(PR4, P[-1]).rank() == 4)
uB = sp.Matrix([1, 2, 0, -1])
zeta = kp(uB, sp.Matrix([1, 0, 0, 0]))
psi = U*zeta
bar = lambda v: v.H*G[0]
etaU = [1, -1, -1, -1]
sb = [sp.nsimplify(sp.simplify((bar(psi)*etaU[mu]*G[mu]*psi)[0])) for mu in range(4)]
ub = [sp.simplify((uB.H*gam[0]*etaU[mu]*gam[mu]*uB)[0]) for mu in range(4)]
R['example psi = U.(u (x) e_1), u = (1, 2, 0, -1): psibar psi = ubar u = 4 and psibar Gamma^mu psi = ubar gamma^mu u'] = (
    sp.simplify((bar(psi)*psi)[0] - 4) == 0 and sp.simplify((uB.H*gam[0]*uB)[0] - 4) == 0 and all(sp.simplify(x - y_) == 0 for x, y_ in zip(sb, ub)))
R['NUM example psi = U.(u (x) e_1), u = (1, 2, 0, -1): psibar psi, psibar Gamma^mu psi (mu = 0..3)'] = [4] + [float(x) for x in sb]
rng = np.random.default_rng(2610)
z1 = rng.integers(-4, 5, 16) + 1j*rng.integers(-4, 5, 16); z2 = rng.integers(-4, 5, 16) + 1j*rng.integers(-4, 5, 16)
Un = np.array(U.evalf().tolist(), dtype=complex); G0n = Gn[0]; g0n = np.array(gam[0].tolist(), dtype=complex)
def bar4(z, a_): return np.conj(z.reshape(4, 4)[:, a_]) @ g0n
lhs = np.conj(Un @ z1) @ G0n @ (Un @ z2); rhs = sum(bar4(z1, a_) @ z2.reshape(4, 4)[:, a_] for a_ in range(4))
lhsm = [np.conj(Un @ z1) @ G0n @ (etaU[mu]*Gn[mu]) @ (Un @ z2) for mu in range(4)]
rhsm = [sum(bar4(z1, a_) @ (etaU[mu]*np.array(gam[mu].tolist(), dtype=complex)) @ z2.reshape(4, 4)[:, a_] for a_ in range(4)) for mu in range(4)]
R['four brane Dirac spinors: psibar psi\' = Sum_a psibar_a psi\'_a and psibar Gamma^mu psi\' = Sum_a psibar_a gamma^mu psi\'_a (random integer spinors)'] = (
    abs(lhs - rhs) < 1e-10 and all(abs(x - y_) < 1e-10 for x, y_ in zip(lhsm, rhsm)))
b0, mps = sp.Rational(13, 10), sp.Rational(7, 10); qv = [2*sp.pi*t for t in (2, -1, sp.Rational(1, 2), 3)]
Xkk = kp(g5, sum((qv[k]*rho[k] for k in range(4)), sp.zeros(4)))/b0
R['Kaluza-Klein operator X = gamma5 (x) (rho.q)/b0: X^2 = -|q|^2/b0^2 1_16 and X anticommutes with gamma_mu (x) 1_4'] = (
    sp.simplify(Xkk*Xkk + sum(t*t for t in qv)/b0**2*I16) == sp.zeros(16) and all(sp.simplify(Xkk*Rb[mu] + Rb[mu]*Xkk) == sp.zeros(16) for mu in range(4)))

# ==================================================================================================== H. matter on the brane
u0_, u1_, u2_, u3_, fW = sp.symbols('u0:4 f_W', positive=True)
sq = 1/sp.sqrt(u0_*u1_*u2_*u3_)                    # sqrt|h| for h = diag(1/u0, -1/u1, -1/u2, -1/u3): h^00 = u0, h^ii = -u_i
t00 = (2/sq)*sp.diff(-fW*sq, u0_); t11 = (2/sq)*(-sp.diff(-fW*sq, u1_))
R['brane tension: T_mu nu = (2/sqrt|h|) dS/dh^mu nu = f_W h_mu nu for S = -Int sqrt|h| f_W; T^0_0 = f_W, T^1_1 = f_W (p = -f_W, w = -1)'] = (
    sp.simplify(t00 - fW/u0_) == 0 and sp.simplify(t11 + fW/u1_) == 0 and sp.simplify(u0_*t00) == fW and sp.simplify(-u1_*t11) == fW)
aS, bS = sp.symbols('a_S b_S', positive=True)
def stress(Wf):
    return [sp.simplify(Wf/(aS**3*bS**4)), sp.simplify(-sp.diff(Wf, aS)/(3*aS**2*bS**4)), sp.simplify(-sp.diff(Wf, bS)/(4*aS**3*bS**3))]
r0 = sp.Symbol('r_0', positive=True)
st_t = stress(aS**3*fW/Lc**4); st_r = stress(aS**3*r0*aS**-4/Lc**4)
th = sp.Symbol('theta', positive=True)
c2 = sp.integrate(sp.cos(th)**4*sp.sin(th)**3, (th, 0, sp.pi/2))/sp.integrate(sp.cos(th)**2*sp.sin(th)**3, (th, 0, sp.pi/2))
R['smeared sources: tension eps = f_W/V_perp, pa = -eps, pb = 0; brane radiation pa = eps/3, pb = 0; bulk radiation <cos^2> = 3/7, <sin^2> = 4/7: pa = pb = eps/7'] = (
    sp.simplify(st_t[0] - fW/(bS**4*Lc**4)) == 0 and sp.simplify(st_t[1] + st_t[0]) == 0 and st_t[2] == 0
    and sp.simplify(st_r[1] - st_r[0]/3) == 0 and st_r[2] == 0 and c2 == sp.Rational(3, 7) and sp.Rational(1, 3)*c2 == sp.Rational(1, 7)
    and sp.Rational(1, 4)*(1 - c2) == sp.Rational(1, 7))
Tstar = (M8**6/MPl)**0.2
m8s, mps_, Ts = sp.symbols('M_8 M_Pl T', positive=True)
R['normalcy temperature: (T^7/M8^6)/(T^2/M_Pl) = 1 at T* = (M8^6/M_Pl)^(1/5) = 13.3 GeV; V_perp T^4 T^3/M_Pl^2 = T^7/M8^6'] = (
    sp.simplify(((Ts**7/m8s**6)/(Ts**2/mps_)).subs(Ts, (m8s**6/mps_)**sp.Rational(1, 5))) == 1 and abs(Tstar - 13.3) < 0.05
    and sp.simplify((mps_**2/m8s**6)*Ts**4*Ts**3/mps_**2 - Ts**7/m8s**6) == 0)
R['NUM T* = (M8^6/M_Pl)^(1/5) (GeV), f_W/rho_Lambda for f_W = M8^4, 1/(2 pi M8) (GeV^-1)'] = (Tstar, M8**4/rhoL, 1/(2*math.pi*M8))
R['the tension is vacuum energy: kappa f_W/V_perp = f_W/M_Pl^2 (V_perp = kappa M_Pl^2); f_W = M8^4 is 3e62 to 5e62 rho_Lambda'] = (
    sp.simplify(kap*fW/(kap*MP**2) - fW/MP**2) == 0 and 3e62 < M8**4/rhoL < 5e62)
yr = sp.Symbol('y', positive=True)
R['Green\'s function of the four transverse directions: (1/y^3) d/dy (y^3 d/dy y^-2) = 0, flux 2 pi^2 y^3 d/dy y^-2 = -4 pi^2; sqrt(kappa f_W) = 1/M8 for f_W = M8^4'] = (
    sp.simplify(sp.diff(yr**3*sp.diff(yr**-2, yr), yr)) == 0 and sp.simplify(2*sp.pi**2*yr**3*sp.diff(yr**-2, yr)) == -4*sp.pi**2
    and sp.simplify(sp.sqrt(m8s**-6*m8s**4) - 1/m8s) == 0)

if __name__ == '__main__':
    for k_, v_ in R.items():
        print('%s : %s' % (k_, v_))
    nb = [k_ for k_, v_ in R.items() if isinstance(v_, bool)]
    nf = sum(1 for k_ in nb if not R[k_])
    print('boolean checks: %d   failed: %d' % (len(nb), nf))
    print('VERDICT:', 'ALL CHECKS PASSED' if nf == 0 else 'SOMETHING FAILED')
    sys.exit(0 if nf == 0 else 1)
