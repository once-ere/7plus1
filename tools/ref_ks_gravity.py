#!/usr/bin/env python3
"""Independent Python (sympy / numpy / scipy / mpmath) check of sections 15-16 of 7+1.nb:
the generalized self-coupling of dirac16complex (all nine ranks of the Clifford algebra, the 9x9 Fierz matrix,
Hartree-Fock as Hartree with Fierz-transformed couplings) and the Kohn-Sham theory of the self-coupled field in
gravitational fields: the covariant nonlinear field equation, the chiral form in curved spacetime, the homogeneous
cosmological Kohn-Sham equations with their energy-momentum tensor, and the self-consistent coupling of the
Kohn-Sham matter to Einstein-Lovelock gravity in the primordial, modern and far-future epochs.

Run from the repository root:   python3 tools/ref_ks_gravity.py
Every line of output is  name : value ; the script exits with status 1 if any check is not True.
"""
import itertools, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from ref71 import G17, eta17

R = {}
G = {A: np.array(G17[A].tolist(), dtype=float) for A in range(9)}
eta = np.diag([1.0] + [-1.0]*7)
Gup = {A: eta[A, A]*G[A] for A in range(8)}
I16 = np.eye(16)
omega = G[8]; chi = 1j*omega

# ======================================================================= 1. the nine channels
def prod(mats):
    out = I16.copy()
    for M in mats: out = out @ M
    return out
Ks = [K for r in range(9) for K in itertools.combinations(range(8), r)]
sig = lambda r: (r*(r - 1)//2) % 2
Th_up = {K: (1j**sig(len(K)))*prod([Gup[A] for A in K]) for K in Ks}            # Theta^K
etaK = {K: float(np.prod([eta[A, A] for A in K])) if K else 1.0 for K in Ks}
Th_dn = {K: etaK[K]*Th_up[K] for K in Ks}                                      # Theta_K
R['256 products in 9 ranks: 1, 8, 28, 56, 70, 56, 28, 8, 1'] = [sum(1 for K in Ks if len(K) == r) for r in range(9)] == [1, 8, 28, 56, 70, 56, 28, 8, 1]
R['every bilinear psibar Theta^K psi is real: Gamma17[0] Theta^K Hermitian (Theta^K = i^((r(r-1)/2) mod 2) Gamma^K)'] = all(
    np.allclose((G[0] @ Th_up[K]).conj().T, G[0] @ Th_up[K]) for K in Ks)
R['Theta_K Theta^K = 1_16 and tr(Theta_K Theta^J) = 16 delta_KJ'] = all(np.allclose(Th_dn[K] @ Th_up[K], I16) for K in Ks) and np.allclose(
    np.array([[np.trace(Th_dn[K] @ Th_up[J]) for J in Ks] for K in Ks]), 16*np.eye(256))
comp = sum(np.einsum('cd,ab->abcd', Th_dn[K], Th_up[K]) for K in Ks)
R['Fierz completeness: sum_K (Theta_K)_cd (Theta^K)_ab = 16 delta_ad delta_bc'] = np.allclose(comp, 16*np.einsum('ad,bc->abcd', I16, I16))
full = tuple(range(8))
R['Hodge duality: Theta^K omega17 = +- Theta^(complement of K) (up to a phase)'] = all(
    np.allclose(abs(np.vdot(Th_up[K] @ omega, Th_up[tuple(sorted(set(full) - set(K)))]))/16, 1.0) for K in Ks)
rng = np.random.default_rng(71)
def rspin(n=16): return rng.normal(size=n) + 1j*rng.normal(size=n)
def Ichan(r, pb, ps):
    return sum(etaK[K]*(pb @ Th_up[K] @ ps)**2 for K in Ks if len(K) == r)
ps, pb = rspin(), rspin()
bil = lambda O: pb @ O @ ps
iS = bil(I16)**2; iV = sum(bil(Gup[A])*bil(G[A]) for A in range(8))
iP = bil(1j*chi)**2; iA = sum(bil(Gup[A] @ chi)*bil(G[A] @ chi) for A in range(8))
R['stage-2 channels: I_0 = (psibar psi)^2, I_1 = (psibar Gamma^A psi)(psibar Gamma_A psi), I_8 = -(psibar i chi psi)^2, I_7 = -(axial)'] = (
    np.isclose(Ichan(0, pb, ps), iS) and np.isclose(Ichan(1, pb, ps), iV) and np.isclose(Ichan(8, pb, ps), -iP) and np.isclose(Ichan(7, pb, ps), -iA))
# local Lorentz invariance of every channel
wgen = rng.normal(size=(8, 8)); wgen = wgen - wgen.T                              # omega_AB antisymmetric
Sgen = sum(wgen[A, B]*(Gup[A] @ Gup[B] - Gup[B] @ Gup[A])/4 for A in range(8) for B in range(8))/2
S = expm(0.3*Sgen); Si = expm(-0.3*Sgen)
R['local Lorentz invariance: I_r(S psi, psibar S^-1) = I_r(psi, psibar), r = 0..8 (random finite boost-rotation)'] = all(
    np.isclose(Ichan(r, pb @ Si, S @ ps), Ichan(r, pb, ps)) for r in range(9)) and np.allclose(Si, G[0] @ S.conj().T @ G[0])

# ======================================================================= 2. the 9x9 Fierz matrix
def fierz_closed(r, s):
    return sp.Rational(1, 16)*sum((-1)**(r*s - l)*sp.binomial(s, l)*sp.binomial(8 - s, r - l) for l in range(0, min(r, s) + 1))
F = sp.Matrix(9, 9, lambda r, s: fierz_closed(r, s))
Fnum = np.zeros((9, 9))
indep = True
for s_ in range(9):
    Js = [J for J in Ks if len(J) == s_]
    vals = []
    for J in Js:
        for r in range(9):
            pass
    for r in range(9):
        vr = [sum(np.trace(Th_dn[K] @ Th_dn[J] @ Th_up[K] @ Th_up[J]) for K in Ks if len(K) == r).real/256 for J in Js[:3] + Js[-2:]]
        Fnum[r, s_] = vr[0]
        indep = indep and np.allclose(vr, vr[0])
R['Fierz matrix F_rs = (1/256) sum_{K in r} tr(Theta_K Theta_J Theta^K Theta^J), independent of the representative J of rank s'] = indep
R['Fierz matrix closed form F_rs = (1/16) sum_l (-1)^(rs-l) C(s,l) C(8-s,r-l)'] = np.allclose(np.array(F.tolist(), dtype=float), Fnum)
R['F^2 = 1 (the Fierz transformation is an involution)'] = F*F == sp.eye(9)
ev = sorted(np.linalg.eigvals(np.array(F.tolist(), dtype=float)).real.round(10))
R['eigenvalues of F: number of +1, number of -1'] = (ev.count(1.0), ev.count(-1.0))
R['Fierz matrix (rows r, columns s) times 16'] = (F*16).tolist()
a1, b1, c1, d1 = rspin(), rspin(), rspin(), rspin()
lhs = [sum((a1 @ Th_dn[K] @ b1)*(c1 @ Th_up[K] @ d1) for K in Ks if len(K) == r) for r in range(9)]
rhs = [sum(float(F[r, s_])*sum((a1 @ Th_dn[J] @ d1)*(c1 @ Th_up[J] @ b1) for J in Ks if len(J) == s_) for s_ in range(9)) for r in range(9)]
R['Fierz identity for commuting spinors: sum_{K in r} (a Theta_K b)(c Theta^K d) = sum_s F_rs sum_{J in s} (a Theta_J d)(c Theta^J b)'] = np.allclose(lhs, rhs)
# Hartree-Fock = Hartree with G = (1 - F^T) g
Rr = rng.normal(size=(16, 16)) + 1j*rng.normal(size=(16, 16))
gr = rng.normal(size=9)
Dch = [sum(np.trace(Th_dn[J] @ Rr)*np.trace(Th_up[J] @ Rr) for J in Ks if len(J) == s_) for s_ in range(9)]
Ex = -sum(gr[r]/2*sum(np.trace(Th_dn[K] @ Rr @ Th_up[K] @ Rr) for K in Ks if len(K) == r) for r in range(9))
EH = sum(gr[s_]/2*Dch[s_] for s_ in range(9))
Gr = (np.eye(9) - np.array(F.tolist(), dtype=float).T) @ gr
R['Hartree-Fock energy E_H + E_x = sum_s (G_s/2) D_s with G = (1 - F^T) g (random density matrix, random couplings)'] = np.isclose(EH + Ex, sum(Gr[s_]/2*Dch[s_] for s_ in range(9)))
rS_, rho_ = sp.symbols('rho_S rho', real=True)
gS_, gV_, gP_, gA_ = sp.symbols('g_S g_V g_P g_A', real=True)
gvec = sp.Matrix([gS_, gV_, 0, 0, 0, 0, 0, -gA_, -gP_])
Gvec = (sp.eye(9) - F.T)*gvec
Ex2 = sp.expand(sum((Gvec[s_] - gvec[s_])/2*(rS_**2 if s_ == 0 else rho_**2 if s_ == 1 else 0) for s_ in range(9)))
Ex_stage2 = -(gS_/32)*(rS_**2 + rho_**2) + (gP_/32)*(rS_**2 - rho_**2) - (gV_/32)*(8*rS_**2 - 6*rho_**2) + (gA_/32)*(8*rS_**2 + 6*rho_**2)
R['uniform matter R = (rho_S + rho Gamma17[0])/16: only D_0 = rho_S^2, D_1 = rho^2; with g_0 = gS, g_1 = gV, g_7 = -gA, g_8 = -gP the exchange energy is that of stage 2 (section 8.7)'] = (
    sp.expand(Ex2 - Ex_stage2) == 0)
Rmat = lambda rs, rh: (rs*I16 + rh*G[0])/16
Dun = [sum(np.trace(Th_dn[J] @ Rmat(0.7, 1.3))*np.trace(Th_up[J] @ Rmat(0.7, 1.3)) for J in Ks if len(J) == s_) for s_ in range(9)]
R['uniform density matrix: D_0 = rho_S^2, D_1 = rho^2, D_2..D_8 = 0'] = np.allclose(Dun, [0.49, 1.69] + [0]*7)
gst2 = sp.Matrix([-100000, 10000, 0, 0, 0, 0, 0, 0, 0])
G2 = (sp.eye(9) - F.T)*gst2
R['stage-2 couplings gS = -100000, gV = 10000: G_0 = -98750, G_1 = 20000 (the HF couplings of section 8.7)'] = (G2[0], G2[1]) == (-98750, 20000)
R['effective couplings G_0, G_1 in terms of all nine g_r'] = (sp.expand(((sp.eye(9) - F.T)*sp.Matrix(sp.symbols('g0:9')))[0]),
                                                              sp.expand(((sp.eye(9) - F.T)*sp.Matrix(sp.symbols('g0:9')))[1]))
gex = sp.Matrix([-100000, 10000, 3000, 0, 0, 0, 0, 0, 0])
Gex = (sp.eye(9) - F.T)*gex
R['example g = (-100000, 10000, 3000, 0, ..., 0): G_0, G_1'] = (Gex[0], Gex[1])

# ======================================================================= 3. the nonlinear field equation (c-number derivative)
psb_s = sp.Matrix([[sp.Symbol('pb%d' % i) for i in range(16)]])
ps_n = [sp.Rational(int(rng.integers(-3, 4)), 1) + sp.I*int(rng.integers(-3, 4)) for _ in range(16)]
ps_v = sp.Matrix(ps_n)
gsym = [sp.Rational(int(rng.integers(-5, 6)), 1) for _ in range(9)]
ThU_s = {K: sp.Matrix(np.round(Th_up[K].real).astype(int)) + sp.I*sp.Matrix(np.round(Th_up[K].imag).astype(int)) for K in Ks}
Lint = -sum(gsym[len(K)]/2*etaK[K]*((psb_s*ThU_s[K]*ps_v)[0, 0])**2 for K in Ks)
fe_rhs = sum((gsym[len(K)]*etaK[K]*(psb_s*ThU_s[K]*ps_v)[0, 0]*(ThU_s[K]*ps_v) for K in Ks), sp.zeros(16, 1))
R['nonlinear field equation: dL_int/d psibar_a = -sum_r g_r sum_K (psibar Theta_K psi)(Theta^K psi)_a'] = all(
    sp.expand(sp.diff(Lint, psb_s[0, a_]) + fe_rhs[a_, 0]) == 0 for a_ in range(16))

# ======================================================================= 4. chirality in curved spacetime
Sig = [(Gup[A] @ Gup[B] - Gup[B] @ Gup[A])/4 for A in range(8) for B in range(A + 1, 8)]
R['the chirality chi17 commutes with every generator Sigma^AB, hence with every spin connection Omega_mu'] = all(np.allclose(chi @ s_, s_ @ chi) for s_ in Sig)
X = (sp.Symbol('x0', positive=True),) + sp.symbols('x1:8', real=True)
t = X[0]
a = sp.Function('a', real=True, positive=True)(t); b = sp.Function('b', real=True, positive=True)(t)
GS = {A: G17[A] for A in range(9)}
GSup = {A: eta17[A, A]*GS[A] for A in range(8)}
Om = {0: sp.zeros(16)}
for i in (1, 2, 3): Om[i] = sp.diff(a, t)/2*GS[i]*GS[0]
for k in (4, 5, 6, 7): Om[k] = sp.diff(b, t)/2*GS[k]*GS[0]
Gamu = {0: GSup[0]}
for i in (1, 2, 3): Gamu[i] = GSup[i]/a
for k in (4, 5, 6, 7): Gamu[k] = GSup[k]/b
from ref71 import tauO
I8 = sp.eye(8)
EL = sp.Matrix.vstack(I8, -sp.I*I8)/sp.sqrt(2); ER = sp.Matrix.vstack(-sp.I*I8, I8)/sp.sqrt(2)
Y = (X[0], X[1], X[5])
phL = sp.Matrix([sp.Function('l%d' % i)(*Y) for i in range(8)]); phR = sp.Matrix([sp.Function('r%d' % i)(*Y) for i in range(8)])
psi = EL*phL + ER*phR
mm = sp.Symbol('m', positive=True)
dirac = sp.I*sum((Gamu[mu]*(psi.diff(X[mu]) + Om[mu]*psi) for mu in range(8)), sp.zeros(16, 1)) - mm*psi
Hc = 3*sp.diff(a, t)/(2*a) + 2*sp.diff(b, t)/b
tdL = lambda f: sum((tauO[i]*f.diff(X[i]) for i in (1, 2, 3)), sp.zeros(8, 1))/a + sum((tauO[k]*f.diff(X[k]) for k in (4, 5, 6, 7)), sp.zeros(8, 1))/b
dag = lambda M: M.conjugate().T
R['curved chiral equations: E_L^dag Gamma0 (i Gamma^mu D_mu - m) psi = (i d0 + i(3Ha/2 + 2Hb) + tauO.d/(a,b)) phi_L - m phi_R'] = (
    (dag(EL)*GS[0]*dirac - (sp.I*phL.diff(t) + sp.I*Hc*phL + tdL(phL) - mm*phR)).applyfunc(sp.simplify) == sp.zeros(8, 1))
R['curved chiral equations: E_R^dag Gamma0 (i Gamma^mu D_mu - m) psi = (i d0 + i(3Ha/2 + 2Hb) - tauO.d/(a,b)) phi_R - m phi_L'] = (
    (dag(ER)*GS[0]*dirac - (sp.I*phR.diff(t) + sp.I*Hc*phR - tdL(phR) - mm*phL)).applyfunc(sp.simplify) == sp.zeros(8, 1))

# ======================================================================= 5. homogeneous Kohn-Sham matter in the cosmological spacetime
xk, wk = np.polynomial.legendre.leggauss(96)
xt, wt = np.polynomial.legendre.leggauss(64)
def moments(M, aa, bb, kF):
    """physical densities of a comoving Fermi ball |k| < kF (8 spin states) in the metric (a, b)"""
    k = kF*(xk + 1)/2; wkk = kF/2*wk
    th = np.pi/4*(xt + 1); wth = np.pi/4*wt
    K_, T_ = np.meshgrid(k, th, indexing='ij')
    W = np.outer(wkk, wth)*K_**6*np.cos(T_)**2*np.sin(T_)**3
    Pa2 = (K_*np.cos(T_)/aa)**2; Pb2 = (K_*np.sin(T_)/bb)**2
    om = np.sqrt(Pa2 + Pb2 + M*M)
    pref = 1/(2*np.pi**4*aa**3*bb**4)
    return dict(rho=pref*W.sum(), rhoS=pref*(W*M/om).sum(), ekin=pref*(W*om).sum(),
                pa=pref*(W*Pa2/om).sum()/3, pb=pref*(W*Pb2/om).sum()/4)
m0 = 1.0; gS1, gV1 = -100000.0, 10000.0
def gap(aa, bb, kF, GS_=gS1):
    f = lambda M: M - m0 - GS_*moments(M, aa, bb, kF)['rhoS']
    return brentq(f, 1e-14, m0, xtol=1e-15, rtol=1e-15)
def stress(aa, bb, kF, GS_=gS1, GV_=gV1):
    M = gap(aa, bb, kF, GS_); mo = moments(M, aa, bb, kF)
    eint = -(M - m0)**2/(2*GS_) + GV_/2*mo['rho']**2
    pint = (M - m0)**2/(2*GS_) + GV_/2*mo['rho']**2
    return dict(M=M, rho=mo['rho'], rhoS=mo['rhoS'], eps=mo['ekin'] + eint, pa=mo['pa'] + pint, pb=mo['pb'] + pint)
kF0 = 0.75663615070149289469
st0 = stress(1.0, 1.0, kF0)
R['flat limit a = b = 1, kF = kF0 (stage-2 saturation): m* = 0.3480973401..., E/A - m = -0.0283214279...'] = (
    abs(st0['M'] - 0.34809734014624326514) < 1e-11 and abs(st0['eps']/st0['rho'] - 1 + 0.028321427964738575854) < 1e-11)
R['flat limit: isotropic pressure pa = pb = rho eps_F - eps (Gibbs-Duhem, eps_F = V0 + sqrt(kF^2 + m*^2))'] = (
    abs(st0['pa'] - st0['pb']) < 1e-15 and abs(st0['pa'] - (st0['rho']*(gV1*st0['rho'] + math.sqrt(kF0**2 + st0['M']**2)) - st0['eps'])) < 1e-15)
R['NUM flat: m*, rho, eps, p'] = (st0['M'], st0['rho'], st0['eps'], st0['pa'])
# conservation along an arbitrary expansion history with fixed comoving occupation
def aofx(x): return 1 + 0.3*x + 0.05*x*x
def bofx(x): return 1/(1 + 0.2*x)
def epsx(x): return stress(aofx(x), bofx(x), kF0)
xc = 0.7; h = 1e-3
e_m2, e_m1, e_p1, e_p2 = (epsx(xc + d)['eps'] for d in (-2*h, -h, h, 2*h))
de = (e_m2 - 8*e_m1 + 8*e_p1 - e_p2)/(12*h)
sc = epsx(xc)
Hax = (0.3 + 0.1*xc)/aofx(xc); Hbx = -0.2/(1 + 0.2*xc)
cons = de + 3*Hax*(sc['eps'] + sc['pa']) + 4*Hbx*(sc['eps'] + sc['pb'])
R['energy-momentum conservation eps\' + 3Ha(eps + pa) + 4Hb(eps + pb) = 0 for the KS matter (arbitrary a, b; relative residual)'] = abs(cons/de) < 1e-8
R['NUM conservation residual / eps\''] = cons/de

# ======================================================================= 6. the three epochs (KS matter on given backgrounds)
pK = (1 + 2*math.sqrt(2))/7; qK = (2 - 3*math.sqrt(2))/14
prim = {x: stress(x**pK, x**qK, kF0) for x in (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0)}
R['primordial (Einstein-Kasner background, comoving kF = kF0 at x0 = 1): m*(x0) for x0 = 1e-4 ... 10'] = tuple(prim[x]['M'] for x in sorted(prim))
R['primordial: rho(x0) x0 = kF0^7/(105 pi^4) (density grows like 1/x0)'] = all(abs(prim[x]['rho']*x - kF0**7/(105*math.pi**4)) < 1e-12 for x in prim)
R['primordial: m* decreases toward 0 as x0 -> 0 (blue-shifted ordinary momenta, attractive scalar coupling)'] = (
    prim[1e-4]['M'] < prim[1e-3]['M'] < prim[1e-2]['M'] < prim[1e-1]['M'] < prim[1.0]['M'] and prim[1e-4]['M'] < 0.05)
R['primordial: anisotropic pressure pa/pb at x0 = 1e-3, 1, 10'] = tuple(prim[x]['pa']/prim[x]['pb'] for x in (1e-3, 1.0, 10.0))
mod = {aa: stress(aa, 1.0, kF0) for aa in (1.0, 2.0, 4.0, 10.0, 100.0)}
R['modern (static auxiliary space b = 1): m*(a) for a = 1, 2, 4, 10, 100'] = tuple(mod[aa]['M'] for aa in sorted(mod))
R['modern: eps a^3 -> constant and pa/eps -> 0 (dirac16complex behaves as pressureless matter in ordinary space)'] = (
    abs(mod[100.0]['eps']*100.0**3/(mod[10.0]['eps']*10.0**3) - 1) < 0.05 and mod[100.0]['pa']/mod[100.0]['eps'] < 1e-3)
R['modern: pb/eps at a = 100 (stress along the auxiliary space, balanced by the stabilizing stress)'] = mod[100.0]['pb']/mod[100.0]['eps']
Hds = 0.1
fut = {x: stress(math.exp(Hds*x), 1.0, kF0) for x in (0.0, 10.0, 20.0, 40.0, 80.0)}
R['far future (de Sitter a = e^{0.1 x0}, b = 1): m*(x0) for x0 = 0, 10, 20, 40, 80 -> m = 1'] = tuple(fut[x]['M'] for x in sorted(fut))
R['far future: rho = rho0 e^{-3 H x0} and 1 - m* -> 0 like rho'] = (abs(fut[40.0]['rho']/fut[0.0]['rho'] - math.exp(-12)) < 1e-12
    and abs((1 - fut[80.0]['M'])/(1 - fut[40.0]['M'])/math.exp(-12) - 1) < 0.05)

# ======================================================================= 7. self-consistent Einstein-Lovelock + Kohn-Sham evolution
Ha, Hb, dHa, dHb = sp.symbols('H_a H_b dH_a dH_b', real=True)
# closed forms of E^(k) (00, 11, 44) for the cosmological metric (from the pair-set Lovelock tensors, see tools/ref_gravity.py)
Ecl = {
 (1, 0): 3*Ha**2 + 12*Ha*Hb + 6*Hb**2,
 (1, 1): 2*dHa + 3*Ha**2 + 4*dHb + 10*Hb**2 + 8*Ha*Hb,
 (1, 4): 3*dHa + 6*Ha**2 + 3*dHb + 6*Hb**2 + 9*Ha*Hb,
 (2, 0): -12*Hb*(4*Ha**3 + 18*Ha**2*Hb + 12*Ha*Hb**2 + Hb**3),
 (3, 0): 360*Ha**2*Hb**3*(4*Ha + 3*Hb)}
from lovelock71 import cosmological_closed_forms
CF = cosmological_closed_forms()
R['the 00 closed forms of E^(1), E^(2), E^(3) (pair-set algorithm, tools/lovelock71.py) agree with section 12.3'] = all(
    sp.expand(CF[(kk, 0)] - Ecl[(kk, 0)]) == 0 for kk in (1, 2, 3)) and all(sp.expand(CF[(1, c)] - Ecl[(1, c)]) == 0 for c in (1, 4))
for kk in (2, 3):
    for c in (1, 4):
        Ecl[(kk, c)] = CF[(kk, c)]
al2, al3, Lam, kap = sp.symbols('alpha2 alpha3 Lambda kappa', real=True)
tot = {c: Ecl[(1, c)] + al2*Ecl[(2, c)] + al3*Ecl[(3, c)] - Lam for c in (0, 1, 4)}
R['Einstein-Lovelock equations are linear in (Ha\', Hb\') (quasi-linear second-order system)'] = all(
    sp.degree(sp.expand(tot[c]), dHa) <= 1 and sp.degree(sp.expand(tot[c]), dHb) <= 1 for c in (1, 4)) and all(
    not tot[0].has(dHa) and not tot[0].has(dHb) for _ in (0,))
lin = {c: sp.lambdify((Ha, Hb, al2, al3, Lam), [sp.expand(tot[c]).coeff(dHa), sp.expand(tot[c]).coeff(dHb),
                                                 sp.expand(tot[c]).subs({dHa: 0, dHb: 0})]) for c in (1, 4)}
con = sp.lambdify((Ha, Hb, al2, al3, Lam), tot[0])
def evolve_EL(alpha2, alpha3, Lmb, kappa, x_i, x_f, Ha_i):
    s0 = stress(1.0, 1.0, kF0)
    cfun = lambda hb: con(Ha_i, hb, alpha2, alpha3, Lmb) - kappa*s0['eps']
    nodes = np.linspace(-1.0, 1.0, 5)                       # the constraint is a polynomial of degree <= 4 in Hb
    coef = np.polyfit(nodes, [cfun(z) for z in nodes], 4)
    rts = [z.real for z in np.roots(coef) if abs(z.imag) < 1e-9]
    z0 = min(rts, key=abs); Hb_i = brentq(cfun, z0 - 1e-3*(1 + abs(z0)), z0 + 1e-3*(1 + abs(z0)), xtol=1e-16)
    def rhs(x, y):
        aa, bb, ha, hb = y
        st = stress(aa, bb, kF0)
        c1 = lin[1](ha, hb, alpha2, alpha3, Lmb); c4 = lin[4](ha, hb, alpha2, alpha3, Lmb)
        Mx = np.array([[c1[0], c1[1]], [c4[0], c4[1]]], dtype=float)
        rv = np.array([-kappa*st['pa'] - c1[2], -kappa*st['pb'] - c4[2]], dtype=float)
        dh = np.linalg.solve(Mx, rv)
        return [ha*aa, hb*bb, dh[0], dh[1]]
    sol = solve_ivp(rhs, (x_i, x_f), [1.0, 1.0, Ha_i, Hb_i], method='DOP853', rtol=1e-10, atol=1e-13, dense_output=True)
    res = []
    for x in np.linspace(x_i, x_f, 6):
        aa, bb, ha, hb = sol.sol(x)
        st = stress(aa, bb, kF0)
        res.append((con(ha, hb, alpha2, alpha3, Lmb) - kappa*st['eps'])/(kappa*st['eps']))
    return Hb_i, sol.y[:, -1], max(abs(r_) for r_ in res)
kap_ = 1000.0
HbE, endE, cE = evolve_EL(0.0, 0.0, 0.0, kap_, 0.0, 20.0, 0.05)
R['Einstein + KS matter (kappa = 1000, Lambda = 0, Ha(0) = 0.05): Hb(0) from the constraint, (a, b, Ha, Hb) at x0 = 20, max relative constraint violation'] = (HbE, tuple(endE), cE)
HbL, endL, cL = evolve_EL(20.0, 500.0, 0.0, kap_, 0.0, 20.0, 0.05)
R['Einstein-Lovelock (alpha2 = 20, alpha3 = 500) + KS matter: Hb(0), (a, b, Ha, Hb) at x0 = 20, max relative constraint violation'] = (HbL, tuple(endL), cL)
R['the constraint Sum alpha_k E^(k)0_0 - Lambda = kappa eps is preserved by the evolution (relative violation < 1e-7) in both cases'] = cE < 1e-7 and cL < 1e-7

if __name__ == '__main__':
    bad = [k for k, v in R.items() if isinstance(v, (bool, np.bool_)) and not v]
    nbool = sum(isinstance(v, (bool, np.bool_)) for v in R.values())
    for k, v in R.items():
        print(f'{k} : {v}')
    print(f'boolean checks: {nbool}   failed: {len(bad)}')
    print('VERDICT:', 'ALL CHECKS PASSED' if not bad else 'FAILED: ' + '; '.join(bad))
    sys.exit(1 if bad else 0)
