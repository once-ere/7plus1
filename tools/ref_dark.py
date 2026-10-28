#!/usr/bin/env python3
"""Independent Python checks of section 24 of 7+1.nb: dirac16complex as dark matter and dark energy.

    python3 tools/ref_dark.py

Computes, without reading the notebook: the required dark-matter abundance from the Planck 2018 data; the collider scale M_D of the
model and the Pauli bound for 16 internal states; the stress tensor of a massive dirac16complex particle (sympy) and its cold-matter
limit; the radion displacement by dark matter; the couplings of the lightest Kaluza-Klein mode, the direct-detection and
self-interaction cross sections, the lifetimes through the Yukawa portal; the abundance from the primordial gas, from radion decay
and from the portals (the Boltzmann factor of production below T_RH, the zero mode, the Kaluza-Klein tower through the Yukawa portal,
the pairs of Kaluza-Klein levels through the Higgs portal, the vector portal); the equation of state of the Kohn-Sham matter of
section 8 (mpmath, numerical quadrature instead of the closed forms of the notebook); the Casimir energy of dirac16complex on the
auxiliary torus (Jacobi's r_4 and the Epstein sum).  The capture gas share and the radion decay rates are taken from
tools/ref_coupled.py, which is run here.  Every line 'key : value' is printed; boolean checks are summarized at the end (exit code 0
if all are True)."""
import math, os, re, subprocess, sys
from collections import OrderedDict
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.special import kv as besselk
from scipy.optimize import brentq
from scipy.integrate import quad, dblquad

HERE = os.path.dirname(os.path.abspath(__file__))
R = OrderedDict()

# ---------------------------------------------------------------------------------------------------- constants
MPl = 2.435e18; M8 = 1.0e4; Vperp = MPl**2/M8**6; lp = Vperp**0.25
hbar = 6.582119569e-25; hbarc = 1.973269804e-14; Mpc_km = 3.0856775814913673e19; Gyr_s = 3.15576e16
GeV_per_K = 8.617333262e-14; Tcmb = 2.7255*GeV_per_K
H0 = 67.4/Mpc_km*hbar
Om = 0.315
Og = (math.pi**2/15)*Tcmb**4/(3*H0**2*MPl**2); Or = Og*(1 + 7/8*(4/11)**(4/3)*3.044); OL = 1 - Om - Or
GF = 1.1663788e-5; vEW = (math.sqrt(2)*GF)**-0.5; mH = 125.20; alpha = 1/137.035999
GeV_g = 1.78266192e-24

# ---------------------------------------------------------------------------------------------------- 24.1 requirements
Och2, h = 0.120, 0.674
Oc = Och2/h**2
rhoc_h2 = 3*(100/Mpc_km*hbar)**2*MPl**2/hbarc**3
gs0 = 2 + 21/4*(4/11)
s0 = 2*math.pi**2/45*gs0*Tcmb**3/hbarc**3
Ym = Och2*rhoc_h2/s0
MD = (2*math.pi)**(2/3)*M8
pauli16 = 0.13*(2/16)**0.25
R['NUM Omega_c, rho_c/h^2 (GeV/cm^3), s0 (1/cm^3), m Y (GeV)'] = (Oc, rhoc_h2, s0, Ym)
R['NUM M_D (GeV) = (2 pi)^(2/3) M8'] = MD
R['NUM Pauli bound for g = 16 internal states, 0.13 keV (2/16)^(1/4) (keV)'] = pauli16
R['required abundance m Y = 4.37e-10 GeV (rho_c/h^2 = 1.053e-5 GeV/cm^3, s0 = 2891/cm^3, to 1e-3)'] = (
    abs(rhoc_h2/1.0534e-5 - 1) < 1e-3 and abs(s0/2891.2 - 1) < 1e-3 and abs(Ym/4.372e-10 - 1) < 1e-3)
m8, rR = sp.symbols('M_8 R', positive=True)
R['M_Pl^2 = R^4 M_D^6 with V = (2 pi R)^4 = M_Pl^2/M8^6 gives M_D = (2 pi)^(2/3) M8 = 34 TeV > 11.2 TeV (strongest collider bound)'] = (
    sp.simplify(((2*sp.pi)**sp.Rational(2, 3)*m8)**6*rR**4 - m8**6*(2*sp.pi*rR)**4) == 0 and MD > 11.2e3)

# ---------------------------------------------------------------------------------------------------- 24.2 cold matter
a, b, m, ka, kb = sp.symbols('a b m k_a k_b', positive=True)
W = sp.sqrt(m**2 + ka**2/a**2 + kb**2/b**2)
eps = W/(a**3*b**4); pa = -sp.diff(W, a)/(3*a**2*b**4); pb = -sp.diff(W, b)/(4*a**3*b**3)
R['massive particle: eps V = omega, pa V = k_a^2/(3 a^2 omega), pb V = k_b^2/(4 b^2 omega)'] = all(sp.simplify(e) == 0 for e in (
    eps*a**3*b**4 - W, pa*a**3*b**4 - ka**2/(3*a**2*W), pb*a**3*b**4 - kb**2/(4*b**2*W)))
R['lightest mode (k_b = 0): pa/eps = k_a^2/(3 (m^2 a^2 + k_a^2)), pb = 0, eps ~ a^-3 at fixed b for k_a = 0'] = (
    sp.simplify((pa/eps).subs(kb, 0) - ka**2/(3*(m**2*a**2 + ka**2))) == 0 and sp.simplify(pb.subs(kb, 0)) == 0
    and sp.simplify(a*sp.diff(sp.log(eps.subs({ka: 0, kb: 0})), a)) == -3)
dbeta = H0**2*Oc/(2*100.0)
R['NUM radion displacement by dark matter today H0^2 Omega_c/(2 m_r^2) (m_r = 10 GeV)'] = dbeta
R['dark matter leaves G4 unchanged: displacement < 1e-80'] = dbeta < 1e-80
om_, ol_, hh, tt, zz = sp.symbols('Omega_m Omega_L H t z', positive=True)
R['regime III: Omega_m a^-3/Omega_L = csch^2(3 H t/2) for a = (Omega_m/Omega_L)^(1/3) sinh^(2/3)(3 H t/2), -> 4 exp(-3 H t)'] = (
    sp.simplify((om_/ol_)/((om_/ol_)**sp.Rational(1, 3)*sp.sinh(sp.Rational(3, 2)*hh*tt)**sp.Rational(2, 3))**3 - 1/sp.sinh(sp.Rational(3, 2)*hh*tt)**2) == 0
    and sp.limit(4*sp.exp(2*zz)/(sp.exp(zz) - sp.exp(-zz))**2, zz, sp.oo) == 4)

# ---------------------------------------------------------------------------------------------------- the coupled results (run tools/ref_coupled.py)
rc = subprocess.run([sys.executable, '-B', os.path.join(HERE, 'ref_coupled.py')], capture_output=True, text=True).stdout
def rc_num(prefix):
    line = next(l for l in rc.splitlines() if l.startswith(prefix))
    s = re.sub(r'np\.float64\(([^()]*)\)', r'\1', line.split(' : ', 1)[1])
    return [float(x) for x in re.findall(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?', s)]
gas = rc_num('NUM after capture: gas fraction')
Gbr = rc_num('NUM radion of mass 10 GeV')[0]
Gbulk = rc_num('NUM radion into bulk dirac16complex pairs')[0]
late = rc_num('NUM coupled, gas-dominated late regime')
R['NUM from tools/ref_coupled.py: gas share at tau = 1e5, 3e5; Gamma_brane, Gamma_bulk (GeV); tau Hb at Ha = 1e-14'] = (gas[-2], gas[-1], Gbr, Gbulk, late[1])

# ---------------------------------------------------------------------------------------------------- 24.3 couplings, cross sections, stability
cP, cy, cV, cr, mPl_ = sp.symbols('c_P c_y c_V c_r M_Pl', positive=True)
vp = mPl_**2/m8**6
R['zero-mode couplings: lambda_P/V = c_P M8/M_Pl^2, y/sqrt(V) = c_y M8/M_Pl, G_P/V = c_V/M_Pl^2, g_r/V = c_r/M_Pl^2'] = all(sp.simplify(e) == 0 for e in (
    (cP/m8**5)/vp - cP*m8/mPl_**2, (cy/m8**2)/sp.sqrt(vp) - cy*m8/mPl_, (cV/m8**6)/vp - cV/mPl_**2, (cr/m8**6)/vp - cr/mPl_**2))
lam1, y1, G1 = M8/MPl**2, M8/MPl, 1/MPl**2
def sigSI(c, mx):
    mu = mx*0.9383/(mx + 0.9383)
    return (c*lam1)**2*0.30**2*0.9383**2*mu**2/(math.pi*mH**4)*hbarc**2
def sigm(c, mx): return (c*G1)**2*mx/(4*math.pi)*hbarc**2/GeV_g
def crmax(mx): return math.sqrt(4*math.pi*0.7*GeV_g/(hbarc**2*mx))/G1
def muP(c): return c*y1*vEW/math.sqrt(2)
def tauG(c, mx): return hbar/(9*alpha*GF**2*mx**5*4*(muP(c)/mx)**2/(1024*math.pi**4))
def tau3(c, mx): return hbar/(GF**2*mx**5*(muP(c)/mx)**2/(96*math.pi**3))
# age of the universe from the brane Friedmann equation (independent integration)
age_s = quad(lambda x: 1/(x*(H0/hbar)*math.sqrt(Or/x**4 + Om/x**3 + OL)), 0, 1, epsabs=0, epsrel=1e-12, limit=400)[0]
mmax = brentq(lambda mx: tau3(1, mx) - age_s, 1e-3, 1)
R['NUM lambda_eff, y_eff, G_r per unit coupling (GeV^-1, 1, GeV^-2)'] = (lam1, y1, G1)
R['NUM sigma_SI (cm^2) at c_P = 1, m = 40 GeV'] = sigSI(1, 40.0)
R['NUM sigma/m (cm^2/g) at c_r = 1, m = 1 GeV; c_r allowed by the Bullet Cluster at m = 1 GeV'] = (sigm(1, 1.0), crmax(1.0))
R['NUM mu_P (GeV) at c_y = 1; lifetimes chi -> nu gamma, chi -> 3 nu (s) at c_y = 1, m = 1 GeV and 10 keV'] = (muP(1), tauG(1, 1.0), tau3(1, 1.0), tauG(1, 1e-5), tau3(1, 1e-5))
R['NUM largest m (GeV) with lifetime chi -> 3 nu > age for c_y = 1'] = mmax
R['dark and collisionless: sigma_SI < 1e-50 x 2.2e-48 cm^2, sigma/m < 1e-70 cm^2/g'] = sigSI(1, 40.0) < 1e-50*2.2e-48 and sigm(1, 1.0) < 1e-70
R['stability: lifetime = age at m = 40 MeV (c_y = 1); lifetime ~ 1/(c_y^2 m^3); Gamma(nu gamma)/Gamma(3 nu) = 27 alpha/(8 pi)'] = (
    0.03 < mmax < 0.05 and abs(tau3(2, 0.01)/tau3(1, 0.01) - 0.25) < 1e-12 and abs(tau3(1, 0.02)/tau3(1, 0.01) - 0.125) < 1e-12
    and abs(tau3(1, 0.01)/tauG(1, 0.01) - 27*alpha/(8*math.pi)) < 1e-12)

# ---------------------------------------------------------------------------------------------------- 24.4 abundance: the solutions
f = gas[-1]
def over(trh): return (f/(1 - f))*trh/(4/3*Ym)
frac_bulk = Gbulk/(Gbulk + Gbr)
R['NUM gas share after capture f, overabundance O at T_RH = 1 GeV and 4 MeV'] = (f, over(1.0), over(0.004))
R['NUM fraction of radion decays into the bulk; brane share Gamma_brane/Gamma_bulk'] = (frac_bulk, Gbr/Gbulk)
R['primordial gas overabundance 1.5e7 T_RH/GeV (6e4 at 4 MeV)'] = abs(over(1.0)/1.51e7 - 1) < 0.01 and 5e4 < over(0.004) < 7e4
R['radion decays: more than 1 - 1e-12 into the bulk'] = 1 - frac_bulk < 1e-12

# ---------------------------------------------------------------------------------------------------- 24.4 abundance: the portals
Tstar = (M8**6/MPl)**0.2
mp.mp.dps = 30
def S(x):
    return float(2/(3*mp.pi)*mp.quad(lambda t: t**3*mp.besselk(1, t), [x, x + 20, mp.inf]))
t_ = sp.Symbol('t', positive=True)
exact_int = sp.integrate(t_**3*sp.besselk(1, t_), (t_, 0, sp.oo))
Sstar = S(mH/Tstar)
R['NUM T* (GeV), S(m_h/T*)'] = (Tstar, Sstar)
R['Boltzmann factor: int_0^oo t^3 K_1(t) dt = 3 pi/2 (sympy), S(0) = 1 (mpmath), S(m_h/T*) < 1e-2, 13 < T* < 14 GeV'] = (
    sp.simplify(exact_int - 3*sp.pi/2) == 0 and abs(S(0) - 1) < 1e-10 and Sstar < 1e-2 and 13 < Tstar < 14)
Cfi = 135/(8*math.pi**3*1.66)/106.75**1.5*1.22091e19
def Yfi(c, nch): return Cfi*nch*((c*y1)**2*mH/(16*math.pi))/mH**2
def Yhp(c): return Cfi*((c*lam1*vEW)**2*mH/(8*math.pi))/mH**2
K = tau3(1, 1.0)
def mstar(nch): return math.sqrt(Yfi(1, nch)*K/(Ym*age_s))
def cneed(mx, nch): return math.sqrt(Ym/(mx*Yfi(1, nch)))
def cmax(mx): return math.sqrt(K/(mx**3*age_s))
R['NUM zero mode: Yukawa-portal Y at c_y = 1 (N = 1, S = 1), Higgs-portal Y at c_P = 1 (S = 1)'] = (Yfi(1, 1), Yhp(1))
R['NUM m* (GeV) for N = 1 and N = 12'] = (mstar(1), mstar(12))
R['NUM c_y needed and allowed at m = 10 keV, 10 MeV'] = (cneed(1e-5, 1), cmax(1e-5), cneed(0.01, 1), cmax(0.01))
R['zero mode, Yukawa portal: abundance and stability compatible only below m* < 5.7 keV (N = 1, 12)'] = (
    mstar(1) < 5.7e-6 and mstar(12) < 5.7e-6 and abs(cneed(mstar(1), 1)/cmax(mstar(1)) - 1) < 1e-9)
# the Kaluza-Klein tower: continuum sums (sympy) and the lattice sums with r_4 (numpy)
u_ = sp.Symbol('u', positive=True)
i1 = sp.integrate(2*sp.pi**2*u_**3*(1 - u_**2)**2, (u_, 0, 1)); i2 = sp.integrate(2*sp.pi**2*u_**4*(1 - u_**2)**2, (u_, 0, 1))
TowN = Vperp*mH**4/(192*math.pi**2); TowE = Vperp*mH**5/(315*math.pi**2)
K4 = 1600000
sig = np.zeros(K4 + 1)
for dd in range(1, K4 + 1):
    sig[dd::dd] += dd
kk = np.arange(1, K4 + 1)
r4 = 8*sig[1:] - np.where(kk % 4 == 0, 32*sig[kk//4], 0.0)
R2 = 10**4; kq = np.arange(1, R2 + 1, dtype=float)
lat = 1 + np.sum(r4[:R2]*(1 - kq/R2)**2); latE = np.sum(r4[:R2]*np.sqrt(kq/R2)*(1 - kq/R2)**2)
R['NUM lattice sums for R = 100: Sum r_4(k)(1 - k/R^2)^2 / (pi^2 R^4/12) - 1, Sum r_4(k) sqrt(k)/R (1 - k/R^2)^2 / (16 pi^2 R^4/315) - 1'] = (
    float(lat/(math.pi**2*R2**2/12) - 1), float(latE/(16*math.pi**2*R2**2/315) - 1))
R['Kaluza-Klein tower: int 2 pi^2 u^3 (1-u^2)^2 = pi^2/12, int 2 pi^2 u^4 (1-u^2)^2 = 16 pi^2/315 (sympy); lattice sums agree to 1e-8 (R = 100)'] = (
    sp.simplify(i1 - sp.pi**2/12) == 0 and sp.simplify(i2 - 16*sp.pi**2/315) == 0
    and abs(lat/(math.pi**2*R2**2/12) - 1) < 1e-8 and abs(latE/(16*math.pi**2*R2**2/315) - 1) < 1e-8)
cytow = math.sqrt(Ym/(Yfi(1, 1)*TowE*Sstar))
R['NUM Kaluza-Klein tower below m_h: number, mass sum (GeV), Y_tower and mass per entropy (GeV) at c_y = 1 (N = 1, S = 1); c_y bound at T*'] = (
    TowN, TowE, Yfi(1, 1)*TowN, Yfi(1, 1)*TowE, cytow)
R['NUM tau(chi_n -> 3 nu) (s) at mu_n = 10 GeV, c_y = 1; largest m (GeV) with tau(chi -> 3 nu) > age for c_y = 1e-5'] = (tau3(1, 10.0), (K/(1e-10*age_s))**(1/3))
R['tower through the Yukawa portal: 7.7e17 states below m_h, mean mass (192/315) m_h, c_y bound < 2e-5 at T*'] = (
    abs(TowN/7.688e17 - 1) < 1e-3 and abs(TowE/(TowN*192/315*mH) - 1) < 1e-12 and cytow < 2e-5)
# the Higgs portal: pairs of Kaluza-Klein levels
fI = lambda u2, u1: u1**3*u2**3*(1 - (u1 + u2)**2)**1.5*math.sqrt(1 - (u1 - u2)**2)
Ihp = dblquad(fI, 0, 1, lambda u1: 0, lambda u1: 1 - u1, epsabs=1e-16, epsrel=1e-12)[0]
IhpE = dblquad(lambda u2, u1: (u1 + u2)*fI(u2, u1), 0, 1, lambda u1: 0, lambda u1: 1 - u1, epsabs=1e-16, epsrel=1e-12)[0]
Pairs = (Vperp*mH**4/(8*math.pi**2))**2
def Rhp(c, nch, x): return Yhp(c)*nch*Pairs*mH*IhpE*S(x)
def cpneed(x): return math.sqrt(Ym/Rhp(1, 4, x))
u1s, u2s = sp.symbols('u1 u2')
R['NUM Higgs portal pairs: I_hp, I_hp^E, number of pairs, mass per entropy (GeV) at c_P = 1 (N = 4, S = 1); c_P needed for S = 1 and at T*'] = (
    Ihp, IhpE, Pairs, Rhp(1, 4, 0), cpneed(0), cpneed(mH/Tstar))
R['NUM Higgs portal pairs at c_P = 1, T_RH = T*: fraction of the dark matter'] = Rhp(1, 4, mH/Tstar)/Ym
R['Higgs portal: zero mode Y < 1e-40; lambda(1, u1^2, u2^2) = (1 - (u1+u2)^2)(1 - (u1-u2)^2); pairs need 10 < c_P < 20 (S = 1), 100 < c_P < 200 (T*)'] = (
    Yhp(1) < 1e-40 and sp.expand((1 - u1s**2 - u2s**2)**2 - 4*u1s**2*u2s**2 - (1 - (u1s + u2s)**2)*(1 - (u1s - u2s)**2)) == 0
    and 10 < cpneed(0) < 20 and 100 < cpneed(mH/Tstar) < 200)
TV = (M8**12/MPl)**(1/11)
R['NUM vector portal: T_V (GeV), (T*/T_V)^11'] = (TV, (Tstar/TV)**11)
R['vector portal: 480 < T_V < 500 GeV, (T*/T_V)^11 < 1e-17'] = 480 < TV < 500 and (Tstar/TV)**11 < 1e-17
R['NUM T_RH (GeV) = 5, 8, 10, T*: S(m_h/T_RH), c_y bound of the tower (N = 1), c_P needed by the Higgs-portal pairs (N = 4)'] = [
    (trh, S(mH/trh), math.sqrt(Ym/(Yfi(1, 1)*TowE*S(mH/trh))), cpneed(mH/trh)) for trh in (5.0, 8.0, 10.0, Tstar)]

# ---------------------------------------------------------------------------------------------------- 24.5 equation of state of Kohn-Sham matter (section 8)
mp.mp.dps = 40
c7 = 1/(15*mp.pi**4)
def rho7(kF): return kF**7/(105*mp.pi**4)
def iS7(ms, kF): return mp.quad(lambda k: k**6*ms/mp.sqrt(k**2 + ms**2), [0, kF])
def iE7(ms, kF): return mp.quad(lambda k: k**6*mp.sqrt(k**2 + ms**2), [0, kF])
gS0, gV0 = -mp.mpf(10)**5, mp.mpf(10)**4
def gap7(kF): return mp.findroot(lambda ms: ms - 1 - gS0*c7*iS7(ms, kF), (mp.mpf('1e-8'), mp.mpf(1)), solver='anderson')
def eps7(kF):
    m1 = gap7(kF); return c7*iE7(m1, kF) - (m1 - 1)**2/(2*gS0) + gV0/2*rho7(kF)**2
def fermi7(kF): return gV0*rho7(kF) + mp.sqrt(kF**2 + gap7(kF)**2)
def pres(kF): return fermi7(kF)*rho7(kF) - eps7(kF)
kF0 = mp.findroot(lambda k: eps7(k)/rho7(k) - fermi7(k), (mp.mpf('0.7'), mp.mpf('0.8')), solver='anderson')
tab = [(mp.mpf(i)/100, pres(mp.mpf(i)/100)) for i in range(5, 86)]
wtab = [(k, pk/eps7(k)) for k, pk in tab]
wmin_k, wmin = min(wtab, key=lambda t: t[1])
R['NUM w_min of Kohn-Sham matter and its kF (units m = 1); saturation kF0'] = (float(wmin), float(wmin_k), float(kF0))
R['Kohn-Sham matter: p > 0 at low density, p < 0 just below kF0, p = 0 at kF0, a spinodal region, w >= -0.06'] = (
    tab[0][1] > 0 and pres(kF0 - mp.mpf(1)/100) < 0 and abs(pres(kF0)) < mp.mpf(10)**-20 and any(tab[i + 1][1] < tab[i][1] for i in range(len(tab) - 1))
    and -0.06 < float(wmin) < -0.04)

# ---------------------------------------------------------------------------------------------------- Casimir energy on the auxiliary torus
Z48 = 8*(1 - 4.0**-3)*float(mp.zeta(4)*mp.zeta(3))
Zshell = float(np.sum(r4[:20000]/np.arange(1, 20001, dtype=float)**4))
R['Epstein sum: Sum over Z^4\\{0} of |n|^-8 = 8 (1 - 4^-3) zeta(4) zeta(3) = shell sum to k = 20000 (tail < 3e-8)'] = abs(Z48 - Zshell) < 3e-8
C0 = math.gamma(4)/math.pi**4*Z48
rhoL = 3*H0**2*MPl**2*OL
def cas(x):
    ks = np.arange(1, 401, dtype=float)
    return 16*np.sum(r4[:400]*2*(x/(2*math.pi*np.sqrt(ks)))**4*besselk(4, x*np.sqrt(ks)))/lp**4
xL = brentq(lambda x: math.log(cas(x)/rhoL), 50, 150)
R['NUM rho_Lambda (GeV^4)'] = rhoL
R['NUM Casimir energy / rho_Lambda: dirac16complex (m = 0), graviton (20 dof)'] = (16*C0/lp**4/rhoL, -20*C0/lp**4/rhoL)
R['NUM m l_perp and m (GeV) at which the dirac16complex Casimir energy equals rho_Lambda'] = (xL, xL/lp)
R['NUM kF, p, w at kF = 0.1, 0.3, 0.5, 0.6, 0.67, 0.7, 0.75 (units m = 1)'] = [
    (float(k), float(pk), float(pk/eps7(k))) for k, pk in tab if int(round(float(k)*100)) in (10, 30, 50, 60, 67, 70, 75)]
R['NUM dirac16complex Casimir energy / rho_Lambda at m l_perp = 1, 10, 30, 60'] = [float(cas(x)/rhoL) for x in (1.0, 10.0, 30.0, 60.0)]
R['Casimir: massless > 1e34 rho_Lambda; massive sum -> massless (m l = 1e-3, to 1e-4); = rho_Lambda at m l = 92.2'] = (
    16*C0/lp**4/rhoL > 1e34 and abs(cas(1e-3)/(16*C0/lp**4) - 1) < 1e-4 and abs(xL - 92.2) < 0.05)

# ---------------------------------------------------------------------------------------------------- output
for k, v in R.items():
    print('%s : %s' % (k, v))
nb = [v for v in R.values() if isinstance(v, (bool, np.bool_))]
nf = sum(1 for v in nb if not v)
print('boolean checks: %d   failed: %d' % (len(nb), nf))
print('VERDICT:', 'ALL CHECKS PASSED' if nf == 0 else 'SOMETHING FAILED')
sys.exit(0 if nf == 0 else 1)
