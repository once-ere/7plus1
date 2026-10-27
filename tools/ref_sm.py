#!/usr/bin/env python3
"""Independent Python (sympy / numpy / scipy) check of sections 17-19 of 7+1.nb:
the Standard Model on the world brane of the 7+1 spacetime (all fields, quantum numbers, anomaly cancellation,
electroweak symmetry breaking, quark and lepton masses, CKM and PMNS mixing), dirac16complex at the brane (the
intertwiner that splits its 16 components into four brane Dirac spinors, the Kaluza-Klein masses on the auxiliary
torus), the portals of dirac16complex to the Standard Model, and neutrino flavour conversion (vacuum, matter, the Sun,
supernova resonances, the cosmic neutrino background, mixing with dirac16complex).

Measured inputs (all with their source in the documentation): PDG 2025 summary tables (m_W, m_Z, m_H, quark masses),
PDG 2024 CKM review (eq. 12.28), NuFit-6.0 (arXiv:2410.05380, normal ordering, without SK atmospheric data),
G_F = 1.1663788e-5 GeV^-2, Planck 2018 cosmological parameters.

Run from the repository root:   python3 tools/ref_sm.py
"""
import itertools, math, os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad
from ref71 import G17, eta17

R = {}
# ======================================================================= 1. fields, quantum numbers, anomalies
# (name, SU(3) dim, SU(2) dim, Y, chirality L/R, generations)
F = [('Q', 3, 2, Fr(1, 6), 'L'), ('u_R', 3, 1, Fr(2, 3), 'R'), ('d_R', 3, 1, Fr(-1, 3), 'R'),
     ('L', 1, 2, Fr(-1, 2), 'L'), ('e_R', 1, 1, Fr(-1), 'R'), ('nu_R', 1, 1, Fr(0), 'R')]
sgn = lambda c: 1 if c == 'L' else -1
an = {}
an['SU(3)^2 U(1)'] = sum(sgn(c)*n2*Y for (nm, n3, n2, Y, c) in F if n3 == 3)
an['SU(2)^2 U(1)'] = sum(sgn(c)*n3*Y for (nm, n3, n2, Y, c) in F if n2 == 2)
an['U(1)^3'] = sum(sgn(c)*n3*n2*Y**3 for (nm, n3, n2, Y, c) in F)
an['gravity^2 U(1)'] = sum(sgn(c)*n3*n2*Y for (nm, n3, n2, Y, c) in F)
an['SU(3)^3'] = sum(sgn(c)*n2 for (nm, n3, n2, Y, c) in F if n3 == 3)
R['anomaly cancellation in each generation (with nu_R): SU(3)^2U(1), SU(2)^2U(1), U(1)^3, grav^2U(1), SU(3)^3 all 0'] = all(v == 0 for v in an.values())
R['Witten SU(2) anomaly: number of doublets per generation (3 colours of Q + 1 L) = 4, even'] = sum(n3 for (nm, n3, n2, Y, c) in F if n2 == 2) % 2 == 0
charges = {nm: [Y + Fr(1, 2), Y - Fr(1, 2)] if n2 == 2 else [Y] for (nm, n3, n2, Y, c) in F}
R['electric charges Q = T3 + Y: Q -> (2/3, -1/3), u_R 2/3, d_R -1/3, L -> (0, -1), e_R -1, nu_R 0'] = charges == {
    'Q': [Fr(2, 3), Fr(-1, 3)], 'u_R': [Fr(2, 3)], 'd_R': [Fr(-1, 3)], 'L': [Fr(0), Fr(-1)], 'e_R': [Fr(-1)], 'nu_R': [Fr(0)]}
# relativistic degrees of freedom (all SM fields massless and in equilibrium)
gb = 8*2 + 3*2 + 1*2 + 4                                 # gluons, W^1..3, B (2 polarizations each), Higgs doublet (4 real)
gf = 3*sum(n3*n2*2*(1 if True else 1) for (nm, n3, n2, Y, c) in F if nm != 'nu_R')   # Weyl fields: 2 (particle+antiparticle) per component
R['effective number of relativistic degrees of freedom g* = 28 + (7/8) 90 = 106.75 (106.75 + 5.25 if nu_R were thermal)'] = (
    gb == 28 and gf == 90 and gb + Fr(7, 8)*gf == Fr(427, 4) and gb + Fr(7, 8)*(gf + 6) == 112)

# ======================================================================= 2. electroweak symmetry breaking
g, gp, v = sp.symbols('g gp v', positive=True)
W = sp.symbols('W1:4', real=True); Bf = sp.Symbol('B', real=True)
tau = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
Phi0 = sp.Matrix([0, v/sp.sqrt(2)])
DPhi = -sp.I*(g*sum((W[i]*tau[i]/2 for i in range(3)), sp.zeros(2)) + gp*sp.Rational(1, 2)*Bf*sp.eye(2))*Phi0
mass_term = sp.expand((DPhi.H*DPhi)[0, 0])                 # = (1/2) V^T M2 V
fields = list(W) + [Bf]
M2 = sp.Matrix(4, 4, lambda i, j: sp.diff(mass_term, fields[i], fields[j]))
ev = M2.eigenvals()
evl = sorted([sp.expand(k_) for k_, mult in ev.items() for _ in range(mult)], key=sp.default_sort_key)
R['gauge-boson mass matrix from |D<Phi>|^2: eigenvalues g^2v^2/4 (twice, W+-), (g^2+g\'^2)v^2/4 (Z), 0 (photon)'] = (
    evl == sorted([sp.expand(x_) for x_ in (g**2*v**2/4, g**2*v**2/4, (g**2 + gp**2)*v**2/4, sp.Integer(0))], key=sp.default_sort_key))
GF = 1.1663788e-5; mW = 80.3692; mZ = 91.1880; mH = 125.20
vv = (math.sqrt(2)*GF)**-0.5
gg = 2*mW/vv; cw = mW/mZ; sw2 = 1 - cw*cw; ggp = gg*math.sqrt(sw2)/cw; ee = gg*math.sqrt(sw2)
lam_h = mH**2/(2*vv**2)
R['NUM v = (sqrt2 G_F)^(-1/2) [GeV], g, g\', sin^2 theta_W (on shell), e, 1/alpha(tree), lambda = m_H^2/(2v^2)'] = (vv, gg, ggp, sw2, ee, 4*math.pi/ee**2, lam_h)
R['tree-level relations: m_W = g v/2, m_Z = sqrt(g^2+g\'^2) v/2, rho = m_W^2/(m_Z^2 cos^2) = 1, m_H^2 = 2 lambda v^2'] = (
    abs(gg*vv/2 - mW) < 1e-12 and abs(math.sqrt(gg**2 + ggp**2)*vv/2 - mZ) < 1e-10 and abs(mW**2/(mZ**2*cw**2) - 1) < 1e-14 and abs(2*lam_h*vv**2 - mH**2) < 1e-9)

# ======================================================================= 3. quark masses, Yukawas, CKM
s12, s13, s23, dq = 0.22501, 0.003732, 0.04183, 1.147
def mix(s12_, s13_, s23_, d_):
    c12, c13, c23 = (math.sqrt(1 - s*s) for s in (s12_, s13_, s23_))
    ed = np.exp(1j*d_)
    return np.array([[c12*c13, s12_*c13, s13_/ed],
                     [-s12_*c23 - c12*s23_*s13_*ed, c12*c23 - s12_*s23_*s13_*ed, s23_*c13],
                     [s12_*s23_ - c12*c23*s13_*ed, -c12*s23_ - s12_*c23*s13_*ed, c23*c13]])
Vckm = mix(s12, s13, s23, dq)
R['CKM unitarity V V^dag = 1'] = np.allclose(Vckm @ Vckm.conj().T, np.eye(3), atol=1e-15)
Jq = (Vckm[0, 1]*Vckm[1, 2]*np.conj(Vckm[0, 2])*np.conj(Vckm[1, 1])).imag
R['NUM CKM |V| (rows ud us ub / cd cs cb / td ts tb)'] = np.abs(Vckm).round(6).tolist()
R['NUM Jarlskog J = Im(V_us V_cb V_ub* V_cs*)'] = Jq
R['J = c12 c23 c13^2 s12 s23 s13 sin(delta) = 3.12e-5 (PDG 2024: 3.12 +0.13 -0.12 e-5)'] = abs(Jq - 3.12e-5) < 0.01e-5
rb = -Vckm[0, 0]*np.conj(Vckm[0, 2])/(Vckm[1, 0]*np.conj(Vckm[1, 2]))
R['NUM Wolfenstein lambda = s12, A = s23/s12^2, rhobar + i etabar = -V_ud V_ub*/(V_cd V_cb*)'] = (s12, s23/s12**2, rb.real, rb.imag)
R['Wolfenstein parameters agree with the PDG fit (A 0.826, rhobar 0.1591, etabar 0.3523) within their quoted uncertainties'] = (
    abs(s23/s12**2 - 0.826) < 0.016 and abs(rb.real - 0.1591) < 0.0094 and abs(rb.imag - 0.3523) < 0.0073)
mq_u = np.array([2.16e-3, 1.2730, 172.56]); mq_d = np.array([4.70e-3, 93.5e-3, 4.183])
Yu = math.sqrt(2)*np.diag(mq_u)/vv; Yd = math.sqrt(2)*Vckm @ np.diag(mq_d)/vv
Uu, su, Vuh = np.linalg.svd(Yu); Ud, sd, Vdh = np.linalg.svd(Yd)
R['Yukawa matrices Y_u = sqrt2 diag(m_u,m_c,m_t)/v, Y_d = sqrt2 V diag(m_d,m_s,m_b)/v: singular values x v/sqrt2 = masses, |U_uL^dag U_dL| = |V|'] = (
    np.allclose(sorted(su*vv/math.sqrt(2)), sorted(mq_u), rtol=1e-12) and np.allclose(sorted(sd*vv/math.sqrt(2)), sorted(mq_d), rtol=1e-12)
    and np.allclose(np.abs((Uu[:, ::-1].conj().T @ Ud[:, ::-1])), np.abs(Vckm), atol=1e-9))
R['NUM top Yukawa y_t = sqrt2 m_t/v'] = math.sqrt(2)*172.56/vv

# ======================================================================= 4. leptons, PMNS, Dirac neutrino masses
sn12, sn23, sn13 = math.sqrt(0.307), math.sqrt(0.561), math.sqrt(0.02195); dl = math.radians(177.0)
dm21, dm31 = 7.49e-5, 2.534e-3                                    # eV^2, normal ordering
Upmns = mix(sn12, sn13, sn23, dl)
R['PMNS unitarity'] = np.allclose(Upmns @ Upmns.conj().T, np.eye(3), atol=1e-15)
Jl = (Upmns[0, 1]*Upmns[1, 2]*np.conj(Upmns[0, 2])*np.conj(Upmns[1, 1])).imag
R['NUM PMNS |U| and J_PMNS'] = (np.abs(Upmns).round(5).tolist(), Jl)
m1 = 0.0; m2 = math.sqrt(dm21); m3 = math.sqrt(dm31)
R['NUM Dirac neutrino masses (lightest m1 = 0, eV) and their sum'] = (m1, m2, m3, m1 + m2 + m3)
ml = np.array([0.51099895e-3, 105.6583755e-3, 1776.93e-3])
Ye = math.sqrt(2)*np.diag(ml)/vv; Ynu = math.sqrt(2)*Upmns @ np.diag([m1, m2, m3])*1e-9/vv
Un_, sn_, _ = np.linalg.svd(Ynu)
ordn = np.argsort(sn_)
R['neutrino Yukawa matrix Y_nu = sqrt2 U diag(m1,m2,m3)/v: masses and |U| recovered (singular value decomposition)'] = (
    np.allclose(sn_[ordn]*vv/math.sqrt(2)*1e9, [m1, m2, m3], atol=1e-15) and np.allclose(np.abs(Un_[:, ordn]), np.abs(Upmns), atol=1e-9))
R['NUM largest neutrino Yukawa sqrt2 m3/v'] = math.sqrt(2)*m3*1e-9/vv

# ======================================================================= 5. neutrino oscillations
HBARC = 1.973269804e-16                                           # GeV m
def phase_factor(L_km, E_GeV):                                   # Delta m^2 [eV^2] L / (2E)  -> radians per eV^2
    return 1e-18*L_km*1e3/(2*E_GeV*HBARC)
def P_vac(U, L, E, anti=False):
    Uu_ = U.conj() if anti else U
    ph = np.exp(-1j*np.array([0.0, dm21, dm31])*phase_factor(L, E))
    A = Uu_ @ np.diag(ph) @ Uu_.conj().T                         # A[b, a] = amplitude a -> b
    return np.abs(A.T)**2                                         # P[a, b]
Pd = P_vac(Upmns, 1300.0, 2.5)
R['vacuum: probabilities P(a -> b) form a doubly stochastic matrix (rows and columns sum to 1)'] = np.allclose(Pd.sum(1), 1) and np.allclose(Pd.sum(0), 1)
def P_formula(U, a, b, L, E):
    dms = {(1, 0): dm21, (2, 0): dm31, (2, 1): dm31 - dm21}
    s = 1.0 if a == b else 0.0
    for (j, i), d in dms.items():
        X = U[a, i]*np.conj(U[b, i])*np.conj(U[a, j])*U[b, j]
        D_ = d*phase_factor(L, E)/2
        s += -4*X.real*math.sin(D_)**2 + 2*X.imag*math.sin(2*D_)
    return s
R['vacuum: amplitude result = delta_ab - 4 sum Re(X) sin^2 Delta + 2 sum Im(X) sin 2Delta, X = U_ai U_bi* U_aj* U_bj'] = all(
    abs(P_formula(Upmns, a, b, 1300.0, 2.5) - Pd[a, b]) < 1e-13 for a in range(3) for b in range(3))
Pbar = P_vac(Upmns, 1300.0, 2.5, anti=True)
R['NUM DUNE vacuum (L = 1300 km, E = 2.5 GeV): P(mu->e), P(mubar->ebar), P(mu->mu)'] = (Pd[1, 0], Pbar[1, 0], Pd[1, 1])
R['CPT: P(a->b) = P(bbar->abar); CP asymmetry P(mu->e) - P(mubar->ebar) proportional to J (vanishes for delta = 0, pi)'] = (
    np.allclose(Pd, Pbar.T) and abs(P_vac(mix(sn12, sn13, sn23, math.pi), 1300.0, 2.5)[1, 0] - P_vac(mix(sn12, sn13, sn23, math.pi), 1300.0, 2.5, True)[1, 0]) < 1e-14)
# matter of constant density: V_CC = sqrt2 G_F n_e
NA = 6.02214076e23; HBARC_EVCM = 1.973269804e-5                   # eV cm
def Vcc(rho_gcc, Ye):                                             # eV
    ne = rho_gcc*Ye*NA                                            # cm^-3
    return math.sqrt(2)*GF*1e-18*ne*HBARC_EVCM**3                 # G_F in eV^-2
R['NUM V_CC = sqrt2 G_F n_e for rho Y_e = 1 g/cm^3 [eV] (textbook 7.63e-14 eV)'] = Vcc(1.0, 1.0)
R['V_CC(rho Y_e = 1 g/cm^3) = 7.63e-14 eV'] = abs(Vcc(1.0, 1.0) - 7.63e-14) < 0.01e-14
def P_matter(U, L, E, rho, Ye, anti=False, extra=None):
    Uu_ = U.conj() if anti else U
    H = Uu_ @ np.diag([0.0, dm21, dm31]) @ Uu_.conj().T/(2*E*1e9)   # eV
    Vm = np.diag([Vcc(rho, Ye), 0, 0])*(-1 if anti else 1)
    if extra is not None: Vm = Vm + np.diag(extra)*(-1 if anti else 1)
    Ht = H + Vm
    w, vecs = np.linalg.eigh(Ht)
    Lev = L*1e3/(HBARC*1e9)                                      # km -> eV^-1
    A = vecs @ np.diag(np.exp(-1j*w*Lev)) @ vecs.conj().T
    return np.abs(A.T)**2
Pm = P_matter(Upmns, 1300.0, 2.5, 2.848, 0.5)
R['matter (constant density): unitary evolution, probabilities doubly stochastic'] = np.allclose(Pm.sum(1), 1) and np.allclose(Pm.sum(0), 1)
R['NUM DUNE in matter (rho = 2.848 g/cm^3, Y_e = 0.5): P(mu->e), P(mubar->ebar)'] = (Pm[1, 0], P_matter(Upmns, 1300.0, 2.5, 2.848, 0.5, anti=True)[1, 0])
R['matter: a flavour-universal potential has no effect (only differences of potentials matter)'] = np.allclose(
    P_matter(Upmns, 1300.0, 2.5, 2.848, 0.5, extra=[1e-13]*3), Pm, atol=1e-12)
# two-flavour resonance
th = math.asin(sn12); dm = dm21
def th_matter(E_MeV, V):
    A_ = 2*E_MeV*1e6*V/dm
    return 0.5*math.atan2(math.sin(2*th), math.cos(2*th) - A_)
def eff_mix_num(E_MeV, V):
    H2 = dm/(4*E_MeV*1e6)*np.array([[-math.cos(2*th), math.sin(2*th)], [math.sin(2*th), math.cos(2*th)]]) + np.diag([V, 0])
    w, vecs = np.linalg.eigh(H2)
    return abs(vecs[0, 1])**2                                      # |<nu_e|nu_2m>|^2 = sin^2 theta_m
Vres_test = 5e-12
Eres = dm*math.cos(2*th)/(2*Vres_test)/1e6
R['two-flavour MSW: sin^2 theta_m from diagonalization = formula; maximal mixing (sin^2 2theta_m = 1) at 2 E V = dm^2 cos 2theta'] = (
    all(abs(eff_mix_num(Ex, Vres_test) - math.sin(th_matter(Ex, Vres_test))**2) < 1e-12 for Ex in (0.5, 2.0, 10.0))
    and abs(math.sin(2*th_matter(Eres, Vres_test))**2 - 1) < 1e-12)
# the Sun: n_e(r) = 245 N_A exp(-10.54 r/R) cm^-3 (Bahcall approximation)
Rsun = 6.957e8                                                    # m
def V_sun(r): return math.sqrt(2)*GF*1e-18*245*NA*math.exp(-10.54*r/Rsun)*HBARC_EVCM**3
def Pee_sun_numeric(E_MeV):
    c13sq = 1 - sn13**2
    k = dm/(4*E_MeV*1e6)
    def rhs(r, y):
        V = V_sun(r)*c13sq
        H2 = np.array([[-k*math.cos(2*th) + V/2, k*math.sin(2*th)], [k*math.sin(2*th), k*math.cos(2*th) - V/2]])
        return -1j*(H2 @ y)*(1/(HBARC*1e9))                      # dr in m, H in eV: phase = H r/(hbar c)
    w0, v0 = np.linalg.eigh(np.array([[-k*math.cos(2*th) + V_sun(0)*c13sq/2, k*math.sin(2*th)], [k*math.sin(2*th), k*math.cos(2*th) - V_sun(0)*c13sq/2]]))
    y0 = np.array([1.0, 0.0], complex)                           # nu_e produced at the centre
    sol = solve_ivp(rhs, (0.0, Rsun), y0, method='DOP853', rtol=1e-9, atol=1e-11)
    y1 = sol.y[:, -1]
    vac = np.array([[math.cos(th), math.sin(th)], [-math.sin(th), math.cos(th)]])    # columns: nu_1, nu_2 in flavour basis (e, x)
    p1 = abs(np.vdot(vac[:, 0], y1))**2; p2 = abs(np.vdot(vac[:, 1], y1))**2
    P2 = p1*math.cos(th)**2 + p2*math.sin(th)**2                 # incoherent propagation to the Earth
    return c13sq**2*P2 + sn13**4
def Pee_sun_adiabatic(E_MeV):
    c13sq = 1 - sn13**2
    thm = 0.5*math.atan2(math.sin(2*th), math.cos(2*th) - 2*E_MeV*1e6*V_sun(0)*c13sq/dm)
    return c13sq**2*0.5*(1 + math.cos(2*th)*math.cos(2*thm)) + sn13**4
sunE = {}
for Ex in (0.3, 1.0, 5.0, 10.0):
    sunE[Ex] = (Pee_sun_numeric(Ex), Pee_sun_adiabatic(Ex))
R['NUM Sun: P_ee (numerical through the solar density profile, adiabatic formula) for E = 0.3, 1, 5, 10 MeV'] = sunE
R['Sun: numerical P_ee = adiabatic MSW formula to 2e-3 at E = 0.3, 1, 5, 10 MeV'] = all(abs(sunE[Ex][0] - sunE[Ex][1]) < 2e-3 for Ex in sunE)
R['Sun: P_ee falls from near the vacuum average c13^4 (1 - sin^2(2 theta12)/2) + s13^4 (within 0.03 at 0.3 MeV) to about 0.31 at 10 MeV'] = (
    abs(sunE[0.3][0] - ((1 - sn13**2)**2*(1 - 0.5*math.sin(2*th)**2) + sn13**4)) < 0.03 and sunE[0.3][0] > sunE[1.0][0] > sunE[5.0][0] > sunE[10.0][0])
# supernova MSW resonance densities (E = 20 MeV, Y_e = 0.5): rho_res = dm^2 cos(2 theta) / (2 E V_CC(1 g/cm^3) Y_e)
def rho_res(dmsq, theta, E_MeV, Ye=0.5): return dmsq*math.cos(2*theta)/(2*E_MeV*1e6*Vcc(1.0, Ye))
R['NUM supernova MSW resonance densities for E = 20 MeV, Y_e = 0.5: H (dm31, theta13), L (dm21, theta12) [g/cm^3]'] = (
    rho_res(dm31, math.asin(sn13), 20.0), rho_res(dm21, th, 20.0))
# cosmic neutrino background
zeta3 = 1.2020569031595942; kB = 8.617333262e-5                    # eV/K
Tg = 2.7255; Tnu = (4/11)**(1/3)*Tg
n_nu = 3*zeta3/(4*math.pi**2)*2*(kB*Tnu/HBARC_EVCM)**3            # per flavour, nu + nubar, one helicity each
R['NUM cosmic neutrino background: T_nu = (4/11)^(1/3) T_gamma [K], number density per flavour [cm^-3]'] = (Tnu, n_nu)
R['CNB: T_nu = 1.945 K and n = 112 cm^-3 per flavour'] = abs(Tnu - 1.945) < 0.001 and abs(n_nu - 112) < 0.5

# ======================================================================= 6. dirac16complex at the brane: the intertwiner U
G = {A: np.array(G17[A].tolist(), dtype=complex) for A in range(8)}
g4 = {0: np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]], complex)}
sig = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
for i in range(3):
    g4[i + 1] = np.block([[np.zeros((2, 2)), sig[i]], [-sig[i], np.zeros((2, 2))]])
g5 = 1j*g4[0] @ g4[1] @ g4[2] @ g4[3]
ek = [g4[0], 1j*g4[1], 1j*g4[2], 1j*g4[3]]                         # Hermitian Euclidean Clifford generators
rho = [1j*e_ for e_ in ek]                                         # rho_k (k = 4..7): anti-Hermitian, {rho_k, rho_l} = -2 delta
I4 = np.eye(4)
R2 = {mu: np.kron(g4[mu], I4) for mu in range(4)}
for k in range(4): R2[k + 4] = np.kron(g5, rho[k])
eta = np.diag([1.0] + [-1.0]*7)
R['brane representation: gamma_mu (x) 1_4 (mu = 0..3) and gamma5 (x) rho_k (k = 4..7) satisfy {R_A, R_B} = 2 eta_AB'] = all(
    np.allclose(R2[A] @ R2[B] + R2[B] @ R2[A], 2*eta[A, B]*np.eye(16)) for A in range(8) for B in range(8))
def prods(Gs):
    out = []
    for r in range(9):
        for K in itertools.combinations(range(8), r):
            M = np.eye(16, dtype=complex)
            for A in K: M = M @ Gs[A]
            out.append(M)
    return out
P1, P2 = prods(G), prods(R2)
U = None
for j in range(256):
    Mseed = np.zeros((16, 16), complex); Mseed[j // 16, j % 16] = 1
    Ut = sum(p1 @ Mseed @ np.linalg.inv(p2) for p1, p2 in zip(P1, P2))
    if np.linalg.norm(Ut) > 1e-6:
        U = Ut/math.sqrt(abs(np.trace(Ut.conj().T @ Ut))/16); break
R['intertwiner U (Schur average over the 256 Clifford products): Gamma17[A] U = U R_A for A = 0..7, U unitary'] = all(
    np.allclose(G[A] @ U, U @ R2[A]) for A in range(8)) and np.allclose(U.conj().T @ U, np.eye(16))
chi16 = 1j*np.linalg.multi_dot([G[A] for A in range(8)])
rho4567 = rho[0] @ rho[1] @ rho[2] @ rho[3]
R['8D chirality at the brane: chi17 = U (gamma5 (x) rho4567) U^dag, rho4567 Hermitian with eigenvalues +1, +1, -1, -1'] = (
    np.allclose(chi16, U @ np.kron(g5, rho4567) @ U.conj().T) and np.allclose(rho4567, rho4567.conj().T)
    and np.allclose(sorted(np.linalg.eigvalsh(rho4567)), [-1, -1, 1, 1]))
rng = np.random.default_rng(1719)
z1 = rng.normal(size=16) + 1j*rng.normal(size=16); z2 = rng.normal(size=16) + 1j*rng.normal(size=16)
Ps1, Ps2 = U @ z1, U @ z2
bar16 = lambda P: P.conj() @ G[0]
bar4 = lambda zz: [zz.reshape(4, 4)[:, a].conj() @ g4[0] for a in range(4)]
R['bilinears at the brane: Psibar Psi\' = sum_a psibar_a psi\'_a and Psibar Gamma^mu Psi\' = sum_a psibar_a gamma^mu psi\'_a (mu = 0..3)'] = (
    np.isclose(bar16(Ps1) @ Ps2, sum(bar4(z1)[a] @ z2.reshape(4, 4)[:, a] for a in range(4)))
    and all(np.isclose(bar16(Ps1) @ (eta[mu, mu]*G[mu]) @ Ps2, sum(bar4(z1)[a] @ (eta[mu, mu]*g4[mu]) @ z2.reshape(4, 4)[:, a] for a in range(4))) for mu in range(4)))
# Kaluza-Klein masses on the torus: (gamma.p (x) 1 + gamma5 (x) rho.q/b0 - m) has zero modes iff p^2 = m^2 + |q|^2/b0^2
mPsi, b0 = 0.7, 1.3; qv = np.array([2.0, -1.0, 0.5, 3.0])*2*math.pi
mn = math.sqrt(mPsi**2 + qv @ qv/b0**2)
Xkk = np.kron(g5, sum(qv[k]*rho[k] for k in range(4)))/b0
def det_at(E_):
    pv = np.array([E_, 0.3, -0.2, 0.4])
    gp_ = sum(eta[mu, mu]*pv[mu]*g4[mu] for mu in range(4))       # gamma^mu p_mu with p_mu = eta_mu mu p^mu
    return np.linalg.det(np.kron(gp_, I4) + Xkk - mPsi*np.eye(16))
Eon = math.sqrt(mn**2 + 0.3**2 + 0.2**2 + 0.4**2)
R['Kaluza-Klein: the 16x16 operator gamma.p (x) 1 + gamma5 (x) rho.q/b0 - m is singular exactly on p0^2 - |p|^2 = m^2 + |q|^2/b0^2'] = (
    abs(det_at(Eon)) < 1e-6*abs(det_at(1.1*Eon)) and abs(det_at(1.1*Eon)) > 1e-3)
R['KK operator: X^2 = -|q|^2/b0^2 and X anticommutes with gamma^mu (x) 1'] = np.allclose(Xkk @ Xkk, -(qv @ qv)/b0**2*np.eye(16)) and all(
    np.allclose(Xkk @ np.kron(g4[mu], I4) + np.kron(g4[mu], I4) @ Xkk, 0) for mu in range(4))

# ======================================================================= 7. portals: neutrino mixing with dirac16complex
# Dirac mass matrix in the basis (nu_L(3), psi_L(4)) x (nu_R(3), psi_R(4)): [[mD, mu],[0, mPsi 1_4]]
mD = Upmns @ np.diag([m1, m2, m3])                                 # eV, charged-lepton basis
mPs = 1.0e9                                                        # dirac16complex zero-mode mass: 1 GeV in eV
muP = np.zeros((3, 4), complex); muP[:, 0] = [3e4, 1e4, -2e4]; muP[:, 1] = [1e4, 2e4, 1e4]   # portal masses mu_ia = y_ia v/(sqrt2 sqrt(V_perp)) [eV]
Mfull = np.block([[mD, muP], [np.zeros((4, 3)), mPs*np.eye(4)]])
UL, sv, VRh = np.linalg.svd(Mfull)
order = np.argsort(sv); sv = sv[order]; UL = UL[:, order]
Nact = UL[:3, :3]                                                  # active flavours x light mass states
eta_nu = np.eye(3) - Nact @ Nact.conj().T
R['portal mixing: light masses (eV) and heavy masses (GeV) of the 7x7 Dirac mass matrix'] = (sv[:3].tolist(), (sv[3:]/1e9).tolist())
R['portal mixing: N N^dag = 1 - mu mu^dag/m^2 to second order (non-unitary light mixing), eta_ee'] = (
    np.allclose(eta_nu, muP @ muP.conj().T/mPs**2, atol=1e-14), eta_nu[0, 0].real)
def P_nonunitary(N_, masses, L, E):
    ph = np.exp(-1j*masses**2*phase_factor(L, E))
    A = N_ @ np.diag(ph) @ N_.conj().T
    return np.abs(A.T)**2
Pnu = P_nonunitary(Nact, sv[:3], 1300.0, 2.5)
R['non-unitary oscillations: sum_b P(a->b) = (N N^dag)_aa^2-type deficit < 1 (zero-distance effect), complemented by the heavy states'] = all(
    Pnu[a].sum() < 1 for a in range(3)) and np.allclose(P_nonunitary(Nact, sv[:3], 0.0, 2.5).diagonal(), np.abs(np.diag(Nact @ Nact.conj().T))**2)
R['portal decoupling: for mu -> 0 the light block becomes the PMNS matrix and the oscillation probabilities those of section 19.1'] = (
    np.allclose(P_nonunitary(Upmns, np.array([m1, m2, m3]), 1300.0, 2.5), P_vac(Upmns, 1300.0, 2.5), atol=1e-14))
# vector portal: potential from a dirac16complex background n_Psi (forward scattering), non-universal couplings
epsP = np.array([1.0, 0.0, 0.0])*Vcc(2.848, 0.5)*0.2                # V_P^alpha = G_P^alpha n_Psi, 20% of V_CC for nu_e only
PmP = P_matter(Upmns, 1300.0, 2.5, 2.848, 0.5, extra=epsP)
R['NUM vector portal (V_P,e = 0.2 V_CC): P(mu->e) at DUNE with and without the dirac16complex potential'] = (PmP[1, 0], Pm[1, 0])

# ======================================================================= 8. the modern epoch on the brane: Friedmann equations
H0 = 67.4; Om = 0.315; Or = 9.14e-5; OL = 1 - Om - Or
Mpc_km = 3.0856775814913673e19; Gyr_s = 3.15576e16
t0 = quad(lambda a_: 1/(a_*math.sqrt(Or/a_**4 + Om/a_**3 + OL)), 0, 1, epsabs=1e-12, epsrel=1e-12)[0]/(H0/Mpc_km)/Gyr_s
R['NUM age of the universe from the brane Friedmann equation (Planck 2018: H0 = 67.4, Omega_m = 0.315) [Gyr]'] = t0
R['age 13.79 Gyr (Planck 2018: 13.787 +- 0.020)'] = abs(t0 - 13.79) < 0.03
k8, rho4, Vp = sp.symbols('kappa8 rho4 V_perp', positive=True)
R['ADD relation: 3 Ha^2 = kappa8 rho4/V_perp = 8 pi G4 rho4 with G4 = G8/V_perp'] = sp.simplify(k8*rho4/Vp - 8*sp.pi*(k8/(8*sp.pi)/Vp)*rho4) == 0
MPl = 2.435e18; M8 = 1.0e4
Vperp = MPl**2/M8**6
R['NUM ADD example (reduced 8D Planck mass 10 TeV): V_perp [GeV^-4], V_perp^(1/4) [m], lowest KK mass 2 pi/V_perp^(1/4) [MeV]'] = (
    Vperp, Vperp**0.25*HBARC, 2*math.pi/Vperp**0.25*1e3)

if __name__ == '__main__':
    bad = [k for k, v in R.items() if isinstance(v, (bool, np.bool_)) and not v]
    nbool = sum(isinstance(v, (bool, np.bool_)) for v in R.values())
    for k, v in R.items():
        print(f'{k} : {v}')
    print(f'boolean checks: {nbool}   failed: {len(bad)}')
    print('VERDICT:', 'ALL CHECKS PASSED' if not bad else 'FAILED: ' + '; '.join(bad))
    sys.exit(1 if bad else 0)
