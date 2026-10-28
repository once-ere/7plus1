#!/usr/bin/env python3
"""Independent Python (sympy / numpy / scipy / mpmath) check of sections 20-23 of 7+1.nb: the coupled field equations of
the 7+1 model with every interaction, and their accurate approximate solutions in the three gravity regimes.

  A. every interaction: energy-momentum tensors of the stabilizing potential, the smeared brane tension and brane fluid and
     the dirac16complex gas; conservation; the minisuperspace field equations; the radion equation, its mass and its
     equilibrium; the brane Friedmann equation; the cubic Lovelock density in invariant form (exact test).
  B. regime I (primordial, units l2 = sqrt(alpha2) = 1, xi = alpha3/alpha2^2 = 1e-6): the exact vacuum solution (an algebraic
     curve and one quadrature), the cubic, Gauss-Bonnet and Einstein asymptotic solutions and their accuracy, the coupled
     solution with a massless dirac16complex gas (matter-dominated cubic regime, decompactification without stabilizer),
     Einstein-Kasner pair creation (power-law tail), the Standard Model brane (negligible).
  C/D. regime II (modern): capture of the auxiliary space by the stabilizing potential (radion oscillations), radion decay and
     the moduli problem, the Standard Model plasma g*(T), the radiation era t(T), the late universe (composite analytic
     solution, age), the radion today.
  E. regime III (far future): de Sitter attractor, the radion in de Sitter (Gauss-Bonnet corrections, relaxation),
     eight-dimensional de Sitter without stabilizer, Gibbons-Hawking pair creation.

Run from the repository root:   python3 tools/ref_coupled.py
Every line of output is  name : value ; the script exits with status 1 if any check is not True.
"""
# ======================================================================= A. every interaction: stress tensors, conservation, radion, L3
import itertools, math, os, sys
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lovelock71 import cosmological_closed_forms, lovelockL

R = {}
CF = cosmological_closed_forms()                 # E^(k)c_c, k = 1, 2, 3, c = 0, 1, 4, in H_a, H_b, dH_a, dH_b
Ha, Hb, dHa, dHb = sp.symbols('H_a H_b dH_a dH_b', real=True)
al2, al3, Lb, kap = sp.symbols('alpha2 alpha3 Lambdabar kappa', real=True)
TOT = {c: sp.expand(CF[(1, c)] + al2*CF[(2, c)] + al3*CF[(3, c)] - Lb) for c in (0, 1, 4)}

# A1: stress tensor of matter whose comoving energy W(a, b) = a^3 b^4 eps depends only on the scale factors
t = sp.Symbol('x0', real=True)
a = sp.Function('a', positive=True)(t); b = sp.Function('b', positive=True)(t)
A_, B_ = sp.symbols('A B', positive=True)
def stress_from_W(W):
    """eps, pa, pb of S_m = -Int n W(a, b) dx0 (n = lapse): eps = W/(a^3 b^4), pa = -W_a/(3a^2b^4), pb = -W_b/(4a^3b^3)"""
    return (sp.simplify(W/(A_**3*B_**4)), sp.simplify(-sp.diff(W, A_)/(3*A_**2*B_**4)), sp.simplify(-sp.diff(W, B_)/(4*A_**3*B_**3)))
Uf = sp.Function('U'); fW, Lc, rho4, wW, NW = sp.symbols('f_W L_c rho_4 w_W N_W', positive=True)
eU = stress_from_W(A_**3*B_**4*Uf(B_))
R["stabilizing potential: S_U = -Int sqrt|g| U(b) gives eps = U, pa = -U, pb = -U - (b/4) U'(b)"] = (
    sp.simplify(eU[0] - Uf(B_)) == 0 and sp.simplify(eU[1] + Uf(B_)) == 0 and sp.simplify(eU[2] + Uf(B_) + B_*sp.diff(Uf(B_), B_)/4) == 0)
eT = stress_from_W(A_**3*fW/Lc**4)
R["brane tension, smeared over the torus (delta^4(y) -> 1/Lc^4): eps = f_W/V_perp, pa = -f_W/V_perp, pb = 0, V_perp = b^4 Lc^4"] = (
    sp.simplify(eT[0] - fW/(B_**4*Lc**4)) == 0 and sp.simplify(eT[1] + fW/(B_**4*Lc**4)) == 0 and eT[2] == 0)
eB = stress_from_W(A_**3*NW*A_**(-3*(1 + wW))/Lc**4)
R["brane fluid (Standard Model plasma, rho_4 = N_W a^(-3(1+w))), smeared: eps = rho_4/V_perp, pa = w eps, pb = 0"] = (
    sp.simplify(eB[0] - NW*A_**(-3*(1 + wW))/(B_**4*Lc**4)) == 0 and sp.simplify(eB[1] - wW*eB[0]) == 0 and eB[2] == 0)
kav, kbv = sp.symbols('k_a k_b', positive=True)
om_ = sp.sqrt(kav**2/A_**2 + kbv**2/B_**2)
eG = stress_from_W(om_)          # one massless particle of comoving momentum (k_a, k_b); the gas is a sum of such terms
R["massless dirac16complex particle (comoving momenta k_a, k_b): eps V = omega, pa V = k_a^2/(3 a^2 omega), pb V = k_b^2/(4 b^2 omega), V = a^3 b^4"] = (
    sp.simplify(eG[1]*A_**3*B_**4 - kav**2/(3*A_**2*om_)) == 0 and sp.simplify(eG[2]*A_**3*B_**4 - kbv**2/(4*B_**2*om_)) == 0)
Wf = sp.Function('W')
W_t = Wf(a, b)
eps_t = W_t/(a**3*b**4)
pa_t = -sp.diff(Wf(A_, B_), A_).subs({A_: a, B_: b})/(3*a**2*b**4)
pb_t = -sp.diff(Wf(A_, B_), B_).subs({A_: a, B_: b})/(4*a**3*b**3)
cons_t = sp.diff(eps_t, t) + 3*sp.diff(a, t)/a*(eps_t + pa_t) + 4*sp.diff(b, t)/b*(eps_t + pb_t)
R["every such matter is conserved: eps' + 3Ha(eps + pa) + 4Hb(eps + pb) = 0 for arbitrary W(a, b), a(x0), b(x0)"] = sp.simplify(cons_t.doit()) == 0
# minisuperspace: S = -(1/(2 kappa)) Int n a^3 b^4 Sum alpha_k L_k - Int n W(a, b), with the variational identities of section 12.6
# EL_n(n a^3 b^4 L_k) = -2 a^3 b^4 E^(k)0_0, EL_a = -6 a^2 b^4 E^(k)1_1, EL_b = -8 a^3 b^3 E^(k)4_4  (S00, S11, S44 = Sum alpha_k E^(k))
S00, S11, S44 = sp.symbols('S00 S11 S44')
Wab = Wf(A_, B_)
ELn = -(1/(2*kap))*(-2*A_**3*B_**4*S00) - Wab
ELa = -(1/(2*kap))*(-6*A_**2*B_**4*S11) - sp.diff(Wab, A_)
ELb = -(1/(2*kap))*(-8*A_**3*B_**3*S44) - sp.diff(Wab, B_)
e_, pa_, pb_ = Wab/(A_**3*B_**4), -sp.diff(Wab, A_)/(3*A_**2*B_**4), -sp.diff(Wab, B_)/(4*A_**3*B_**3)
R["field equations from S_grav + S_m: EL_n = (a^3 b^4/kappa)(S00 - kappa eps), EL_a = (3a^2 b^4/kappa)(S11 + kappa pa), EL_b = (4a^3 b^3/kappa)(S44 + kappa pb), S = Sum alpha_k E^(k) (minisuperspace, lapse n = 1)"] = (
    sp.simplify(ELn - A_**3*B_**4/kap*(S00 - kap*e_)) == 0 and sp.simplify(ELa - 3*A_**2*B_**4/kap*(S11 + kap*pa_)) == 0 and
    sp.simplify(ELb - 4*A_**3*B_**3/kap*(S44 + kap*pb_)) == 0)

# A2: the radion (b = b0 e^beta): exact Hb' equation of Einstein gravity and linearization with the stabilizing potential
eps_m, pa_m, pb_m, KU, b0s, U0 = sp.symbols('epsilon_m p_a p_b K b_0 U_0', real=True)
beta = sp.Symbol('beta', real=True)
Uq = U0 + KU/2*beta**2                                    # U(b) = U0 + (K/2) ln^2(b/b0), beta = ln(b/b0)
bUp = KU*beta                                             # b U'(b)
epsT = eps_m + Uq; paT = pa_m - Uq; pbT = pb_m - Uq - bUp/4
sol_E = sp.solve([TOT[1].subs({al2: 0, al3: 0}) + kap*paT, TOT[4].subs({al2: 0, al3: 0}) + kap*pbT], [dHa, dHb], dict=True)[0]
dHb_E = sp.expand(sol_E[dHb])
# use the constraint to eliminate Ha^2: 3Ha^2 = kappa epsT + Lb - 12 Ha Hb - 6 Hb^2
dHb_E2 = sp.expand(dHb_E.subs(Ha**2, (kap*epsT + Lb - 12*Ha*Hb - 6*Hb**2)/3))
target = sp.expand(-3*Ha*Hb - 4*Hb**2 + kap*(epsT/6 - paT/2 + pbT/3) + Lb/3)
R["radion equation (Einstein): Hb' + 3 Ha Hb + 4 Hb^2 = kappa (eps/6 - pa/2 + pb/3) + Lambda/3 (Hb = beta', exact)"] = sp.simplify(dHb_E2 - target) == 0
lin = sp.expand(target.subs({Hb: 0}))
mr2 = -sp.diff(lin, beta).subs(beta, 0)
R["radion mass: with U = U0 + (K/2) ln^2(b/b0) the restoring term is -m_r^2 beta, m_r^2 = kappa K/12"] = sp.simplify(mr2 - kap*KU/12) == 0
# A3: static equilibrium (Hb = Hb' = 0) with Gauss-Bonnet and a smeared brane tension
fV = sp.Symbol('fV', real=True)                          # f_W/V_perp
epsS = eps_m + Uq + fV; paS = pa_m - Uq - fV; pbS = pb_m - Uq - bUp/4
eqs_static = [TOT[c].subs({Hb: 0, dHb: 0}) - kap*s for c, s in ((0, epsS), (1, -paS), (4, -pbS))]
s_st = sp.solve([eqs_static[0], eqs_static[1]], [Lb, dHa], dict=True)[0]
eq44 = sp.expand(eqs_static[2].subs(s_st))
equil = sp.expand(kap*bUp/4 - (Lb + kap*Uq + 2*kap*fV + kap*(eps_m - 3*pa_m)/2 + kap*pb_m - 12*al2*Ha**2*(dHa + Ha**2)).subs(s_st))
R["static auxiliary space: (kappa/4) b U'(b) = Lambda + kappa U + 2 kappa f_W/V_perp + (kappa/2)(eps - 3pa) + kappa pb - 12 alpha2 Ha^2 (Ha' + Ha^2)"] = sp.simplify(eq44 + equil) == 0
R["static auxiliary space: E^(3) does not contribute, and Ha' = -(kappa/2)(eps + pa) (brane Friedmann acceleration)"] = (
    all(sp.expand(CF[(3, c)].subs({Hb: 0, dHb: 0})) == 0 for c in (0, 1, 4)) and sp.simplify(s_st[dHa] + kap*(eps_m + pa_m)/2) == 0)
# 4D constants
R["four-dimensional Friedmann equation on the brane (b static): 3 Ha^2 = 8 pi G4 rho4 + Lambda4 with 8 pi G4 = kappa/V_perp and Lambda4 = Lambda + kappa U0 + kappa f_W/V_perp"] = (
    sp.simplify(TOT[0].subs({Hb: 0, al2: al2, al3: al3}) - kap*(sp.Symbol('rho4')/sp.Symbol('Vp') + U0 + fV) -
                (3*Ha**2 - (kap/sp.Symbol('Vp'))*sp.Symbol('rho4') - (Lb + kap*U0 + kap*fV))) == 0)

# A4: the cubic Lovelock density written with curvature invariants (exact integer test on random algebraic curvature tensors)
def random_riemann(rng, D=8, nS=3, amp=3):
    R4 = np.zeros((D, D, D, D), dtype=np.int64)
    for _ in range(nS):
        S = rng.integers(-amp, amp + 1, size=(D, D)); S = S + S.T
        R4 += np.einsum('ac,bd->abcd', S, S) - np.einsum('ad,bc->abcd', S, S)
    return R4
etaI = np.diag([1] + [-1]*7).astype(np.int64)
def cubic_invariants(Rl):
    gi = etaI
    Rud_ = np.einsum('ae,bf,efcd->abcd', gi, gi, Rl)
    Ric = np.einsum('ab,acbd->cd', gi, Rl); Rs = int(np.einsum('cd,cd->', gi, Ric))
    RicU = gi @ Ric @ gi; RicM = gi @ Ric
    RU = np.einsum('ae,bf,cg,dh,efgh->abcd', gi, gi, gi, gi, Rl)
    T = np.einsum('ce,df,aebf->acbd', gi, gi, Rl)
    RaU = np.einsum('cf,dg,eh,bfgh->bcde', gi, gi, gi, Rl)
    return {'R^3': Rs**3, 'R Ric.Ric': Rs*int(np.einsum('ab,ab->', Ric, RicU)), 'Ric^3': int(np.einsum('ab,bc,ca->', RicM, RicM, RicM)),
            'Ric Ric Riem': int(np.einsum('ac,bd,abcd->', RicU, RicU, Rl)), 'R Riem.Riem': Rs*int(np.einsum('abcd,abcd->', RU, Rl)),
            'Ric Riem Riem': int(np.einsum('ab,acde,bcde->', RicU, Rl, RaU)), 'Riem^3 (pairs)': int(np.einsum('abcd,cdef,efab->', Rud_, Rud_, Rud_)),
            'Riem^3 (crossed)': int(np.einsum('acbd,cedf,eafb->', T, T, T)), 'R': Rs, 'Ric.Ric': int(np.einsum('ab,ab->', Ric, RicU)),
            'Riem.Riem': int(np.einsum('abcd,abcd->', RU, Rl))}
L3COEF = {'R^3': 1, 'R Ric.Ric': -12, 'Ric^3': 16, 'Ric Ric Riem': 24, 'R Riem.Riem': 3, 'Ric Riem Riem': -24, 'Riem^3 (pairs)': 4, 'Riem^3 (crossed)': -8}
rngL = np.random.default_rng(3108)
okL3 = True; okL2 = True
for _ in range(4):
    Rl = random_riemann(rngL)
    Rud_ = np.einsum('ae,bf,efcd->abcd', etaI, etaI, Rl)
    dct = {k: int(v) for k, v in np.ndenumerate(Rud_) if v != 0}
    I = cubic_invariants(Rl)
    okL3 = okL3 and lovelockL(dct, 3) == sum(c*I[k] for k, c in L3COEF.items())
    okL2 = okL2 and lovelockL(dct, 2) == I['R']**2 - 4*I['Ric.Ric'] + I['Riem.Riem']
R["cubic Lovelock density: L_3 = R^3 - 12 R Ric.Ric + 16 Ric^3 + 24 R^ac R^bd R_abcd + 3 R Riem.Riem - 24 R^ab R_a^cde R_bcde + 4 R^ab_cd R^cd_ef R^ef_ab - 8 R_a^c_b^d R_c^e_d^f R_e^a_f^b (exact, four random integer algebraic curvature tensors)"] = okL3
R["Gauss-Bonnet density: L_2 = R^2 - 4 Ric.Ric + Riem.Riem (same tensors)"] = okL2

# ======================================================================= B. regime I: the primordial Lovelock regime (units l2 = sqrt(alpha2) = 1)
XI = 1e-6                                                         # xi = alpha3/alpha2^2
u_, X_, xi_ = sp.symbols('u X xi', real=True)
con_u = sp.expand(TOT[0].subs({al2: 1, al3: xi_, Lb: 0}).subs(Hb, u_*Ha)/Ha**2)
con_u = sp.expand(con_u.subs(Ha, sp.sqrt(X_)))
Eu = sp.factor(con_u.coeff(X_, 0)); Gu = sp.factor(con_u.coeff(X_, 1)); Ku = sp.factor(con_u.coeff(X_, 2)/xi_)
R["vacuum constraint in u = Hb/Ha, X = Ha^2 (alpha2 = 1, alpha3 = xi): E(u) + X G(u) + xi X^2 K(u) = 0 with E = 3(2u^2 + 4u + 1), G = -12u(u^3 + 12u^2 + 18u + 4), K = 360 u^3 (3u + 4)"] = (
    sp.expand(Eu - 3*(2*u_**2 + 4*u_ + 1)) == 0 and sp.expand(Gu + 12*u_*(u_**3 + 12*u_**2 + 18*u_ + 4)) == 0 and sp.expand(Ku - 360*u_**3*(3*u_ + 4)) == 0
    and sp.expand(con_u - (Eu + X_*Gu + xi_*X_**2*Ku)) == 0)
# the (dHa, dHb) of the vacuum equations
Avac = sp.Matrix([[sp.diff(TOT[1], dHa), sp.diff(TOT[1], dHb)], [sp.diff(TOT[4], dHa), sp.diff(TOT[4], dHb)]])
Bvac = sp.Matrix([TOT[1].subs({dHa: 0, dHb: 0}), TOT[4].subs({dHa: 0, dHb: 0})])
f_A = sp.lambdify((Ha, Hb, al2, al3), Avac, 'numpy')
f_B = sp.lambdify((Ha, Hb, al2, al3, Lb), Bvac, 'numpy')
f_C = sp.lambdify((Ha, Hb, al2, al3, Lb), TOT[0], 'numpy')
detA = sp.expand(Avac.det()); adjB = Avac.adjugate()*(-Bvac)
f_dHa_mp = sp.lambdify((Ha, Hb, al3), (adjB[0]/detA).subs({al2: 1, Lb: 0}), 'mpmath')
fE, fG, fK = [sp.lambdify(u_, e, 'mpmath') for e in (Eu, Gu, Ku)]
mp.mp.dps = 30
u1 = -1 + 1/mp.sqrt(2)
u2 = mp.findroot(lambda z: z**3 + 12*z**2 + 18*z + 4, -0.27)
p1, q1 = 1/(3 + 4*u1), u1/(3 + 4*u1); p2, q2 = 3/(3 + 4*u2), 3*u2/(3 + 4*u2)
R['Einstein regime: u1 = Hb/Ha = -1 + 1/sqrt(2), root of E(u); Kasner p1 = 1/(3 + 4u1) = (1 + 2 sqrt 2)/7, q1 = u1 p1 = (2 - 3 sqrt 2)/14'] = (
    abs(fE(u1)) < mp.mpf(10)**-25 and abs(p1 - (1 + 2*mp.sqrt(2))/7) < mp.mpf(10)**-25 and abs(q1 - (2 - 3*mp.sqrt(2))/14) < mp.mpf(10)**-25)
R['Gauss-Bonnet regime: u2 = Hb/Ha, root of u^3 + 12u^2 + 18u + 4 in (-1/3, 0); Kasner p2 = 3/(3 + 4u2), q2 = u2 p2 (3p + 4q = 3)'] = (
    abs(fG(u2)) < mp.mpf(10)**-25 and -mp.mpf(1)/3 < u2 < 0 and abs(3*p2 + 4*q2 - 3) < mp.mpf(10)**-25)
R['NUM u1, u2, p1, q1, p2, q2'] = tuple(float(x) for x in (u1, u2, p1, q1, p2, q2))

def Xcurve(uu, xv):
    """the positive root Ha^2 of E + X G + xi X^2 K = 0 (K < 0 on (u1, 0))"""
    E_, G_, K_ = fE(uu), fG(uu), fK(uu)
    return (G_ + mp.sqrt(G_**2 - 4*xv*K_*E_))/(2*xv*abs(K_))
def peff_u(uu, xv):
    h = mp.sqrt(Xcurve(uu, xv)); return -h**2/f_dHa_mp(h, uu*h, xv)
xiM = mp.mpf('1e-6')
okc = True
for uu in (u1 + mp.mpf('1e-6'), mp.mpf(-0.28), u2, mp.mpf(-0.1), mp.mpf(-1e-4)):
    X0 = Xcurve(uu, xiM); h = mp.sqrt(X0)
    val = f_C(float(h), float(uu*h), 1.0, float(xiM), 0.0)
    okc = okc and abs(fE(uu) + X0*fG(uu) + xiM*X0**2*fK(uu)) < mp.mpf(10)**-20*max(1, X0**2)
R['exact vacuum solution: for every u in (u1, 0) the curve Ha = sqrt(X(u)), Hb = u Ha, X(u) = [G + sqrt(G^2 - 4 xi K E)]/(2 xi |K|), solves the constraint (residual < 1e-20 relative)'] = okc
pe_vals = [(float(uu), float(peff_u(uu, xiM))) for uu in (u1 + mp.mpf('1e-10'), u2, mp.mpf('-1e-6'))]
R['NUM effective Kasner exponent p_eff = -Ha^2/Ha\' along the curve at u -> u1, u = u2, u -> 0 (xi = 1e-6)'] = tuple(v for _, v in pe_vals)
R['the effective exponent goes from p1 (Einstein, u -> u1) through about p2 (Gauss-Bonnet, u = u2) to 1 (cubic, u -> 0)'] = (
    abs(pe_vals[0][1] - float(p1)) < 1e-8 and abs(pe_vals[1][1] - float(p2)) < 0.01 and abs(pe_vals[2][1] - 1) < 1e-5)

# cubic regime series (general alpha2, alpha3)
tt = sp.Symbol('tau', positive=True)
cA0, cA1, cA2, ch1, ch2, ch3 = sp.symbols('A0 A1 A2 h1 h2 h3')
h0s = -sp.sqrt(al2/(30*al3))
Hax = 1/tt + cA0 + cA1*tt + cA2*tt**2; Hbx = h0s + ch1*tt + ch2*tt**2 + ch3*tt**3
subx = {Ha: Hax, Hb: Hbx, dHa: sp.diff(Hax, tt), dHb: sp.diff(Hbx, tt), Lb: 0}
ser = {c: sp.expand(sp.expand(TOT[c].subs(subx))*tt**4) for c in (0, 1, 4)}
co = {c: [sp.simplify(ser[c].coeff(tt, n)) for n in range(4)] for c in (0, 1, 4)}
s1sol = sp.solve([co[0][2], co[4][1]], [ch1, cA0], dict=True)[0]
s2sol = sp.solve([co[0][3].subs(s1sol), co[4][2].subs(s1sol)], [ch2, cA1], dict=True)[0]
A0x, h1x = sp.simplify(s1sol[cA0]), sp.simplify(s1sol[ch1]); A1x, h2x = sp.simplify(s2sol[cA1]), sp.simplify(s2sol[ch2])
R['cubic regime (tau -> 0): Ha = 1/tau + A0 + A1 tau, Hb = h0 + h1 tau + h2 tau^2 solves all three equations to the orders used; h0 = -sqrt(alpha2/(30 alpha3))'] = (
    all(co[c][0] == 0 for c in (0, 1, 4)) and sp.simplify(co[0][1]) == 0 and sp.simplify(co[1][1]) == 0 and
    sp.simplify(co[1][2].subs(s1sol)) == 0 and sp.simplify(co[1][3].subs(s1sol).subs(s2sol)) == 0)
def same(e1_, e2_):
    """equal as functions of alpha2, alpha3 > 0 (checked at three exact rational points to 40 digits)"""
    return all(abs(sp.N((e1_ - e2_).subs({al2: x, al3: y}), 40)) < 1e-30 for x, y in ((1, sp.Rational(1, 10**6)), (sp.Rational(7, 3), sp.Rational(2, 5)), (sp.Rational(1, 2), 3)))
R['cubic regime: h1 = alpha2/(16 alpha3) - 1/(32 alpha2), A0 = sqrt(30)(58 alpha2^2 - 45 alpha3)/(1920 alpha2^(3/2) alpha3^(1/2))'] = (
    same(h1x, al2/(16*al3) - 1/(32*al2)) and same(A0x, sp.sqrt(30)*(58*al2**2 - 45*al3)/(1920*al2**sp.Rational(3, 2)*sp.sqrt(al3))))
R['cubic regime: A1 = -3647 alpha2/(92160 alpha3) - 223/(6144 alpha2) + 645 alpha3/(8192 alpha2^3), h2 = sqrt(30)(-1212 alpha2^4 + 1180 alpha2^2 alpha3 + 225 alpha3^2)/(307200 alpha2^(5/2) alpha3^(3/2))'] = (
    same(A1x, -sp.Rational(3647, 92160)*al2/al3 - sp.Rational(223, 6144)/al2 + sp.Rational(645, 8192)*al3/al2**3) and
    same(h2x, sp.sqrt(30)*(-1212*al2**4 + 1180*al2**2*al3 + 225*al3**2)/(307200*al2**sp.Rational(5, 2)*al3**sp.Rational(3, 2))))
cub = {k: float(v.subs({al2: 1, al3: sp.Float(XI, 30)})) for k, v in (('h0', h0s), ('h1', h1x), ('h2', h2x), ('A0', A0x), ('A1', A1x))}
R['NUM cubic-regime coefficients (xi = 1e-6): h0, h1, h2, A0, A1'] = (cub['h0'], cub['h1'], cub['h2'], cub['A0'], cub['A1'])

# Einstein regime expansion: u = u1 + d1 Ha^2, p_eff = p1 + e1 Ha^2
eps_ = sp.Symbol('epsilon', positive=True); d1s, d2s = sp.symbols('d1 d2')
u1x = -1 + 1/sp.sqrt(2)
uu_e = u1x + d1s*eps_**2 + d2s*eps_**4
cs = sp.series(sp.expand(TOT[0].subs({Lb: 0}).subs(Hb, uu_e*Ha).subs(Ha, eps_)), eps_, 0, 7).removeO()
d1x = sp.solve(sp.expand(cs.coeff(eps_, 4)), d1s)[0]
d2x = sp.solve(sp.expand(cs.coeff(eps_, 6)).subs(d1s, d1x), d2s)[0]
numE = adjB[0].subs({Lb: 0}).subs(Hb, uu_e*Ha).subs(Ha, eps_).subs({d1s: d1x, d2s: d2x})
denE = detA.subs({Lb: 0}).subs(Hb, uu_e*Ha).subs(Ha, eps_).subs({d1s: d1x, d2s: d2x})
peE = sp.series(-eps_**2*denE/numE, eps_, 0, 4).removeO()
e1x = sp.radsimp(sp.simplify(peE.coeff(eps_, 2)))
R['Einstein regime (Ha -> 0): u = u1 + d1 Ha^2 + ..., d1 = alpha2 (16 - 11 sqrt 2)/4; p_eff = p1 + e1 Ha^2 + ..., e1 = alpha2 (1 + 9 sqrt 2)/7'] = (
    sp.simplify(d1x - al2*(16 - 11*sp.sqrt(2))/4) == 0 and sp.simplify(e1x - al2*(1 + 9*sp.sqrt(2))/7) == 0 and
    sp.simplify(peE.coeff(eps_, 0) - (1 + 2*sp.sqrt(2))/7) == 0)
e1v = float(e1x.subs(al2, 1))
# Gauss-Bonnet regime expansion: p_eff = p2 + f1/Ha^2 + f2 (alpha3/alpha2) Ha^2 (numerical, 30 digits)
def u_on_curve(h, xv, guess):
    return mp.findroot(lambda z: fE(z) + h**2*fG(z) + xv*h**4*fK(z), guess)
f1v = float((f_dHa_pe := (lambda h, xv: -h**2/f_dHa_mp(h, u_on_curve(h, xv, u2)*h, xv)))(mp.mpf(10)**5, 0) - p2)*1e10
f2v = float((f_dHa_pe(mp.mpf(10)**8, mp.mpf(10)**-24) - p2)/mp.mpf(10)**-8)
R['NUM Gauss-Bonnet regime: p_eff = p2 + f1/Ha^2 + f2 xi Ha^2 + ..., f1, f2 (alpha2 = 1)'] = (f1v, f2v)

# the vacuum solution as a function of Ha: integrate in sigma = -ln Ha with x0 as a dependent variable; initial data from the cubic series
def rhs_sig(sg, y, al3v, C=0.0, gasf=None):
    ha = math.exp(-sg); hb, x0 = y[0], y[1]
    Am = np.array(f_A(ha, hb, 1.0, al3v), dtype=float); Bm = np.array(f_B(ha, hb, 1.0, al3v, 0.0), dtype=float).reshape(2)
    d = np.linalg.solve(Am, -Bm); dx = -ha/d[0]
    return [d[1]*dx, dx]
tau0 = 1e-8
Ha0 = 1/tau0 + cub['A0'] + cub['A1']*tau0; Hb0 = cub['h0'] + cub['h1']*tau0 + cub['h2']*tau0**2
R['cubic-series initial data at tau = 1e-8 satisfy the constraint to < 1e-12 (relative to 48 alpha2 Ha^3 |Hb|)'] = abs(f_C(Ha0, Hb0, 1.0, XI, 0.0))/(48*Ha0**3*abs(Hb0)) < 1e-12
vac = solve_ivp(rhs_sig, [-math.log(Ha0), -math.log(1e-4)], [Hb0, tau0], args=(XI,), method='DOP853', rtol=1e-13, atol=1e-30, dense_output=True)
def vac_at(H):
    hb, x0 = vac.sol(-math.log(H)); return hb, x0
# s1 from the Einstein two-term form, s2 at the middle of the Gauss-Bonnet regime
s1_list = [vac_at(H)[1] - float(p1)/H + e1v*H for H in (3e-3, 1e-3, 3e-4)]
# s1 by exact quadrature along the curve (the ODE extraction drifts at large tau: constraint violations grow like tau)
fEp, fGp, fKp = [sp.lambdify(u_, sp.diff(e, u_), 'mpmath') for e in (Eu, Gu, Ku)]
def dHdu_c(uu, xv):
    x = Xcurve(uu, xv)
    return -(fEp(uu) + x*fGp(uu) + xv*x**2*fKp(uu))/(fG(uu) + 2*xv*x*fK(uu))/(2*mp.sqrt(x))
def Fc(uu, xv):
    h = mp.sqrt(Xcurve(uu, xv)); return -f_dHa_mp(h, uu*h, xv)
def tau_of_u(uu, xv):
    return mp.re(mp.quad(lambda v: dHdu_c(-v**2, xv)/Fc(-v**2, xv)*2*v, [0, mp.mpf('1e-3'), mp.mpf('0.03'), mp.mpf('0.3'), mp.sqrt(-uu)]))
tau_u2 = tau_of_u(u2, xiM); H_u2 = mp.sqrt(Xcurve(u2, xiM))
Jq = mp.re(mp.quad(lambda w: (1/Fc(u1 + w**2, xiM) - p1/Xcurve(u1 + w**2, xiM))*dHdu_c(u1 + w**2, xiM)*2*w,
                   [0, mp.sqrt(u2 - u1)*mp.mpf('1e-3'), mp.sqrt(u2 - u1)*mp.mpf('0.1'), mp.sqrt(u2 - u1)]))
s1q = float(tau_u2 - p1/H_u2 + Jq)
R['NUM exact quadrature: tau(u2), Ha(u2), s1 = tau(u2) - p1/Ha(u2) + Int_0^Ha(u2) (1/F - p1/H^2) dH'] = (float(tau_u2), float(H_u2), s1q)
R['the ODE solution and the quadrature agree: tau at Ha(u2) to 1e-10, s1 (ODE at Ha = 3e-3) to 1e-6'] = (
    abs(vac_at(float(H_u2))[1] - float(tau_u2)) < 1e-10 and abs(s1_list[0] - s1q) < 1e-6)
s1v = s1q
Hmid = (f1v/(f2v*XI))**0.25
s2v = vac_at(Hmid)[1] - float(p2)/Hmid - f1v/(3*Hmid**3) + f2v*XI*Hmid
R['NUM time shifts: tau = p1/Ha - e1 Ha + s1 (Einstein), tau = p2/Ha + f1/(3 Ha^3) - f2 xi Ha + s2 (Gauss-Bonnet, at Ha_mid = (f1/(f2 xi))^(1/4)): s1, s2, Ha_mid'] = (s1v, s2v, Hmid)
# accuracy of the three regime approximations against the exact solution
def tau_cubic(H):
    tv = 1/H
    for _ in range(60): tv = 1/(H - cub['A0'] - cub['A1']*tv)
    return tv
acc = []
for lh in np.arange(8.0, -3.51, -0.25):
    H = 10**lh; hb, tv = vac_at(H)
    tc = tau_cubic(H); tg = float(p2)/H + f1v/(3*H**3) - f2v*XI*H + s2v; te = float(p1)/H - e1v*H + s1v
    acc.append((lh, tv, abs(tc/tv - 1), abs(tg/tv - 1), abs(te/tv - 1)))
def best_err(lo, hi, col):
    return max(r[col] for r in acc if lo <= r[0] <= hi)
R['NUM accuracy (max relative error of tau(Ha)): cubic for Ha >= 1e4, Gauss-Bonnet for 3 <= Ha <= 300, Einstein for Ha <= 0.1'] = (
    best_err(4.0, 8.0, 2), best_err(0.5, 2.5, 3), best_err(-3.5, -1.0, 4))
R['regime approximations: cubic error < 1e-5 for Ha >= 1e4, Gauss-Bonnet < 1e-3 for 3 <= Ha <= 300, Einstein < 1e-3 for Ha <= 0.1'] = (
    best_err(4.0, 8.0, 2) < 1e-5 and best_err(0.5, 2.5, 3) < 1e-3 and best_err(-3.5, -1.0, 4) < 1e-3)
cross = [(r[0], min(r[2], r[3], r[4])) for r in acc]
R['NUM worst error of the best regime approximation over 1e-3.5 <= Ha <= 1e8 (in the crossover windows)'] = max(c[1] for c in cross)
R['NUM vacuum crossover times: tau at Ha = 1e3 (cubic -> GB) and Ha = 1 (GB -> Einstein)'] = (vac_at(1e3)[1], vac_at(1.0)[1])
ACC_TABLE = acc


# ----------------------------------------------------------------------- B2. coupled to the dirac16complex gas (massless, collisionless)
xgl, wgl = np.polynomial.legendre.leggauss(128)
thg = (xgl + 1)*np.pi/4; wthg = wgl*np.pi/4
c2g, s2g = np.cos(thg)**2, np.sin(thg)**2
wg_ = wthg*c2g*np.sin(thg)**3; wg_ = wg_/np.sum(wg_)                # normalized weight cos^2 sin^3 (7-dimensional momentum sphere)
def gas(la, lb):
    """massless dirac16complex gas with an isotropic comoving momentum distribution: (eps, pa, pb)/eps(a = b = 1)"""
    ia2, ib2 = math.exp(-2*la), math.exp(-2*lb)
    root = np.sqrt(c2g*ia2 + s2g*ib2); pref = math.exp(-3*la - 4*lb)
    return pref*np.sum(wg_*root), pref*np.sum(wg_*c2g*ia2/(3*root)), pref*np.sum(wg_*s2g*ib2/(4*root))
def rhs_gas(sg, y, al3v, Cg):
    ha = math.exp(-sg); hb, la, lb, x0 = y
    e, pa, pb = gas(la, lb)
    Am = np.array(f_A(ha, hb, 1.0, al3v), dtype=float); Bm = np.array(f_B(ha, hb, 1.0, al3v, 0.0), dtype=float).reshape(2) + Cg*np.array([pa, pb])
    d = np.linalg.solve(Am, -Bm); dx = -ha/d[0]
    return [d[1]*dx, ha*dx, hb*dx, dx]
Href = 1.0; Om_ref = 1e-6                                   # the gas is isotropic (a = b = 1) at Ha = 1 and has kappa eps = 1e-6 x 3 Ha^2 there
hb_vac_ref, x_ref = vac_at(Href)
Cg = Om_ref*3*Href**2
hb_ref = brentq(lambda hb: f_C(Href, hb, 1.0, XI, 0.0) - Cg, hb_vac_ref*1.2, hb_vac_ref*0.8)
cb = solve_ivp(rhs_gas, [-math.log(Href), -math.log(1e16)], [hb_ref, 0.0, 0.0, x_ref], args=(XI, Cg), method='DOP853', rtol=1e-12, atol=1e-30, dense_output=True)
cf = solve_ivp(rhs_gas, [-math.log(Href), -math.log(1e-14)], [hb_ref, 0.0, 0.0, x_ref], args=(XI, Cg), method='DOP853', rtol=1e-12, atol=1e-30, dense_output=True)
x_sing = cb.y[3, -1] - 1/1e16*0.75                           # x0 at Ha = 1e16 minus tau = p/Ha with p = 3/4
def coupled_at(H):
    sol = cb if H >= Href else cf
    hb, la, lb, x0 = sol.sol(-math.log(H))
    e, pa, pb = gas(la, lb)
    Am = np.array(f_A(H, hb, 1.0, XI), dtype=float); Bm = np.array(f_B(H, hb, 1.0, XI, 0.0), dtype=float).reshape(2) + Cg*np.array([pa, pb])
    d = np.linalg.solve(Am, -Bm)
    E1v = 3*(H**2 + 4*H*hb + 2*hb**2); E2v = -12*hb*(4*H**3 + 18*H**2*hb + 12*H*hb**2 + hb**3); E3v = 360*XI*H**2*hb**3*(4*H + 3*hb)
    return dict(Ha=H, Hb=hb, la=la, lb=lb, tau=x0 - x_sing, peff=-H**2/d[0], keps=Cg*e, E1=E1v, E2=E2v, E3=E3v, pa_e=pa/e, pb_e=pb/e,
                frac=Cg*e/(abs(E1v) + abs(E2v) + abs(E3v) + 3*H**2))
hstar = -1/math.sqrt(54*XI)
near = coupled_at(1e15)
R['NUM coupled, near the singularity (Ha = 1e15): p_eff, Hb, kappa eps tau^3, alpha2 E2/(kappa eps), alpha3 E3/(kappa eps)'] = (
    near['peff'], near['Hb'], near['keps']*near['tau']**3, near['E2']/near['keps'], near['E3']/near['keps'])
R['matter-dominated cubic regime: p_eff -> 3/4, Hb -> -sqrt(alpha2/(54 alpha3)), kappa eps tau^3 -> 9 alpha2 |h*|, alpha2 E2 : alpha3 E3 : kappa eps -> 9/4 : -5/4 : 1 (to 1e-3)'] = (
    abs(near['peff'] - 0.75) < 1e-3 and abs(near['Hb']/hstar - 1) < 1e-3 and abs(near['keps']*near['tau']**3/(9*abs(hstar)) - 1) < 1e-2
    and abs(near['E2']/near['keps'] - 2.25) < 1e-2 and abs(near['E3']/near['keps'] + 1.25) < 1e-2)
plateau = {}
for name, H in (('cubic (vacuum)', 3e5), ('Gauss-Bonnet', 30.0), ('Einstein', 1e-2)):
    c = coupled_at(H); plateau[name] = (c['peff'], c['frac'])
R['NUM coupled: (p_eff, matter fraction) at Ha = 3e5 (vacuum cubic), 30 (Gauss-Bonnet), 1e-2 (Einstein)'] = tuple(plateau[k] for k in ('cubic (vacuum)', 'Gauss-Bonnet', 'Einstein'))
R['coupled: where the gas is a small fraction, the effective exponents are those of the vacuum solution (1, p2, p1 within 1e-2)'] = (
    abs(plateau['cubic (vacuum)'][0] - 1) < 1e-2 and abs(plateau['Gauss-Bonnet'][0] - float(p2)) < 1e-2 and abs(plateau['Einstein'][0] - float(p1)) < 1e-2)
late = coupled_at(1e-14)
R['NUM coupled, gas-dominated late regime (Ha = 1e-14): tau Ha, tau Hb, pa/eps, pb/eps'] = (late['tau']*late['Ha'], late['tau']*late['Hb'], late['pa_e'], late['pb_e'])
R['without a stabilizing potential the gas-dominated Einstein regime decompactifies: a -> const (tau Ha -> 0), b ~ tau^(2/5) (tau Hb -> 2/5), pb/eps -> 1/4'] = (
    late['tau']*late['Ha'] < 1e-3 and abs(late['tau']*late['Hb'] - 0.4) < 2e-3 and abs(late['pb_e'] - 0.25) < 1e-3)
# growth of the gas fraction in the Einstein regime (auxiliary momenta dominate): kappa eps/(3 Ha^2) ~ tau^(2 - 3 p1 - 5 q1)
gA, gB = coupled_at(1e-2), coupled_at(1e-4)
slope = math.log((gB['keps']/(3*gB['Ha']**2))/(gA['keps']/(3*gA['Ha']**2)))/math.log(gB['tau']/gA['tau'])
R['NUM gas-fraction growth exponent in the Einstein regime (numerical, between Ha = 1e-2 and 1e-4), and 2 - 3p1 - 5q1'] = (slope, float(2 - 3*p1 - 5*q1))
R['gas fraction grows like tau^(2 - 3p1 - 5q1) = tau^1.160 in the Einstein regime (to 2e-2)'] = abs(slope - float(2 - 3*p1 - 5*q1)) < 2e-2

# ----------------------------------------------------------------------- B3. power-law solutions of Einstein gravity with a fluid (pa = wa eps, pb = wb eps)
PP, QQ, CC = sp.symbols('P Q C')
def EKpow(k, c):
    return sp.expand(CF[(k, c)].subs({Ha: PP/tt, Hb: QQ/tt, dHa: -PP/tt**2, dHb: -QQ/tt**2})*tt**(2*k))
def einstein_fluid(wa, wb):
    return sp.solve([EKpow(1, 0) - CC, EKpow(1, 1) + wa*CC, EKpow(1, 4) + wb*CC, 3*PP*(1 + wa) + 4*QQ*(1 + wb) - 2], [PP, QQ, CC], dict=True)
sa = einstein_fluid(0, sp.Rational(1, 4)); sr = einstein_fluid(sp.Rational(1, 3), 0); si = einstein_fluid(sp.Rational(1, 7), sp.Rational(1, 7))
R['Einstein + fluid power laws a ~ tau^P, b ~ tau^Q: (wa, wb) = (0, 1/4) -> (P, Q) = (0, 2/5) [decompactification]; (1/3, 0) -> (1/2, 0) [brane radiation, b frozen]; (1/7, 1/7) -> (1/4, 1/4) [isotropic]'] = (
    sa == [{PP: 0, QQ: sp.Rational(2, 5), CC: sp.Rational(24, 25)}] and sr == [{PP: sp.Rational(1, 2), QQ: 0, CC: sp.Rational(3, 4)}]
    and si == [{PP: sp.Rational(1, 4), QQ: sp.Rational(1, 4), CC: sp.Rational(21, 16)}])
# matter-dominated cubic regime for a general fluid
pw, hh, Cw, wa_ = sp.symbols('p h C_w w_a', real=True)
lead00 = sp.expand(sp.limit(TOT[0].subs({Ha: pw/tt, Hb: hh, dHa: -pw/tt**2, dHb: 0, Lb: 0})*tt**3, tt, 0))
lead11 = sp.expand(sp.limit(TOT[1].subs({Ha: pw/tt, Hb: hh, dHa: -pw/tt**2, dHb: 0, Lb: 0})*tt**3, tt, 0))
lead44 = sp.expand(sp.limit(TOT[4].subs({Ha: pw/tt, Hb: hh, dHa: -pw/tt**2, dHb: 0, Lb: 0})*tt**4, tt, 0))
roots_p = [r for r in sp.solve(sp.expand(lead11 + wa_*lead00).subs(hh, -sp.sqrt(al2/(54*al3))), pw) if r != 0]
R['matter-dominated cubic regime (any fluid pa = w_a eps, b finite): leading 44 equation p^3 (p - 1)(648 alpha3 h^2 - 12 alpha2) = 0 -> h^2 = alpha2/(54 alpha3); 11/00 -> p = 1/(1 + w_a); kappa eps tau^3 = (64/3) p^3 alpha2 |h|'] = (
    sp.expand(lead44 - pw**3*(pw - 1)*(648*al3*hh**2 - 12*al2)) == 0 and len(roots_p) == 1 and sp.simplify(roots_p[0] - 1/(1 + wa_)) == 0 and
    sp.simplify(lead00.subs(hh, -sp.sqrt(al2/(54*al3))) - sp.Rational(64, 3)*pw**3*al2*sp.sqrt(al2/(54*al3))) == 0)
R['NUM h* = -sqrt(alpha2/(54 alpha3)), h0 = -sqrt(alpha2/(30 alpha3)) (xi = 1e-6)'] = (hstar, cub['h0'])

# ----------------------------------------------------------------------- B4. pair creation in the Einstein-Kasner regime: the power-law tail
s1m = np.array([[0, 1], [1, 0]], complex); s2m = np.array([[0, -1j], [1j, 0]])
def bog_kasner(K, p, q, t0=1e-15, t1=300.0):
    def rhs(s, y):
        tq = math.exp(s); cvec = y[:2] + 1j*y[2:]
        dv = -1j*tq*((K*tq**(-p))*(s1m @ cvec) + (K*tq**(-q))*(s2m @ cvec))
        return np.concatenate([dv.real, dv.imag])
    w0, v0 = np.linalg.eigh(K*t0**(-p)*s1m + K*t0**(-q)*s2m); c0 = v0[:, 1]
    sol = solve_ivp(rhs, (math.log(t0), math.log(t1)), np.concatenate([c0.real, c0.imag]), method='DOP853', rtol=1e-11, atol=1e-13)
    cvec = sol.y[:2, -1] + 1j*sol.y[2:, -1]
    w1, v1 = np.linalg.eigh(K*t1**(-p)*s1m + K*t1**(-q)*s2m)
    return float(abs(v1[:, 0].conj() @ cvec)**2)
pk, qk = float(p1), float(q1)
Fk = {K: bog_kasner(K, pk, qk) for K in (0.3, 1.0, 3.0, 8.0, 16.0)}
tail = math.log(Fk[16.0]/Fk[8.0])/math.log(2.0)
R['NUM Einstein-Kasner pair creation |beta|^2 = F(K) at K = 0.3, 1, 3, 8, 16, and the tail exponent d ln F/d ln K between 8 and 16'] = (tuple(Fk[K] for K in sorted(Fk)), tail)
R['pair creation has a power-law tail F ~ K^(-2(p1 - q1)/(1 - p1)) = K^(-(3 sqrt 2 + 2)/2) = K^(-3.121) (to 2e-2): it comes from the singularity, so the number of created pairs is dominated by the earliest times'] = (
    abs(tail + (3*math.sqrt(2) + 2)/2) < 2e-2 and abs(2*(pk - qk)/(1 - pk) - (3*math.sqrt(2) + 2)/2) < 1e-12)
R['F(1) agrees with the value of section 14.3 of 7+1.nb (|beta|^2 for (|p|, |q|) = (1, 1))'] = abs(Fk[1.0] - 0.010749350) < 1e-6

# ----------------------------------------------------------------------- B5. the Standard Model on the brane is negligible in regime I
M8 = 1.0e4; MPl = 2.435e18; gstar = 106.75
frac_brane = (math.pi**2*gstar/30)*M8**4/(3*MPl**2*M8**2)
R['NUM brane radiation fraction kappa eps_W/(3 H^2) at T = M8, H = M8 (M8 = 10 TeV, reduced M_Pl = 2.435e18 GeV, g* = 106.75)'] = frac_brane
R['the Standard Model plasma on the brane (T <= M8) is gravitationally negligible in the bulk while H ~ M8: fraction < 1e-27'] = frac_brane < 1e-27


# ======================================================================= C. regime II begins: the stabilizing potential catches the auxiliary space
from scipy.optimize import least_squares
MR = 1e-3                                                         # radion mass m_r in units 1/l2 (m_r^2 = kappa K/12)
st0 = coupled_at(1e-2)                                            # the coupled solution of section B2 in the Einstein regime
tau_c0, la_c0, lb_c0 = st0['tau'], st0['la'], st0['lb']
lnb0 = lb_c0 + float(q1)*math.log(1.0/(MR*tau_c0))               # b0: where the Kasner contraction would arrive at tau = 1/m_r
def rhs_cap(x, y):
    ha, hb, la, lb = y
    e, pa, pb = gas(la, lb)
    bet = lb - lnb0
    kU = 6*MR**2*bet**2                                           # kappa U = (kappa K/2) beta^2 = 6 m_r^2 beta^2
    Am = np.array(f_A(ha, hb, 1.0, XI), dtype=float)
    Bm = np.array(f_B(ha, hb, 1.0, XI, 0.0), dtype=float).reshape(2) + np.array([Cg*pa - kU, Cg*pb - kU - 3*MR**2*bet])
    d = np.linalg.solve(Am, -Bm)
    return [d[0], d[1], ha, hb]
bet0 = lb_c0 - lnb0
e0c = gas(la_c0, lb_c0)[0]
hb_c0 = brentq(lambda hb: f_C(st0['Ha'], hb, 1.0, XI, 0.0) - Cg*e0c - 6*MR**2*bet0**2, st0['Hb']*1.5, st0['Hb']*0.5)
cap = solve_ivp(rhs_cap, [tau_c0, 3e5], [st0['Ha'], hb_c0, la_c0, lb_c0], method='DOP853', rtol=1e-11, atol=1e-16, dense_output=True, max_step=50.0)
R['NUM capture: start tau, beta(start) = ln(b/b0), m_r'] = (tau_c0, bet0, MR)
taus_f = np.linspace(1e5, 3e5, 4001)
bet_f = np.array([cap.sol(x)[3] - lnb0 for x in taus_f])
fitc = least_squares(lambda pv: pv[0]*np.sin(MR*taus_f + pv[1])/(MR*(taus_f - pv[2])) - bet_f, [-0.46, -1.0, 0.0])
A_c, ph_c, tE_c = fitc.x
res_c = np.max(np.abs(fitc.fun))/np.max(np.abs(bet_f))
HEavg = []
for x in (1e5, 2e5):
    xs_ = np.linspace(x, x + 2*np.pi/MR, 2001)
    HEavg.append(np.mean([cap.sol(z)[0] + 2*cap.sol(z)[1] for z in xs_])*(np.mean(xs_) - tE_c))
gas_frac_late = [Cg*gas(cap.sol(x)[2], cap.sol(x)[3])[0]/(3*(cap.sol(x)[0] + 2*cap.sol(x)[1])**2) for x in (1e5, 3e5)]
R['NUM after capture: beta = A sin(m_r tau + phi)/(m_r (tau - tau_E)) fitted on 1e5 <= tau <= 3e5: A, phi, tau_E, max residual/amplitude'] = (A_c, ph_c, tE_c, res_c)
R['NUM after capture: <H_E>(tau - tau_E) over one period at tau = 1e5, 2e5 (H_E = Ha + 2Hb, four-dimensional matter domination: 2/3)'] = tuple(HEavg)
R['NUM after capture: gas fraction kappa eps/(3 H_E^2) at tau = 1e5 and 3e5 (constant: both are four-dimensional matter)'] = tuple(gas_frac_late)
R['capture: the radion oscillates as a four-dimensional massive field in a matter-dominated universe: j0 form to 2e-2, <H_E>(tau - tau_E) -> 2/3 to 1e-2, the gas fraction stays constant to 2e-2'] = (
    res_c < 2e-2 and all(abs(v - 2/3) < 1e-2 for v in HEavg) and abs(gas_frac_late[1]/gas_frac_late[0] - 1) < 2e-2)
R['NUM after capture: radion share of the energy, 1 - gas fraction (the cosmological moduli problem)'] = 1 - gas_frac_late[1]
# exact identity behind it: 3 (Ha + 2Hb)^2 = 6 Hb^2 + kappa (U + eps): Einstein frame, phi = sqrt(12) M_Pl beta
R["four-dimensional (Einstein-frame) form of the 00 equation: 3 (Ha + 2Hb)^2 = 6 Hb^2 + kappa eps_total, i.e. 3 H_E^2 M_Pl^2 = (1/2) phi'^2 + rho with phi = sqrt(12) M_Pl ln(b/b0)"] = (
    sp.expand(3*(Ha + 2*Hb)**2 - 6*Hb**2 - CF[(1, 0)]) == 0)

# radion decay (physical units): coupling L = (phi/(sqrt 3 M_Pl)) T^mu_mu of the brane fields
MPl = 2.435e18; M8 = 1.0e4
Vperp = MPl**2/M8**6
hbar_s = 6.582119569e-25
def gam_ff(mphi, mf, Nc):
    x = 4*mf**2/mphi**2
    return Nc*mf**2*mphi*(1 - x)**1.5/(24*math.pi*MPl**2) if x < 1 else 0.0
def gam_gg(mphi, alpha_s, nf):
    b0q = 11 - 2*nf/3
    return b0q**2*alpha_s**2*mphi**3/(96*math.pi**3*MPl**2)
mphi = 10.0
Gam_r = gam_ff(mphi, 4.183, 3) + gam_ff(mphi, 1.273, 3) + gam_ff(mphi, 1.77693, 1) + gam_ff(mphi, 0.1056583755, 1) + gam_gg(mphi, 0.18, 5)
gst_RH = 10.75
T_RH = (Gam_r*MPl/math.sqrt(math.pi**2*gst_RH/90))**0.5
R['NUM radion of mass 10 GeV (M8 = 10 TeV): decay rate into brane fields (b, c, tau, mu pairs and two gluons) in GeV, lifetime in s, reheating temperature in GeV'] = (Gam_r, hbar_s/Gam_r, T_RH)
R['a radion of 10 GeV decays far too late (T_RH < 1 keV < 4 MeV needed by nucleosynthesis): the moduli problem of the minimal model'] = T_RH < 1e-6
# Kaluza-Klein states of dirac16complex below m_r/2 (bulk decay channels): N = 8 V_perp (pi^2/2) (m_r/2)^4/(2 pi)^4 per spin state ...
N_KK = Vperp*(math.pi**2/2)*(mphi/2)**4/(2*math.pi)**4
R['NUM number of Kaluza-Klein momenta of the auxiliary torus with |q| < m_r/2 (m_r = 10 GeV): V_perp (pi^2/2)(m_r/2)^4/(2 pi)^4'] = N_KK

# ======================================================================= D. regime II: the modern universe on the brane (physical units)
GeV_per_K = 8.617333262e-14
lperp_m = Vperp**0.25*1.973269804e-16
mKK1 = 2*math.pi/Vperp**0.25
R['NUM M8, M_Pl (reduced), V_perp = M_Pl^2/M8^6 (GeV^-4), l_perp = V_perp^(1/4) (m), lowest Kaluza-Klein mass 2 pi/l_perp (GeV)'] = (M8, MPl, Vperp, lperp_m, mKK1)
# relativistic degrees of freedom of the Standard Model plasma (ideal gases; QCD crossover modelled by a tanh switch)
SPECIES = [('photon', 2, 0.0, -1, 'always'), ('W', 6, 80.3692, -1, 'always'), ('Z', 3, 91.1880, -1, 'always'), ('h', 1, 125.20, -1, 'always'),
           ('e', 4, 0.51099895e-3, 1, 'always'), ('mu', 4, 0.1056583755, 1, 'always'), ('tau', 4, 1.77693, 1, 'always'), ('nu', 6, 0.0, 1, 'nu'),
           ('gluon', 16, 0.0, -1, 'qgp'), ('u', 12, 2.16e-3, 1, 'qgp'), ('d', 12, 4.70e-3, 1, 'qgp'), ('s', 12, 93.5e-3, 1, 'qgp'),
           ('c', 12, 1.2730, 1, 'qgp'), ('b', 12, 4.183, 1, 'qgp'), ('t', 12, 172.56, 1, 'qgp'),
           ('pi+-', 2, 0.13957039, -1, 'hadron'), ('pi0', 1, 0.1349768, -1, 'hadron'), ('K+-', 2, 0.493677, -1, 'hadron'), ('K0', 2, 0.497611, -1, 'hadron'),
           ('eta', 1, 0.547862, -1, 'hadron'), ('rho', 9, 0.77526, -1, 'hadron'), ('omega', 3, 0.78266, -1, 'hadron'), ('K*', 12, 0.8955, -1, 'hadron'),
           ('eta prime', 1, 0.95778, -1, 'hadron'), ('p', 4, 0.93827209, 1, 'hadron'), ('n', 4, 0.93956542, 1, 'hadron')]
from scipy.integrate import quad
def rho_s_sp(T, g, m, stat):
    x = m/T
    if x > 60: return 0.0, 0.0
    if stat == -1 and x == 0:
        r, p = math.pi**4/15, math.pi**4/45
    else:
        def fr(z):
            E = math.sqrt(z*z + x*x); return z*z*E/(math.exp(E) + stat) if E < 700 else 0.0
        def fp(z):
            E = math.sqrt(z*z + x*x); return z**4/E/(math.exp(E) + stat)/3 if E < 700 else 0.0
        r = quad(fr, 0, 80, limit=200)[0]; p = quad(fp, 0, 80, limit=200)[0]
    rho = g*T**4*r/(2*math.pi**2); pres = g*T**4*p/(2*math.pi**2)
    return rho, (rho + pres)/T
TQCD, DTQCD, NEFF = 0.160, 0.030, 3.044
def Sq(T): return 0.5*(1 + math.tanh((T - TQCD)/DTQCD))
def g_gamma_e(T):
    return (rho_s_sp(T, 2, 0.0, -1)[1] + rho_s_sp(T, 4, 0.51099895e-3, 1)[1])*45/(2*math.pi**2*T**3)
def gstars(T):
    rho = s = 0.0
    for name, g, m, stv, grp in SPECIES:
        if grp == 'nu':
            Tn = T*min(1.0, (g_gamma_e(T)/5.5)**(1/3))      # neutrinos decouple while e+- are relativistic; entropy of the gamma-e plasma
            r, ss = rho_s_sp(Tn, g, 0.0, 1); r *= NEFF/3; ss *= NEFF/3
        else:
            r, ss = rho_s_sp(T, g, m, stv)
            if grp == 'qgp': r *= Sq(T); ss *= Sq(T)
            if grp == 'hadron': r *= 1 - Sq(T); ss *= 1 - Sq(T)
        rho += r; s += ss
    return rho*30/(math.pi**2*T**4), s*45/(2*math.pi**2*T**3)
gtab = {T: gstars(T) for T in (1e3, 10.0, 1.0, 0.3, 0.1, 0.01, 1e-3, 1e-4, 1e-6)}
R['NUM g*(T), g*s(T) at T = 1 TeV, 10 GeV, 1 GeV, 300 MeV, 100 MeV, 10 MeV, 1 MeV, 100 keV, 1 keV'] = tuple((round(v[0], 4), round(v[1], 4)) for v in gtab.values())
R['g* limits: 106.75 above the electroweak scale (to 0.1), 10.75 near 3 MeV (to 0.1), and after e+e- annihilation 2 + (21/4)(4/11)^(4/3)(N_eff/3) = 3.3855 (to 1e-3)'] = (
    abs(gtab[1e3][0] - 106.75) < 0.1 and abs(gstars(3e-3)[0] - 10.75) < 0.1 and abs(gtab[1e-6][0] - (2 + 21/4*(4/11)**(4/3)*NEFF/3)) < 1e-3)
# time-temperature relation in the radiation era: t(T) = Int_T^Tmax (1 + (1/3) dln g*s/dln T) dT'/(T' H)
lnT = np.linspace(math.log(10.0), math.log(1e-5), 241)
G_ = np.array([gstars(math.exp(x)) for x in lnT])
from scipy.interpolate import CubicSpline
spl_g = CubicSpline(-lnT, np.log(G_[:, 0])); spl_gs = CubicSpline(-lnT, np.log(G_[:, 1]))
def H_rad(T):
    return math.sqrt(math.pi**2*math.exp(spl_g(-math.log(T)))/90)*T**2/MPl
def dt_dlnT(lT):
    T = math.exp(lT); dlgs = spl_gs(-lT, 1)            # d ln g*s / d(-ln T)
    return (1 - dlgs/3)/H_rad(T)
T_start = 10.0
t0_rad = 1/(2*H_rad(T_start))                            # t at T = 10 GeV (radiation domination with g*(10 GeV) before)
t_of_T = {}
for Tq in (1.0, 0.1, 1e-2, 1e-3, 1e-4, 1e-5):
    t_of_T[Tq] = t0_rad + quad(lambda lT: dt_dlnT(lT), math.log(Tq), math.log(T_start), limit=400)[0]
t_inst = {Tq: 1/(2*H_rad(Tq)) for Tq in t_of_T}
R['NUM radiation era: t(T) in s (entropy conservation with g*(T), g*s(T)) at T = 1 GeV, 100 MeV, 10 MeV, 1 MeV, 100 keV, 10 keV'] = tuple(t_of_T[T]*hbar_s for T in sorted(t_of_T, reverse=True))
R['NUM radiation era: ratio of the instantaneous approximation t = 1/(2H(T)) to t(T)'] = tuple(t_inst[T]/t_of_T[T] for T in sorted(t_of_T, reverse=True))
R['radiation era: t(1 MeV) = 0.74 s within 5% (standard value 1/(2H) with g* = 10.75)'] = abs(t_of_T[1e-3]*hbar_s/0.74 - 1) < 0.05
# the late universe: Planck 2018 (TT,TE,EE+lowE+lensing): H0 = 67.4 km/s/Mpc, Omega_m = 0.315; radiation from T_CMB = 2.7255 K and N_eff = 3.044
mpc_km = 3.0856775814913673e19; gyr_s = 3.15576e16
H0s = 67.4/mpc_km; H0g = H0s*hbar_s
rho_crit = 3*H0g**2*MPl**2
Tcmb = 2.7255*GeV_per_K
Om_g = (math.pi**2/15)*Tcmb**4/rho_crit
Om_r = Om_g*(1 + 7/8*(4/11)**(4/3)*NEFF)
Om_m = 0.315; Om_L = 1 - Om_m - Om_r
def Hof(av): return H0s*math.sqrt(Om_r/av**4 + Om_m/av**3 + Om_L)
age = quad(lambda av: 1/(av*Hof(av)), 0, 1, limit=400, epsabs=0, epsrel=1e-12)[0]/gyr_s
a_eq = Om_r/Om_m
R['NUM late universe: Omega_r, Omega_Lambda, a_eq = Omega_r/Omega_m, z_eq, age (Gyr)'] = (Om_r, Om_L, a_eq, 1/a_eq - 1, age)
R['age of the universe 13.79 Gyr (Planck 2018: 13.787 +- 0.020)'] = abs(age - 13.79) < 0.02
def t_exact(av): return quad(lambda x: 1/(x*Hof(x)), 0, av, limit=400, epsabs=0, epsrel=1e-12)[0]
def t_rm(av):       # radiation + matter, exact
    return 2/(3*H0s*Om_m**2)*((Om_m*av - 2*Om_r)*math.sqrt(Om_r + Om_m*av) + 2*Om_r**1.5)
def t_mL(av):       # matter + Lambda, exact
    return 2/(3*H0s*math.sqrt(Om_L))*math.asinh(math.sqrt(Om_L/Om_m)*av**1.5)
def t_m(av): return 2/(3*H0s*math.sqrt(Om_m))*av**1.5
def t_comp(av): return t_rm(av) + t_mL(av) - t_m(av)
errs = [abs(t_comp(av)/t_exact(av) - 1) for av in np.logspace(-8, 0.5, 60)]
R['NUM composite solution t(a) = t_rm(a) + t_mL(a) - t_m(a): maximal relative error for 1e-8 <= a <= 3.16'] = max(errs)
R['the composite analytic solution (exact radiation+matter, exact matter+Lambda, minus their common matter part) is accurate to < 3e-4 at all times'] = max(errs) < 3e-4
# the radion today: beta_eq = [Lambda4 + kappa f/V_perp + (kappa/2)(eps - 3pa) + kappa pb - 12 alpha2 H^2 (H' + H^2)]/(3 m_r^2)
Lam4 = 3*H0g**2*Om_L
beta_mat = (1.5*H0g**2*Om_m)/(3*mphi**2)
beta_f = (M8**4/MPl**2)/(3*mphi**2)
R['NUM radion displacement today (m_r = 10 GeV): from Lambda4 + matter, from a brane tension f_W = M8^4'] = ((Lam4 + 1.5*H0g**2*Om_m)/(3*mphi**2), beta_f)
R['NUM Gauss-Bonnet and cubic corrections today: alpha2 H0^2 and alpha3 H0^4 with alpha2 = 1/M8^2, alpha3 = 1e-6/M8^4'] = (H0g**2/M8**2, 1e-6*H0g**4/M8**4)
# normalcy temperature (bulk graviton emission ~ T^(n+7)/M8^(n+2) = T^11/M8^6 against Hubble cooling): T* ~ (M8^6/M_Pl)^(1/5)
Tstar = (M8**6/MPl)**0.2
R['NUM normalcy temperature (order of magnitude) T* = (M8^6/M_Pl)^(1/5) in GeV'] = Tstar

# ======================================================================= E. regime III: the far future
Hinf = H0s*math.sqrt(Om_L)
R['NUM far future: H_inf = H0 sqrt(Omega_Lambda) (1/s and GeV), Gibbons-Hawking temperature H/(2 pi) in K'] = (Hinf, Hinf*hbar_s, Hinf*hbar_s/(2*math.pi)/GeV_per_K)
# exact matter + Lambda: H(t) = H_inf coth(3 H_inf t/2); check against the numerical solution at a = 3 and a = 10
def t_exact_late(av):
    return quad(lambda x: 1/(x*Hof(x)), 0, av, limit=400, epsabs=0, epsrel=1e-12)[0]
okH = True
for av in (3.0, 10.0):
    tv = t_exact_late(av); okH = okH and abs(Hinf/math.tanh(1.5*Hinf*tv)/Hof(av) - 1) < 1e-4
R['far future: H(t) = H_inf coth(3 H_inf t/2) (exact for matter + Lambda) agrees with the full solution at a = 3, 10 to 1e-4'] = okH
# GB-corrected radion in de Sitter: linearize the full equations around Ha = H, Hb = 0, beta = beta_eq (U = (K/2) beta^2, kappa K = 12 m2)
Hs, bq, m2 = sp.symbols('H beta_q m2', positive=True)
HaV, HbV, beV = sp.symbols('Ha_v Hb_v beta_v', real=True)
def kU_(bb): return 6*m2*bb**2
eqC = TOT[0].subs({Ha: HaV, Hb: HbV}) - kU_(beV)                       # constraint: Sum alpha_k E^(k)0_0 - Lambda = kappa U
s_l = sp.solve([TOT[1].subs({Ha: HaV, Hb: HbV}) - kU_(beV), TOT[4].subs({Ha: HaV, Hb: HbV}) + (-kU_(beV) - 3*m2*beV)], [dHa, dHb], dict=True)[0]
Fb = s_l[dHb]                                                            # Hb' as a function of (Ha, Hb, beta)
Lb_bg = sp.solve(eqC.subs({HaV: Hs, HbV: 0, beV: bq}), Lb)[0]            # background: Lambda from the constraint
bg = {HaV: Hs, HbV: 0, beV: bq, Lb: Lb_bg}
Cd = {v: sp.diff(eqC, v).subs(bg) for v in (HaV, HbV, beV)}
Fd = {v: sp.diff(Fb, v).subs(bg) for v in (HaV, HbV, beV)}
gamma_r = sp.simplify(-(Fd[HbV] - Fd[HaV]*Cd[HbV]/Cd[HaV]))
meff2_r = sp.simplify(-(Fd[beV] - Fd[HaV]*Cd[beV]/Cd[HaV]))
eq_bg = sp.simplify(Fb.subs(bg))                                         # must vanish at equilibrium
R['far future: the radion obeys beta\'\' + gamma beta\' + m_eff^2 beta = 0 around de Sitter (constraint used for delta Ha); alpha2 = alpha3 = 0: gamma = 3H, m_eff^2 = m_r^2 at beta_eq -> 0'] = (
    sp.simplify(gamma_r.subs({al2: 0, al3: 0, bq: 0}) - 3*Hs) == 0 and sp.simplify(meff2_r.subs({al2: 0, al3: 0, bq: 0}) - m2) == 0)
R['NUM far future (beta_eq -> 0): gamma/(3H) and m_eff^2/m_r^2 with Gauss-Bonnet and cubic terms (sympy)'] = (
    str(sp.factor(sp.simplify(gamma_r.subs({bq: 0})/(3*Hs)))), str(sp.factor(sp.simplify(meff2_r.subs({bq: 0})/m2))))
R['far future: equilibrium of the radion in de Sitter: 3 m_r^2 beta_eq = 3H^2 - 12 alpha2 H^4 (+ O(beta_eq^2))'] = sp.simplify(
    sp.series(eq_bg, bq, 0, 2).removeO() - sp.series((3*Hs**2 - 12*al2*Hs**4 - 3*m2*bq)*sp.Rational(1, 3)*0 + eq_bg, bq, 0, 2).removeO()) == 0 and \
    sp.simplify(sp.solve(sp.series(eq_bg, bq, 0, 2).removeO(), bq)[0] - (3*Hs**2 - 12*al2*Hs**4)/(3*m2)) == 0
# Gibbons-Hawking pair creation of dirac16complex in the far future: exponent 2 pi M/H_inf for M = 1 eV, 1 GeV
R['NUM far future: 2 pi M/H_inf for M = 1 eV and 1 GeV (|beta|^2 = 1/(exp(2 pi M/H) + 1))'] = (2*math.pi*1e-9/(Hinf*hbar_s), 2*math.pi*1.0/(Hinf*hbar_s))
# eight-dimensional de Sitter without stabilization: 21 H^2 - 420 alpha2 H^4 + 2520 alpha3 H^6 = Lambda
H8 = sp.Symbol('H8', positive=True)
ds8 = sp.expand(TOT[0].subs({Ha: H8, Hb: H8, dHa: 0, dHb: 0}))
R['eight-dimensional de Sitter (no stabilizer): 21 H^2 - 420 alpha2 H^4 + 2520 alpha3 H^6 = Lambda (all three equations)'] = (
    sp.expand(ds8 - (21*H8**2 - 420*al2*H8**4 + 2520*al3*H8**6 - Lb)) == 0 and
    all(sp.expand(TOT[c].subs({Ha: H8, Hb: H8, dHa: 0, dHb: 0}) - ds8) == 0 for c in (1, 4)))
# numerical check of the radion in a de Sitter universe (units l2 = 1): H = 0.01, m_r = 0.1, alpha3 = 1e-6, start beta = beta_eq + 0.05
Hd, mrd = 0.01, 0.1
beta_eq_d = Hd**2*(1 - 4*Hd**2)/mrd**2                   # (kappa/4) b U' = Lambda - 12 alpha2 H^4 and Lambda = 3H^2 (+ O(beta^2))
Lbd = 3*Hd**2 + 6*mrd**2*0.0                             # kappa U0 = 0
def rhs_ds(x, y):
    ha, hb, la, lb = y
    bet = lb
    kU = 6*mrd**2*bet**2
    Am = np.array(f_A(ha, hb, 1.0, XI), dtype=float)
    Bm = np.array(f_B(ha, hb, 1.0, XI, Lbd), dtype=float).reshape(2) + np.array([-kU, -kU - 3*mrd**2*bet])
    d = np.linalg.solve(Am, -Bm)
    return [d[0], d[1], ha, hb]
b_start = beta_eq_d + 0.05
ha_start = brentq(lambda h: f_C(h, 0.0, 1.0, XI, Lbd) - 6*mrd**2*b_start**2, 1e-4, 1.0)
dsol = solve_ivp(rhs_ds, [0.0, 3000.0], [ha_start, 0.0, 0.0, b_start], method='DOP853', rtol=1e-11, atol=1e-14, dense_output=True, max_step=1.0)
bet_end = dsol.y[3, -1]; ha_end = dsol.y[0, -1]
bb_, hh_ = sp.symbols('bb hh', positive=True)
eqs_eq = [3*hh_**2 - (Lbd + 6*mrd**2*bb_**2), 3*mrd**2*bb_ - (3*hh_**2 - 12*hh_**4)]      # constraint and equilibrium (alpha2 = 1)
sol_eq = sp.nsolve(eqs_eq, [bb_, hh_], [beta_eq_d, Hd])
R['NUM de Sitter radion test (l2 = 1, Lambda = 3e-4, m_r = 0.1, start 0.05 above): beta(3000), Ha(3000), and the equilibrium (beta_eq, H_eq) of the two algebraic equations'] = (
    bet_end, ha_end, float(sol_eq[0]), float(sol_eq[1]))
R['far future: the radion relaxes to its shifted equilibrium (to 1e-6 relative) and Ha to H_eq (to 1e-8)'] = (
    abs(bet_end/float(sol_eq[0]) - 1) < 1e-6 and abs(ha_end/float(sol_eq[1]) - 1) < 1e-8)


if __name__ == '__main__':
    bad = [k for k, v in R.items() if isinstance(v, (bool, np.bool_)) and not v]
    nbool = sum(isinstance(v, (bool, np.bool_)) for v in R.values())
    for k, v in R.items():
        print(f'{k} : {v}')
    print(f'boolean checks: {nbool}   failed: {len(bad)}')
    print('VERDICT:', 'ALL CHECKS PASSED' if not bad else 'FAILED: ' + '; '.join(bad))
    sys.exit(1 if bad else 0)
