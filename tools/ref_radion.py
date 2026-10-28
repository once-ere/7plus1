"""Independent checks of section 26 of 7+1.nb (the radion: the size of the four extra dimensions, on the 4-torus and on the 4-sphere).

Shares no code with the notebook: its own curvature routine for diagonal metrics (sympy), its own Euler-Lagrange operator, its own lattice
sieve (numpy), and high-precision numerical checks (mpmath) where the notebook is exact.  Prints one line per result, 'name : value'; the
lines whose names begin with NUM are compared with the numbers printed by section 26 (tools/verify_7p1_nb.py); the others are boolean
checks.  Exit code 0 if every boolean check is True.

Curvature conventions of 7+1.nb (mostly minus): R^r_smn = d_m G^r_ns - d_n G^r_ms + G^r_ml G^l_ns - G^r_nl G^l_ms, R_sn = R^r_srn,
R = g^sn R_sn, G^A_B = R^A_B - delta R/2; field equations G^A_B - Lambda delta^A_B = kappa T^A_B."""
import math, sys, random
import numpy as np
import sympy as sp
import mpmath as mp

R = {}
I = sp.I


# ---------------------------------------------------------------------------------------------------- curvature of diagonal metrics
def christoffel(gd, xs):
    n = len(xs)
    G = [[[sp.S(0)]*n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                if a == b == c:
                    G[a][b][c] = sp.diff(gd[a], xs[a])/(2*gd[a])
                elif a == b:
                    G[a][b][c] = sp.diff(gd[a], xs[c])/(2*gd[a])
                elif a == c:
                    G[a][b][c] = sp.diff(gd[a], xs[b])/(2*gd[a])
                elif b == c:
                    G[a][b][c] = -sp.diff(gd[b], xs[a])/(2*gd[a])
    return G


def ricci(gd, xs):
    n = len(xs); G = christoffel(gd, xs)
    return [[sp.simplify(sum(sp.diff(G[r][m][s], xs[r]) - sp.diff(G[r][r][s], xs[m])
                             + sum(G[r][r][l]*G[l][m][s] - G[r][m][l]*G[l][r][s] for l in range(n)) for r in range(n)))
             for m in range(n)] for s in range(n)]


def scalar(gd, xs):
    Ric = ricci(gd, xs)
    return sp.simplify(sum(Ric[i][i]/gd[i] for i in range(len(xs))))


def einstein(gd, xs):
    Ric = ricci(gd, xs)
    Rs = sp.simplify(sum(Ric[i][i]/gd[i] for i in range(len(xs))))
    return [sp.simplify(Ric[i][i]/gd[i] - Rs/2) for i in range(len(xs))]


def EL(L, q, t):
    """Euler-Lagrange expression of L(q, q', q'') for the function q(t)"""
    q0, q1, q2 = sp.symbols('_q0 _q1 _q2')
    Ls = L.subs(sp.Derivative(q, (t, 2)), q2).subs(sp.Derivative(q, t), q1).subs(q, q0)
    back = {q0: q, q1: sp.diff(q, t), q2: sp.diff(q, t, 2)}
    return sp.simplify(sp.diff(Ls, q0).subs(back) - sp.diff(sp.diff(Ls, q1).subs(back), t) + sp.diff(sp.diff(Ls, q2).subs(back), t, 2))


t = sp.Symbol('x0'); X = [t] + list(sp.symbols('x1:8'))
x4, x5, x6, x7 = X[4:]
kap, Lam, r0, b0, Lc = sp.symbols('kappa Lambda r_0 b_0 L_c', positive=True)
N, a, B, ph = (sp.Function(f)(t) for f in ('N', 'a', 'B', 'phi'))
U = sp.Function('U')
gam = [1, sp.sin(x4)**2, sp.sin(x4)**2*sp.sin(x5)**2, sp.sin(x4)**2*sp.sin(x5)**2*sp.sin(x6)**2]

# ==================================================================================================== A. definition and reduction
Ric0 = ricci([1, -1, -1, -1] + [-r0**2*g for g in gam], X)
R['the 4-sphere of radius r0: Ricci_ij = 3 gamma_ij, off-diagonal 0, R = -12/r0^2 in the signature of 7+1.nb; volume of the unit S^4 = 8 pi^2/3'] = (
    all(sp.simplify(Ric0[i][i] - 3*gam[i - 4]) == 0 for i in range(4, 8))
    and all(sp.simplify(Ric0[i][j]) == 0 for i in range(8) for j in range(8) if i != j)
    and sp.simplify(scalar([1, -1, -1, -1] + [-r0**2*g for g in gam], X) + 12/r0**2) == 0
    and sp.simplify(sp.integrate(sp.sin(x4)**3*sp.sin(x5)**2*sp.sin(x6), (x4, 0, sp.pi), (x5, 0, sp.pi), (x6, 0, sp.pi), (x7, 0, 2*sp.pi)) - sp.Rational(8, 3)*sp.pi**2) == 0)


def reduction_ok(gint, rint, V0):
    gd = [N**2*sp.exp(-4*B), -a**2*sp.exp(-4*B), -a**2*sp.exp(-4*B), -a**2*sp.exp(-4*B)] + [-sp.exp(2*B)*g for g in gint]
    Rs = scalar(gd, X); MP2 = V0/kap
    LJ = -(1/(2*kap))*N*a**3*sp.exp(-4*B)*V0*(Rs + 2*Lam) - N*a**3*sp.exp(-4*B)*V0*U(B)
    RE = scalar([N**2, -a**2, -a**2, -a**2], X[:4])
    V = sp.exp(-4*B)*(MP2*Lam + V0*U(B)) - MP2/2*sp.exp(-4*B)*rint
    LE = -(MP2/2)*N*a**3*RE + 6*MP2*a**3*sp.diff(B, t)**2/N - N*a**3*V
    D_ = sp.simplify(LJ - LE)
    return not any(Rs.has(x) for x in X[4:]) and all(EL(D_, q, t) == 0 for q in (N, a, B))


R['radion on T^4: the reduced action equals -(M_Pl^2/2) R_E + 6 M_Pl^2 (d beta)^2 - e^(-4 beta)(M_Pl^2 Lambda + V0 U) up to a total derivative (Euler-Lagrange for N, a, beta), M_Pl^2 = V0/kappa'] = (
    reduction_ok([b0**2]*4, 0, b0**4*Lc**4))
R['radion on S^4: the same with the potential -(M_Pl^2/2) e^(-4 beta) 12/(r0^2 e^(2 beta)) from the curvature of the sphere, V0 = (8 pi^2/3) r0^4'] = (
    reduction_ok([r0**2*g for g in gam], 12/(r0**2*sp.exp(2*B)), r0**4*sp.Rational(8, 3)*sp.pi**2))


def kinetic_ok(n, beta, coef):
    xs = [t] + list(sp.symbols('y1:%d' % (4 + n)))
    gd = [N**2*sp.exp(-n*beta*ph), -a**2*sp.exp(-n*beta*ph), -a**2*sp.exp(-n*beta*ph), -a**2*sp.exp(-n*beta*ph)] + [-sp.exp(2*beta*ph)]*n
    LJ = -(1/(2*kap))*N*a**3*sp.exp(-n*beta*ph)*scalar(gd, xs)
    RE = scalar([N**2, -a**2, -a**2, -a**2], xs[:4])
    LE = -(1/(2*kap))*N*a**3*RE + (1/(2*kap))*coef*a**3*sp.diff(ph, t)**2/N
    D_ = sp.simplify(LJ - LE)
    return all(EL(D_, q, t) == 0 for q in (N, a, ph))


R["n extra dimensions (T^n, n = 1..4): kinetic term (M_Pl^2/2)(n(n+2)/2)(dB)^2; the note's beta_n = sqrt(2/(n(n+2))) gives (M_Pl^2/2)(d phi)^2, not (M_Pl^2/2)(1/2)(d phi)^2; beta_n = 1/sqrt(n(n+2)) gives the 1/2"] = all(
    kinetic_ok(n, 1, sp.Rational(n*(n + 2), 2)) and kinetic_ok(n, sp.sqrt(sp.Rational(2, n*(n + 2))), 1)
    and not kinetic_ok(n, sp.sqrt(sp.Rational(2, n*(n + 2))), sp.Rational(1, 2)) and kinetic_ok(n, 1/sp.sqrt(n*(n + 2)), sp.Rational(1, 2)) for n in range(1, 5))
al, be_ = -1/sp.sqrt(3), 1/sp.sqrt(12)
R["the note's eight-dimensional alpha = -1/sqrt(3), beta = 1/sqrt(12): 2 alpha + 4 beta = 0, 4 beta = 2/sqrt(3), sqrt(12) beta = 1"] = (
    sp.simplify(2*al + 4*be_) == 0 and sp.simplify(4*be_ - 2/sp.sqrt(3)) == 0 and sp.simplify(sp.sqrt(12)*be_ - 1) == 0)
def kinetic_coef(n):
    """the coefficient c of (1/(2 kappa)) a^3 phi'^2/N that makes the reduction exact, solved from the Euler-Lagrange expression for N"""
    xs = [t] + list(sp.symbols('y1:%d' % (4 + n))); cc = sp.Symbol('cc')
    gd = [N**2*sp.exp(-n*ph), -a**2*sp.exp(-n*ph), -a**2*sp.exp(-n*ph), -a**2*sp.exp(-n*ph)] + [-sp.exp(2*ph)]*n
    LJ = -(1/(2*kap))*N*a**3*sp.exp(-n*ph)*scalar(gd, xs)
    RE = scalar([N**2, -a**2, -a**2, -a**2], xs[:4])
    LE = -(1/(2*kap))*N*a**3*RE + (1/(2*kap))*cc*a**3*sp.diff(ph, t)**2/N
    return sp.simplify(sp.solve(EL(sp.expand(LJ - LE), N, t), cc)[0])


kc = [kinetic_coef(n) for n in range(1, 5)]
R['the kinetic coefficient solved from the reduction (n = 1..4) is n(n+2)/2'] = all(sp.simplify(kc[n - 1] - sp.Rational(n*(n + 2), 2)) == 0 for n in range(1, 5))
R['NUM general n: kinetic coefficient n(n+2)/2 for n = 1..4'] = tuple(float(c) for c in kc)
R["NUM the note's beta_n^2 = 2/(n(n+2)) and the canonical 1/(n(n+2)), n = 4"] = (float(1/kc[3]), float(1/(2*kc[3])))

# ==================================================================================================== B. potential and mass
aa, bb = sp.Function('a')(t), sp.Function('b')(t)
Ha, Hb, dHa, dHb = sp.symbols('Ha Hb dHa dHb')


def toH(e):
    e = e.subs(sp.Derivative(aa, (t, 2)), (dHa + Ha**2)*aa).subs(sp.Derivative(bb, (t, 2)), (dHb + Hb**2)*bb)
    return sp.expand(sp.simplify(e.subs(sp.Derivative(aa, t), Ha*aa).subs(sp.Derivative(bb, t), Hb*bb)))


ET = einstein([1, -aa**2, -aa**2, -aa**2] + [-bb**2]*4, X)
ES = einstein([1, -aa**2, -aa**2, -aa**2] + [-bb**2*g for g in gam], X)
RT = scalar([1, -aa**2, -aa**2, -aa**2] + [-bb**2]*4, X); RS_ = scalar([1, -aa**2, -aa**2, -aa**2] + [-bb**2*g for g in gam], X)
R['the cosmological metric with the sphere of radius b: G^A_A = torus + (6, 6, 6, 6, 3, 3, 3, 3)/b^2, R = torus - 12/b^2'] = (
    all(sp.simplify(toH(ES[i]) - toH(ET[i]) - c/bb**2) == 0 for i, c in enumerate((6, 6, 6, 6, 3, 3, 3, 3)))
    and sp.simplify(toH(RS_) - toH(RT) + 12/bb**2) == 0)
eps, pa, pb, bs = sp.symbols('epsilon p_a p_b b_s', positive=True)
E00, E11, E44 = (toH(ES[i]).subs(bb, bs) for i in (0, 1, 4))
sol = sp.solve([sp.Eq(E11 - Lam, -kap*pa), sp.Eq(E44 - Lam, -kap*pb)], [dHa, dHb], dict=True)[0]
rad = sp.expand(sol[dHb]).subs(Ha**2, sp.expand((kap*eps + Lam - 12*Ha*Hb - 6*Hb**2 - 6/bs**2)/3))
R["radion equation on S^4: Hb' + 3 Ha Hb + 4 Hb^2 = kappa (eps/6 - pa/2 + pb/3) + Lambda/3 - 3/b^2"] = (
    sp.simplify(E00 - (3*(Ha**2 + 4*Ha*Hb + 2*Hb**2) + 6/bs**2)) == 0
    and sp.simplify(rad - (-3*Ha*Hb - 4*Hb**2 + kap*(eps/6 - pa/2 + pb/3) + Lam/3 - 3/bs**2)) == 0)

bet, U0, U1, KU, MP2s, fx, q2 = sp.symbols('beta U_0 U_1 K_U M2 f q2', real=True)
Uq = U0 + U1*bet + KU/2*bet**2
def Vpot(curv, flux): return MP2s*(sp.exp(-4*bet)*(Lam + kap*Uq) - curv*6*sp.exp(-6*bet)/r0**2 + flux*sp.exp(-12*bet))
def RHS(e, p_a, p_b, lam): return kap*(e/6 - p_a/2 + p_b/3) + lam/3 - 3/(r0**2*sp.exp(2*bet))
V0_ = Vpot(0, 0)
sT = sp.solve([V0_.subs(bet, 0), sp.diff(V0_, bet).subs(bet, 0)], [Lam, U1], dict=True)
R['radion on T^4 with U: flat vacuum Lambda + kappa U0 = 0, U1 = 0; m_r^2 = kappa K_U/12'] = (
    len(sT) == 1 and sp.simplify(sT[0][Lam] + kap*U0) == 0 and sT[0][U1] == 0
    and sp.simplify((sp.diff(V0_, bet, 2).subs(bet, 0)/(12*MP2s)).subs(sT[0]) - kap*KU/12) == 0)
V1_ = Vpot(1, 0)
sS = sp.solve([V1_.subs(bet, 0), sp.diff(V1_, bet).subs(bet, 0)], [Lam, U1], dict=True)[0]
Ub = Uq.subs(sS); R2 = RHS(Ub, -Ub, -Ub - sp.diff(Ub, bet)/4, sS[Lam])
mS = sp.simplify((sp.diff(V1_, bet, 2).subs(bet, 0)/(12*MP2s)).subs(sS))
R['radion on S^4 with U: Lambda + kappa U0 = 6/r0^2, kappa U1 = -12/r0^2, m_r^2 = kappa K_U/12 - 2/r0^2, from V and from the 8D radion equation'] = (
    sp.simplify(sS[Lam] + kap*U0 - 6/r0**2) == 0 and sp.simplify(kap*sS[U1] + 12/r0**2) == 0 and sp.simplify(mS - (kap*KU/12 - 2/r0**2)) == 0
    and sp.simplify(R2.subs(bet, 0)) == 0 and sp.simplify(-sp.diff(R2, bet).subs(bet, 0) - (kap*KU/12 - 2/r0**2)) == 0)
VL = V1_.subs({U0: 0, U1: 0, KU: 0})
LamL = sp.solve(sp.diff(VL, bet).subs(bet, 0), Lam)[0]
HE2 = sp.simplify(VL.subs(bet, 0).subs(Lam, LamL)/(3*MP2s)); m2L = sp.simplify(sp.diff(VL, bet, 2).subs(bet, 0).subs(Lam, LamL)/(12*MP2s))
dS = [sp.simplify(e.subs({Ha: 1/r0, Hb: 0, dHa: 0, dHb: 0, bs: r0}) - 9/r0**2) for e in (E00, E11, E44)]
R['S^4 with Lambda only: dS4 x S^4 (r0^2 = 9/Lambda, H^2 = 1/r0^2) solves the 00, 11, 44 equations; maximum of V with m^2 = -6 H^2, also from the 8D radion equation'] = (
    sp.simplify(LamL*r0**2) == 9 and sp.simplify(HE2*r0**2) == 1 and sp.simplify(m2L/HE2) == -6 and dS == [0, 0, 0]
    and sp.simplify(RHS(0, 0, 0, 9/r0**2).subs(bet, 0)) == 0 and sp.simplify(-sp.diff(RHS(0, 0, 0, 9/r0**2), bet).subs(bet, 0)*r0**2) == -6)
VF = Vpot(1, fx).subs({U0: 0, U1: 0, KU: 0})
sF = sp.solve([VF.subs(bet, 0), sp.diff(VF, bet).subs(bet, 0)], [Lam, fx], dict=True)[0]
# the flux from its comoving energy W(a, b) = a^3 b^4 q^2/(2 b^8): eps = W/(a^3 b^4), pa = -dW/da/(3 a^2 b^4), pb = -dW/db/(4 a^3 b^3)
As, Bs = sp.symbols('A_s B_s', positive=True)
W = As**3*Bs**4*q2/(2*Bs**8)
st = [sp.simplify(W/(As**3*Bs**4)), sp.simplify(-sp.diff(W, As)/(3*As**2*Bs**4)), sp.simplify(-sp.diff(W, Bs)/(4*As**3*Bs**3))]
rhoF = sF[fx]/kap*sp.exp(-8*bet)
RF = RHS(rhoF, -rhoF, rhoF, sF[Lam])
mF = sp.simplify(sp.diff(VF, bet, 2).subs(bet, 0).subs(sF)*r0**2/(12*MP2s))
R['S^4 with a 4-form flux: eps = rho_F, pa = -rho_F, pb = +rho_F; Minkowski needs Lambda = 9/(2 r0^2), kappa rho_F = 3/(2 r0^2); m^2 = 6/r0^2 from V and from the 8D radion equation'] = (
    [sp.simplify(st[0] - q2/(2*Bs**8)), sp.simplify(st[1] + q2/(2*Bs**8)), sp.simplify(st[2] - q2/(2*Bs**8))] == [0, 0, 0]
    and sp.simplify(sF[Lam]*r0**2) == sp.Rational(9, 2) and sp.simplify(sF[fx]*r0**2) == sp.Rational(3, 2) and mF == 6
    and sp.simplify(RF.subs(bet, 0)) == 0 and sp.simplify(-sp.diff(RF, bet).subs(bet, 0)*r0**2) == 6)
# the note's Goldberger-Wise potential: static internal profile Phi(x4, x5)
Phi = sp.Function('Phi')(x4, x5); mm = sp.Symbol('m', positive=True)
gdS = [N**2*sp.exp(-4*B), -a**2*sp.exp(-4*B), -a**2*sp.exp(-4*B), -a**2*sp.exp(-4*B)] + [-r0**2*sp.exp(2*B)*g for g in gam]
sqG = N*a**3*sp.exp(-4*B)*r0**4*sp.sin(x4)**3*sp.sin(x5)**2*sp.sin(x6)
lhs = sqG*(sp.Rational(1, 2)*sum(sp.diff(Phi, X[i])**2/gdS[i] for i in range(4, 8)) - mm**2*Phi**2/2)
rhs = -N*a**3*sp.exp(-4*B)*r0**4*sp.sin(x4)**3*sp.sin(x5)**2*sp.sin(x6)*(sp.exp(-2*B)/2*sum(sp.diff(Phi, X[i])**2/(r0**2*gam[i - 4]) for i in range(4, 8)) + mm**2*Phi**2/2)
R["the note's Goldberger-Wise potential: for a static Phi(y) the Einstein-frame energy carries e^(-4 beta) (gradient part e^(-6 beta)), not e^(+4 beta)"] = sp.simplify(lhs - rhs) == 0
R['NUM sphere with U: Lambda + kappa U0, kappa U1, m_r^2 - kappa K_U/12 (units 1/r0^2)'] = (
    float(sp.simplify((sS[Lam] + kap*U0)*r0**2)), float(sp.simplify(kap*sS[U1]*r0**2)), float(sp.simplify((mS - kap*KU/12)*r0**2)))
R['NUM sphere with Lambda only: Lambda r0^2, H^2 r0^2, m^2/H^2'] = (float(LamL*r0**2), float(HE2*r0**2), float(m2L/HE2))
R['NUM sphere with a 4-form flux: Lambda r0^2, kappa rho_F r0^2, m^2 r0^2'] = (float(sF[Lam]*r0**2), float(sF[fx]*r0**2), float(mF))

# ==================================================================================================== C. the radion and the brane
Ms, nn = sp.symbols('M n', positive=True)
pp = sp.series(-Ms*sp.exp(-nn*bet/2), bet, 0, 2).removeO()
g2 = ((nn/2)/sp.sqrt(nn*(nn + 2)/2))**2
R['radion coupling to a mass on the brane: -M e^(-n beta/2) = -M + (n/2) M beta; coupling sqrt(n/(2(n+2)))/M_Pl = 1/(sqrt(3) M_Pl) for n = 4; fifth force 2 (coupling M_Pl)^2 = n/(n+2) = 2/3'] = (
    sp.simplify(pp - (-Ms + nn*Ms*bet/2)) == 0 and sp.simplify(g2 - nn/(2*(nn + 2))) == 0 and g2.subs(nn, 4) == sp.Rational(1, 3)
    and sp.simplify(2*g2 - nn/(nn + 2)) == 0)
hbarc = 1.973269804e-16        # GeV m
MPl, M8 = 2.435e18, 1.0e4
vEW = (math.sqrt(2)*1.1663788e-5)**-0.5
EW = (hbarc/38.6e-6*1e9, hbarc/10.0, vEW/(math.sqrt(3)*MPl))
R['Eot-Wash: m_r > hbar c/38.6 micrometres = 5.1 meV; range of a 10 GeV radion 2.0e-17 m; gamma = v/(sqrt(3) M_Pl) = 5.8e-17'] = (
    abs(EW[0]/5.112e-3 - 1) < 1e-3 and abs(EW[1]/1.9733e-17 - 1) < 1e-4 and abs(EW[2]/5.838e-17 - 1) < 1e-3)
R['NUM Eot-Wash bound on the radion mass (eV), range of a 10 GeV radion (m), gamma = v/(sqrt(3) M_Pl)'] = EW

# ==================================================================================================== D. S^4 and T^4 side by side
th = sp.Symbol('theta')
def dl(l): return (2*l + 3)*(l + 1)*(l + 2)//6
ok = True
for l in range(7):
    C = sp.gegenbauer(l, sp.Rational(3, 2), sp.cos(th)); C1 = sp.gegenbauer(l, sp.Rational(3, 2), 1)
    ok = ok and sp.simplify(sp.diff(sp.sin(th)**3*sp.diff(C, th), th)/sp.sin(th)**3 + l*(l + 3)*C) == 0 and C1 == (l + 1)*(l + 2)//2
    ok = ok and sp.simplify(2*sp.pi**2*sp.integrate(sp.expand((C/C1)**2*sp.sin(th)**3), (th, 0, sp.pi)) - sp.Rational(8, 3)*sp.pi**2/dl(l)) == 0
R['zonal harmonics on S^4 (l = 0..6): eigenvalue -l(l+3), C_l(1) = (l+1)(l+2)/2, Int (C_l/C_l(1))^2 = V/d_l; d_l = binomial(l+4,4) - binomial(l+2,4) (l <= 20)'] = (
    ok and all(dl(l) == math.comb(l + 4, 4) - (math.comb(l + 2, 4) if l >= 2 else 0) for l in range(21)))
def NS4(x):
    s, l = 0, 0
    while l*(l + 3) <= x*x:
        s += dl(l); l += 1
    return s
weyl = tuple(NS4(x)/(x**4/12) - 1 for x in (10, 100, 1000))
R["Weyl's law on S^4: N(x)/(x^4/12) - 1 = -1/x^2 at x = 10, 100, 1000"] = all(abs(w + 1/x**2) < 1e-9 for w, x in zip(weyl, (10, 100, 1000)))
R['NUM S^4 multiplicities d_l, l = 0..8'] = tuple(dl(l) for l in range(9))
R['NUM S^4 Weyl: N(x)/(x^4/12) - 1 at x = mu r = 10, 100, 1000'] = weyl
# the embedding of S^4 in R^5, the north pole, SO(4)
emb = [r0*sp.sin(x4)*sp.sin(x5)*sp.sin(x6)*sp.cos(x7), r0*sp.sin(x4)*sp.sin(x5)*sp.sin(x6)*sp.sin(x7), r0*sp.sin(x4)*sp.sin(x5)*sp.cos(x6),
       r0*sp.sin(x4)*sp.cos(x5), r0*sp.cos(x4)]
J = sp.Matrix([[sp.diff(e, x) for x in (x4, x5, x6, x7)] for e in emb])
Npole = [e.subs(x4, 0) for e in emb]
kept = [(i, j) for i in range(5) for j in range(i + 1, 5) if all(v == 0 for v in [Npole[i]*(1 if k == j else 0) - Npole[j]*(1 if k == i else 0) for k in range(5)])]
R['the brane at the north pole of S^4: the embedding induces r0^2 gamma; north pole (0,0,0,0,r0); 6 of the 10 rotations of R^5 fix it (SO(4)), 4 move it'] = (
    sp.simplify(J.T*J - r0**2*sp.diag(*gam)) == sp.zeros(4) and Npole == [0, 0, 0, 0, r0] and len(kept) == 6)
# Killing spinors of S^4 (numerical, 30 digits, random points)
s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -I], [I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]]); Z2 = sp.zeros(2)
def blk(A_, B_, C_, D_): return sp.Matrix(sp.BlockMatrix([[A_, B_], [C_, D_]]))
g0 = sp.diag(1, 1, -1, -1); gis = [blk(Z2, s, -s, Z2) for s in (s1, s2, s3)]
rhoK = [I*g0, -gis[0], -gis[1], -gis[2]]
gt = [I*m for m in rhoK]
ths = [x4, x5, x6, x7]
hf = [r0, r0*sp.sin(x4), r0*sp.sin(x4)*sp.sin(x5), r0*sp.sin(x4)*sp.sin(x5)*sp.sin(x6)]
def om(mu, i, j):
    v = 0
    if mu == i: v += sp.diff(hf[i], ths[j])/hf[j]
    if mu == j: v -= sp.diff(hf[j], ths[i])/hf[i]
    return v
Sc = [sum((om(mu, i, j)*gt[i]*gt[j]/4 for i in range(4) for j in range(4) if i != j), sp.zeros(4)) for mu in range(4)]
E4 = sp.eye(4)
def KS(s):
    return ((sp.cos(x4/2)*E4 + s*I*sp.sin(x4/2)*gt[0])*(sp.cos(x5/2)*E4 + sp.sin(x5/2)*gt[0]*gt[1])
            *(sp.cos(x6/2)*E4 + sp.sin(x6/2)*gt[1]*gt[2])*(sp.cos(x7/2)*E4 + sp.sin(x7/2)*gt[2]*gt[3]))
random.seed(26)
pts = [{x4: sp.Rational(random.randint(200, 2900), 1000), x5: sp.Rational(random.randint(200, 2900), 1000), x6: sp.Rational(random.randint(200, 2900), 1000),
        x7: sp.Rational(random.randint(0, 6200), 1000), r0: sp.Rational(random.randint(500, 2000), 1000)} for _ in range(4)]   # exact rational points
worst = 0.0
for s in (1, -1):
    M = KS(s)
    for mu in range(4):
        E = sp.diff(M, ths[mu]) + Sc[mu]*M - s*I/(2*r0)*hf[mu]*gt[mu]*M
        for p in pts:
            worst = max(worst, max(abs(complex(v)) for v in E.subs(p).evalf(30)))
    Dm = sum((gt[a_]*(sp.diff(M, ths[a_]) + Sc[a_]*M)/hf[a_] for a_ in range(4)), sp.zeros(4)) - s*2*I/r0*M
    for p in pts:
        worst = max(worst, max(abs(complex(v)) for v in Dm.subs(p).evalf(30)))
herm = all(sp.simplify(g - g.H) == sp.zeros(4) for g in gt) and all(sp.simplify(gt[i]*gt[j] + gt[j]*gt[i] - 2*(1 if i == j else 0)*E4) == sp.zeros(4) for i in range(4) for j in range(4))
R['Killing spinors of S^4: nabla eps_s = s (i/(2 r0)) gt eps_s and D eps_s = (2 s i/r0) eps_s (s = +-1; residual < 1e-25 at four random points, 30 digits); Friedrich bound 4/r0^2'] = (
    herm and worst < 1e-25 and sp.simplify(4*(12/r0**2)/(4*3) - (2/r0)**2) == 0)
spw = sum(8*math.comb(l + 3, 3) for l in range(999))/(1000**4/3) - 1
R["Dirac spectrum on S^4: Weyl's law for spinors L^4/3 versus the degeneracies 2 x 4 binomial(l+3,3) (2.0e-3 at L = 1000)"] = abs(spw) < 5e-3
R['NUM Dirac on S^4: degeneracy 4 binomial(l + 3, 3) of +i(l + 2)/r, l = 0..5; Weyl ratio - 1 at L = 1000'] = tuple([4*math.comb(l + 3, 3) for l in range(6)] + [spw])
# numbers today, with the volume of the torus
Vp = MPl**2/M8**6; lp = Vp**0.25
rS = (3*Vp/(8*math.pi**2))**0.25
K = int(math.floor((10.0*lp/(4*math.pi))**2))
sig = np.zeros(K + 1, dtype=np.int64)
for d in range(1, K + 1):
    sig[d::d] += d
r4 = 8*sig.copy(); r4[4::4] -= 32*sig[1:K//4 + 1]
lattice = 1 + int(r4[1:].sum())
Nmu5 = Vp*5.0**4/(32*math.pi**2)
today = (rS, rS*hbarc, 2*math.pi/lp*1e3, 2/rS*1e3, 2/rS*1e3, lattice, NS4(5*rS), Nmu5)
R['the sphere with the volume of the torus today: r = 688.9 GeV^-1; lowest levels 4.03 MeV (torus) and 2.90 MeV (sphere); dirac16complex on the sphere 2.90 MeV; levels below 5 GeV agree with Weyl'] = (
    abs(rS/688.94 - 1) < 1e-4 and abs(today[3]/2.903 - 1) < 1e-3 and abs(today[6]/Nmu5 - 1) < 1e-3 and abs(lattice/Nmu5 - 1) < 1e-4)
R['NUM today: r_S4 (GeV^-1, m), lowest masses (MeV): torus 2 pi/l_perp, sphere scalar 2/r, sphere dirac16complex 2/r; levels below 5 GeV: torus lattice, sphere exact, Weyl'] = today

# ==================================================================================================== E. flat versus warped
ty, y = sp.symbols('t y', real=True); xs5 = [ty] + list(sp.symbols('u1:4')) + [y]
A = sp.Function('A')(y); k = sp.Symbol('k', positive=True)
G5 = einstein([sp.exp(-2*A), -sp.exp(-2*A), -sp.exp(-2*A), -sp.exp(-2*A), -1], xs5)
c5 = sp.expand(G5[0]).coeff(sp.Derivative(A, (y, 2)))                  # coefficient of A'' in G^0_0
jump = sp.diff(k*y, y) - sp.diff(-k*y, y)                                # A' jumps from -k to +k at the UV brane (A = k|y|)
fUV, fIR = sp.solve(sp.Eq(c5*jump, kap*sp.Symbol('f')), sp.Symbol('f'))[0], sp.solve(sp.Eq(c5*(-jump), kap*sp.Symbol('f')), sp.Symbol('f'))[0]
R['RS1 (5D, mostly minus): G^mu_nu = 3 A\'\' - 6 A\'^2, G^y_y = -6 A\'^2; bulk Lambda = -6 k^2; tensions f_UV = +6k/kappa, f_IR = -6k/kappa'] = (
    all(sp.simplify(G5[i] - (3*sp.diff(A, y, 2) - 6*sp.diff(A, y)**2)) == 0 for i in range(4)) and sp.simplify(G5[4] + 6*sp.diff(A, y)**2) == 0
    and all(sp.simplify(g.subs(sp.Derivative(A, (y, 2)), 0).subs(sp.Derivative(A, y), k) + 6*k**2) == 0 for g in G5)
    and c5 == 3 and jump == 2*k and sp.simplify(fUV - 6*k/kap) == 0 and sp.simplify(fIR + 6*k/kap) == 0)
g5 = [sp.exp(-2*A), -sp.exp(-2*A), -sp.exp(-2*A), -sp.exp(-2*A)]
Kmix = [sp.simplify(sp.diff(gg, y)/(2*gg)) for gg in g5]
Kjump = [sp.simplify(K_.subs(sp.Derivative(A, y), k) - K_.subs(sp.Derivative(A, y), -k)) for K_ in Kmix]
R['Israel junction at the UV brane: K^mu_nu = -A\' delta, [K^mu_nu] - delta [K] = 6k = 3 x (jump of A\'), the delta part of G^mu_nu; in the signature (-,+,+,+,+) [K - g K] = -kappa S'] = (
    all(sp.simplify(K_ + sp.diff(A, y)) == 0 for K_ in Kmix) and all(sp.simplify(kj - sum(Kjump) - c5*jump) == 0 for kj in Kjump) and sp.simplify(c5*jump - 6*k) == 0)
Rr = sp.Symbol('R_rs', positive=True)
R['RS1: M_Pl^2 = (1/kappa) Int e^(-2k|y|) dy over [-pi R, pi R] = (1 - e^(-2 pi k R))/(kappa k)'] = (
    sp.simplify(2*sp.integrate(sp.exp(-2*k*y), (y, 0, sp.pi*Rr))/kap - (1 - sp.exp(-2*sp.pi*k*Rr))/(kap*k)) == 0)
T = sp.Function('T')(ty)
R5 = scalar([sp.exp(-2*k*T*y), -sp.exp(-2*k*T*y), -sp.exp(-2*k*T*y), -sp.exp(-2*k*T*y), -T**2], xs5)
Ts, T1, T2 = sp.symbols('T_s T_1 T_2', positive=True)
Lb = 2*sp.integrate(sp.expand((-(1/(2*kap))*sp.exp(-4*k*T*y)*T*R5).subs(sp.Derivative(T, (ty, 2)), T2).subs(sp.Derivative(T, ty), T1).subs(T, Ts)), (y, 0, sp.pi), conds='none')
c2 = sp.simplify(sp.diff(Lb, T2)); L1 = sp.expand(sp.simplify(Lb - c2*T2 - sp.diff(c2, Ts)*T1**2))
kin = sp.simplify(sp.diff(L1, T1, 2)/2)
R['RS1 radion: kinetic term (3/(kappa k)) (d e^(-k pi T))^2, canonical mu = sqrt(6/(kappa k)) e^(-k pi T)'] = (
    not c2.has(T1) and not c2.has(T2) and sp.simplify(kin - 3/(kap*k)*(k*sp.pi*sp.exp(-k*sp.pi*Ts))**2) == 0 and sp.simplify(sp.diff(L1, T1).subs(T1, 0)) == 0)
xs8w = [ty] + list(sp.symbols('u1:4')) + [y] + list(sp.symbols('z1:4'))
G8w = einstein([sp.exp(-2*A), -sp.exp(-2*A), -sp.exp(-2*A), -sp.exp(-2*A), -1, -1, -1, -1], xs8w)
A1, A2 = sp.symbols('A1 A2')
solw = sp.solve([sp.Eq(3*A2 - 6*A1**2, -6*A1**2), sp.Eq(-6*A1**2, 4*A2 - 10*A1**2)], [A1, A2], dict=True)
R['8D warped along one extra direction: G^y_y = -6 A\'^2, G^z_z = 4 A\'\' - 10 A\'^2; with only Lambda A\'\' = 0 and A\' = 0'] = (
    all(sp.simplify(G8w[i] - (3*sp.diff(A, y, 2) - 6*sp.diff(A, y)**2)) == 0 for i in range(4)) and sp.simplify(G8w[4] + 6*sp.diff(A, y)**2) == 0
    and all(sp.simplify(G8w[i] - (4*sp.diff(A, y, 2) - 10*sp.diff(A, y)**2)) == 0 for i in range(5, 8)) and all(s_[A1] == 0 and s_[A2] == 0 for s_ in solw))
R["NUM RS1: coefficient of A'' in G^0_0, f_UV kappa/k, f_IR kappa/k, f^2 kappa k (canonical radion)"] = (
    float(c5), float(sp.simplify(fUV*kap/k)), float(sp.simplify(fIR*kap/k)), float(sp.simplify(2*kin/(k*sp.pi*sp.exp(-k*sp.pi*Ts))**2*kap*k)))
R['NUM strength of the radion coupling: warped with <mu> = 3 TeV over 7+1, i.e. sqrt(3) M_Pl/(3 TeV)'] = math.sqrt(3)*MPl/3000.0
# Goldberger-Wise in RS1: the bulk equation, its solution, the potential, the minimum and the radion mass
mg, C1, C2, nu = sp.symbols('m_g C_1 C_2 nu', positive=True)
Ph = sp.Function('Phi_y')(y)
Lgw = sp.exp(-4*k*y)*(-sp.Rational(1, 2)*sp.diff(Ph, y)**2 - sp.Rational(1, 2)*mg**2*Ph**2)
p0, p1 = sp.symbols('_p0 _p1')
Ls = Lgw.subs(sp.Derivative(Ph, y), p1).subs(Ph, p0)
ELgw = sp.diff(Ls, p0).subs({p0: Ph, p1: sp.diff(Ph, y)}) - sp.diff(sp.diff(Ls, p1).subs({p0: Ph, p1: sp.diff(Ph, y)}), y)
sol_gw = sp.exp(2*k*y)*(C1*sp.exp(nu*k*y) + C2*sp.exp(-nu*k*y)); nuv = sp.sqrt(4 + mg**2/k**2)
mu_, fm, vA, vB, ep = sp.symbols('mu f_mu v_A v_B epsilon', positive=True)
Vgw = 4*k*(mu_/fm)**4*(vB - vA*(mu_/fm)**ep)**2
mu0 = fm*(vB/vA)**(1/ep)
V2gw = sp.simplify(sp.powsimp(sp.diff(Vgw, mu_, 2).subs(mu_, mu0), force=True))
mgw = sp.simplify((V2gw/(k**2*ep**2*vB**2*(mu0/fm)**2)).subs(fm, sp.sqrt(6/(kap*k)))/kap)
mgw_note = sp.simplify((V2gw/(k**2*ep**2*vB**2*(mu0/fm)**2)).subs(fm, sp.sqrt(12/(kap*k)))/kap)
R["Goldberger-Wise: Phi'' - 4k Phi' - m^2 Phi = 0 from the bulk action, solved by e^(2ky)(C1 e^(nu k y) + C2 e^(-nu k y)), nu = sqrt(4 + m^2/k^2) = 2 + eps; V'(<mu>) = 0 at (<mu>/f)^eps = v_B/v_A, V''(<mu>) = 8 k eps^2 v_B^2 <mu>^2/f^4; m_r^2 = (4/3) kappa k^2 eps^2 v_B^2 (<mu>/f)^2 (note: 2/3)"] = (
    sp.simplify(ELgw - sp.exp(-4*k*y)*(sp.diff(Ph, y, 2) - 4*k*sp.diff(Ph, y) - mg**2*Ph)) == 0
    and sp.simplify((sp.diff(sol_gw, y, 2) - 4*k*sp.diff(sol_gw, y) - mg**2*sol_gw).subs(nu, nuv)) == 0
    and sp.simplify(sp.series(nuv, mg, 0, 3).removeO() - (2 + mg**2/(4*k**2))) == 0
    and sp.simplify(sp.powsimp(sp.diff(Vgw, mu_).subs(mu_, mu0), force=True)) == 0
    and sp.simplify(V2gw - 8*k*ep**2*vB**2*mu0**2/fm**4) == 0
    and sp.simplify((sp.expand_log(-sp.log(mu0/fm), force=True) - (4*k**2/mg**2)*(sp.log(vA) - sp.log(vB))).subs(ep, mg**2/(4*k**2))) == 0
    and mgw == sp.Rational(4, 3) and mgw_note == sp.Rational(2, 3))
R["NUM Goldberger-Wise: m_r^2/(kappa k^2 eps^2 v_B^2 (<mu>/f)^2) with f^2 = 6/(kappa k) and with the note's f^2 = 12/(kappa k)"] = (float(mgw), float(mgw_note))

# ==================================================================================================== F. Higgs-radion mixing on the brane
x1_, x2_, x3_ = X[1], X[2], X[3]
Om = sp.Function('Omega')(t, x1_)
R4 = scalar([Om**2, -Om**2, -Om**2, -Om**2], [t, x1_, x2_, x3_])
c4 = sp.simplify(Om**4*R4/(Om*(sp.diff(Om, t, 2) - sp.diff(Om, x1_, 2))))
R['conformally flat h = Omega^2 eta: sqrt|h| R(h) = -6 Omega box Omega, box = d0^2 - d1^2'] = (c4 == -6)


def EL2(L, q):
    """Euler-Lagrange expression of L(q, q_0, q_1, q_00, q_01, q_11) for the field q(x0, x1)"""
    names = [(), (t,), (x1_,), (t, t), (t, x1_), (x1_, x1_)]
    syms = sp.symbols('_e0:6')
    ders = [q] + [sp.diff(q, *nm) for nm in names[1:]]
    Ls = L
    order = [3, 4, 5, 1, 2, 0]                                    # the highest derivatives first
    for i in order:
        Ls = Ls.subs(ders[i], syms[i])
    back = {syms[i]: ders[i] for i in range(6)}
    out = sp.diff(Ls, syms[0]).subs(back)
    out -= sp.diff(sp.diff(Ls, syms[1]).subs(back), t) + sp.diff(sp.diff(Ls, syms[2]).subs(back), x1_)
    out += sp.diff(sp.diff(Ls, syms[3]).subs(back), t, 2) + sp.diff(sp.diff(Ls, syms[4]).subs(back), t, x1_) + sp.diff(sp.diff(Ls, syms[5]).subs(back), x1_, 2)
    return sp.simplify(out)


rf, hf_ = sp.Function('r')(t, x1_), sp.Function('h')(t, x1_)
o1, o2, epz, vv, lamH, xi, mr0, mh0 = sp.symbols('o_1 o_2 epsilon_z v lambda_H xi m_r0 m_h0', real=True)
OmE = 1 + o1*epz*rf + o2*epz**2*rf**2; PhE = vv + epz*hf_
Rh = scalar([OmE**2, -OmE**2, -OmE**2, -OmE**2], [t, x1_, x2_, x3_])
Lbr = OmE**4*(sp.Rational(1, 2)*OmE**-2*(sp.diff(PhE, t)**2 - sp.diff(PhE, x1_)**2) - lamH/4*(PhE**2 - vv**2)**2 + xi/2*Rh*PhE**2) \
    + epz**2*(sp.Rational(1, 2)*(sp.diff(rf, t)**2 - sp.diff(rf, x1_)**2) - sp.Rational(1, 2)*mr0**2*rf**2)
HL1 = sp.simplify(sp.diff(Lbr, epz).subs(epz, 0)); HL2 = sp.simplify((sp.diff(Lbr, epz, 2)/2).subs(epz, 0))
gm = o1*vv
Lq = (sp.Rational(1, 2)*(sp.diff(hf_, t)**2 - sp.diff(hf_, x1_)**2) + sp.Rational(1, 2)*(1 + 6*xi*gm**2)*(sp.diff(rf, t)**2 - sp.diff(rf, x1_)**2)
      + 6*xi*gm*(sp.diff(hf_, t)*sp.diff(rf, t) - sp.diff(hf_, x1_)*sp.diff(rf, x1_)) - sp.Rational(1, 2)*(2*lamH*vv**2)*hf_**2 - sp.Rational(1, 2)*mr0**2*rf**2)
ELr2 = sp.expand(EL2(HL2, rf))
Khr = -ELr2.coeff(sp.Derivative(hf_, (t, 2))); Krr = -ELr2.coeff(sp.Derivative(rf, (t, 2)))
R['the Higgs boson on the brane: linear terms are total derivatives; to second order L = (1/2)(dh)^2 + (1/2)(1 + 6 xi gamma^2)(dr)^2 + 6 xi gamma dh.dr - (1/2) m_h0^2 h^2 - (1/2) m_r0^2 r^2 up to total derivatives; o_2 drops out'] = (
    EL2(HL1, rf) == 0 and EL2(HL1, hf_) == 0 and EL2(HL2 - Lq, rf) == 0 and EL2(HL2 - Lq, hf_) == 0
    and sp.simplify(Khr - 6*xi*gm) == 0 and sp.simplify(Krr - (1 + 6*xi*gm**2)) == 0 and not ELr2.has(o2) and not sp.expand(EL2(HL2, hf_)).has(o2))
g_, Zs, th = sp.symbols('gamma Z theta', positive=True)
Km = sp.Matrix([[1, 6*xi*g_], [6*xi*g_, 1 + 6*xi*g_**2]]); Knote = sp.Matrix([[1, 6*xi*g_], [6*xi*g_, 1]])
Sm = sp.Matrix([[1, -6*xi*g_/Zs], [0, 1/Zs]]); Z2 = 1 + 6*xi*g_**2*(1 - 6*xi)
Mp = Sm.T*sp.diag(mh0**2, mr0**2)*Sm
Q = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])
off = sp.simplify(sp.expand_trig((Q.T*Mp*Q)[0, 1]) - (-(Mp[0, 0] - Mp[1, 1])*sp.sin(2*th)/2 + Mp[0, 1]*sp.cos(2*th)).expand(trig=True))
off_inv = sp.simplify(sp.expand_trig((Q*Mp*Q.T)[0, 1]) - ((Mp[0, 0] - Mp[1, 1])*sp.sin(2*th)/2 + Mp[0, 1]*sp.cos(2*th)).expand(trig=True))
R["the shift h = h' - (6 xi gamma/Z) r', r = r'/Z is canonical iff Z^2 = 1 + 6 xi gamma^2 (1 - 6 xi) (note's literal kinetic terms: 1 - 36 xi^2 gamma^2); tan 2 theta = 12 gamma xi Z m_h0^2/(m_r0^2 - m_h0^2 (Z^2 - 36 xi^2 gamma^2)) diagonalizes with (h', r') = Q(theta)(h_m, r_m); the inverse matrix needs -theta; C_h = cos - (6 xi gamma/Z) sin, C_r = sin/Z"] = (
    sp.simplify((Sm.T*Km*Sm).subs(Zs, sp.sqrt(Z2)) - sp.eye(2)) == sp.zeros(2)
    and sp.simplify((Sm.T*Knote*Sm)[1, 1]*Zs**2 - (1 - 36*xi**2*g_**2)) == 0 and sp.simplify((Sm.T*Knote*Sm)[0, 1]) == 0
    and off == 0 and off_inv == 0
    and sp.simplify(2*Mp[0, 1]/(Mp[0, 0] - Mp[1, 1]) - 12*g_*xi*Zs*mh0**2/(mr0**2 - mh0**2*(Zs**2 - 36*xi**2*g_**2))) == 0
    and sp.simplify((Sm*Q)[:, 0] - sp.Matrix([sp.cos(th) - 6*xi*g_*sp.sin(th)/Zs, sp.sin(th)/Zs])) == sp.zeros(2, 1))
mp.mp.dps = 60
def hr_num(gam, xi_, mr, mh):
    z = mp.sqrt(1 + 6*xi_*gam**2*(1 - 6*xi_))
    tt = mp.atan(12*gam*xi_*z*mh**2/(mr**2 - mh**2*(z**2 - 36*xi_**2*gam**2)))/2
    ch = mp.cos(tt) - 6*xi_*gam*mp.sin(tt)/z; cr = mp.sin(tt)/z
    return (z - 1, tt, (ch + cr*gam*mp.mpf(-11)/3/(mp.mpf(-13)/2))**2 - 1, (ch + cr*gam*7/(mp.mpf(-4)/3))**2 - 1)
hr71 = hr_num(mp.mpf(EW[2]), 1, mp.mpf(10), mp.mpf(125.20))
mp.mp.dps = 15
R['Higgs-radion mixing in 7+1 (gamma = v/(sqrt(3) M_Pl), xi = 1, m_r0 = 10 GeV): |Z - 1| < 1e-30, |theta| = 3.5e-16, rates change by < 1e-30; warped gamma = v/(3 TeV) = 0.082'] = (
    abs(hr71[0]) < 1e-30 and abs(abs(hr71[1])/3.53e-16 - 1) < 1e-2 and abs(hr71[2]) < 1e-30 and abs(hr71[3]) < 1e-30 and abs(vEW/3000/0.0821 - 1) < 1e-2)
R['NUM Higgs-radion mixing: c in sqrt|h| R(h) = c Omega box Omega, K_hr/(xi gamma), (K_rr - 1)/(xi gamma^2)'] = (
    float(c4), float(sp.simplify(Khr/(xi*gm))), float(sp.simplify((Krr - 1)/(xi*gm**2))))
R['NUM Higgs-radion mixing in 7+1 (gamma = v/(sqrt(3) M_Pl), xi = 1, m_r0 = 10 GeV): Z - 1, theta, Gamma(h -> gamma gamma)/SM - 1, Gamma(h -> g g)/SM - 1; warped gamma = v/(3 TeV)'] = (
    tuple(float(x) for x in hr71) + (vEW/3000.0,))

# ---------------------------------------------------------------------------------------------------- output
nb = 0; nf = 0
for key, v in R.items():
    print('%s : %r' % (key, v))
    if isinstance(v, bool):
        nb += 1; nf += (not v)
print('boolean checks: %d   failed: %d' % (nb, nf))
print('VERDICT: ' + ('ALL CHECKS PASSED' if nf == 0 else 'FAILURES'))
sys.exit(0 if nf == 0 else 1)
