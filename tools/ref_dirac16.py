#!/usr/bin/env python3
"""Independent Python (sympy / mpmath) check of sections 7-9 of 7+1.nb:
the Dirac equation of the 16-component field dirac16complex on 7+1 spacetime with the REAL 16x16
matrices Gamma17, its chiral 8-component form, its quantization, and the 16-component and chiral
Kohn-Sham theories of the self-coupled field.

Run from the repository root:   python3 tools/ref_dirac16.py
Every line of output is  name : value ; the script exits with status 1 if any check is not True.
"""
import itertools, math, sys
import sympy as sp
import mpmath as mp
from ref71 import *                      # G17 (real Gamma17), tauO, eta17, blk, products, ...

R = {}                                   # name -> result (True/False or a number)
I8, I16, Z8 = sp.eye(8), sp.eye(16), sp.zeros(8)
G = {A: G17[A] for A in range(9)}        # G[0..7]: Gamma17[A] (lower index); G[8] = omega = G0 G1 ... G7
eta = eta17                              # diag(+1,-1,...,-1): simply the squares of the octonion basis (e_0^2 = +1, e_i^2 = -1)
Gup = {A: eta[A, A]*G[A] for A in range(8)}            # Gamma^A = eta^{AA} Gamma_A
beta = G[0]
alpha = {i: G[0]*Gup[i] for i in range(1, 8)}          # alpha^i = Gamma^0 Gamma^i
omega = G[8]
chi = sp.I*omega                                       # chirality
PL, PR = (I16 - chi)/2, (I16 + chi)/2
EL = sp.Matrix.vstack(I8, -sp.I*I8)/sp.sqrt(2)         # 16x8 isometry onto the chi = -1 subspace
ER = sp.Matrix.vstack(-sp.I*I8, I8)/sp.sqrt(2)         # 16x8 isometry onto the chi = +1 subspace
W = sp.Matrix.hstack(EL, ER)
sig = {i: sp.I*tauO[i] for i in range(1, 8)}           # 8x8 Hermitian sigma_i = i tauO_i
dag = lambda M: M.conjugate().T
zero = lambda M: sp.simplify(M) == sp.zeros(*M.shape)

# ---------------------------------------------------------------- 1. the real gamma matrices
R['Gamma17 real, entries 0,+-1'] = all(set(G[A]) <= {0, 1, -1} for A in range(8))
R['{Gamma_A,Gamma_B} = 2 eta17_AB'] = all(G[A]*G[B] + G[B]*G[A] == 2*eta[A, B]*I16 for A in range(8) for B in range(8))
R['Gamma0 symmetric, Gamma_i antisymmetric'] = G[0].T == G[0] and all(G[i].T == -G[i] for i in range(1, 8))
R['Gamma^A dagger = Gamma0 Gamma^A Gamma0'] = all(dag(Gup[A]) == G[0]*Gup[A]*G[0] for A in range(8))
R['alpha^i, beta real symmetric'] = beta.T == beta and all(alpha[i].T == alpha[i] and set(alpha[i]) <= {0, 1, -1} for i in alpha)
R['{alpha_i,alpha_j}=2delta, {alpha_i,beta}=0, beta^2=1'] = (all(alpha[i]*alpha[j] + alpha[j]*alpha[i] == 2*(i == j)*I16 for i in alpha for j in alpha)
    and all(alpha[i]*beta + beta*alpha[i] == sp.zeros(16) for i in alpha) and beta*beta == I16)
R['block forms alpha^i = [[0,tauO_i],[-tauO_i,0]], beta = [[0,1],[1,0]]'] = (beta == blk(Z8, I8, I8, Z8)
    and all(alpha[i] == blk(Z8, tauO[i], -tauO[i], Z8) for i in alpha))
R['omega = Gamma17[8] = [[0,-1],[1,0]], omega^2 = -1'] = omega == blk(Z8, -I8, I8, Z8) and omega*omega == -I16

# ---------------------------------------------------------------- 2. the Dirac equation in position space
X = sp.symbols('x0:8', real=True)
m = sp.symbols('m', positive=True)
psi = sp.Matrix([sp.Function('psi%d' % a)(*X) for a in range(16)])       # 16 arbitrary component functions
Dop = lambda f: sum((Gup[A]*f.diff(X[A]) for A in range(8)), sp.zeros(f.shape[0], 1))
box = lambda f: sum((eta[A, A]*f.diff(X[A], 2) for A in range(8)), sp.zeros(f.shape[0], 1))
dirac = lambda f, mm: sp.I*Dop(f) - mm*f
R['(i Dslash + m)(i Dslash - m) psi = -(box + m^2) psi'] = zero(sp.I*Dop(dirac(psi, m)) + m*dirac(psi, m) + box(psi) + m**2*psi)
hpsi = lambda f, mm: sum((-sp.I*alpha[i]*f.diff(X[i]) for i in range(1, 8)), sp.zeros(16, 1)) + beta*mm*f
R['Gamma0 (i d0 - h) psi = (i Dslash - m) psi'] = zero(G[0]*(sp.I*psi.diff(X[0]) - hpsi(psi, m)) - dirac(psi, m))

# ---------------------------------------------------------------- 3. plane waves and the real spinors u_s, v_s
P = sp.symbols('p1:8', real=True)
psq = sum(q**2 for q in P)
En = sp.sqrt(m**2 + psq)
tauP = sum((P[i-1]*tauO[i] for i in range(1, 8)), sp.zeros(8))
Mm = m*I8 + tauP                                   # M(p) = m + p^i tauO_i
hp = lambda Q, mm: sum((Q[i-1]*alpha[i] for i in range(1, 8)), sp.zeros(16)) + mm*beta
R['h(p) = [[0,M],[M^T,0]]'] = hp(P, m) == blk(Z8, Mm, Mm.T, Z8)
R['h(p)^2 = E^2 1'] = zero((hp(P, m)*hp(P, m)).expand() - (m**2 + psq)*I16)
R['M M^T = M^T M = E^2 1_8'] = zero((Mm*Mm.T).expand() - (m**2 + psq)*I8) and zero((Mm.T*Mm).expand() - (m**2 + psq)*I8)
e8 = lambda s: I8[:, s]
u = lambda s, Mx, Ex: sp.Matrix.vstack(Ex*e8(s), Mx.T*e8(s))/sp.sqrt(Ex)
v = lambda s, Mx, Ex: sp.Matrix.vstack(Ex*e8(s), -Mx*e8(s))/sp.sqrt(Ex)
Mtil = m*I8 - tauP                                 # M(-p) = M(p)^T
pslash = lambda E0, Q: E0*G[0] + sum((Q[i-1]*G[i] for i in range(1, 8)), sp.zeros(16))
Pm = [-q for q in P]
R['h(p) u_s = E u_s,  h(-p) v_s = -E v_s'] = all(zero(hp(P, m)*u(s, Mm, En) - En*u(s, Mm, En)) and zero(hp(Pm, m)*v(s, Mm, En) + En*v(s, Mm, En)) for s in range(8))
R['(pslash - m) u_s = 0, (pslash + m) v_s = 0'] = all(zero((pslash(En, P) - m*I16)*u(s, Mm, En)) and zero((pslash(En, P) + m*I16)*v(s, Mm, En)) for s in range(8))
US = sp.Matrix.hstack(*[u(s, Mm, En) for s in range(8)]); VS = sp.Matrix.hstack(*[v(s, Mm, En) for s in range(8)])
VSt = sp.Matrix.hstack(*[v(s, Mtil, En) for s in range(8)])          # v_s(p~), p~ = (E, -p)
R['u_s^T u_r = 2E delta, v_s^T v_r = 2E delta'] = zero((US.T*US).expand() - 2*En*I8) and zero((VS.T*VS).expand() - 2*En*I8)
R['u_s(p)^T v_r(p~) = 0'] = zero((US.T*VSt).expand())
R['sum_s u_s ubar_s = pslash + m, sum_s v_s vbar_s = pslash - m'] = zero(US*US.T*G[0] - pslash(En, P) - m*I16) and zero(VS*VS.T*G[0] - pslash(En, P) + m*I16)
R['completeness sum u u^T + sum v(p~) v(p~)^T = 2E 1_16'] = zero(US*US.T + VSt*VSt.T - 2*En*I16)
Up = sp.Matrix.hstack(US, VSt)/sp.sqrt(2*En)
R['U_p orthogonal and U_p^T h(p) U_p = diag(E 1_8, -E 1_8)'] = zero(Up.T*Up - I16) and zero(Up.T*hp(P, m)*Up - sp.diag(En*I8, -En*I8))
R['charge conjugation: (pslash + m) omega u_s = 0'] = all(zero((pslash(En, P) + m*I16)*omega*u(s, Mm, En)) for s in range(8))
R['omega u_s = sum_r v_r (-M^T/E)_rs'] = zero(omega*US - VS*(-Mm.T/En))
# Lorentz covariance for the boost of 7+1.nb section 5.2 at beta = 2 log 2 (cosh(beta/2) = 5/4, sinh(beta/2) = 3/4)
ch, sh = sp.Rational(5, 4), sp.Rational(3, 4)
Sb, Sbi = ch*I16 + sh*G[0]*G[1], ch*I16 - sh*G[0]*G[1]
Lam = sp.Matrix(8, 8, lambda B, K: sp.nsimplify((eta[B, B]*G[B]*Sb*G[K]*Sbi).trace()/16))
R['S Gamma_K S^-1 = sum_B Lambda_BK Gamma_B (boost)'] = all(Sb*G[K]*Sbi == sum((Lam[B, K]*G[B] for B in range(8)), sp.zeros(16)) for K in range(8))
p0 = [2, 2, 4, 0, 0, 0, 0]; M0 = I8 + sum((p0[i-1]*tauO[i] for i in range(1, 8)), sp.zeros(8)); E0 = 5
pp = Lam*sp.Matrix([E0] + p0)
R['Lorentz covariance: (pslash(Lambda p) - m) S u_s(p) = 0'] = all(zero((pslash(pp[0], list(pp[1:])) - I16)*Sb*u(s, M0, E0)) for s in range(8))

# ---------------------------------------------------------------- 4. chirality: no real chirality in 7+1, complex projectors
import numpy as np
Gn = [np.array(G[A].tolist(), dtype=float) for A in range(8)]
sysA = np.vstack([np.kron(np.eye(16), g.T) + np.kron(g, np.eye(16)) for g in Gn])   # J G_A + G_A J = 0, J real 16x16
R['real matrices anticommuting with all Gamma17[A]: dimension 1 (rank 255)'] = int(np.linalg.matrix_rank(sysA)) == 255
R['omega anticommutes with all Gamma17[A]'] = all(omega*G[A] + G[A]*omega == sp.zeros(16) for A in range(8))
R['chi = i omega: Hermitian, chi^2 = 1, {chi,Gamma_A} = 0'] = dag(chi) == chi and chi*chi == I16 and all(chi*G[A] + G[A]*chi == sp.zeros(16) for A in range(8))
Sig = lambda A, B: (G[A]*G[B] - G[B]*G[A])/4
R['[chi, Sigma_AB] = 0 (chirality is Lorentz invariant)'] = all(chi*Sig(A, B) == Sig(A, B)*chi for A in range(8) for B in range(A + 1, 8))
R['P_L, P_R not real'] = any(sp.im(x) != 0 for x in PL) and any(sp.im(x) != 0 for x in PR)
R['P_L+P_R = 1, P_L^2 = P_L, P_R^2 = P_R, P_L P_R = 0, rank 8 each'] = (PL + PR == I16 and PL*PL == PL and PR*PR == PR
    and PL*PR == sp.zeros(16) and PL.rank() == 8 and PR.rank() == 8)
R['Gamma_A P_L = P_R Gamma_A,  P_L^dagger Gamma0 = Gamma0 P_R'] = all(G[A]*PL == PR*G[A] for A in range(8)) and dag(PL)*G[0] == G[0]*PR
R['E_L, E_R isometries: E^dag E = 1_8, E_L^dag E_R = 0, E_L E_L^dag = P_L, E_R E_R^dag = P_R'] = (
    dag(EL)*EL == I8 and dag(ER)*ER == I8 and dag(EL)*ER == sp.zeros(8) and EL*dag(EL) == PL and ER*dag(ER) == PR)
R['chi E_L = -E_L, chi E_R = E_R, W unitary'] = chi*EL == -EL and chi*ER == ER and dag(W)*W == I16
sigP = sum((P[i-1]*sig[i] for i in range(1, 8)), sp.zeros(8))
R['W^dag h(p) W = [[-sigma.p, m],[m, sigma.p]]'] = zero(dag(W)*hp(P, m)*W - blk(-sigP, m*I8, m*I8, sigP))
R['sigma_i Hermitian, {sigma_i,sigma_j} = 2 delta'] = all(dag(sig[i]) == sig[i] for i in sig) and all(sig[i]*sig[j] + sig[j]*sig[i] == 2*(i == j)*I8 for i in sig for j in sig)
R['chiral blocks: E_L^dag beta E_L = 0, E_L^dag beta E_R = 1, E_L^dag alpha E_L = -sigma, E_R^dag alpha E_R = sigma, E_L^dag alpha E_R = 0'] = (
    dag(EL)*beta*EL == sp.zeros(8) and dag(EL)*beta*ER == I8 and all(dag(EL)*alpha[i]*EL == -sig[i] and dag(ER)*alpha[i]*ER == sig[i]
    and dag(EL)*alpha[i]*ER == sp.zeros(8) for i in alpha))
# the two coupled 8-component equations, for 16 arbitrary component functions
phiL, phiR = dag(EL)*psi, dag(ER)*psi
tdL = lambda f: sum((tauO[i]*f.diff(X[i]) for i in range(1, 8)), sp.zeros(8, 1))
R['E_L^dag Gamma0 (i Dslash - m) psi = (i d0 + tauO.d) phi_L - m phi_R'] = zero(dag(EL)*G[0]*dirac(psi, m) - (sp.I*phiL.diff(X[0]) + tdL(phiL) - m*phiR))
R['E_R^dag Gamma0 (i Dslash - m) psi = (i d0 - tauO.d) phi_R - m phi_L'] = zero(dag(ER)*G[0]*dirac(psi, m) - (sp.I*phiR.diff(X[0]) - tdL(phiR) - m*phiL))
f8 = sp.Matrix([sp.Function('f%d' % a)(*X) for a in range(8)])
R['(i d0 - tauO.d)(i d0 + tauO.d) f = -box f (8 components)'] = zero((sp.I*(sp.I*f8.diff(X[0]) + tdL(f8)).diff(X[0]) - tdL(sp.I*f8.diff(X[0]) + tdL(f8))) + box(f8))
R['psi = E_L phi_L + E_R phi_R'] = zero(EL*phiL + ER*phiR - psi)
# bilinears in chiral form
R['chiral form of the bilinears (beta i chi, chi, alpha chi)'] = (dag(EL)*beta*sp.I*chi*ER == sp.I*I8 and dag(ER)*beta*sp.I*chi*EL == -sp.I*I8
    and dag(EL)*chi*EL == -I8 and dag(ER)*chi*ER == I8 and all(dag(EL)*alpha[i]*chi*EL == sig[i] and dag(ER)*alpha[i]*chi*ER == sig[i] for i in alpha))

# ---------------------------------------------------------------- 5. quantization: canonical anticommutators from the mode algebra
t = X[0]; xs7 = X[1:]; ys7 = sp.symbols('y1:8', real=True); Vb = sp.symbols('V', positive=True)
mom = [(2, 2, 4, 0, 0, 0, 0), (-2, -2, -4, 0, 0, 0, 0), (0, 0, 0, 2, 2, 0, 0), (0, 0, 0, -2, -2, 0, 0), (0, 0, 0, 0, 0, 0, 0)]
Emo = [sp.sqrt(1 + sum(c*c for c in q)) for q in mom]                 # m = 1: E = 5, 5, 3, 3, 1
Mmo = [I8 + sum((q[i-1]*tauO[i] for i in range(1, 8)), sp.zeros(8)) for q in mom]
dot = lambda q, x: sum(qi*xi for qi, xi in zip(q, x))
def field(x):   # dict operator -> 16-vector coefficient;   ('b',n,s) annihilates a particle, ('dd',n,s) creates an antiparticle
    out = {}
    for n, q in enumerate(mom):
        for s in range(8):
            out[('b', n, s)] = u(s, Mmo[n], Emo[n])*sp.exp(-sp.I*Emo[n]*t + sp.I*dot(q, x))/sp.sqrt(2*Emo[n]*Vb)
            out[('dd', n, s)] = v(s, Mmo[n], Emo[n])*sp.exp(sp.I*Emo[n]*t - sp.I*dot(q, x))/sp.sqrt(2*Emo[n]*Vb)
    return out
adjop = {'b': 'bd', 'dd': 'd', 'bd': 'b', 'd': 'dd'}
def acomm(F1, F2dag):  # {F1_a, F2dag_b} for F1 = sum op*c, F2dag = sum adj(op)*conj(c);  {b,bd} = {dd,d} = 1 for equal labels
    out = sp.zeros(16)
    for key, c1 in F1.items():
        k2 = (adjop[key[0]],) + key[1:]
        if k2 in F2dag: out += c1*F2dag[k2].T
    return out
F1 = field(xs7)
F2dag = {(adjop[k[0]],) + k[1:]: c.conjugate() for k, c in field(ys7).items()}
CAR = acomm(F1, F2dag)
target = sum((sp.exp(sp.I*dot(q, [a - b for a, b in zip(xs7, ys7)])) for q in mom), 0)/Vb*I16
R['{psi_a(t,x), psi_b^dag(t,y)} = delta_ab (1/V) sum_q e^{iq(x-y)} (5 momenta, abstract mode algebra)'] = zero(sp.expand(CAR - target))
F2 = field(ys7)
R['{psi_a(t,x), psi_b(t,y)} = 0'] = all(adjop[k[0]] + '' not in [kk[0] for kk in F2] for k in F1)   # no adjoint pairs among b, dd

# ---------------------------------------------------------------- 6. Grassmann algebra
def parity(seq):
    s, inv = list(seq), 0
    for i in range(len(s)):
        for j in range(i + 1, len(s)):
            if s[i] > s[j]: inv += 1
    return -1 if inv % 2 else 1
def gmul(x, y):
    out = {}
    for s1, c1 in x.items():
        for s2, c2 in y.items():
            if set(s1) & set(s2): continue
            key = tuple(sorted(s1 + s2)); out[key] = out.get(key, 0) + parity(s1 + s2)*c1*c2
    return {k: c for k, c in out.items() if c != 0}
def gadd(*xs):
    out = {}
    for x in xs:
        for s, c in x.items(): out[s] = out.get(s, 0) + c
    return out
def gscal(a, x): return {s: a*c for s, c in x.items()}
def gdL(x, k):   # left derivative d/dtheta_k
    out = {}
    for s, c in x.items():
        if k in s:
            j = s.index(k); key = s[:j] + s[j+1:]; out[key] = out.get(key, 0) + (-1)**j*c
    return out
def gsimp(x): return {s: sp.simplify(c) for s, c in x.items() if sp.simplify(c) != 0}
# (a) a complex Grassmann field built from c-number plane-wave solutions: psi = sum_k theta_k U_k(x)
U1 = u(0, Mmo[0], Emo[0])*sp.exp(-sp.I*(Emo[0]*t - dot(mom[0], xs7)))
U2 = v(2, Mmo[2], Emo[2])*sp.exp(sp.I*(Emo[2]*t - dot(mom[2], xs7)))
U3 = u(4, Mmo[4], Emo[4])*sp.exp(-sp.I*Emo[4]*t)
Us = [U1, U2, U3]
R['the three U_k solve the Dirac equation (m = 1)'] = all(zero(dirac(Uk, 1)) for Uk in Us)
psiG = [{(k,): Us[k][a] for k in range(3)} for a in range(16)]                       # theta_k = generator k
psibarG = [{(3 + k,): (Us[k].conjugate().T*G[0])[a] for k in range(3)} for a in range(16)]   # thetabar_k = generator 3+k
R['Grassmann components anticommute: psi_a psi_b + psi_b psi_a = 0, psi_a^2 = 0'] = all(
    gsimp(gadd(gmul(psiG[a], psiG[b]), gmul(psiG[b], psiG[a]))) == {} for a in range(16) for b in range(a, 16))
def bil(Mx):   # psibar Mx psi
    return gadd(*[gscal(Mx[a, b], gmul(psibarG[a], psiG[b])) for a in range(16) for b in range(16) if Mx[a, b] != 0])
jG = {A: bil(Gup[A]) for A in range(8)}
div = {}
for A in range(8):
    for s, c in jG[A].items(): div[s] = div.get(s, 0) + sp.diff(c, X[A])
R['current conservation d_A (psibar Gamma^A psi) = 0 (Grassmann field)'] = all(sp.simplify(c) == 0 for c in div.values())
# (b) the nonlinear field equation of the self-coupled field, by the left derivative with respect to psibar_a
gS, gP, gV, gA = sp.symbols('g_S g_P g_V g_A', real=True)
psiP = [{(a,): 1} for a in range(16)]; psibP = [{(16 + a,): 1} for a in range(16)]
def bilP(Mx): return gadd(*[gscal(Mx[a, b], gmul(psibP[a], psiP[b])) for a in range(16) for b in range(16) if Mx[a, b] != 0])
iChi = sp.I*chi
Ops = {'S': (gS, [(I16, I16)]), 'P': (gP, [(iChi, iChi)]), 'V': (gV, [(Gup[A], G[A]) for A in range(8)]), 'A': (gA, [(Gup[A]*chi, G[A]*chi) for A in range(8)])}
Lint = {}
for key, (g, pairs) in Ops.items():
    for O1, O2 in pairs: Lint = gadd(Lint, gscal(-g/2, gmul(bilP(O1), bilP(O2))))
ok = True
for a in range(16):
    lhs = gdL(Lint, 16 + a)
    rhs = {}
    for key, (g, pairs) in Ops.items():
        for O1, O2 in pairs:   # -g (psibar O2 psi)(O1 psi)_a
            O1psi = gadd(*[gscal(O1[a, b], psiP[b]) for b in range(16) if O1[a, b] != 0])
            rhs = gadd(rhs, gscal(-g, gmul(bilP(O2), O1psi)))
    if gsimp(gadd(lhs, gscal(-1, rhs))) != {}: ok = False; break
R['dL_int/dpsibar_a = -[gS(psibar psi) psi + gP (psibar i chi psi) i chi psi + gV (psibar G_A psi) G^A psi + gA (psibar G_A chi psi) G^A chi psi]_a'] = ok

# ---------------------------------------------------------------- 7. Kohn-Sham: Hermitian KS Hamiltonian, chiral blocks, potentials = derivatives of E_H
Ssym, Psym = sp.symbols('S P', real=True)
Vs = sp.symbols('V0:8', real=True); As = sp.symbols('A0:8', real=True)          # V_A, A_A with LOWER index
hKS = lambda Q: (sum((alpha[i]*(Q[i-1] + Vs[i]) for i in range(1, 8)), sp.zeros(16)) + beta*(m + Ssym) + Vs[0]*I16
                 + beta*iChi*Psym + chi*As[0] + sum((alpha[i]*chi*As[i] for i in range(1, 8)), sp.zeros(16)))
R['h_KS Hermitian for real S, P, V_A, A_A'] = zero(dag(hKS(P)) - hKS(P))
epsy = sp.symbols('epsilon', real=True)
R['Gamma0 [pslash(eps,p) - m - S - i chi P - Gamma^A (V_A + chi A_A)] = eps - h_KS(p)  (plane wave e^{-i(eps t - p.x)})'] = zero(
    G[0]*(pslash(epsy, P) - (m + Ssym)*I16 - iChi*Psym - sum((Gup[A]*(Vs[A]*I16 + chi*As[A]) for A in range(8)), sp.zeros(16)))
    - (epsy*I16 - hKS(P)))
VL = [Vs[A] - As[A] for A in range(8)]; VR = [Vs[A] + As[A] for A in range(8)]
Mc = m + Ssym + sp.I*Psym
sigQ = lambda Q: sum((Q[i-1]*sig[i] for i in range(1, 8)), sp.zeros(8))
hchi = blk(VL[0]*I8 - sigQ([P[i-1] + VL[i] for i in range(1, 8)]), Mc*I8, sp.conjugate(Mc)*I8, VR[0]*I8 + sigQ([P[i-1] + VR[i] for i in range(1, 8)]))
R['W^dag h_KS W = [[V^L_0 - sigma.(p+V^L), M],[M^*, V^R_0 + sigma.(p+V^R)]], V^L = V - A, V^R = V + A, M = m + S + iP'] = zero(dag(W)*hKS(P)*W - hchi)
# 16-component Hartree energy and its derivative (z = orbital components, zb = their conjugates, independent symbols)
zs16, zbs16 = sp.symbols('z0:16'), sp.symbols('zb0:16')
Zc, ZBc = sp.Matrix(zs16), sp.Matrix(zbs16)
dens = lambda Mx, Zv=Zc, ZBv=ZBc: sp.expand((ZBv.T*Mx*Zv)[0, 0])          # zb^T Mx z   (Mx = Gamma0 O  ->  psibar O psi)
def hartree(Zv, ZBv):
    rS_, rP_ = dens(beta, Zv, ZBv), dens(beta*iChi, Zv, ZBv)
    j_ = [dens(beta*Gup[A], Zv, ZBv) for A in range(8)]; j5_ = [dens(beta*Gup[A]*chi, Zv, ZBv) for A in range(8)]
    EH_ = sp.expand(gS/2*rS_**2 + gP/2*rP_**2 + gV/2*sum(eta[A, A]*j_[A]**2 for A in range(8)) + gA/2*sum(eta[A, A]*j5_[A]**2 for A in range(8)))
    return EH_, rS_, rP_, j_, j5_
EH, rS, rP, jv, j5 = hartree(Zc, ZBc)
Sv, Pv = gS*rS, gP*rP
Vl = [gV*eta[A, A]*jv[A] for A in range(8)]; Al = [gA*eta[A, A]*j5[A] for A in range(8)]     # V_A, A_A (lower index)
vKS = (beta*Sv + beta*iChi*Pv + sum((G[0]*Gup[A]*Vl[A] for A in range(8)), sp.zeros(16))
       + sum((G[0]*Gup[A]*chi*Al[A] for A in range(8)), sp.zeros(16)))
R['Hartree potentials: dE_H/dphi*_a = (v_KS phi)_a with S = gS rhoS, P = gP rhoP, V_A = gV j_A, A_A = gA j5_A'] = all(
    sp.expand(sp.diff(EH, zbs16[a]) - (vKS*Zc)[a]) == 0 for a in range(16))
R['densities real: Gamma0 O Hermitian for O = 1, i chi, Gamma^A, Gamma^A chi'] = all(dag(G[0]*O) == G[0]*O for O in [I16, iChi] + [Gup[A] for A in range(8)] + [Gup[A]*chi for A in range(8)])
# chiral form of the Hartree energy and of its derivatives
zL8, zLb8, zR8, zRb8 = sp.symbols('zl0:8'), sp.symbols('zlb0:8'), sp.symbols('zr0:8'), sp.symbols('zrb0:8')
ZL, ZLb, ZR, ZRb = map(sp.Matrix, (zL8, zLb8, zR8, zRb8))
psiC = EL*ZL + ER*ZR
psiCbT = ZLb.T*dag(EL) + ZRb.T*dag(ER)                       # psi^dagger as a row
densC = lambda Mx: sp.expand((psiCbT*Mx*psiC)[0, 0])
kap, kaps = sp.expand((ZLb.T*ZR)[0, 0]), sp.expand((ZRb.T*ZL)[0, 0])     # kappa = phi_L^dag phi_R, kappa^* = phi_R^dag phi_L
sbarA = [I8] + [-sig[i] for i in range(1, 8)]; sA = [I8] + [sig[i] for i in range(1, 8)]
jL = [sp.expand((ZLb.T*sbarA[A]*ZL)[0, 0]) for A in range(8)]; jR = [sp.expand((ZRb.T*sA[A]*ZR)[0, 0]) for A in range(8)]
R['chiral densities: rhoS = kappa + kappa^*, rhoP = i(kappa - kappa^*), j = jL + jR, j5 = jR - jL'] = (
    sp.expand(densC(beta) - kap - kaps) == 0 and sp.expand(densC(beta*iChi) - sp.I*(kap - kaps)) == 0
    and all(sp.expand(densC(beta*Gup[A]) - jL[A] - jR[A]) == 0 and sp.expand(densC(beta*Gup[A]*chi) - jR[A] + jL[A]) == 0 for A in range(8)))
EHc16 = sp.expand(gS/2*densC(beta)**2 + gP/2*densC(beta*iChi)**2 + gV/2*sum(eta[A, A]*densC(beta*Gup[A])**2 for A in range(8))
                  + gA/2*sum(eta[A, A]*densC(beta*Gup[A]*chi)**2 for A in range(8)))
EHchi = sp.expand(gS/2*(kap + kaps)**2 - gP/2*(kap - kaps)**2 + (gV + gA)/2*sum(eta[A, A]*(jL[A]**2 + jR[A]**2) for A in range(8))
                  + (gV - gA)*sum(eta[A, A]*jL[A]*jR[A] for A in range(8)))
R['E_H = gS/2 (k+k*)^2 - gP/2 (k-k*)^2 + (gV+gA)/2 (jL.jL + jR.jR) + (gV-gA) jL.jR'] = sp.expand(EHc16 - EHchi) == 0
VLl = [(gV + gA)*eta[A, A]*jL[A] + (gV - gA)*eta[A, A]*jR[A] for A in range(8)]      # V^L_A (lower)
VRl = [(gV + gA)*eta[A, A]*jR[A] + (gV - gA)*eta[A, A]*jL[A] for A in range(8)]
dEdk, dEdks = gS*(kap + kaps) - gP*(kap - kaps), gS*(kap + kaps) + gP*(kap - kaps)
R['V^L_A = V_A - A_A, V^R_A = V_A + A_A (Hartree, chiral currents)'] = all(
    sp.expand(VLl[A] - (gV*eta[A, A]*densC(beta*Gup[A]) - gA*eta[A, A]*densC(beta*Gup[A]*chi))) == 0
    and sp.expand(VRl[A] - (gV*eta[A, A]*densC(beta*Gup[A]) + gA*eta[A, A]*densC(beta*Gup[A]*chi))) == 0 for A in range(8))
R['M - m = dE_H/dkappa = gS(k+k*) - gP(k-k*) = S + iP'] = sp.expand(dEdk - (gS*densC(beta) + sp.I*gP*densC(beta*iChi))) == 0
lhsL = [sp.diff(EHchi, zLb8[a]) for a in range(8)]
rhsL = sum((VLl[A]*sbarA[A]*ZL for A in range(8)), sp.zeros(8, 1)) + dEdk*ZR
lhsR = [sp.diff(EHchi, zRb8[a]) for a in range(8)]
rhsR = sum((VRl[A]*sA[A]*ZR for A in range(8)), sp.zeros(8, 1)) + dEdks*ZL
R['dE_H/dphi_L* = (V^L_0 - sigma_i V^L_i) phi_L + (M - m) phi_R ; dE_H/dphi_R* = (V^R_0 + sigma_i V^R_i) phi_R + (M - m)^* phi_L'] = (
    all(sp.expand(lhsL[a] - rhsL[a]) == 0 for a in range(8)) and all(sp.expand(lhsR[a] - rhsR[a]) == 0 for a in range(8)))
# chiral U(1)_A:  psi -> e^{i theta chi} psi = (cos theta - sin theta omega) psi  (a REAL matrix), cos theta = 3/5, sin theta = 4/5
cth, sth = sp.Rational(3, 5), sp.Rational(4, 5)
Uth = cth*I16 + sp.I*sth*chi
e2 = (cth + sp.I*sth)**2
R['U_theta = cos(theta) - sin(theta) omega is real orthogonal; U E_L = e^{-i theta} E_L, U E_R = e^{i theta} E_R'] = (
    Uth == cth*I16 - sth*omega and Uth.T*Uth == I16 and zero(Uth*EL - (cth - sp.I*sth)*EL) and zero(Uth*ER - (cth + sp.I*sth)*ER))
R['U^dag beta U = beta (cos 2theta + i sin 2theta chi);  vector and axial currents invariant'] = (
    zero(dag(Uth)*beta*Uth - beta*(sp.re(e2)*I16 + sp.I*sp.im(e2)*chi))
    and all(zero(dag(Uth)*beta*Gup[A]*Uth - beta*Gup[A]) and zero(dag(Uth)*beta*Gup[A]*chi*Uth - beta*Gup[A]*chi) for A in range(8)))
EHrot = hartree(Uth*Zc, (Uth.conjugate())*ZBc)[0]
kap16, kaps16 = dens(EL*dag(ER)), dens(ER*dag(EL))
dEH = sp.expand(EHrot - EH)
R['E_H(e^{i theta chi} psi) - E_H(psi) = (gS-gP)/2 [(e^{4i theta}-1) k^2 + (e^{-4i theta}-1) k*^2]  ->  U(1)_A symmetric iff gP = gS'] = (
    sp.expand(dEH - (gS - gP)/2*((e2**2 - 1)*kap16**2 + (sp.conjugate(e2)**2 - 1)*kaps16**2)) == 0 and sp.expand(dEH.subs(gP, gS)) == 0 and dEH != 0)

# ---------------------------------------------------------------- 8. exchange (Fock) energy of the contact couplings; Fierz completeness
B256 = [Mx for S_, Mx in products(G)]                                            # the 256 ordered products Gamma_K
Bn = np.array([np.array(Mx.tolist(), dtype=np.int64) for Mx in B256])          # entries 0, +-1
R['256 products Gamma_K are orthogonal: Gamma_K^-1 = Gamma_K^T'] = bool(all((Bn[k].T @ Bn[k] == np.eye(16, dtype=np.int64)).all() for k in range(256)))
R['tr(Gamma_K^-1 Gamma_L) = 16 delta_KL (256 x 256)'] = bool((np.einsum('kab,lab->kl', Bn, Bn) == 16*np.eye(256, dtype=np.int64)).all())
C4 = np.einsum('kab,kdc->abcd', Bn, Bn)                                          # sum_K (Gamma_K)_ab (Gamma_K^-1)_cd
D4 = 16*np.einsum('ad,cb->abcd', np.eye(16, dtype=np.int64), np.eye(16, dtype=np.int64))
R['Fierz completeness: sum_K (Gamma_K)_ab (Gamma_K^-1)_cd = 16 delta_ad delta_cb'] = bool((C4 == D4).all())
chan = {'S': (gS, [(I16, I16)]), 'P': (gP, [(iChi, iChi)]), 'V': (gV, [(Gup[A], G[A]) for A in range(8)]),
        'A': (gA, [(Gup[A]*chi, G[A]*chi) for A in range(8)])}
Ex = lambda Rm: sp.expand(sum((-g/2*(O1*Rm*O2*Rm).trace() for g, prs in chan.values() for O1, O2 in prs), 0))
Sigx = lambda Rm: sp.expand(sum((-g*O2*Rm*O1 for g, prs in chan.values() for O1, O2 in prs), sp.zeros(16)))
EHR = lambda Rm: sp.expand(sum((g/2*(O1*Rm).trace()*(O2*Rm).trace() for g, prs in chan.values() for O1, O2 in prs), 0))
import random
random.seed(71)
R0 = sp.Matrix(16, 16, lambda a, b: random.randint(-3, 3) + sp.I*random.randint(-3, 3))
D0 = sp.Matrix(16, 16, lambda a, b: random.randint(-3, 3) + sp.I*random.randint(-3, 3))
tt = sp.symbols('t')
R['delta E_x = tr[Sigma_x delta R] (random R, dR), Sigma_x = -gS R - gP (i chi) R (i chi) - gV G_A R G^A - gA G_A chi R G^A chi'] = (
    sp.expand(sp.diff(Ex(R0 + tt*D0), tt).subs(tt, 0) - (Sigx(R0)*D0).trace()) == 0)
R['E_x = (1/2) tr[Sigma_x R]'] = sp.expand(Ex(R0) - (Sigx(R0)*R0).trace()/2) == 0
R['E_H = (gS/2)(tr R)^2 + ... = sum_c (g_c/2) tr(O_c R) tr(O^c R) and dE_H = tr[Sigma_H dR]'] = sp.expand(
    sp.diff(EHR(R0 + tt*D0), tt).subs(tt, 0) - sum(((g*(O2*R0).trace()*O1*D0).trace() for g, prs in chan.values() for O1, O2 in prs), 0)) == 0
# Fierz form of the exchange terms:  tr[O1 R O2 R] = (1/16) sum_K tr(Gamma_K^-1 R) tr(O1 Gamma_K O2 R)
R0n = np.array(R0.tolist(), dtype=complex)
def fierz_ok(O1, O2):
    O1n, O2n = np.array(O1.tolist(), dtype=complex), np.array(O2.tolist(), dtype=complex)
    lhs = np.trace(O1n @ R0n @ O2n @ R0n)
    rhs = sum(np.trace(Bn[k].T @ R0n)*np.trace(O1n @ Bn[k] @ O2n @ R0n) for k in range(256))/16
    return abs(lhs - rhs) < 1e-9
R['Fierz: tr[O1 R O2 R] = (1/16) sum_K tr(Gamma_K^-1 R) tr(O1 Gamma_K O2 R) for all channel vertices'] = all(
    fierz_ok(O1, O2) for g, prs in chan.values() for O1, O2 in prs)
# uniform static matter at rest:  R = (rhoS + rho Gamma_0)/16
rSu, rhou = sp.symbols('rho_S rho', real=True)
Ru = (rSu*I16 + rhou*G[0])/16
R['uniform R: tr R = rho_S, tr(Gamma^0 R) = rho'] = sp.simplify(Ru.trace() - rSu) == 0 and sp.simplify((Gup[0]*Ru).trace() - rhou) == 0
Rp = (US*US.T*G[0])/(2*En); Rmp = (US.subs({q: -q for q in P}, simultaneous=True)*US.subs({q: -q for q in P}, simultaneous=True).T*G[0])/(2*En)
R['(1/2) sum_s [u_s(p) ubar_s(p) + u_s(-p) ubar_s(-p)]/(2E) = (E Gamma_0 + m)/(2E)'] = zero((Rp + Rmp)/2 - (En*G[0] + m*I16)/(2*En))
ExU = Ex(Ru)
ExU_closed = -gS/32*(rSu**2 + rhou**2) + gP/32*(rSu**2 - rhou**2) - gV/32*(8*rSu**2 - 6*rhou**2) + gA/32*(8*rSu**2 + 6*rhou**2)
R['E_x(uniform) = -gS/32 (rS^2+rho^2) + gP/32 (rS^2-rho^2) - gV/32 (8rS^2-6rho^2) + gA/32 (8rS^2+6rho^2)'] = sp.expand(ExU - ExU_closed) == 0
Sx_c = rSu/16*(-gS + gP - 8*gV + 8*gA); Vx_c = rhou/16*(-gS - gP + 6*gV + 6*gA)
R['Sigma_x(uniform) = S_x + V_x Gamma_0, S_x = rS/16(-gS+gP-8gV+8gA), V_x = rho/16(-gS-gP+6gV+6gA)'] = zero(Sigx(Ru) - Sx_c*I16 - Vx_c*G[0])
R['S_x = dE_x/drho_S, V_x = dE_x/drho'] = sp.expand(sp.diff(ExU_closed, rSu) - Sx_c) == 0 and sp.expand(sp.diff(ExU_closed, rhou) - Vx_c) == 0
GSe = gS + (-gS + gP - 8*gV + 8*gA)/16; GVe = gV + (-gS - gP + 6*gV + 6*gA)/16
R['uniform Hartree-Fock = Hartree with G_S = gS + (-gS+gP-8gV+8gA)/16, G_V = gV + (-gS-gP+6gV+6gA)/16'] = sp.expand(
    gS/2*rSu**2 + gV/2*rhou**2 + ExU_closed - GSe/2*rSu**2 - GVe/2*rhou**2) == 0
R['uniform matter: Hartree terms of P and A vanish (rho_P = 0, j5 = 0, j^i = 0)'] = (sp.simplify((G[0]*iChi*Ru).trace()) == 0
    and all(sp.simplify((Gup[A]*chi*Ru).trace()) == 0 for A in range(8)) and all(sp.simplify((Gup[i]*Ru).trace()) == 0 for i in range(1, 8)))
# Wick's theorem (Hartree + exchange) checked in the 3-particle Fock space of the 16 modes psi_a at one point:
#   <Phi| :(psi^dag G0 O psi)(psi^dag G0 O' psi): |Phi> = tr(O R) tr(O' R) - tr(O R O' R),  R = sum_k phi_k phibar_k
basis3 = list(itertools.combinations(range(16), 3)); idx3 = {S_: n for n, S_ in enumerate(basis3)}
def cdc(a, b):                    # matrix of c_a^dag c_b on the 3-particle sector
    Mx = np.zeros((len(basis3), len(basis3)), dtype=complex)
    for n, S_ in enumerate(basis3):
        if b not in S_: continue
        sgn = (-1)**S_.index(b); T_ = [s for s in S_ if s != b]
        if a in T_: continue
        sgn *= (-1)**sum(1 for s in T_ if s < a); T2 = tuple(sorted(T_ + [a]))
        Mx[idx3[T2], n] += sgn
    return Mx
CDC = {(a, b): cdc(a, b) for a in range(16) for b in range(16)}
Bop = lambda Mx: sum(Mx[a, b]*CDC[(a, b)] for a in range(16) for b in range(16) if Mx[a, b] != 0)
rng = np.random.default_rng(71)
Q_, _ = np.linalg.qr(rng.normal(size=(16, 3)) + 1j*rng.normal(size=(16, 3)))       # three orthonormal orbitals phi_k (columns)
Phi = np.array([np.linalg.det(Q_[list(S_), :].T) for S_ in basis3])                  # amplitudes det[phi_k(s_j)]
Rn = sum(np.outer(Q_[:, k], Q_[:, k].conj()) for k in range(3)) @ np.array(G[0].tolist(), dtype=complex)   # R = sum phi phi^dag Gamma0
okW = abs(np.vdot(Phi, Phi) - 1) < 1e-12
for g, prs in chan.values():
    for O1, O2 in prs:
        O1n, O2n = np.array(O1.tolist(), dtype=complex), np.array(O2.tolist(), dtype=complex)
        M1, M2 = np.array((G[0]*O1).tolist(), dtype=complex), np.array((G[0]*O2).tolist(), dtype=complex)
        lhs = np.vdot(Phi, Bop(M1) @ (Bop(M2) @ Phi)) - np.vdot(Phi, Bop(M1 @ M2) @ Phi)       # normal ordering removes B_{M1 M2}
        rhs = np.trace(O1n @ Rn)*np.trace(O2n @ Rn) - np.trace(O1n @ Rn @ O2n @ Rn)
        okW = okW and abs(lhs - rhs) < 1e-10
R['Wick (3 particles, 16 modes): <:(psibar O psi)(psibar O\' psi):> = tr(O R) tr(O\' R) - tr(O R O\' R) for all 18 vertex pairs'] = bool(okW)

# ---------------------------------------------------------------- 9. uniform matter, no-sea Hartree and Hartree-Fock (numbers, m = 1)
mp.mp.dps = 30
pi_ = mp.pi
c7 = 1/(15*pi_**4)                                     # 8 int d^7k/(2pi)^7 f(|k|) = c7 int_0^inf k^6 f(k) dk
rhoF = lambda k: k**7/(105*pi_**4)
rhoSF = lambda ms, k: c7*mp.quad(lambda q: q**6*ms/mp.sqrt(q*q + ms*ms), [0, k])
ekinF = lambda ms, k: c7*mp.quad(lambda q: q**6*mp.sqrt(q*q + ms*ms), [0, k])
m1 = mp.mpf(1)
R['c7 = 8 |S^6|/(2pi)^7 = 1/(15 pi^4) and rho = c7 kF^7/7 = kF^7/(105 pi^4)'] = (abs(8*(2*pi_**mp.mpf(3.5)/mp.gamma(mp.mpf(3.5)))/(2*pi_)**7 - c7) < mp.mpf(10)**-28
    and abs(c7*mp.quad(lambda q: q**6, [0, mp.mpf('0.7')]) - rhoF(mp.mpf('0.7'))) < mp.mpf(10)**-28)
def gapF(GS, k):                                       # m* = m + GS rhoS(m*)  (unique root in (0, m) for GS < 0)
    return mp.findroot(lambda ms: ms - m1 - GS*rhoSF(ms, k), (mp.mpf(10)**-8, m1), solver='anderson')
def epsF(k, GS, GV):                                   # energy density of the self-consistent state
    ms = gapF(GS, k)
    return ekinF(ms, k) - (ms - m1)**2/(2*GS) + GV/2*rhoF(k)**2
def fermiF(k, GS, GV):                                 # Fermi energy  V0 + sqrt(kF^2 + m*^2)
    return GV*rhoF(k) + mp.sqrt(k*k + gapF(GS, k)**2)
gS1, gV1 = mp.mpf(-100000), mp.mpf(10000)
kF0 = mp.findroot(lambda k: epsF(k, gS1, gV1)/rhoF(k) - fermiF(k, gS1, gV1), mp.mpf('0.75'))
ms0 = gapF(gS1, kF0); rho0 = rhoF(kF0); rS0 = rhoSF(ms0, kF0); V00 = gV1*rho0; eps0 = epsF(kF0, gS1, gV1)
R['NUM Hartree saturation kF0'] = kF0; R['NUM Hartree rho0'] = rho0; R['NUM Hartree m*'] = ms0; R['NUM Hartree rhoS'] = rS0
R['NUM Hartree V0'] = V00; R['NUM Hartree E/A - m'] = eps0/rho0 - 1; R['NUM Hartree eps_F'] = fermiF(kF0, gS1, gV1)
R['Hartree gap: m* - m - gS rhoS(m*) = 0'] = abs(ms0 - 1 - gS1*rS0) < mp.mpf(10)**-25
R['m* is a stationary point of eps(m*) = ekin(m*) - (m*-m)^2/(2gS) + gV rho^2/2'] = abs(mp.diff(lambda x: ekinF(x, kF0) - (x - 1)**2/(2*gS1), ms0)) < mp.mpf(10)**-20
sumeps = V00*rho0 + ekinF(ms0, kF0)                    # sum of occupied KS eigenvalues  (per volume)
R['total energy = sum eps_k - E_H (double counting) = ekin - (m*-m) rhoS + E_H'] = (abs(eps0 - (sumeps - gS1/2*rS0**2 - gV1/2*rho0**2)) < mp.mpf(10)**-25
    and abs(eps0 - (ekinF(ms0, kF0) - (ms0 - 1)*rS0 + gS1/2*rS0**2 + gV1/2*rho0**2)) < mp.mpf(10)**-25)
hh = mp.mpf(10)**-8
dEdrho = lambda k, GS, GV: (epsF(k + hh, GS, GV) - epsF(k - hh, GS, GV))/(rhoF(k + hh) - rhoF(k - hh))
R['Hugenholtz-Van Hove at saturation: E/A = eps_F = d eps/d rho'] = (abs(eps0/rho0 - fermiF(kF0, gS1, gV1)) < mp.mpf(10)**-25
    and abs(dEdrho(kF0, gS1, gV1) - fermiF(kF0, gS1, gV1)) < mp.mpf(10)**-12)
R['d(eps/rho)/dkF = 0 at kF0 (minimum: E/A(0.95 kF0), E/A(1.05 kF0) larger)'] = (abs((epsF(kF0 + hh, gS1, gV1)/rhoF(kF0 + hh) - epsF(kF0 - hh, gS1, gV1)/rhoF(kF0 - hh))/(2*hh)) < mp.mpf(10)**-12
    and epsF(kF0*mp.mpf('0.95'), gS1, gV1)/rhoF(kF0*mp.mpf('0.95')) > eps0/rho0 and epsF(kF0*mp.mpf('1.05'), gS1, gV1)/rhoF(kF0*mp.mpf('1.05')) > eps0/rho0)
k1 = mp.mpf(1)
R['Hugenholtz-Van Hove mu = d eps/d rho = eps_F at kF = 1 (not saturated)'] = abs(dEdrho(k1, gS1, gV1) - fermiF(k1, gS1, gV1)) < mp.mpf(10)**-12
R['NUM Hartree kF=1: m*'] = gapF(gS1, k1); R['NUM Hartree kF=1: E/A - m'] = epsF(k1, gS1, gV1)/rhoF(k1) - 1; R['NUM Hartree kF=1: eps_F'] = fermiF(k1, gS1, gV1)
# Hartree-Fock (gP = gA = 0): G_S = (15/16) gS - gV/2,  G_V = (11/8) gV - gS/16, at kF = 3/4
GS1 = gS1 + (-gS1 - 8*gV1)/16; GV1 = gV1 + (-gS1 + 6*gV1)/16
R['NUM HF couplings G_S, G_V'] = (GS1, GV1)
k34 = mp.mpf(3)/4
msHF = gapF(GS1, k34); rhoHF = rhoF(k34); rSHF = rhoSF(msHF, k34)
SxHF, VxHF = rSHF/16*(-gS1 - 8*gV1), rhoHF/16*(-gS1 + 6*gV1)
R['HF self-consistency with the exchange self-energy: m* - m = gS rhoS + S_x, V0 = gV rho + V_x'] = (abs(msHF - 1 - (gS1*rSHF + SxHF)) < mp.mpf(10)**-25
    and abs(GV1*rhoHF - (gV1*rhoHF + VxHF)) < mp.mpf(10)**-25)
ExHF = -gS1/32*(rSHF**2 + rhoHF**2) - gV1/32*(8*rSHF**2 - 6*rhoHF**2)
R['HF energy: ekin - (m*-m) rhoS + E_H + E_x = eps(G_S, G_V)'] = abs(ekinF(msHF, k34) - (msHF - 1)*rSHF + gS1/2*rSHF**2 + gV1/2*rhoHF**2 + ExHF - epsF(k34, GS1, GV1)) < mp.mpf(10)**-25
R['HF Hugenholtz-Van Hove: d eps/d rho = eps_F at kF = 3/4'] = abs(dEdrho(k34, GS1, GV1) - fermiF(k34, GS1, GV1)) < mp.mpf(10)**-12
R['NUM HF kF=3/4: m*'] = msHF; R['NUM HF kF=3/4: V0'] = GV1*rhoHF; R['NUM HF kF=3/4: E/A - m'] = epsF(k34, GS1, GV1)/rhoHF - 1
R['NUM Hartree kF=3/4: m*'] = gapF(gS1, k34); R['NUM Hartree kF=3/4: E/A - m'] = epsF(k34, gS1, gV1)/rhoHF - 1

# ---------------------------------------------------------------- 10. chiral KS with the Dirac sea: NJL vacuum (m = 0 and m = 1/100), cutoff Lambda = 1
Lam = mp.mpf(1)
Jn = lambda x: mp.quad(lambda q: q**6/mp.sqrt(q*q + x*x), [0, Lam])/(120*pi_**4)      # int_{|p|<Lambda} d^7p/(2pi)^7 1/E
Kn = lambda x: mp.quad(lambda q: q**6*mp.sqrt(q*q + x*x), [0, Lam])/(120*pi_**4)      # int_{|p|<Lambda} d^7p/(2pi)^7 E
gc = -90*pi_**4/Lam**6
gN = mp.mpf(6)/5*gc
R['J(0) = Lambda^6/(720 pi^4), g_c = -1/(8 J(0)) = -90 pi^4/Lambda^6'] = abs(Jn(0) - Lam**6/(720*pi_**4)) < mp.mpf(10)**-28 and abs(gc + 1/(8*Jn(0))) < mp.mpf(10)**-22
Etr = lambda Mz, mm: -8*Kn(abs(Mz)) + 8*Jn(abs(Mz))*mp.re((Mz - mm)*mp.conj(Mz)) + 32*gN*abs(Mz)**2*Jn(abs(Mz))**2
x0 = mp.findroot(lambda x: 1 + 8*gN*Jn(x), mp.mpf('0.4'))
R['NUM NJL g_c'] = gc; R['NUM NJL g = 6/5 g_c'] = gN; R['NUM NJL |M| (m = 0)'] = x0
R['NUM NJL kappa = -4 M J(|M|) (M real)'] = -4*x0*Jn(x0); R['NUM NJL E(|M|) - E(0)'] = Etr(x0, 0) - Etr(0, 0)
R['NJL m = 0: gap 1 + 8 g J(|M|) = 0 solved, M = 2 g kappa^*, stationary and lower than M = 0'] = (abs(1 + 8*gN*Jn(x0)) < mp.mpf(10)**-25
    and abs(x0 - 2*gN*(-4*x0*Jn(x0))) < mp.mpf(10)**-25 and abs(mp.diff(lambda x: Etr(x, 0), x0)) < mp.mpf(10)**-20 and Etr(x0, 0) < Etr(0, 0))
R['NJL m = 0: E(|M| e^{i theta}) independent of theta (degenerate chiral vacua)'] = all(abs(Etr(x0*mp.expjpi(th), 0) - Etr(x0, 0)) < mp.mpf(10)**-25 for th in [mp.mpf(1)/7, mp.mpf(1)/2, mp.mpf(5)/3])
mq = mp.mpf(1)/100
xm = mp.findroot(lambda x: x - mq + 8*gN*Jn(x)*x, mp.mpf('0.4'))
R['NUM NJL M (m = 1/100)'] = xm
R['NJL m = 1/100: M = m - 8 g J(M) M real, stationary; E(M e^{i theta}) > E(M) for theta != 0'] = (abs(mp.diff(lambda x: Etr(x, mq), xm)) < mp.mpf(10)**-20
    and all(Etr(xm*mp.expjpi(th), mq) > Etr(xm, mq) for th in [mp.mpf(1)/50, mp.mpf(1)/2, mp.mpf(1)]))
# one momentum of the sea, numerically: h_chi = [[-sigma.p, M],[M^*, sigma.p]], sea projector, kappa_p = -4 M^*/E,  sum eps_sea = -8E
pv = [0.3, -0.2, 0.1, 0.25, 0.0, -0.15, 0.05]
Mz = complex(x0)*np.exp(0.7j)
sp_ = sum(pv[i-1]*np.array(sig[i].tolist(), dtype=complex) for i in range(1, 8))
hc = np.block([[-sp_, Mz*np.eye(8)], [np.conj(Mz)*np.eye(8), sp_]])
ev_, evec = np.linalg.eigh(hc)
Ep = np.sqrt(sum(q*q for q in pv) + abs(Mz)**2)
Pim = evec[:, :8] @ evec[:, :8].conj().T
R['sea at one momentum: eigenvalues -E (x8), +E (x8); kappa_p = tr[(Pi_-)_{RL}] = -4 M^*/E; tr[Pi_- h_0] = -8E + 8 Re[(M-m) M^*]/E'] = bool(
    np.allclose(ev_, [-Ep]*8 + [Ep]*8) and abs(np.trace(Pim[8:, :8]) + 4*np.conj(Mz)/Ep) < 1e-12
    and abs(np.trace(Pim @ np.block([[-sp_, np.zeros((8, 8))], [np.zeros((8, 8)), sp_]])) - (-8*Ep + 8*abs(Mz)**2/Ep)) < 1e-12)

if __name__ == '__main__':
    bad = [k for k, v in R.items() if isinstance(v, (bool, np.bool_)) and not v]
    nbool = sum(isinstance(v, (bool, np.bool_)) for v in R.values())
    for k, v in R.items():
        if isinstance(v, tuple): v = ', '.join(mp.nstr(x, 20) if isinstance(x, mp.mpf) else str(x) for x in v)
        elif isinstance(v, mp.mpf): v = mp.nstr(v, 20)
        print(f'{k} : {v}')
    print(f'boolean checks: {nbool}   failed: {len(bad)}')
    print('VERDICT:', 'ALL CHECKS PASSED' if not bad else 'FAILED: ' + '; '.join(bad))
    sys.exit(1 if bad else 0)
