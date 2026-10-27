"""The curvature and Lovelock engine of tools/ref_gravity.py and tools/ref_ks_gravity.py (7+1 spacetime, D = 8).

geometry(g)        inverse metric, Christoffel symbols chr[l][m][n] = {l; m n}, Riemann R^r_{s m n} (Landau-Lifshitz convention)
ricci(Riem, gi)    Ricci tensor R_{s n} = R^r_{s r n}, scalar curvature, mixed Einstein tensor G^m_n
rud(Riem, gi)      R^{c d}_{a b} = g^{d s} R^c_{s a b}
lovelockE(Rud, k, mu, nu), lovelockL(Rud, k)   Lovelock tensor E^(k)mu_nu and density L_k by the pair-set algorithm
cosmological_closed_forms()   E^(k) (00, 11, 44), k = 1, 2, 3, of ds^2 = dx0^2 - a^2 dx_(1..3)^2 - b^2 dx_(4..7)^2 in Ha, Hb, Ha', Hb'
"""
import itertools, math
import sympy as sp

D = 8
X = (sp.Symbol('x0', positive=True),) + sp.symbols('x1:8', real=True)

def geometry(g, simp=sp.simplify):
    """inverse metric, Christoffel symbols chr[l][m][n], Riemann tensor R^r_{s m n} (dict of non-zero entries)"""
    gi = g.inv().applyfunc(simp)
    chr_ = [[[simp(sum(gi[l, r]*(sp.diff(g[r, n], X[m]) + sp.diff(g[r, m], X[n]) - sp.diff(g[m, n], X[r]))
                       for r in range(D))/2) for n in range(D)] for m in range(D)] for l in range(D)]
    Riem = {}
    for r in range(D):
        for s in range(D):
            for m in range(D):
                for n in range(m + 1, D):
                    v = simp(sp.diff(chr_[r][s][n], X[m]) - sp.diff(chr_[r][s][m], X[n])
                             + sum(chr_[r][m][l]*chr_[l][s][n] - chr_[r][n][l]*chr_[l][s][m] for l in range(D)))
                    if v != 0:
                        Riem[(r, s, m, n)] = v; Riem[(r, s, n, m)] = -v
    return gi, chr_, Riem

def ricci(Riem, gi, simp=sp.simplify):
    Ric = sp.Matrix(D, D, lambda s, n: simp(sum(Riem.get((r, s, r, n), 0) for r in range(D))))
    Rs = simp(sum(gi[s, n]*Ric[s, n] for s in range(D) for n in range(D)))
    Gmix = (gi*Ric - sp.eye(D)*Rs/2).applyfunc(simp)
    return Ric, Rs, Gmix

def rud(Riem, gi, simp=sp.simplify):
    """R^{c d}_{a b} = g^{d s} R^c_{s a b}"""
    out = {}
    for (c, s, a_, b_), v in Riem.items():
        for d in range(D):
            if gi[d, s] != 0:
                out[(c, d, a_, b_)] = out.get((c, d, a_, b_), 0) + gi[d, s]*v
    out = {k: simp(v) for k, v in out.items()}
    return {k: v for k, v in out.items() if v != 0}

# ======================================================================= Lovelock tensors (pair-set algorithm)
def perm_sign(up, dn):
    """delta^{up}_{dn} for dn a permutation of up"""
    pos = {v: i for i, v in enumerate(up)}
    p = [pos[v] for v in dn]; sgn = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; sgn = -sgn
    return sgn

def pair_sets(idx, k):
    """all sets of k disjoint pairs (i < j) taken from the sorted list idx"""
    if k == 0:
        yield (); return
    if len(idx) < 2*k:
        return
    f = idx[0]
    for j in range(1, len(idx)):
        for ps in pair_sets(idx[1:j] + idx[j + 1:], k - 1):
            yield ((f, idx[j]),) + ps
    for ps in pair_sets(idx[1:], k):
        yield ps

def ordered_seqs(idx, k):
    for ps in pair_sets(sorted(idx), k):
        for perm in itertools.permutations(ps):
            yield perm

def lovelockE(Rud, k, mu, nu):
    """E^(k) mu_nu = -(1/2^(k+1)) delta^{mu a1 b1 .. ak bk}_{nu c1 d1 .. ck dk} R^{c1 d1}_{a1 b1} ... R^{ck dk}_{ak bk}
       = -(2^(k-1) k!) sum_{upper pair-sets U} sum_{lower ordered pair-sequences L} sign * prod_i R^{L_i}_{U_i}"""
    if k == 0:
        return -sp.Rational(1, 2) if mu == nu else 0
    tot = 0
    for U in pair_sets([i for i in range(D) if i != mu], k):
        up = (mu,) + tuple(x for pr in U for x in pr)
        if nu not in up: continue
        for L in ordered_seqs(sorted(set(up) - {nu}), k):
            dn = (nu,) + tuple(x for pr in L for x in pr)
            term = perm_sign(up, dn)
            for i in range(k):
                v = Rud.get((L[i][0], L[i][1], U[i][0], U[i][1]), 0)
                if v == 0:
                    term = 0; break
                term = term*v
            tot += term
    return -(2**(k - 1))*math.factorial(k)*tot

def lovelockL(Rud, k):
    """L_k = 2^-k delta^{a1 b1 ..}_{c1 d1 ..} R^{c1 d1}_{a1 b1} ... = 2^k k! sum_U sum_L sign prod R"""
    if k == 0:
        return 1
    tot = 0
    for U in pair_sets(list(range(D)), k):
        up = tuple(x for pr in U for x in pr)
        for L in ordered_seqs(list(up), k):
            dn = tuple(x for pr in L for x in pr)
            term = perm_sign(up, dn)
            for i in range(k):
                v = Rud.get((L[i][0], L[i][1], U[i][0], U[i][1]), 0)
                if v == 0:
                    term = 0; break
                term = term*v
            tot += term
    return 2**k*math.factorial(k)*tot


def cosmological_closed_forms():
    t = X[0]
    a = sp.Function('a', real=True, positive=True)(t); b = sp.Function('b', real=True, positive=True)(t)
    g = sp.diag(1, -a**2, -a**2, -a**2, -b**2, -b**2, -b**2, -b**2)
    gi, chr_, Riem = geometry(g)
    Rd = rud(Riem, gi)
    Ha, Hb, dHa, dHb = sp.symbols('H_a H_b dH_a dH_b', real=True)
    def hub(e):
        e = e.subs({sp.Derivative(a, (t, 2)): (dHa + Ha**2)*a, sp.Derivative(b, (t, 2)): (dHb + Hb**2)*b})
        e = e.subs({sp.Derivative(a, t): Ha*a, sp.Derivative(b, t): Hb*b})
        return sp.expand(sp.simplify(e))
    return {(k, c): hub(lovelockE(Rd, k, c, c)) for k in (1, 2, 3) for c in (0, 1, 4)}
