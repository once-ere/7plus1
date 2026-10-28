"""Numbers shared by the documents RADION*, the radion parts of WORLD_BRANE* and AUXILIARY_TORUS*, and the README parts of section 26 of
7+1.nb (the radion: the size of the four extra dimensions, on the 4-torus and on the 4-sphere): the output keys of section 26 of the
evaluated notebook and of tools/ref_radion.py, and the comparison of the two programs.  Nothing that can be computed is typed by hand."""
from doc_stage3_common import fmt, nums_in

NB = {  # short name -> output key of section 26 of 7+1.nb
    'kin': 'general n: kinetic coefficient n(n+2)/2 for n = 1..4',
    'beta2': "the note's beta_n^2 = 2/(n(n+2)) and the canonical 1/(n(n+2)), n = 4",
    'sphU': 'sphere with U: Lambda + kappa U0, kappa U1, m_r^2 - kappa K_U/12 (units 1/r0^2)',
    'sphL': 'sphere with Lambda only: Lambda r0^2, H^2 r0^2, m^2/H^2',
    'sphF': 'sphere with a 4-form flux: Lambda r0^2, kappa rho_F r0^2, m^2 r0^2',
    'ew': 'Eot-Wash bound on the radion mass (eV), range of a 10 GeV radion (m), gamma = v/(Sqrt[3] M_Pl)',
    'dl': 'S^4 multiplicities d_l, l = 0..8',
    'weyl': 'S^4 Weyl: N(x)/(x^4/12) - 1 at x = mu r = 10, 100, 1000',
    'dirac': 'Dirac on S^4: degeneracy 4 Binomial[l + 3, 3] of +i(l + 2)/r, l = 0..5; Weyl ratio - 1 at L = 1000',
    'today': 'today: r_S4 (GeV^-1, m), lowest masses (MeV): torus 2 Pi/l_perp, sphere scalar 2/r, sphere dirac16complex 2/r; levels below '
             '5 GeV: torus lattice, sphere exact, Weyl',
    'rs1': "RS1: coefficient of A'' in G^0_0, f_UV kappa/k, f_IR kappa/k, f^2 kappa k (canonical radion)",
    'strength': 'strength of the radion coupling: warped with <mu> = 3 TeV over 7+1, i.e. Sqrt[3] M_Pl/(3 TeV)',
    'gw': "Goldberger-Wise: m_r^2/(kappa k^2 eps^2 v_B^2 (<mu>/f)^2) with f^2 = 6/(kappa k) and with the note's f^2 = 12/(kappa k)",
    'hrk': 'Higgs-radion mixing: c in Sqrt|h| R(h) = c Omega box Omega, K_hr/(xi gamma), (K_rr - 1)/(xi gamma^2)',
    'hr71': 'Higgs-radion mixing in 7+1 (gamma = v/(Sqrt[3] M_Pl), xi = 1, m_r0 = 10 GeV): Z - 1, theta, Gamma(h -> gamma gamma)/SM - 1, '
            'Gamma(h -> g g)/SM - 1; warped gamma = v/(3 TeV)'}
PY = {  # short name -> output key of tools/ref_radion.py
    'kin': 'NUM general n: kinetic coefficient n(n+2)/2 for n = 1..4',
    'beta2': "NUM the note's beta_n^2 = 2/(n(n+2)) and the canonical 1/(n(n+2)), n = 4",
    'sphU': 'NUM sphere with U: Lambda + kappa U0, kappa U1, m_r^2 - kappa K_U/12 (units 1/r0^2)',
    'sphL': 'NUM sphere with Lambda only: Lambda r0^2, H^2 r0^2, m^2/H^2',
    'sphF': 'NUM sphere with a 4-form flux: Lambda r0^2, kappa rho_F r0^2, m^2 r0^2',
    'ew': 'NUM Eot-Wash bound on the radion mass (eV), range of a 10 GeV radion (m), gamma = v/(sqrt(3) M_Pl)',
    'dl': 'NUM S^4 multiplicities d_l, l = 0..8',
    'weyl': 'NUM S^4 Weyl: N(x)/(x^4/12) - 1 at x = mu r = 10, 100, 1000',
    'dirac': 'NUM Dirac on S^4: degeneracy 4 binomial(l + 3, 3) of +i(l + 2)/r, l = 0..5; Weyl ratio - 1 at L = 1000',
    'today': 'NUM today: r_S4 (GeV^-1, m), lowest masses (MeV): torus 2 pi/l_perp, sphere scalar 2/r, sphere dirac16complex 2/r; levels '
             'below 5 GeV: torus lattice, sphere exact, Weyl',
    'rs1': "NUM RS1: coefficient of A'' in G^0_0, f_UV kappa/k, f_IR kappa/k, f^2 kappa k (canonical radion)",
    'strength': 'NUM strength of the radion coupling: warped with <mu> = 3 TeV over 7+1, i.e. sqrt(3) M_Pl/(3 TeV)',
    'gw': "NUM Goldberger-Wise: m_r^2/(kappa k^2 eps^2 v_B^2 (<mu>/f)^2) with f^2 = 6/(kappa k) and with the note's f^2 = 12/(kappa k)",
    'hrk': 'NUM Higgs-radion mixing: c in sqrt|h| R(h) = c Omega box Omega, K_hr/(xi gamma), (K_rr - 1)/(xi gamma^2)',
    'hr71': 'NUM Higgs-radion mixing in 7+1 (gamma = v/(sqrt(3) M_Pl), xi = 1, m_r0 = 10 GeV): Z - 1, theta, Gamma(h -> gamma gamma)/SM - 1, '
            'Gamma(h -> g g)/SM - 1; warped gamma = v/(3 TeV)'}

# the groups compared between the two programs: (label, key, indices or None for all, tolerance)
GROUPS = [
    ('kinetic coefficient $n(n+2)/2$, $n=1,\\dots,4$ (solved from the reduction)', 'kin', None, 0.0),
    ('$\\beta_n^2$ of the note and the canonical $\\beta_n^2$, $n=4$', 'beta2', None, 1e-15),
    ('$S^4$ with $U$: $(\\Lambda+\\kappa U_0)r_0^2$, $\\kappa U_1r_0^2$, $(m_r^2-\\kappa K_U/12)r_0^2$', 'sphU', None, 0.0),
    ('$S^4$ with $\\Lambda$ only: $\\Lambda r_0^2$, $H^2r_0^2$, $m^2/H^2$', 'sphL', None, 0.0),
    ('$S^4$ with a 4-form flux: $\\Lambda r_0^2$, $\\kappa\\rho_Fr_0^2$, $m^2r_0^2$', 'sphF', None, 1e-15),
    ('E\\"ot-Wash bound, range of a 10 GeV radion, $\\gamma=v/(\\sqrt3M_{\\rm Pl})$', 'ew', None, 1e-10),
    ('multiplicities $d_l$ on $S^4$, $l=0,\\dots,8$', 'dl', None, 0.0),
    ("Weyl's law on $S^4$: $N/(x^4/12)-1$ at $x=10,100,1000$", 'weyl', None, 1e-9),
    ('Dirac degeneracies on $S^4$ and the spinor Weyl ratio', 'dirac', None, 1e-9),
    ('today: $r_{S^4}$, lowest masses, numbers of levels below 5 GeV', 'today', None, 1e-12),
    ('RS1: $3$, $f_{\\rm UV}\\kappa/k$, $f_{\\rm IR}\\kappa/k$, $f^2\\kappa k$', 'rs1', None, 0.0),
    ('strength of the warped coupling over the flat one', 'strength', None, 1e-12),
    ('Goldberger--Wise mass coefficient: $\\frac43$ and the note\'s $\\frac23$', 'gw', None, 1e-15),
    ('Higgs--radion: $c=-6$, $K_{hr}/(\\xi\\gamma)$, $(K_{rr}-1)/(\\xi\\gamma^2)$', 'hrk', None, 0.0),
    ('Higgs--radion mixing in 7+1: $Z-1$, $\\theta$, rate changes; warped $\\gamma$', 'hr71', None, 1e-12)]


def nbv(D, k):
    return nums_in(D['nums'].get('rd', {}).get(NB[k], ''))


def pyv(D, k):
    return nums_in(D['ref'].get('rd', ('', '', {}))[2].get(PY[k], ''))


def max_rel(a, b):
    if not a or not b or len(a) != len(b) or any(x is None for x in list(a) + list(b)):
        return None
    return max(abs(x - y)/max(abs(x), abs(y), 1e-300) for x, y in zip(a, b))


def compare(D):
    """(label, notebook values, Python values, tolerance, max relative difference) for every group printed by both programs"""
    out = []
    for lab, key, idx, tol in GROUPS:
        a, b = nbv(D, key), pyv(D, key)
        if idx is not None:
            a = [a[i] for i in idx if i < len(a)]; b = [b[i] for i in idx if i < len(b)]
        out.append((lab, a, b, tol, max_rel(a, b)))
    return out


def worst(D):
    return max([r for *_, r in compare(D) if r is not None], default=None)


def counts(D):
    """(checks of section 26, of them True, boolean checks of tools/ref_radion.py, of them True)"""
    c = D['checks'].get('rd', [])
    r = [(k, v) for k, v in D.get('ref', {}).get('rd', ('', '', {}))[2].items() if v in ('True', 'False')]
    return len(c), sum(1 for x in c if x[1] == 'True'), len(r), sum(1 for k, v in r if v == 'True')


def g(v, i, n=6, default='?'):
    try:
        return fmt(v[i], n)
    except (IndexError, TypeError):
        return default


def val(D, key, i, n=6, default='?'):
    """the i-th number of the notebook group key, formatted with n significant digits"""
    return g(nbv(D, key), i, n, default)


def ival(D, key, i, default='?'):
    v = nbv(D, key)
    try:
        return '%d' % round(v[i])
    except (IndexError, TypeError):
        return default
