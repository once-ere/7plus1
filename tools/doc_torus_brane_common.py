"""Numbers shared by the documents AUXILIARY_TORUS*, WORLD_BRANE* and the README parts of section 25 of 7+1.nb (the auxiliary torus and
the world brane): the output keys of section 25 of the evaluated notebook and of tools/ref_torus_brane.py, the comparison of the two
programs, and the numbers of the earlier sections (17-24) that the documents quote.  Nothing that can be computed is typed by hand."""
from doc_stage3_common import fmt, nums_in, mnum

NB = {  # short name -> output key of section 25 of 7+1.nb
    'orbits': 'orbit sizes of the shells k = 1..12',
    'r4': 'r_4(k), k = 0..16',
    'NR': 'N(R) for R = 10, 30, 100, 300, 1000',
    'dev': 'N(R)/(Pi^2 R^4/2) - 1 for R = 10, 30, 100, 300, 1000',
    'today': 'V_perp (GeV^-4), l_perp (GeV^-1), l_perp (m), l_perp (fm), 2 Pi/l_perp (MeV)',
    'scales': 'l_perp/proton radius, Bohr radius/l_perp, R = l_perp/(2 Pi) (GeV^-1, m), M_D = (2 Pi)^(2/3) M8 (GeV)',
    'nmu': 'N(mu) = V_perp mu^4/(32 Pi^2) for mu = 5 GeV, m_h, M8; lattice count below 5 GeV',
    'levels': 'first levels: k, r_4(k), Dirac fermions 4 r_4(k), mass 2 Pi Sqrt[k]/l_perp (MeV), k = 1..8',
    'pois': 'images S(r), Kaluza-Klein (4 Pi^2/(3 r)) Sum Exp[-2 Pi r |m|], relative difference, at r = 1/4, 1/2, 1, 2',
    'z52': 'Z(5/2) = 7 zeta(5/2) zeta(3/2), shell sum with tail',
    'nearfar': 'near field (r^5 S - 1)/(r^5 Z(5/2)) at r = 1/20, 1/10; far field (3 r S/(4 Pi^2) - 1)/(8 Exp[-2 Pi r]) at r = 1, 2, 3',
    'chr': 'Christoffel {4 00} at the brane (inhomogeneous metric, x0 = 1/2)',
    'gperm': 'Gamma17[A] as signed permutations: column j has the entry sign(c) in row |c|, c listed for A = 0..7',
    'reU': 'Re(2 Sqrt[2] U)', 'imU': 'Im(2 Sqrt[2] U)',
    'spinor': 'example psi = U.(u (x) e_1), u = (1, 2, 0, -1): psibar psi, psibar Gamma^mu psi (mu = 0..3)',
    'tstar': 'T* = (M8^6/M_Pl)^(1/5) (GeV), f_W/rho_Lambda for f_W = M8^4, 1/(2 Pi M8) (GeV^-1)'}
PY = {  # short name -> output key of tools/ref_torus_brane.py
    'orbits': 'NUM orbit sizes of the shells k = 1..12',
    'r4': 'NUM r_4(k), k = 0..16',
    'NR': 'NUM N(R) for R = 10, 30, 100, 300, 1000',
    'dev': 'NUM N(R)/(pi^2 R^4/2) - 1 for R = 10, 30, 100, 300, 1000',
    'today': 'NUM V_perp (GeV^-4), l_perp (GeV^-1), l_perp (m), l_perp (fm), 2 pi/l_perp (MeV)',
    'scales': 'NUM l_perp/proton radius, Bohr radius/l_perp, R = l_perp/(2 pi) (GeV^-1, m), M_D = (2 pi)^(2/3) M8 (GeV)',
    'nmu': 'NUM N(mu) = V_perp mu^4/(32 pi^2) for mu = 5 GeV, m_h, M8; lattice count below 5 GeV',
    'levels': 'NUM first levels: k, r_4(k), Dirac fermions 4 r_4(k), mass 2 pi sqrt(k)/l_perp (MeV), k = 1..8',
    'pois': 'NUM images S(r), Kaluza-Klein (4 pi^2/(3 r)) Sum exp(-2 pi r |m|), relative difference, at r = 1/4, 1/2, 1, 2',
    'z52': 'NUM Z(5/2) = 7 zeta(5/2) zeta(3/2), shell sum with tail',
    'nearfar': 'NUM near field (r^5 S - 1)/(r^5 Z(5/2)) at r = 1/20, 1/10; far field (3 r S/(4 pi^2) - 1)/(8 exp(-2 pi r)) at r = 1, 2, 3',
    'chr': 'NUM Christoffel {4 00} at the brane (inhomogeneous metric, x0 = 1/2)',
    'gperm': 'NUM Gamma17[A] as signed permutations: column j has the entry sign(c) in row |c|, c listed for A = 0..7',
    'reU': 'NUM Re(2 sqrt(2) U)', 'imU': 'NUM Im(2 sqrt(2) U)',
    'spinor': 'NUM example psi = U.(u (x) e_1), u = (1, 2, 0, -1): psibar psi, psibar Gamma^mu psi (mu = 0..3)',
    'tstar': 'NUM T* = (M8^6/M_Pl)^(1/5) (GeV), f_W/rho_Lambda for f_W = M8^4, 1/(2 pi M8) (GeV^-1)',
    'quat': 'NUM quaternion products e1 e2, e2 e3, e3 e1: coefficients of e3, e1, e2'}

# the groups compared between the two programs: (label, key, indices or None for all, tolerance)
GROUPS = [
    ('orbit sizes of the shells $k=1,\\dots,12$', 'orbits', None, 0.0),
    ('$r_4(k)$, $k=0,\\dots,16$', 'r4', None, 0.0),
    ('Weyl: $N(R)$, $R=10,30,100,300,1000$', 'NR', None, 0.0),
    ('Weyl: $N(R)/(\\pi^2R^4/2)-1$', 'dev', None, 1e-9),
    ('$V_\\perp$, $\\ell_\\perp$ (GeV$^{-1}$, m, fm), $2\\pi/\\ell_\\perp$', 'today', None, 1e-12),
    ('$\\ell_\\perp/r_p$, $a_0/\\ell_\\perp$, $R=\\ell_\\perp/2\\pi$, $M_D$', 'scales', None, 1e-12),
    ('$N(\\mu)$ at 5 GeV, $m_h$, $M_8$; lattice count below 5 GeV', 'nmu', None, 1e-12),
    ('the first levels ($k$, $r_4$, $4r_4$, mass)', 'levels', None, 1e-12),
    ('image sums $S(r)$ and Kaluza–Klein sums, $r=\\frac14,\\frac12,1,2$', 'pois', (0, 1, 3, 4, 6, 7, 9, 10), 1e-9),
    ('$Z(5/2)$: closed form, shell sum', 'z52', None, 1e-9),
    ('near-field and far-field ratios', 'nearfar', None, 1e-6),
    ('Christoffel symbol $\\{4\\,00\\}$ at the brane', 'chr', None, 1e-15),
    ('$\\Gamma_A$ as signed permutations', 'gperm', None, 0.0),
    ('$\\mathrm{Re}(2\\sqrt2\\,U)$', 'reU', None, 0.0),
    ('$\\mathrm{Im}(2\\sqrt2\\,U)$', 'imU', None, 0.0),
    ('example spinor: $\\bar\\psi\\psi$, $\\bar\\psi\\Gamma^\\mu\\psi$', 'spinor', None, 1e-15),
    ('$T_*$, $f_W/\\rho_\\Lambda$, $1/(2\\pi M_8)$', 'tstar', None, 1e-12)]


def nbv(D, k):
    return nums_in(D['nums'].get('tb', {}).get(NB[k], ''))


def pyv(D, k):
    return nums_in(D['ref'].get('tb', ('', '', {}))[2].get(PY[k], ''))


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
    c = D['checks'].get('tb', [])
    r = [(k, v) for k, v in D.get('ref', {}).get('tb', ('', '', {}))[2].items() if v in ('True', 'False')]
    return len(c), sum(1 for x in c if x[1] == 'True'), len(r), sum(1 for k, v in r if v == 'True')


def g(v, i, n=6, default='?'):
    try:
        return fmt(v[i], n)
    except (IndexError, TypeError):
        return default


def matrix(v, rows, cols):
    """a flat list of numbers as a rows x cols list of lists, or None"""
    if not v or len(v) != rows*cols:
        return None
    return [v[r*cols:(r + 1)*cols] for r in range(rows)]


def U8(D):
    """2 sqrt(2) U from the notebook (Re and Im parts) as a 16 x 16 list of complex numbers, or None"""
    re, im = matrix(nbv(D, 'reU'), 16, 16), matrix(nbv(D, 'imU'), 16, 16)
    if re is None or im is None:
        return None
    return [[complex(re[i][k], im[i][k]) for k in range(16)] for i in range(16)]


def cnum_tex(z):
    z = complex(z)
    if z == 0: return '0'
    if z.imag == 0: return '%d' % round(z.real)
    if z.real == 0:
        s = round(z.imag)
        return {1: 'i', -1: '-i'}.get(s, '%di' % s)
    return '%d%+di' % (round(z.real), round(z.imag))


def U8_tex(D):
    M = U8(D)
    if M is None:
        return None
    rows = [' & '.join(cnum_tex(z) for z in row) for row in M]
    return r'\left(\begin{smallmatrix}' + r' \\ '.join(rows) + r'\end{smallmatrix}\right)'


def gamma_codes(D):
    v = nbv(D, 'gperm')
    return matrix([int(x) for x in v], 8, 16) if v and len(v) == 128 else None


# ------------------------------------------------------------------------------------------------ numbers of the earlier sections
def sm(D, key):
    return D['nums'].get('sm', {}).get(key, '')


def cp(D, key):
    import doc_coupled
    return doc_coupled.nbv(D, key)


def dk(D, key):
    import doc_dark
    return doc_dark.nbv(D, key)


def grav(D, key):
    return D['nums'].get('grav', {}).get(key, '')


def first_num(s):
    v = nums_in(s or '')
    return v[0] if v else None
