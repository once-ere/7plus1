"""Content of COUPLED_SOLUTIONS (Markdown + LaTeX): the coupled field equations of the 7+1 model with every interaction in
equation form, and their accurate approximate solutions in the three gravity regimes — primordial (with its cubic,
Gauss-Bonnet and Einstein sub-regimes), modern, far future (sections 20-23 of 7+1.nb).  Every equation used is written here;
every number is parsed from the evaluated notebook or from the output of tools/ref_coupled.py, or computed here."""
import math
from docgen import Doc
from doc_stage3_common import fmt, nums_in, break_sum
from doc_common import signature_block

# ------------------------------------------------------------------------------------------------ numbers of both programs
NB = {  # short name -> output key of sections 20-23 of 7+1.nb
    'u12': 'u1, u2', 'pq': 'p1, q1, p2, q2', 'peff': 'p_eff at u -> u1, u2, u -> 0',
    'cub': 'cubic: h0, h1, h2, A0, A1 (xi = 10^-6)', 'de': 'Einstein: d1, e1 (alpha2 = 1)', 'f12': 'Gauss-Bonnet: f1, f2',
    'quad': 'tau(u2), Ha(u2), s1, s2, Ha_mid',
    'acc': 'max errors: cubic (Ha >= 10^4), GB (3 .. 300), Einstein (<= 0.1); worst of the best approximation in the crossovers',
    'cross': 'tau at Ha = 10^3 (cubic -> GB) and Ha = 1 (GB -> Einstein)',
    'near': 'near the singularity (Ha = 10^15): p_eff, Hb, kappa eps tau^3, alpha2 E2/(kappa eps), alpha3 E3/(kappa eps)',
    'plat': '(p_eff, gas fraction) at Ha = 3 10^5, 30, 10^-2', 'late': 'late (Ha = 10^-14): tau Ha, tau Hb, pa/eps, pb/eps',
    'grow': 'Einstein-regime growth exponent: numerical, 2 - 3 p1 - 5 q1', 'Fk': 'F(1), F(8), F(16)',
    'tail': 'tail exponent (K = 8 .. 16), -(3 Sqrt[2] + 2)/2', 'brane': 'brane fraction at T = H = M8',
    'cap0': 'start tau, beta(start)', 'capfit': 'fit A, vartheta, tau_E; max residual/amplitude',
    'capHE': '<H_E>(tau - tau_E) at tau = 10^5, 2 10^5', 'capgas': 'gas share at tau = 10^5, 3 10^5', 'caprad': 'radion share 1 - gas share',
    'decay': 'Gamma (GeV), lifetime (s), T_RH (GeV)', 'nkk': 'Kaluza-Klein momenta below m_r/2',
    'bulk': 'radion into bulk dirac16complex pairs (m << m_r/2): Gamma_bulk lattice, continuum (GeV), lifetime (s), Gamma_bulk/Gamma_brane',
    'vperp': 'V_perp (GeV^-4), l_perp (m), 2 Pi/l_perp (GeV)', 'tstar': 'normalcy temperature (M8^6/M_Pl)^(1/5) (GeV)',
    'gstar': 'g*, g*s at 1 TeV, 10 GeV, 1 GeV, 300 MeV, 100 MeV, 10 MeV, 1 MeV, 100 keV, 1 keV',
    'tT': 't(T) in s at 1 GeV, 100 MeV, 10 MeV, 1 MeV, 100 keV, 10 keV', 'tinst': '1/(2 Ha(T)) / t(T)',
    'late4': 'Omega_r, Omega_Lambda, z_eq, age (Gyr)', 'comp': 'max error of the composite solution',
    'radnow': 'radion displacement today (m_r = 10 GeV): Lambda4 + matter, tension f_W = M8^4', 'gbnow': 'alpha2 H0^2, alpha3 H0^4',
    'hinf': 'H_inf (1/s, GeV), T_GH (K)', 'gh': '2 Pi M/H_inf for M = 1 eV, 1 GeV',
    'ds': 'de Sitter test: beta(3000), Ha(3000), beta_eq, H_eq',
    'shares': 'shares (Einstein, Gauss-Bonnet, cubic, gas) at Ha = 3 10^5, 30, 10^-2, 10^15'}
PY = {  # short name -> output key of tools/ref_coupled.py
    'u12pq': 'NUM u1, u2, p1, q1, p2, q2',
    'peff': "NUM effective Kasner exponent p_eff = -Ha^2/Ha' along the curve at u -> u1, u = u2, u -> 0 (xi = 1e-6)",
    'cub': 'NUM cubic-regime coefficients (xi = 1e-6): h0, h1, h2, A0, A1',
    'f12': 'NUM Gauss-Bonnet regime: p_eff = p2 + f1/Ha^2 + f2 xi Ha^2 + ..., f1, f2 (alpha2 = 1)',
    'quad': 'NUM exact quadrature: tau(u2), Ha(u2), s1 = tau(u2) - p1/Ha(u2) + Int_0^Ha(u2) (1/F - p1/H^2) dH',
    'shift': 'NUM time shifts: tau = p1/Ha - e1 Ha + s1 (Einstein), tau = p2/Ha + f1/(3 Ha^3) - f2 xi Ha + s2 (Gauss-Bonnet, at Ha_mid = (f1/(f2 xi))^(1/4)): s1, s2, Ha_mid',
    'acc': 'NUM accuracy (max relative error of tau(Ha)): cubic for Ha >= 1e4, Gauss-Bonnet for 3 <= Ha <= 300, Einstein for Ha <= 0.1',
    'worst': 'NUM worst error of the best regime approximation over 1e-3.5 <= Ha <= 1e8 (in the crossover windows)',
    'cross': 'NUM vacuum crossover times: tau at Ha = 1e3 (cubic -> GB) and Ha = 1 (GB -> Einstein)',
    'near': 'NUM coupled, near the singularity (Ha = 1e15): p_eff, Hb, kappa eps tau^3, alpha2 E2/(kappa eps), alpha3 E3/(kappa eps)',
    'plat': 'NUM coupled: (p_eff, gas fraction kappa eps/(S1 + S2 + S3 + kappa eps)) at Ha = 3e5 (vacuum cubic), 30 (Gauss-Bonnet), 1e-2 (Einstein)',
    'shares': 'NUM shares (Einstein S1, Gauss-Bonnet S2, cubic S3, gas kappa eps) of the 00 equation at Ha = 3e5, 30, 1e-2, 1e15',
    'late': 'NUM coupled, gas-dominated late regime (Ha = 1e-14): tau Ha, tau Hb, pa/eps, pb/eps',
    'grow': 'NUM gas-fraction growth exponent in the Einstein regime (numerical, between Ha = 1e-2 and 1e-4), and 2 - 3p1 - 5q1',
    'hst': 'NUM h* = -sqrt(alpha2/(54 alpha3)), h0 = -sqrt(alpha2/(30 alpha3)) (xi = 1e-6)',
    'Fk': 'NUM Einstein-Kasner pair creation |beta|^2 = F(K) at K = 0.3, 1, 3, 8, 16, and the tail exponent d ln F/d ln K between 8 and 16',
    'brane': 'NUM brane radiation fraction kappa eps_W/(3 H^2) at T = M8, H = M8 (M8 = 10 TeV, reduced M_Pl = 2.435e18 GeV, g* = 106.75)',
    'cap0': 'NUM capture: start tau, beta(start) = ln(b/b0), m_r',
    'capfit': 'NUM after capture: tau_E from <H_E> over the period after tau = 2.9e5, then beta = A sin(m_r tau + vartheta)/(m_r (tau - tau_E)) fitted on 1e5 <= tau <= 3e5: A, vartheta, tau_E, max residual/amplitude',
    'capHE': 'NUM after capture: <H_E>(tau - tau_E) over one period at tau = 1e5, 2e5 (H_E = Ha + 2Hb, four-dimensional matter domination: 2/3)',
    'capgas': 'NUM after capture: gas fraction kappa eps/(3 H_E^2) at tau = 1e5 and 3e5 (constant: both are four-dimensional matter)',
    'caprad': 'NUM after capture: radion share of the energy, 1 - gas fraction (the cosmological moduli problem)',
    'decay': 'NUM radion of mass 10 GeV (M8 = 10 TeV): decay rate into brane fields (b, c, tau, mu pairs and two gluons) in GeV, lifetime in s, reheating temperature in GeV',
    'nkk': 'NUM number of Kaluza-Klein momenta of the auxiliary torus with |q| < m_r/2 (m_r = 10 GeV): V_perp (pi^2/2)(m_r/2)^4/(2 pi)^4',
    'bulk': 'NUM radion into bulk dirac16complex pairs (m << m_r/2): Gamma_bulk lattice, continuum (GeV), lifetime (s), Gamma_bulk/Gamma_brane',
    'consts': 'NUM M8, M_Pl (reduced), V_perp = M_Pl^2/M8^6 (GeV^-4), l_perp = V_perp^(1/4) (m), lowest Kaluza-Klein mass 2 pi/l_perp (GeV)',
    'gstar': 'NUM g*(T), g*s(T) at T = 1 TeV, 10 GeV, 1 GeV, 300 MeV, 100 MeV, 10 MeV, 1 MeV, 100 keV, 1 keV',
    'tT': 'NUM radiation era: t(T) in s (entropy conservation with g*(T), g*s(T)) at T = 1 GeV, 100 MeV, 10 MeV, 1 MeV, 100 keV, 10 keV',
    'tinst': 'NUM radiation era: ratio of the instantaneous approximation t = 1/(2H(T)) to t(T)',
    'thalf': 'NUM radiation era: relative change of t(10 keV) when the grid is halved (error estimate)',
    'late4': 'NUM late universe: Omega_r, Omega_Lambda, a_eq = Omega_r/Omega_m, z_eq, age (Gyr)',
    'comp': 'NUM composite solution t(a) = t_rm(a) + t_mL(a) - t_m(a): maximal relative error for a = 10^(-8), 10^(-7.75), ..., 10^0.5',
    'radnow': 'NUM radion displacement today (m_r = 10 GeV): from Lambda4 + matter, from a brane tension f_W = M8^4',
    'gbnow': 'NUM Gauss-Bonnet and cubic corrections today: alpha2 H0^2 and alpha3 H0^4 with alpha2 = 1/M8^2, alpha3 = 1e-6/M8^4',
    'tstar': 'NUM normalcy temperature (order of magnitude) T* = (M8^6/M_Pl)^(1/5) in GeV',
    'hinf': 'NUM far future: H_inf = H0 sqrt(Omega_Lambda) (1/s and GeV), Gibbons-Hawking temperature H/(2 pi) in K',
    'gh': 'NUM far future: 2 pi M/H_inf for M = 1 eV and 1 GeV (|beta|^2 = 1/(exp(2 pi M/H) + 1))',
    'ds': 'NUM de Sitter radion test (l2 = 1, Lambda = 3e-4, m_r = 0.1, start 0.05 above): beta(3000), Ha(3000), and the equilibrium (beta_eq, H_eq) of the two algebraic equations'}


def nbv(D, k):
    return nums_in(D['nums'].get('cp', {}).get(NB[k], ''))


def pyv(D, k):
    return nums_in(D['ref'].get('cp', ('', '', {}))[2].get(PY[k], ''))


def _sel(v, idx):
    try:
        return [v[i] for i in idx]
    except (IndexError, TypeError):
        return []


def compare(D):
    """(label, notebook values, Python values, tolerance, max relative difference) for every number printed by both programs"""
    rows = [
        ('Kasner constants $u_1,u_2,p_1,q_1,p_2,q_2$', nbv(D, 'u12') + nbv(D, 'pq'), pyv(D, 'u12pq'), 1e-14),
        ('$p_{\\rm eff}$ at $u\\to u_1$, $u=u_2$, $u\\to0$', nbv(D, 'peff'), pyv(D, 'peff'), 1e-10),
        ('cubic-regime coefficients $h_0,h_1,h_2,A_0,A_1$', nbv(D, 'cub'), pyv(D, 'cub'), 1e-13),
        ('Gauss–Bonnet coefficients $f_1,f_2$', nbv(D, 'f12'), pyv(D, 'f12'), 1e-9),
        ('$\\tau(u_2)$, $H_a(u_2)$, $s_1$ (quadrature)', nbv(D, 'quad')[:3], pyv(D, 'quad'), 1e-12),
        ('$s_2$, $H_{a,\\rm mid}$', nbv(D, 'quad')[3:5], pyv(D, 'shift')[1:3], 1e-8),
        ('maximal errors of the regime formulas', nbv(D, 'acc'), pyv(D, 'acc') + pyv(D, 'worst'), 1e-3),
        ('crossover times $\\tau(H_a=10^3)$, $\\tau(H_a=1)$', nbv(D, 'cross'), pyv(D, 'cross'), 1e-9),
        ('coupled, near the singularity', nbv(D, 'near'), pyv(D, 'near'), 1e-8),
        ('coupled: $(p_{\\rm eff}$, gas fraction$)$ in the three Lovelock regimes', nbv(D, 'plat'), pyv(D, 'plat'), 1e-6),
        ('coupled, gas-dominated late regime', nbv(D, 'late'), pyv(D, 'late'), 2e-3),
        ('growth exponent of the gas fraction (Einstein regime)', nbv(D, 'grow'), pyv(D, 'grow'), 1e-8),
        ('pair creation $F(1),F(8),F(16)$ and the tail exponent', nbv(D, 'Fk') + nbv(D, 'tail')[:1], _sel(pyv(D, 'Fk'), (1, 3, 4, 5)), 1e-6),
        ('brane plasma fraction at $T=H=M_8$', nbv(D, 'brane'), pyv(D, 'brane'), 1e-12),
        ('capture: start $\\tau$, $\\beta$', nbv(D, 'cap0'), pyv(D, 'cap0')[:2], 1e-8),
        ('capture: $A,\\vartheta,\\tau_E$, residual', nbv(D, 'capfit'), pyv(D, 'capfit'), 1e-6),
        ('capture: $\\langle H_E\\rangle(\\tau-\\tau_E)$', nbv(D, 'capHE'), pyv(D, 'capHE'), 1e-8),
        ('capture: gas share, radion share', nbv(D, 'capgas') + nbv(D, 'caprad'), pyv(D, 'capgas') + pyv(D, 'caprad'), 1e-8),
        ('radion decay: $\\Gamma$, lifetime, $T_{\\rm RH}$', nbv(D, 'decay'), pyv(D, 'decay'), 1e-12),
        ('Kaluza–Klein momenta below $m_r/2$', nbv(D, 'nkk'), pyv(D, 'nkk'), 1e-12),
        ('radion decay into the bulk: lattice sum, continuum limit, lifetime, ratio', nbv(D, 'bulk'), pyv(D, 'bulk'), 1e-9),
        ('$V_\\perp$, $\\ell_\\perp$, $2\\pi/\\ell_\\perp$', nbv(D, 'vperp'), pyv(D, 'consts')[2:5], 1e-12),
        ('normalcy temperature $T_*$', nbv(D, 'tstar'), pyv(D, 'tstar'), 1e-12),
        ('$g_*(T)$, $g_{*s}(T)$', nbv(D, 'gstar'), pyv(D, 'gstar'), 1e-6),
        ('$t(T)$', nbv(D, 'tT'), pyv(D, 'tT'), 1e-6),
        ('$1/(2H_a)$ divided by $t(T)$', nbv(D, 'tinst'), pyv(D, 'tinst'), 1e-6),
        ('$\\Omega_r$, $\\Omega_\\Lambda$, $z_{\\rm eq}$, age', nbv(D, 'late4'), _sel(pyv(D, 'late4'), (0, 1, 3, 4)), 1e-9),
        ('maximal error of the composite late-time solution', nbv(D, 'comp'), pyv(D, 'comp'), 1e-6),
        ('radion displacement today', nbv(D, 'radnow'), pyv(D, 'radnow'), 1e-12),
        ('$\\alpha_2H_0^2$, $\\alpha_3H_0^4$', nbv(D, 'gbnow'), pyv(D, 'gbnow'), 1e-12),
        ('$H_\\infty$, Gibbons–Hawking temperature', nbv(D, 'hinf'), pyv(D, 'hinf'), 1e-12),
        ('$2\\pi M/H_\\infty$', nbv(D, 'gh'), pyv(D, 'gh'), 1e-12),
        ('de Sitter radion test', nbv(D, 'ds'), pyv(D, 'ds'), 1e-9),
        ('shares of the four terms of the 00 equation (section 9)', nbv(D, 'shares'), pyv(D, 'shares'), 1e-6)]
    out = []
    for lab, a, b, tol in rows:
        out.append((lab, a, b, tol, max_rel(a, b)))
    return out


def max_rel(a, b):
    if not a or not b or len(a) != len(b) or any(x is None for x in list(a) + list(b)):
        return None
    return max(abs(x - y)/max(abs(x), abs(y), 1e-300) for x, y in zip(a, b))


def g(v, i, n=6, default='?'):
    try:
        return fmt(v[i], n)
    except (IndexError, TypeError):
        return default


# ================================================================================================================ document
def build(D, today):
    ck = D['checks'].get('cp', []); ref = D['ref'].get('cp', ('', '', {}))[2]
    nck = len(ck); nok = sum(1 for c in ck if c[1] == 'True')
    rb = [(k, v) for k, v in ref.items() if v in ('True', 'False')]
    cmp_ = compare(D)
    worst = max([r for *_, r in cmp_ if r is not None], default=None)
    u12, pq, peff, cub, de, f12 = nbv(D, 'u12'), nbv(D, 'pq'), nbv(D, 'peff'), nbv(D, 'cub'), nbv(D, 'de'), nbv(D, 'f12')
    quad, acc, cross = nbv(D, 'quad'), nbv(D, 'acc'), nbv(D, 'cross')
    near, plat, late, grow = nbv(D, 'near'), nbv(D, 'plat'), nbv(D, 'late'), nbv(D, 'grow')
    Fk, tail, brane = nbv(D, 'Fk'), nbv(D, 'tail'), nbv(D, 'brane')
    cap0, capfit, capHE, capgas, caprad = nbv(D, 'cap0'), nbv(D, 'capfit'), nbv(D, 'capHE'), nbv(D, 'capgas'), nbv(D, 'caprad')
    decay, nkk, vperp, tstar = nbv(D, 'decay'), nbv(D, 'nkk'), nbv(D, 'vperp'), nbv(D, 'tstar')
    gst, tT, tinst, late4, comp = nbv(D, 'gstar'), nbv(D, 'tT'), nbv(D, 'tinst'), nbv(D, 'late4'), nbv(D, 'comp')
    radnow, gbnow, hinf, gh, ds = nbv(D, 'radnow'), nbv(D, 'gbnow'), nbv(D, 'hinf'), nbv(D, 'gh'), nbv(D, 'ds')
    thalf = pyv(D, 'thalf')
    hst = pyv(D, 'hst')

    d = Doc('The coupled field equations of the 7+1 model and their solutions in the three gravity regimes',
            'Every interaction in equation form; accurate approximate solutions in the primordial regime (cubic, Gauss–Bonnet and '
            'Einstein sub-regimes), the modern regime and the far future, checked against high-accuracy numerical solutions — '
            'sections 20–23 of the notebook 7+1.nb',
            ['Repository: github.com/once-ere/7plus1 — prepared for Nash',
             'Notebook: 7+1.nb (Wolfram Language 15.0.1), evaluated ' + D['log'].get('finished', ''),
             'Document generated by tools/make_docs.py on ' + today],
            'The 7+1 model has gravity with all three Lovelock terms, the 16-component field dirac16complex with its mass, its nine '
            'self-couplings and its spin-connection coupling to gravity, a stabilizing potential of the auxiliary space, and the Standard '
            'Model on a world brane with three portals to dirac16complex. This document writes every one of these interactions as an '
            'equation, derives the coupled field equations and their homogeneous cosmological form, and solves them in the three gravity '
            'regimes: (I) the PRIMORDIAL regime near $x_0=0$, where the cubic, then the Gauss–Bonnet, then the Einstein term dominates; '
            '(II) the MODERN regime, where the auxiliary space is held by its potential and gravity on the brane is four-dimensional; '
            '(III) the FAR FUTURE, de Sitter expansion. Each regime has closed-form solutions (exact where they exist) whose accuracy is '
            'measured against the numerical solution of the full equations. Sections 20–23 of 7+1.nb contain %d checks (%d True); the '
            'independent Python program tools/ref_coupled.py contains %d (%d True); the %d groups of numbers printed by both agree to a '
            'largest relative difference of %s.' % (nck, nok, len(rb), sum(1 for k, v in rb if v == 'True'), len(cmp_),
                                                    fmt(worst, 2) if worst is not None else '?'))

    # ==============================================================================================
    d.h1('1. Results at a glance')
    d.ul(['**Every interaction** of the model is a term of one action, $S=S_0+S_1+S_2+S_3+S_\\psi+S_m+S_{\\rm self}+S_U+S_f+S_{\\rm SM}+S_P$ '
          '(section 3). The cubic Lovelock density in invariant form is '
          '$L_3=R^3-12RR_{\\mu\\nu}R^{\\mu\\nu}+16R^\\mu{}_\\nu R^\\nu{}_\\rho R^\\rho{}_\\mu+24R^{\\mu\\rho}R^{\\nu\\sigma}R_{\\mu\\nu\\rho\\sigma}'
          '+3RR_{\\mu\\nu\\rho\\sigma}R^{\\mu\\nu\\rho\\sigma}-24R^{\\mu\\nu}R_\\mu{}^{\\rho\\sigma\\tau}R_{\\nu\\rho\\sigma\\tau}'
          '+4R^{\\mu\\nu}{}_{\\rho\\sigma}R^{\\rho\\sigma}{}_{\\lambda\\tau}R^{\\lambda\\tau}{}_{\\mu\\nu}-8R_\\mu{}^\\rho{}_\\nu{}^\\sigma R_\\rho{}^\\lambda{}_\\sigma{}^\\tau R_\\lambda{}^\\mu{}_\\tau{}^\\nu$ (checked exactly).',
          '**The radion** $\\beta=\\ln(b/b_0)$ obeys $\\beta\'\'+3H_a\\beta\'+4\\beta\'^2=\\kappa(\\varepsilon/6-p_a/2+p_b/3)+\\Lambda/3$; with '
          '$U=U_0+\\tfrac{K_U}{2}\\beta^2$ its mass is $m_r^2=\\kappa K_U/12$; on the brane $8\\pi G_4=\\kappa/V_\\perp$, '
          '$\\Lambda_4=\\Lambda+\\kappa U_0+\\kappa f_W/V_\\perp$ (section 5).',
          '**Regime I, vacuum: exact solution.** With $u=H_b/H_a$ the constraint is $P_E(u)+H_a^2P_G(u)+\\xi H_a^4P_C(u)=0$ — quadratic in '
          '$H_a^2$ — so the vacuum solution is an explicit curve plus one quadrature. Its limits are the cubic regime '
          '($a\\propto\\tau$, $H_b\\to-\\sqrt{\\alpha_2/(30\\alpha_3)}$), the Gauss–Bonnet Kasner regime ($p_2=%s$, $q_2=%s$) and the Einstein '
          'Kasner regime ($p_1=%s$, $q_1=%s$); their asymptotic formulas are accurate to %s, %s and %s inside the regimes (section 6).'
          % (g(pq, 2, 6), g(pq, 3, 6), g(pq, 0, 6), g(pq, 1, 6), g(acc, 0, 2), g(acc, 1, 2), g(acc, 2, 2)),
          '**Regime I with dirac16complex:** near $x_0=0$ the gas wins: $a\\propto\\tau^{3/4}$, $H_b\\to h_*=-\\sqrt{\\alpha_2/(54\\alpha_3)}$ '
          '(for any fluid $p_a=w_a\\varepsilon$: $a\\propto\\tau^{1/(1+w_a)}$); later the gas dominates the Einstein regime and, WITHOUT a '
          'stabilizing potential, the auxiliary space re-expands, $b\\propto\\tau^{2/5}$ (section 6.4).',
          '**Regime II:** the potential $U(b)$ catches the auxiliary space; the radion then oscillates, '
          '$\\beta\\approx A\\sin(m_r\\tau+\\vartheta)/(m_r(\\tau-\\tau_E))$, and carries %s of the energy — the cosmological moduli problem of '
          'the minimal model: a 10 GeV radion decays in %s s, but into the Kaluza–Klein modes of the bulk; through brane fields alone it would '
          'reheat the brane to only %s GeV. After a later reheating, the brane obeys the Friedmann equations '
          'with $G_4=G_8/V_\\perp$; the composite analytic solution of the late universe is accurate to %s and gives the age %s Gyr (section 7).'
          % (g(caprad, 0, 3), g(nbv(D, 'bulk'), 2, 2), g(decay, 2, 2), g(comp, 0, 2), g(late4, 3, 5)),
          '**Regime III:** de Sitter with $H_\\infty=H_0\\Omega_\\Lambda^{1/2}$, $H_a=H_\\infty\\coth(3H_\\infty t/2)$; the radion is stable '
          '(friction $\\gamma$, mass $m_{\\rm eff}$ given in closed form with the Gauss–Bonnet and cubic terms) and settles at '
          '$\\beta_{\\rm eq}=H^2(1-4\\alpha_2H^2)/m_r^2$ (section 8).',
          '**The size of every interaction** in the three regimes is tabulated in section 9.'])

    # ==============================================================================================
    d.h1('2. Conventions and objects used')
    d.p('Units $\\hbar=c=1$; coordinates $X=(x_0,x_1,\\dots,x_7)$: $x_0$ is time (the direction of the octonion unit $e_0=1$), '
        '$x_1,x_2,x_3$ ordinary space, $x_4,\\dots,x_7$ the four auxiliary space-like coordinates. Greek indices $\\mu,\\nu,\\dots=0,\\dots,7$ are '
        'coordinate indices, capital Latin $A,B,\\dots=0,\\dots,7$ frame indices; repeated indices are summed. The frame metric is '
        '$\\eta=\\mathrm{diag}(+1,-1,\\dots,-1)$, simply the squares of the octonion basis.')
    signature_block(d, eta=r'\eta', eta71=r'\eta_{71}', gam=r'\Gamma')
    d.p('The real 16×16 Dirac matrices, the Dirac conjugate and the chirality:')
    d.align([r'\Gamma_0&=\begin{pmatrix}0&1_8\\1_8&0\end{pmatrix},\quad \Gamma_i=\begin{pmatrix}\tau^O_i&0\\0&-\tau^O_i\end{pmatrix},\quad \Gamma_A\Gamma_B+\Gamma_B\Gamma_A=2\eta_{AB}1_{16},',
             r'\Gamma^A&=\eta^{AA}\Gamma_A,\quad \bar\psi=\psi^\dagger\Gamma_0,\quad \chi=i\,\Gamma_0\Gamma_1\cdots\Gamma_7,'])
    d.p('with real antisymmetric 8×8 matrices $\\tau^O_i$ ($(\\tau^O_i)^2=-1_8$, mutually anticommuting). Curvature (Landau–Lifshitz):')
    d.align([r'\{{}^{\lambda}_{\mu\nu}\}&=\tfrac12g^{\lambda\rho}(\partial_\mu g_{\rho\nu}+\partial_\nu g_{\rho\mu}-\partial_\rho g_{\mu\nu}),\qquad '
             r'R^\rho{}_{\sigma\mu\nu}=\partial_\mu\{{}^{\rho}_{\sigma\nu}\}-\partial_\nu\{{}^{\rho}_{\sigma\mu}\}+\{{}^{\rho}_{\mu\lambda}\}\{{}^{\lambda}_{\sigma\nu}\}-\{{}^{\rho}_{\nu\lambda}\}\{{}^{\lambda}_{\sigma\mu}\},',
             r'R_{\sigma\nu}&=R^\rho{}_{\sigma\rho\nu},\qquad R=g^{\sigma\nu}R_{\sigma\nu},\qquad G^\mu{}_\nu=R^\mu{}_\nu-\tfrac12\delta^\mu_\nu R .'])
    d.p('The generalized Kronecker delta $\\delta^{\\mu_1\\cdots\\mu_n}_{\\nu_1\\cdots\\nu_n}=\\det[\\delta^{\\mu_i}_{\\nu_j}]$ gives the Lovelock '
        'densities and tensors ($k=0,1,2,3$; $L_0=1$, $E^{(0)\\mu}{}_\\nu=-\\delta^\\mu_\\nu/2$, $L_1=R$, $E^{(1)}=G$; $E^{(4)}\\equiv0$ in eight dimensions):')
    d.align([r'L_k&=\frac1{2^k}\,\delta^{\alpha_1\beta_1\cdots\alpha_k\beta_k}_{\gamma_1\delta_1\cdots\gamma_k\delta_k}\,R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}\cdots R^{\gamma_k\delta_k}{}_{\alpha_k\beta_k},',
             r'E^{(k)\mu}{}_\nu&=-\frac1{2^{k+1}}\,\delta^{\mu\alpha_1\beta_1\cdots\alpha_k\beta_k}_{\nu\gamma_1\delta_1\cdots\gamma_k\delta_k}\,R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}\cdots R^{\gamma_k\delta_k}{}_{\alpha_k\beta_k},'
             r'\qquad \delta\!\left(\sqrt{|g|}\,L_k\right)=\sqrt{|g|}\,E^{(k)}_{\mu\nu}\,\delta g^{\mu\nu}.'])
    d.p('Frames, frame matrices and the spin connection (the coupling of dirac16complex to gravity): with $g_{\\mu\\nu}=\\eta_{AB}e^A{}_\\mu e^B{}_\\nu$,')
    d.align([r'\Gamma^\mu(X)&=e_A{}^\mu\Gamma^A,\qquad \omega_\mu{}^A{}_B=e^A{}_\nu\big(\partial_\mu e_B{}^\nu+\{{}^{\nu}_{\mu\lambda}\}e_B{}^\lambda\big),',
             r'\Omega_\mu&=\tfrac14\eta_{AC}\,\omega_\mu{}^C{}_B\,\Gamma^A\Gamma^B,\qquad D_\mu\psi=\partial_\mu\psi+\Omega_\mu\psi .'])
    d.p('The homogeneous cosmological metric of the three regimes, its Hubble rates (a prime is $d/dx_0$) and the time since the '
        'singularity $\\tau=x_0-x_s$:')
    d.eq(r"ds^2=n^2dx_0^2-a^2(dx_1^2+dx_2^2+dx_3^2)-b^2(dx_4^2+\cdots+dx_7^2),\qquad n=1,\qquad H_a=\frac{a'}{a},\qquad H_b=\frac{b'}{b}.")
    d.p('The couplings: $\\kappa=8\\pi G_8=M_8^{-6}$ ($M_8$ the reduced eight-dimensional Planck mass), $\\alpha_1=1$, $\\alpha_2$ (Gauss–Bonnet, '
        'dimension length²) and $\\alpha_3$ (cubic, length⁴), $\\alpha_0=2\\Lambda$ with the bulk cosmological constant $\\Lambda$. The '
        'Gauss–Bonnet length $\\ell_2=\\sqrt{\\alpha_2}$ and the ratio $\\xi=\\alpha_3/\\alpha_2^2$ set the primordial regime. The auxiliary space '
        'is a torus $T^4$ of coordinate period $L_c$, proper volume $V_\\perp=b^4L_c^4$; the Standard Model lives on the world brane '
        '$W=\\{x_4=\\dots=x_7=0\\}$ with induced metric $h_{\\mu\\nu}=g_{\\mu\\nu}(x,y=0)$ ($y=(x_4,\\dots,x_7)$).')
    d.p('**The three gravity regimes.** The regimes are the three epochs of the model in which gravity behaves in three qualitatively '
        'different ways: (I) the primordial regime near $x_0=0$, governed by the higher Lovelock terms (inside it the cubic, the '
        'Gauss–Bonnet and the Einstein term take over one after another — three sub-regimes); (II) the modern regime, in which the '
        'auxiliary space is held by its stabilizing potential and gravity on the brane is four-dimensional Einstein gravity; (III) the far '
        'future, in which the four-dimensional cosmological constant dominates (de Sitter).')

    # ==============================================================================================
    d.h1('3. Every interaction in equation form (section 20.1)')
    d.p('The action of the model is the sum of the following terms ($\\sqrt{|g|}\\,d^8X$ is the volume element of the bulk, '
        '$\\sqrt{|h|}\\,d^4x$ that of the brane; a brane integral is $\\int_Wd^4x\\sqrt{|h|}(\\dots)=\\int d^8X\\sqrt{|g|}\\,\\delta^4(y)\\,b^{-4}(\\dots)$):')
    d.align([r'S_0&=-\frac1{2\kappa}\int d^8X\sqrt{|g|}\,2\Lambda &&\text{(cosmological term, }k=0),',
             r'S_1&=-\frac1{2\kappa}\int d^8X\sqrt{|g|}\,R &&\text{(Einstein–Hilbert, }k=1),',
             r'S_2&=-\frac{\alpha_2}{2\kappa}\int d^8X\sqrt{|g|}\,\big(R^2-4R_{\mu\nu}R^{\mu\nu}+R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}\big) &&\text{(Gauss–Bonnet)},',
             r'S_3&=-\frac{\alpha_3}{2\kappa}\int d^8X\sqrt{|g|}\,L_3 &&\text{(cubic Lovelock)},',
             r'S_\psi&=\int d^8X\sqrt{|g|}\;\bar\psi\,i\Gamma^\mu(X)\big(\partial_\mu+\Omega_\mu\big)\psi &&\text{(dirac16complex, gravity)},',
             r'S_m&=-\int d^8X\sqrt{|g|}\;m\,\bar\psi\psi &&\text{(mass)},',
             r'S_{\rm self}&=-\frac12\sum_{r=0}^{8}g_r\int d^8X\sqrt{|g|}\sum_{|K|=r}\eta_K\big(\bar\psi\Theta^K\psi\big)^2 &&\text{(self-couplings)},',
             r'S_U&=-\int d^8X\sqrt{|g|}\;U(b) &&\text{(stabilizing potential)},',
             r'S_f&=-\int_Wd^4x\sqrt{|h|}\;f_W &&\text{(brane tension)},',
             r'S_{\rm SM}&=\int_Wd^4x\sqrt{|h|}\;\mathcal L_{\rm SM} &&\text{(Standard Model)},',
             r'S_P&=\int_Wd^4x\sqrt{|h|}\;\Big[-\sum_{i,s}y_{is}\,\bar L_i\tilde\Phi P_R\psi_s+\text{h.c.} &&\text{(portals)}',
             r'&\qquad-G_P^f\,(\bar\psi\Gamma^\mu\psi)|_W(\bar f\gamma_\mu f)-\lambda_P\,\Phi^\dagger\Phi\,(\bar\psi\psi)|_W\Big].'])
    d.p('Here the cubic Lovelock density in invariant form is (checked EXACTLY in section 20.1 of 7+1.nb against the generalized-delta '
        'definition on four random integer tensors with all algebraic symmetries of the Riemann tensor, $R_{abcd}=S_{ac}S_{bd}-S_{ad}S_{bc}$ '
        'summed over three random symmetric integer matrices $S$):')
    d.align([r'L_3&=R^3-12R\,R_{\mu\nu}R^{\mu\nu}+16R^\mu{}_\nu R^\nu{}_\rho R^\rho{}_\mu+24R^{\mu\rho}R^{\nu\sigma}R_{\mu\nu\rho\sigma}+3R\,R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}',
             r'&\quad-24R^{\mu\nu}R_\mu{}^{\rho\sigma\tau}R_{\nu\rho\sigma\tau}+4R^{\mu\nu}{}_{\rho\sigma}R^{\rho\sigma}{}_{\lambda\tau}R^{\lambda\tau}{}_{\mu\nu}-8R_\mu{}^\rho{}_\nu{}^\sigma R_\rho{}^\lambda{}_\sigma{}^\tau R_\lambda{}^\mu{}_\tau{}^\nu',
             r'&=\frac18\,\delta^{\alpha_1\beta_1\alpha_2\beta_2\alpha_3\beta_3}_{\gamma_1\delta_1\gamma_2\delta_2\gamma_3\delta_3}R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}R^{\gamma_2\delta_2}{}_{\alpha_2\beta_2}R^{\gamma_3\delta_3}{}_{\alpha_3\beta_3}.'])
    d.p('The self-couplings use, for every strictly increasing multi-index $K=(A_1<\\dots<A_r)$ of rank $r$ (256 in all), the Clifford products')
    d.eq(r'\Theta^K=i^{\sigma_r}\,\Gamma^{A_1}\cdots\Gamma^{A_r},\quad \sigma_r=\tfrac{r(r-1)}2\bmod2,\qquad \eta_K=\eta_{A_1A_1}\cdots\eta_{A_rA_r},\qquad \Theta^{()}=1_{16},')
    d.p('so that every bilinear $\\bar\\psi\\Theta^K\\psi$ is real. Expanded about flat space, $g_{\\mu\\nu}=\\eta_{\\mu\\nu}+h_{\\mu\\nu}$, the '
        'gravitational couplings of every field are the vertex $-\\tfrac12h_{\\mu\\nu}T^{\\mu\\nu}+O(h^2)$ with its energy–momentum tensor '
        '$T^{\\mu\\nu}$ (section 4). The Standard Model Lagrangian on the brane (frame of $h$, covariant derivative '
        '$D_\\mu=\\nabla_\\mu-ig_sG^a_\\mu T^a-igW^I_\\mu\\tau^I/2-ig\'YB_\\mu$) is')
    d.align([r'\mathcal L_{\rm SM}&=-\tfrac14G^a_{\mu\nu}G^{a\mu\nu}-\tfrac14W^I_{\mu\nu}W^{I\mu\nu}-\tfrac14B_{\mu\nu}B^{\mu\nu}+(D_\mu\Phi)^\dagger(D^\mu\Phi)-\mu^2\Phi^\dagger\Phi-\lambda(\Phi^\dagger\Phi)^2',
             r'&\quad+\sum_f\bar f\,i\gamma^\mu D_\mu f-\big[\bar Q_LY_u\tilde\Phi u_R+\bar Q_LY_d\Phi d_R+\bar L_LY_e\Phi e_R+\bar L_LY_\nu\tilde\Phi\nu_R+\text{h.c.}\big],',
             r'G^a_{\mu\nu}&=\partial_\mu G^a_\nu-\partial_\nu G^a_\mu+g_sf^{abc}G^b_\mu G^c_\nu,\qquad W^I_{\mu\nu}=\partial_\mu W^I_\nu-\partial_\nu W^I_\mu+g\,\epsilon^{IJK}W^J_\mu W^K_\nu,\qquad B_{\mu\nu}=\partial_\mu B_\nu-\partial_\nu B_\mu,'])
    d.p('which contains the gauge–fermion vertices $g_s\\bar q\\gamma^\\mu T^aqG^a_\\mu$, $g\\bar f_L\\gamma^\\mu(\\tau^I/2)f_LW^I_\\mu$, '
        '$g\'Y\\bar f\\gamma^\\mu fB_\\mu$, the cubic and quartic gauge self-couplings in $-\\tfrac14F^2$, the Higgs–gauge couplings in '
        '$|D\\Phi|^2$, the Higgs self-couplings $-\\mu^2|\\Phi|^2-\\lambda|\\Phi|^4$ and the Yukawa couplings ($\\tilde\\Phi=i\\sigma_2\\Phi^*$). '
        'In the portal term $L_i$ are the lepton doublets, $\\psi_s$ the brane-right-handed parts of dirac16complex at the brane, $f$ any '
        'Standard Model fermion. Finally, the stabilizing potential used in sections 6–8 is')
    d.eq(r'U(b)=U_0+\frac{K_U}{2}\Big(\ln\frac{b}{b_0}\Big)^2,\qquad b^8=\det(-g_{kl})\ \ (k,l=4,\dots,7),')
    d.p('whose minimum $b_0$ fixes the size of the auxiliary space. One more interaction follows from these terms rather than being '
        'added: the RADION–BRANE coupling. The brane fields see $h=(b_0/b)^4g_E$ with the four-dimensional Einstein-frame metric $g_E$, '
        'so that with the canonical radion $\\phi=\\sqrt{12}\\,M_{\\rm Pl}\\,\\ln(b/b_0)$')
    d.align([r'\mathcal L_{\phi W}&=\frac{\phi}{\sqrt3\,M_{\rm Pl}}\,T^\mu{}_\mu\big|_{\rm brane},\qquad T^\mu{}_\mu\supset m_f\bar ff,',
             r'T^\mu{}_\mu&\supset-\frac{b_{\rm QCD}\,\alpha_s}{8\pi}G^a_{\mu\nu}G^{a\mu\nu}\quad(\text{trace anomaly},\ b_{\rm QCD}=11-\tfrac23n_f).'])

    # ==============================================================================================
    d.h1('4. The coupled field equations and the energy–momentum tensors (section 20.2)')
    d.p('Varying the action with respect to $g^{\\mu\\nu}$, $\\bar\\psi$ and the Standard Model fields gives the coupled field equations:')
    d.align([r'\sum_{k=0}^{3}\alpha_kE^{(k)\mu}{}_\nu&=\kappa\Big[T_\psi+T_U+\frac{\delta^4(y)}{b^4}\big(T_f+T_{\rm SM}+T_P\big)\Big]^\mu{}_\nu,',
             r'\big(i\Gamma^\mu D_\mu-m\big)\psi&-\sum_{r=0}^{8}g_r\sum_{|K|=r}\big(\bar\psi\Theta_K\psi\big)\Theta^K\psi',
             r'&=\frac{\delta^4(y)}{b^4}\Big[\sum_{i,s}y^*_{is}\,(\tilde\Phi^\dagger L_i)\,\chi_s+G_P^f\,(\bar f\gamma_\mu f)\,\Gamma^\mu\psi+\lambda_P\,\Phi^\dagger\Phi\,\psi\Big],',
             r'\nabla_\mu F^{a\mu\nu}&=j^{a\nu}_{\rm SM}+j^{a\nu}_P,',
             r'0&=D_\mu D^\mu\Phi+\mu^2\Phi+2\lambda(\Phi^\dagger\Phi)\Phi+(\text{Yukawa})+\lambda_P(\bar\psi\psi)|_W\Phi,',
             r'i\gamma^\mu D_\mu f&=(\text{Yukawa masses})+G_P^f(\bar\psi\Gamma_\mu\psi)|_W\,\gamma^\mu f,'])
    d.p('with $\\Theta_K=\\eta_K\\Theta^K$, $\\chi_s$ the portal spinor of the brane-right-handed component $s$, and $j_P$ the portal '
        'currents. Every energy–momentum tensor is $T_{\\mu\\nu}=(2/\\sqrt{|g|})\\,\\delta S/\\delta g^{\\mu\\nu}$, i.e. '
        '$\\delta S_{\\rm matter}=\\tfrac12\\int\\sqrt{|g|}\\,T_{\\mu\\nu}\\delta g^{\\mu\\nu}$, with $T^0{}_0=\\varepsilon$. For the homogeneous metric '
        'every source is $T^\\mu{}_\\nu=\\mathrm{diag}(\\varepsilon,-p_a,-p_a,-p_a,-p_b,-p_b,-p_b,-p_b)$. For matter whose comoving energy '
        'depends only on the scale factors at fixed occupation numbers, $W(a,b)=a^3b^4\\varepsilon$, the action is $-\\int nW\\,dx_0$ and')
    d.eq(r'\varepsilon=\frac{W}{a^3b^4},\qquad p_a=-\frac{\partial_aW}{3a^2b^4},\qquad p_b=-\frac{\partial_bW}{4a^3b^3},\qquad '
         r"\varepsilon'+3H_a(\varepsilon+p_a)+4H_b(\varepsilon+p_b)=0 .")
    d.p('The sources of the model in this form (all checked in section 20.2 of 7+1.nb):')
    d.table(['source', '$W(a,b)$', '$\\varepsilon$', '$p_a$', '$p_b$'],
            [['stabilizing potential', '$a^3b^4U(b)$', '$U$', '$-U$', "$-U-\\tfrac b4U'(b)$"],
             ['brane tension (smeared, $\\delta^4(y)\\to1/L_c^4$)', '$a^3f_W/L_c^4$', '$f_W/V_\\perp$', '$-f_W/V_\\perp$', '$0$'],
             ['brane fluid, $\\rho_4\\propto a^{-3(1+w)}$', '$a^3\\rho_4/L_c^4$', '$\\rho_4/V_\\perp$', '$w\\,\\varepsilon$', '$0$'],
             ['one massless dirac16complex particle', '$\\omega=\\sqrt{k_a^2/a^2+k_b^2/b^2}$', '$\\omega/(a^3b^4)$',
              '$k_a^2/(3a^2\\omega\\,a^3b^4)$', '$k_b^2/(4b^2\\omega\\,a^3b^4)$']], widths=[0.27, 0.22, 0.13, 0.19, 0.19], small=True)
    d.p('The minisuperspace variation of $S=-\\frac1{2\\kappa}\\int n\\,a^3b^4\\sum_k\\alpha_kL_k\\,dx_0-\\int nW\\,dx_0$ (with the identities '
        '$\\mathrm{EL}_n(na^3b^4L_k)=-2a^3b^4E^{(k)0}{}_0$, $\\mathrm{EL}_a=-6a^2b^4E^{(k)1}{}_1$, $\\mathrm{EL}_b=-8a^3b^3E^{(k)4}{}_4$ at $n=1$) gives exactly')
    d.eq(r'\sum_k\alpha_kE^{(k)0}{}_0=\kappa\varepsilon,\qquad \sum_k\alpha_kE^{(k)1}{}_1=-\kappa p_a,\qquad \sum_k\alpha_kE^{(k)4}{}_4=-\kappa p_b .')
    d.p('The dirac16complex gas used in the primordial regime is a sum of such particles: a collisionless gas with an isotropic comoving '
        'momentum distribution (the massless limit of the Kohn–Sham matter, valid while $m\\ll H_a$ and the self-couplings are negligible); '
        'with $k_a=k\\cos\\theta$, $k_b=k\\sin\\theta$ and the weight $\\cos^2\\theta\\sin^3\\theta$ of the seven-dimensional momentum sphere,')
    d.align([r'\kappa\varepsilon&=\frac{C_g}{a^3b^4}\big\langle\omega\big\rangle,\qquad \kappa p_a=\frac{C_g}{a^3b^4}\Big\langle\frac{\cos^2\theta}{3a^2\omega}\Big\rangle,\qquad '
             r'\kappa p_b=\frac{C_g}{a^3b^4}\Big\langle\frac{\sin^2\theta}{4b^2\omega}\Big\rangle,\qquad \omega=\sqrt{\frac{\cos^2\theta}{a^2}+\frac{\sin^2\theta}{b^2}},',
             r'\langle f\rangle&=\frac{\int_0^{\pi/2}f(\theta)\cos^2\theta\sin^3\theta\,d\theta}{\int_0^{\pi/2}\cos^2\theta\sin^3\theta\,d\theta},'])
    d.p('evaluated with 128-point Gauss–Legendre quadrature in $\\theta$. The normalization $C_g=3\\times10^{-6}$ (units $\\ell_2=1$, section 6) '
        'makes $\\kappa\\varepsilon=10^{-6}\\times3H_a^2$ where $a=b=1$ and $H_a=1$.')

    # ==============================================================================================
    d.h1('5. The homogeneous equations, the radion and the four-dimensional constants (section 20.3)')
    d.p('With the closed forms of the Lovelock tensors of the cosmological metric ($H_a\'=dH_a/dx_0$, $H_b\'=dH_b/dx_0$):')
    cf = D.get('closed_forms_tex', {})
    names = {0: '0}{}_0', 1: '1}{}_1', 4: '4}{}_4'}
    for k in (1, 2, 3):
        lines = []
        for c in (0, 1, 4):
            lines += break_sum('E^{(%d)%s' % (k, names[c]), cf.get((k, c), '?'), 95, '.' if (k == 3 and c == 4) else ',')
        d.align(lines)
    d.p('the three equations of section 4 are ordinary differential equations; the 11 and 44 equations are linear in $(H_a\',H_b\')$ and the 00 '
        'equation is a constraint, preserved by the evolution. Written out, the constraint is')
    d.eq(r'3(H_a^2+4H_aH_b+2H_b^2)-12\alpha_2H_b(4H_a^3+18H_a^2H_b+12H_aH_b^2+H_b^3)+360\alpha_3H_a^2H_b^3(4H_a+3H_b)-\Lambda=\kappa\varepsilon .')
    d.p('**The radion** is $\\beta=\\ln(b/b_0)$, $H_b=\\beta\'$. In Einstein gravity ($\\alpha_2=\\alpha_3=0$) the 11 and 44 equations, with the '
        'constraint used to remove $H_a^2$, give exactly')
    d.eq(r"\beta''+3H_a\beta'+4\beta'^2=\kappa\Big(\frac{\varepsilon}{6}-\frac{p_a}{2}+\frac{p_b}{3}\Big)+\frac{\Lambda}{3}.")
    d.p('With $U=U_0+\\tfrac{K_U}2\\beta^2$ ($bU\'=K_U\\beta$) the potential enters the right-hand side through $\\kappa(U/3-bU\'/12)$, i.e. as '
        'the restoring term $-(\\kappa K_U/12)\\beta$:')
    d.eq(r'm_r^2=\frac{\kappa K_U}{12}.')
    d.p('**Static auxiliary space** ($H_b=H_b\'=0$): the cubic term drops out ($E^{(3)}=0$), the 00 and 11 equations are the Friedmann '
        'equations of the brane and the 44 equation fixes where the radion sits:')
    d.align([r"3H_a^2&=8\pi G_4\,\rho_4+\Lambda_4,\qquad H_a'=-\frac\kappa2(\varepsilon+p_a),\qquad 8\pi G_4=\frac{\kappa}{V_\perp},\qquad \Lambda_4=\Lambda+\kappa U_0+\frac{\kappa f_W}{V_\perp},",
             r"\frac\kappa4\,b\,U'(b)&=\Lambda+\kappa U+\frac{2\kappa f_W}{V_\perp}+\frac\kappa2(\varepsilon-3p_a)+\kappa p_b-12\alpha_2H_a^2(H_a'+H_a^2),"])
    d.p('where $\\varepsilon$, $p_a$, $p_b$ are now the matter sources besides $U$ and the tension. With $U$ as above the equilibrium is '
        '$\\beta_{\\rm eq}=[\\dots]/(3m_r^2)$, the bracket being the right-hand side of the last equation. **Einstein frame.** With '
        '$H_E=H_a+2H_b$ (the expansion rate of $a_E=ab^2$) the constraint of Einstein gravity is exactly')
    d.eq(r"3H_E^2=6\beta'^2+\kappa\varepsilon\qquad\Longleftrightarrow\qquad 3M_{\rm Pl}^2H_E^2=\tfrac12\phi'^2+\rho,\qquad \phi=\sqrt{12}\,M_{\rm Pl}\,\beta,")
    d.p('so the canonical radion is $\\phi=\\sqrt{12}M_{\\rm Pl}\\beta$ and its coupling to the brane is the $\\mathcal L_{\\phi W}$ of section 3.')

    # ==============================================================================================
    d.h1('6. Regime I: the primordial Lovelock regime (section 21)')
    d.p('Units: $\\ell_2=\\sqrt{\\alpha_2}=1$; then the vacuum equations contain the single parameter $\\xi=\\alpha_3/\\alpha_2^2$, here $\\xi=10^{-6}$ '
        '(this separates the three Lovelock sub-regimes clearly; the crossover times below are given for any $\\xi$). $\\Lambda$ is negligible '
        '($\\Lambda\\ell_2^2\\approx(H_0/M_8)^2\\approx10^{-92}$, section 9) and set to 0. If $\\ell_2=1/M_8$ with $M_8=10$ TeV, the unit of time is '
        '$6.6\\times10^{-29}$ s.')
    d.h2('6.1 The exact vacuum solution: an algebraic curve and one quadrature')
    d.p('With $u=H_b/H_a$ and $X=H_a^2$ the vacuum constraint divided by $H_a^2$ is')
    d.align([r'0&=P_E(u)+X\,P_G(u)+\xi X^2P_C(u),',
             r'P_E&=3(2u^2+4u+1),\quad P_G=-12u(u^3+12u^2+18u+4),\quad P_C=360u^3(3u+4).'])
    d.p('On $u\\in(u_1,0)$, $u_1=-1+1/\\sqrt2$, one has $P_E>0$ and $P_C<0$, so the quadratic has exactly one positive root and the vacuum '
        'solution is the explicit curve')
    d.eq(r'H_a(u)=\sqrt{X(u)},\qquad H_b(u)=u\sqrt{X(u)},\qquad X(u)=\frac{P_G+\sqrt{P_G^2-4\xi P_CP_E}}{2\xi\,|P_C|},')
    d.p('along which the solution moves ($H_a\'$ from the 11 and 44 equations); the time is one quadrature,')
    d.eq(r"\tau(H_a)=\int_{H_a}^{\infty}\frac{dH}{F(H)},\qquad F=-H_a'>0,\qquad p_{\rm eff}=-\frac{H_a^2}{H_a'}\ \ (=p\ \text{for}\ a\propto\tau^p).")
    d.p('The two Kasner solutions ($a\\propto\\tau^p$, $b\\propto\\tau^q$) at the ends of the Lovelock hierarchy are')
    d.align([r'\text{Einstein:}\quad &P_E(u_1)=0,\quad u_1=-1+\tfrac1{\sqrt2},\quad p_1=\frac1{3+4u_1}=\frac{1+2\sqrt2}7,',
             r'&q_1=u_1p_1=\frac{2-3\sqrt2}{14}\quad(3p_1+4q_1=1),',
             r'\text{Gauss–Bonnet:}\quad &u_2^3+12u_2^2+18u_2+4=0,\ u_2\in(-\tfrac13,0),\quad p_2=\frac3{3+4u_2},',
             r'&q_2=u_2p_2\quad(3p_2+4q_2=3).'])
    d.p('Numerically $u_1=%s$, $u_2=%s$, $p_1=%s$, $q_1=%s$, $p_2=%s$, $q_2=%s$. The effective exponent along the curve is %s at $u\\to u_1$, '
        '%s at $u=u_2$ and %s at $u\\to0$ ($\\xi=10^{-6}$): it runs from the Einstein value through about the Gauss–Bonnet value to 1.'
        % (g(u12, 0, 10), g(u12, 1, 10), g(pq, 0, 10), g(pq, 1, 10), g(pq, 2, 10), g(pq, 3, 10), g(peff, 0, 8), g(peff, 1, 8), g(peff, 2, 8)))
    d.h2('6.2 The accurate approximate solutions of the three Lovelock sub-regimes')
    d.p('**Cubic regime** ($\\tau\\ll\\tau_{32}$). The series $H_a=1/\\tau+A_0+A_1\\tau+\\dots$, $H_b=h_0+h_1\\tau+h_2\\tau^2+\\dots$ solves the '
        'three equations order by order with (general $\\alpha_2,\\alpha_3$)')
    d.align([r'h_0&=-\sqrt{\frac{\alpha_2}{30\alpha_3}},\qquad h_1=\frac{\alpha_2}{16\alpha_3}-\frac1{32\alpha_2},\qquad A_0=\frac{\sqrt{30}\,(58\alpha_2^2-45\alpha_3)}{1920\,\alpha_2^{3/2}\alpha_3^{1/2}},',
             r'A_1&=-\frac{3647\,\alpha_2}{92160\,\alpha_3}-\frac{223}{6144\,\alpha_2}+\frac{645\,\alpha_3}{8192\,\alpha_2^3},\qquad '
             r'h_2=\frac{\sqrt{30}\,(-1212\alpha_2^4+1180\alpha_2^2\alpha_3+225\alpha_3^2)}{307200\,\alpha_2^{5/2}\alpha_3^{3/2}},'])
    d.p('i.e. for $\\xi=10^{-6}$: $h_0=%s$, $h_1=%s$, $h_2=%s$, $A_0=%s$, $A_1=%s$. Ordinary space expands like $a\\propto\\tau$ while $b$ '
        'stays finite and shrinks at the finite rate $|h_0|$; the cubic and Gauss–Bonnet terms balance ($\\alpha_3E^{(3)}\\approx-\\alpha_2E^{(2)}$).'
        % (g(cub, 0, 10), g(cub, 1, 10), g(cub, 2, 10), g(cub, 3, 10), g(cub, 4, 10)))
    d.p('**Gauss–Bonnet regime** ($\\tau_{32}\\ll\\tau\\ll\\ell_2$). Expanding the curve about $u_2$: $p_{\\rm eff}=p_2+f_1/H_a^2+f_2\\,\\xi H_a^2+\\dots$ '
        'with $f_1=%s$, $f_2=%s$ (units $\\alpha_2=1$), so that' % (g(f12, 0, 10), g(f12, 1, 10)))
    d.eq(r'\tau=\frac{p_2}{H_a}+\frac{f_1}{3H_a^3}-f_2\,\xi H_a+s_2,\qquad a\propto(\tau-s_2)^{p_2},\qquad b\propto(\tau-s_2)^{q_2},\qquad '
         r'H_{a,\rm mid}=\Big(\frac{f_1}{f_2\,\xi}\Big)^{1/4},')
    d.p('with the time shift $s_2=%s$ fixed at $H_{a,\\rm mid}=%s$.' % (g(quad, 3, 8), g(quad, 4, 8)))
    d.p('**Einstein regime** ($\\tau\\gg\\ell_2$). Expanding about $u_1$: $u=u_1+d_1H_a^2+\\dots$, $d_1=\\alpha_2(16-11\\sqrt2)/4$ (%s), and '
        '$p_{\\rm eff}=p_1+e_1H_a^2+\\dots$, $e_1=\\alpha_2(1+9\\sqrt2)/7$ (%s), so that' % (g(de, 0, 10), g(de, 1, 10)))
    d.eq(r'\tau=\frac{p_1}{H_a}-e_1H_a+s_1,\qquad H_a=\frac{p_1}{\tau-s_1}-\frac{e_1p_1^2}{(\tau-s_1)^3}+\dots,\qquad a\propto(\tau-s_1)^{p_1},\qquad b\propto(\tau-s_1)^{q_1},')
    d.p('with $s_1=\\tau(u_2)-p_1/H_a(u_2)+\\int_0^{H_a(u_2)}\\big(1/F-p_1/H^2\\big)dH=%s$ from the exact quadrature ($\\tau(u_2)=%s$, '
        '$H_a(u_2)=%s$).' % (g(quad, 2, 10), g(quad, 0, 10), g(quad, 1, 10)))
    d.p('**Crossovers and accuracy.** The cubic → Gauss–Bonnet crossover is at $\\tau_{32}\\approx|q_2|/|h_0|=|q_2|\\sqrt{30\\xi}\\,\\ell_2$, the '
        'Gauss–Bonnet → Einstein crossover at $\\tau_{21}\\approx\\ell_2$; on the exact solution $\\tau(H_a=10^3)=%s$ and $\\tau(H_a=1)=%s$. '
        'Against the exact solution (integrated from cubic-series data at $\\tau=10^{-8}$ in the variable $\\sigma=-\\ln H_a$, and checked against '
        'the quadrature), the largest relative errors of $\\tau(H_a)$ are:' % (g(cross, 0, 6), g(cross, 1, 6)))
    d.table(['approximation', 'window', 'largest relative error'],
            [['cubic series ($A_0$, $A_1$; $h_0$, $h_1$, $h_2$)', '$H_a\\ge10^4$', g(acc, 0, 2)],
             ['Gauss–Bonnet Kasner with corrections and $s_2$', '$3\\le H_a\\le300$', g(acc, 1, 2)],
             ['Einstein Kasner with correction and $s_1$', '$H_a\\le0.1$', g(acc, 2, 2)],
             ['best of the three, in the crossover windows', '$10^{-3.5}\\le H_a\\le10^8$', g(acc, 3, 2)]], widths=[0.45, 0.25, 0.30], small=True)
    d.p('Inside each regime the formulas are accurate to better than $10^{-3}$; in the two crossover windows the exact curve-plus-quadrature '
        'solution of section 6.1 is used.')
    d.h2('6.3 Coupled to dirac16complex: the matter-dominated cubic regime')
    d.p('Couple the gas of section 4 (isotropic, $a=b=1$, at $H_a=1$, where $\\kappa\\varepsilon=10^{-6}\\times3H_a^2$). Near $\\tau=0$ any fluid with '
        '$p_a=w_a\\varepsilon$ ($w_a<1$) wins over the vacuum terms. With $b$ finite and $H_a=p/\\tau$ the leading orders of the three equations are')
    d.align([r'\tau^3\big(\textstyle\sum\alpha_kE^{(k)0}{}_0\big)&\to48\,h\,p^3(30\alpha_3h^2-\alpha_2),\qquad \tau^3\big(\textstyle\sum\alpha_kE^{(k)1}{}_1\big)\to32\,h\,p^2(p-1)(18\alpha_3h^2-\alpha_2),',
             r'\tau^4\big(\textstyle\sum\alpha_kE^{(k)4}{}_4\big)&\to p^3(p-1)(648\alpha_3h^2-12\alpha_2),'])
    d.p('so that ($p_b$ is of lower order)')
    d.eq(r'h_*=-\sqrt{\frac{\alpha_2}{54\alpha_3}},\qquad p=\frac1{1+w_a},\qquad \kappa\varepsilon\,\tau^3\to\frac{64}3\,p^3\alpha_2|h_*|,\qquad '
         r'\alpha_2E^{(2)}:\alpha_3E^{(3)}:\kappa\varepsilon\to\frac94:-\frac54:1 .')
    d.p('For the gas $w_a\\to1/3$ as $a/b\\to0$: $a\\propto\\tau^{3/4}$ and $\\kappa\\varepsilon\\tau^3\\to9\\alpha_2|h_*|$. The coupled numerical solution '
        'gives at $H_a=10^{15}$: $p_{\\rm eff}=%s$, $H_b=%s$ ($h_*=%s$), $\\kappa\\varepsilon\\tau^3=%s$ ($9|h_*|=%s$), $\\alpha_2E^{(2)}/\\kappa\\varepsilon=%s$, '
        '$\\alpha_3E^{(3)}/\\kappa\\varepsilon=%s$.' % (g(near, 0, 8), g(near, 1, 8), g(hst, 0, 8), g(near, 2, 8),
                                                     (fmt(9*abs(hst[0]), 8) if hst else '?'), g(near, 3, 7), g(near, 4, 7)))
    d.p('**The size of each term.** In a Kasner regime the monomials of the dominant curvature term cancel each other (that is what '
        'makes $P_E(u_1)=0$ or $P_G(u_2)=0$), so its net value would hide it. The SIZE of a curvature term of the 00 equation is therefore '
        'the sum of the absolute values of its monomials in $H_a$ and $H_b$, and the GAS FRACTION is $\\kappa\\varepsilon$ divided by the sum of '
        'all sizes:')
    d.align([r'S_1&=3\big(H_a^2+4H_a|H_b|+2H_b^2\big),\qquad S_2=12\alpha_2|H_b|\big(4H_a^3+18H_a^2|H_b|+12H_aH_b^2+|H_b|^3\big),',
             r'S_3&=360\alpha_3H_a^2|H_b|^3\big(4H_a+3|H_b|\big),\qquad f_{\rm gas}=\frac{\kappa\varepsilon}{S_1+S_2+S_3+\kappa\varepsilon}\quad(H_a>0).'])
    d.p('Moving forward in time the gas fraction DECREASES in the cubic regime ($\\varepsilon\\propto\\tau^{-4}$ against $\\tau^{-3}$) and in the '
        'Gauss–Bonnet regime ($\\varepsilon\\propto\\tau^{-4(p_2+q_2)}$ against $\\tau^{-4}$), so the vacuum regimes of section 6.2 follow: '
        '$(p_{\\rm eff},f_{\\rm gas})$ = (%s, %s) at $H_a=3\\times10^5$, (%s, %s) at $H_a=30$, (%s, %s) at $H_a=10^{-2}$. In the Einstein regime '
        'the auxiliary momenta are blue-shifted ($p_b\\to\\varepsilon/4$, $p_a\\to0$) and $\\kappa\\varepsilon/(3H_a^2)$ — and with it $f_{\\rm gas}$, since there '
        '$S_1\\approx3(1+4|u_1|+2u_1^2)H_a^2$ — GROWS like $\\tau^{2-3p_1-5q_1}=\\tau^{%s}$ (measured between $H_a=10^{-2}$ and $10^{-4}$: %s).' % (g(plat, 0, 6), g(plat, 1, 3), g(plat, 2, 6), g(plat, 3, 3), g(plat, 4, 6), g(plat, 5, 3),
                             g(grow, 1, 6), g(grow, 0, 6)))
    d.h2('6.4 Without a stabilizing potential the auxiliary space re-expands')
    d.p('Einstein gravity with a fluid $p_a=w_a\\varepsilon$, $p_b=w_b\\varepsilon$ has the power-law solutions $a\\propto\\tau^P$, $b\\propto\\tau^Q$ with '
        '$3P(1+w_a)+4Q(1+w_b)=2$ and')
    d.align([r'\tau^2E^{(1)0}{}_0&=3(P^2+4PQ+2Q^2)=\mathcal C,\qquad \tau^2E^{(1)1}{}_1=3P^2+8PQ-2P+10Q^2-4Q=-w_a\mathcal C,',
             r'\tau^2E^{(1)4}{}_4&=3(2P^2+3PQ-P+2Q^2-Q)=-w_b\mathcal C\qquad(\mathcal C=\kappa\varepsilon\tau^2):'])
    d.table(['fluid $(w_a,w_b)$', '$(P,Q)$', '$\\mathcal C=\\kappa\\varepsilon\\tau^2$', 'meaning'],
            [['$(0,\\tfrac14)$ (auxiliary radiation)', '$(0,\\tfrac25)$', '$\\tfrac{24}{25}$', 'ordinary space stops, the auxiliary space expands'],
             ['$(\\tfrac13,0)$ (brane radiation)', '$(\\tfrac12,0)$', '$\\tfrac34$', 'ordinary radiation era, $b$ frozen'],
             ['$(\\tfrac17,\\tfrac17)$ (isotropic radiation)', '$(\\tfrac14,\\tfrac14)$', '$\\tfrac{21}{16}$', 'isotropic expansion']],
            widths=[0.3, 0.14, 0.18, 0.38], small=True)
    d.p('The coupled solution confirms the first line: at $H_a=10^{-14}$, $\\tau H_a=%s$, $\\tau H_b=%s$, $p_a/\\varepsilon=%s$, $p_b/\\varepsilon=%s$. '
        'Once the gas dominates, the 7+1 universe DECOMPACTIFIES unless a stabilizing potential holds $b$ (section 7.1). The Standard Model plasma '
        'cannot do it: confined to the brane at $T\\le M_8$, its smeared density is a fraction '
        '$(\\pi^2g_*/30)M_8^4/(3M_{\\rm Pl}^2M_8^2)=%s$ of the curvature terms at $H=M_8$ ($M_{\\rm Pl}=2.435\\times10^{18}$ GeV, $g_*=106.75$).'
        % (g(late, 0, 4), g(late, 1, 6), g(late, 2, 3), g(late, 3, 6), g(brane, 0, 3)))
    d.h2('6.5 Pair creation of dirac16complex in the primordial regime')
    d.p('A massless mode on the Einstein–Kasner background obeys $i\\,dc/d\\tau=(K\\tau^{-p_1}\\sigma_1+K\\tau^{-q_1}\\sigma_2)c$ in the scaling '
        'variable $K$ (Pauli matrices $\\sigma_{1,2}$); the Bogoliubov coefficient $F(K)=|\\beta|^2$ is the probability that a positive-energy state '
        'at early times ends in a negative-energy state at late times. $F(1)=%s$, $F(8)=%s$, $F(16)=%s$: for large $K$ it falls as a POWER,'
        % (g(Fk, 0, 8), g(Fk, 1, 6), g(Fk, 2, 6)))
    d.eq(r'F(K)\propto K^{-2(p_1-q_1)/(1-p_1)}=K^{-(3\sqrt2+2)/2}=K^{-3.1213\ldots}\qquad(\text{measured between }K=8\text{ and }16:\ %s),'
         % g(tail, 0, 6))
    d.p('because near $x_0=0$ the phase grows only like $K\\tau^{1-p_1}$ while the direction of the mode Hamiltonian turns like $\\tau^{p_1-q_1}$: '
        'the non-adiabaticity comes from the singularity itself. The modes that become non-adiabatic at the time $\\tau$ have physical momenta '
        'of order $1/\\tau$, so their physical number density is of order $\\tau^{-7}$ (seven space dimensions); per comoving volume, with '
        '$a^3b^4=\\tau^{3p_1+4q_1}=\\tau$, the number of pairs created per logarithmic interval of time grows like')
    d.eq(r'\frac{dN_{\rm comoving}}{d\ln\tau}\propto\tau^{-7}\,\tau^{3p_1+4q_1}=\tau^{-6}\qquad(\tau\to0):')
    d.p('the production is dominated by the earliest, Planckian times, where it gives a gas fraction of order one — the gas of section 6.3.')

    # ==============================================================================================
    d.h1('7. Regime II: the modern regime (section 22)')
    d.h2('7.1 The stabilizing potential catches the auxiliary space')
    d.p('Add $U=\\tfrac{K_U}2\\beta^2$ ($U_0=0$, $\\kappa K_U=12m_r^2$, $m_r=10^{-3}/\\ell_2$) to the coupled solution of section 6.3 at $H_a=10^{-2}$ '
        '($\\tau=%s$, $\\beta=%s$), with $b_0$ where the Kasner contraction would arrive at $\\tau=1/m_r$. When the expansion rate falls below '
        '$m_r$ the radion is caught and oscillates; its energy then behaves as four-dimensional matter, and the radion equation of section 5 in '
        'a matter-dominated universe, $\\beta\'\'+\\frac2\\tau\\beta\'+m_r^2\\beta=0$, has the exact solution that is the accurate approximation here:'
        % (g(cap0, 0, 6), g(cap0, 1, 6)))
    d.eq(r'\beta\approx\frac{A\sin(m_r\tau+\vartheta)}{m_r(\tau-\tau_E)},\qquad H_E=H_a+2H_b\approx\frac{2}{3(\tau-\tau_E)},\qquad\text{corrections }O\big(1/(m_r\tau)\big).')
    d.p('With $\\tau_E$ from the average of $H_E$ over the last period ($\\tau_E=%s$) and $A,\\vartheta$ by linear least squares on '
        '$10^5\\le\\tau\\le3\\times10^5$: $A=%s$, $\\vartheta=%s$, largest deviation %s of the amplitude; '
        '$\\langle H_E\\rangle(\\tau-\\tau_E)=%s$ and %s at $\\tau=10^5$ and $2\\times10^5$ (exact value 2/3). The gas share stays constant (%s, %s): '
        'with frozen auxiliary momenta (Kaluza–Klein masses) the gas is matter too. The radion carries %s of the energy: the cosmological '
        'MODULI PROBLEM of the minimal model.' % (g(capfit, 2, 7), g(capfit, 0, 7), g(capfit, 1, 7), g(capfit, 3, 2), g(capHE, 0, 6),
                                                 g(capHE, 1, 6), g(capgas, 0, 4), g(capgas, 1, 4), g(caprad, 0, 4)))
    d.h2('7.2 Radion decay and the Kaluza–Klein channels')
    d.p('In physical units ($M_8=10$ TeV, $M_{\\rm Pl}=2.435\\times10^{18}$ GeV, $V_\\perp=M_{\\rm Pl}^2/M_8^6=%s$ GeV$^{-4}$, $\\ell_\\perp=V_\\perp^{1/4}=%s$ m, '
        'lowest Kaluza–Klein mass $2\\pi/\\ell_\\perp=%s$ GeV) the coupling $\\mathcal L_{\\phi W}$ of section 3 gives'
        % (g(vperp, 0, 6), g(vperp, 1, 6), g(vperp, 2, 6)))
    d.eq(r'\Gamma(\phi\to f\bar f)=\frac{N_c\,m_f^2\,m_r}{24\pi M_{\rm Pl}^2}\Big(1-\frac{4m_f^2}{m_r^2}\Big)^{3/2},\qquad '
         r'\Gamma(\phi\to gg)=\frac{b_{\rm QCD}^2\alpha_s^2\,m_r^3}{96\pi^3M_{\rm Pl}^2},\qquad T_{\rm RH}=\Big[\frac{\Gamma M_{\rm Pl}}{(\pi^2g_*/90)^{1/2}}\Big]^{1/2}.')
    d.p('For $m_r=10$ GeV ($b,c,\\tau,\\mu$ pairs and gluons, $\\alpha_s=0.18$, $n_f=5$, $g_*=10.75$) the brane channels give $\\Gamma_{\\rm brane}=%s$ GeV, '
        'a lifetime of %s s, and would reheat the brane only to $T_{\\rm RH}=%s$ GeV — far below the 4 MeV that nucleosynthesis needs.'
        % (g(decay, 0, 4), g(decay, 1, 4), g(decay, 2, 4)))
    bulk = nbv(D, 'bulk')
    d.p('**Decay into the bulk.** The radion also couples to every Kaluza–Klein mode of the bulk fields. A mode of dirac16complex with '
        'Kaluza–Klein momentum $\\mu_n=2\\pi|n|/\\ell_\\perp$ ($n\\in\\mathbb Z^4$) has, in the four-dimensional Einstein frame ($h=e^{-4\\beta}g_E$), the mass')
    d.eq(r'M_E=e^{-2\beta}\sqrt{m^2+\mu_n^2e^{-2\beta}},\qquad g_n=-\frac{\partial M_E/\partial\beta}{\sqrt{12}M_{\rm Pl}}\Big|_{\beta=0}=\frac{2M_n+\mu_n^2/M_n}{\sqrt{12}\,M_{\rm Pl}},\qquad M_n=\sqrt{m^2+\mu_n^2},')
    d.p('(a brane fermion has the coupling $m_f/(\\sqrt3M_{\\rm Pl})$ of $\\mathcal L_{\\phi W}$, the case $\\mu_n=0$). A scalar coupled as $-g\\phi\\bar\\chi\\chi$ '
        'decays at the rate $g^2m_r(1-4M^2/m_r^2)^{3/2}/(8\\pi)$; for $m\\ll\\mu_n$ each Kaluza–Klein level (four brane Dirac spinors) gives '
        '$\\Gamma_n=3\\mu_n^2m_r(1-4\\mu_n^2/m_r^2)^{3/2}/(8\\pi M_{\\rm Pl}^2)$. The sum over the %s momenta with $2\\mu_n<m_r$, done shell by shell '
        'with the number $r_4(k)=8\\sigma(k)-32\\sigma(k/4)$ of lattice vectors with $|n|^2=k$ ($\\sigma$ the divisor sum; $\\sigma(k/4)=0$ unless '
        '$4\\mid k$; checked against a direct count for $k\\le60$), equals its continuum limit:' % g(nkk, 0, 4))
    d.eq(r'\Gamma_{\rm bulk}=\sum_n\Gamma_n\approx\frac{3m_r}{8\pi M_{\rm Pl}^2}\,\frac{V_\perp}{8\pi^2}\Big(\frac{m_r}2\Big)^6\frac{B(3,\tfrac52)}2,\qquad '
         r'B(3,\tfrac52)=\frac{\Gamma(3)\Gamma(\tfrac52)}{\Gamma(\tfrac{11}2)}=\frac{16}{315},')
    d.p('$\\Gamma_{\\rm bulk}=%s$ GeV (continuum %s GeV), a lifetime of %s s: %s times faster than the brane channels (the Kaluza–Klein '
        'gravitons add channels of the same strength; for $m>m_r/2$ the dirac16complex channels close). The radion decays into the BULK: its '
        'energy goes to Kaluza–Klein modes, which return it to the brane only through couplings suppressed by the volume (a Kaluza–Klein '
        'graviton of mass $M$ decays into brane fields at a rate of order $M^3/M_{\\rm Pl}^2$). This is the cosmological moduli problem of '
        'the minimal model in the form it takes with large extra dimensions; the model needs the radion and bulk energy left by regime I to be '
        'removed before nucleosynthesis (for example by a later period of inflation on the brane, which is not part of the model). The modern '
        'solutions below start after such a reheating, below the normalcy temperature $T_*\\approx(M_8^6/M_{\\rm Pl})^{1/5}=%s$ GeV (above it the '
        'brane plasma would lose its energy to bulk gravitons faster than it cools by expansion; an order-of-magnitude estimate).'
        % (g(bulk, 0, 6), g(bulk, 1, 6), g(bulk, 2, 3), g(bulk, 3, 2), g(tstar, 0, 4)))
    d.h2('7.3 The Standard Model plasma: $g_*(T)$, $g_{*s}(T)$ and the time–temperature relation')
    d.p('With $b$ static the brane obeys $3M_{\\rm Pl}^2H_a^2=\\rho_4+M_{\\rm Pl}^2\\Lambda_4$ and $H_a\'=-(\\rho_4+p_4)/(2M_{\\rm Pl}^2)$. In the radiation '
        'era $\\rho_4=(\\pi^2/30)g_*T^4$ and the entropy $s=(2\\pi^2/45)g_{*s}T^3$ is conserved, $sa^3=$ const:')
    d.align([r'H_a&=\Big(\frac{\pi^2g_*}{90}\Big)^{1/2}\frac{T^2}{M_{\rm Pl}},\qquad \ln a=-\ln T-\tfrac13\ln g_{*s}(T)+\text{const},',
             r't(T)&=t(T_1)+\int\frac{d\ln a}{H_a}=t(T_1)+\int_T^{T_1}\frac{d\tilde T}{\tilde T}\,\frac{1+\frac13\,d\ln g_{*s}/d\ln\tilde T}{H_a(\tilde T)}.'])
    d.p('Each species $i$ ($x=m_i/T$, degeneracy $g_i$; bosons $+$, fermions with alternating signs) contributes')
    d.align([r'\rho_i&=\frac{g_iT^4}{2\pi^2}\sum_{n\ge1}(\pm1)^{n+1}\frac{3(nx)^2K_2(nx)+(nx)^3K_1(nx)}{n^4},',
             r'p_i&=\frac{g_iT^4}{2\pi^2}\sum_{n\ge1}(\pm1)^{n+1}\frac{(nx)^2K_2(nx)}{n^4},\qquad s_i=\frac{\rho_i+p_i}T,'])
    d.p('$g_*=30\\rho/(\\pi^2T^4)$, $g_{*s}=45s/(2\\pi^2T^3)$, for all Standard Model particles with their masses; the QCD crossover is modelled by '
        'the switch $S(T)=\\tfrac12[1+\\tanh((T-0.160\\,\\mathrm{GeV})/0.030\\,\\mathrm{GeV})]$ between quarks and gluons and the lightest hadrons '
        '($\\pi,K,\\eta,\\rho,\\omega,K^*,\\eta\',p,n$) — this model differs from the lattice equation of state by up to about 10%% near the '
        'crossover; the neutrinos decouple while $e^\\pm$ are relativistic, $T_\\nu=T\\,(g_{\\gamma e}(T)/\\tfrac{11}2)^{1/3}$, and are scaled by '
        '$N_{\\rm eff}/3=3.044/3$. The integral is evaluated on 481 points equally spaced in $\\ln T$ (10 GeV to 10 keV), starting from '
        '$t(T_1)=1/(2H_a(T_1))$ at $T_1=10$ GeV; on each step $1/H_a$ is taken to be an exponential function of $\\ln a$ (exact where $g_*$ and '
        '$g_{*s}$ are constant), so a step adds $\\Delta\\ln a\\,(1/H_{i+1}-1/H_i)/\\ln(H_i/H_{i+1})$ (halving the grid changes $t(10\\,\\mathrm{keV})$ by %s).'
        % (g(thalf, 0, 2)))
    trow_T = ['1 TeV', '10 GeV', '1 GeV', '300 MeV', '100 MeV', '10 MeV', '1 MeV', '100 keV', '1 keV']
    if len(gst) == 18:
        d.table(['$T$'] + trow_T, [['$g_*$'] + [fmt(gst[2*i], 6) for i in range(9)], ['$g_{*s}$'] + [fmt(gst[2*i + 1], 6) for i in range(9)]],
                widths=[0.1] + [0.1]*9, small=True)
    if len(tT) == 6 and len(tinst) == 6:
        d.table(['$T$', '1 GeV', '100 MeV', '10 MeV', '1 MeV', '100 keV', '10 keV'],
                [['$t$ (s)'] + [fmt(v, 6) for v in tT], ['$1/(2H_a)$ divided by $t$'] + [fmt(v, 6) for v in tinst]],
                widths=[0.22] + [0.13]*6, small=True)
    if len(tinst) == 6:
        d.p('The instantaneous approximation $t\\approx1/(2H_a)$ is good to about 1%% where $g_*$ is nearly constant; where $g_*$ changes quickly it '
            'is off by %s%% at 100 MeV (QCD crossover) and by %s%% at 100 keV ($e^\\pm$ annihilation).' % (fmt(100*abs(1 - tinst[1]), 2), fmt(100*abs(1 - tinst[4]), 2)))
    d.h2('7.4 The late universe: an accurate analytic solution and the age')
    d.p('Below about 1 MeV the brane contains photons and neutrinos ($\\Omega_r$ from $T_{\\rm CMB}=2.7255$ K and $N_{\\rm eff}=3.044$: '
        '$\\Omega_r=\\Omega_\\gamma[1+\\tfrac78(\\tfrac4{11})^{4/3}N_{\\rm eff}]$, $\\Omega_\\gamma=(\\pi^2/15)T_{\\rm CMB}^4/(3H_0^2M_{\\rm Pl}^2)$), matter '
        '($\\Omega_m=0.315$; its cold dark matter would be dirac16complex if dirac16complex is the dark matter) and $\\Lambda_4$ '
        '($\\Omega_\\Lambda=1-\\Omega_m-\\Omega_r$), with $H_0=67.4$ km/s/Mpc (Planck 2018). The Friedmann equation '
        '$H_a=H_0(\\Omega_ra^{-4}+\\Omega_ma^{-3}+\\Omega_\\Lambda)^{1/2}$ has exact two-component solutions, and their composite')
    d.align([r't_{rm}(a)&=\frac{2}{3H_0\Omega_m^2}\Big[(\Omega_ma-2\Omega_r)\sqrt{\Omega_r+\Omega_ma}+2\Omega_r^{3/2}\Big],\qquad '
             r't_{m\Lambda}(a)=\frac{2}{3H_0\Omega_\Lambda^{1/2}}\,\mathrm{arcsinh}\Big(\sqrt{\tfrac{\Omega_\Lambda}{\Omega_m}}\,a^{3/2}\Big),',
             r't(a)&\approx t_{rm}(a)+t_{m\Lambda}(a)-t_m(a),\qquad t_m(a)=\frac{2}{3H_0\Omega_m^{1/2}}\,a^{3/2},'])
    d.p('is accurate to %s at every time (the error is the product of the small radiation and $\\Lambda$ corrections). Numbers: $\\Omega_r=%s$, '
        '$\\Omega_\\Lambda=%s$, $z_{\\rm eq}=\\Omega_m/\\Omega_r-1=%s$, age $=\\int_0^1da/(aH_a)=%s$ Gyr (Planck 2018: $13.787\\pm0.020$).'
        % (g(comp, 0, 2), g(late4, 0, 4), g(late4, 1, 6), g(late4, 2, 6), g(late4, 3, 6)))
    d.p('**The radion today** sits at $\\beta_{\\rm eq}=[\\Lambda_4+\\kappa f_W/V_\\perp+\\tfrac\\kappa2(\\varepsilon-3p_a)+\\kappa p_b-12\\alpha_2H_a^2(H_a\'+H_a^2)]/(3m_r^2)$: '
        'for $m_r=10$ GeV, $%s$ without brane tension and $%s$ with $f_W=M_8^4$ — so $G_4=G_8/V_\\perp$ is constant. The Gauss–Bonnet and cubic '
        'corrections are $\\alpha_2H_0^2=%s$ and $\\alpha_3H_0^4=%s$ ($\\alpha_2=1/M_8^2$, $\\alpha_3=10^{-6}/M_8^4$).'
        % (g(radnow, 0, 3), g(radnow, 1, 3), g(gbnow, 0, 3), g(gbnow, 1, 3)))

    # ==============================================================================================
    d.h1('8. Regime III: the far future (section 23)')
    d.p('Once matter has diluted, $\\Lambda_4$ dominates: $H_\\infty=H_0\\Omega_\\Lambda^{1/2}=%s$ s$^{-1}$ (%s GeV), Gibbons–Hawking temperature '
        '$H_\\infty/(2\\pi)=%s$ K. For matter $+\\Lambda$ the solution is exact,' % (g(hinf, 0, 5), g(hinf, 1, 5), g(hinf, 2, 4)))
    d.eq(r'a(t)=\Big(\frac{\Omega_m}{\Omega_\Lambda}\Big)^{1/3}\sinh^{2/3}\Big(\frac32H_\infty t\Big),\qquad H_a(t)=H_\infty\coth\Big(\frac32H_\infty t\Big)\to H_\infty+2H_\infty e^{-3H_\infty t},')
    d.p('and agrees with the full solution (radiation included) to $10^{-4}$ already at $a=3$. Pairs of dirac16complex are created with the '
        'Fermi–Dirac occupation $|\\beta|^2=1/(e^{2\\pi M/H}+1)$; $2\\pi M/H_\\infty=%s$ for $M=1$ eV: none are.' % g(gh, 0, 3))
    d.p('**The radion in de Sitter.** Linearizing the full Einstein–Lovelock equations (with $U=\\tfrac{K_U}2\\beta^2$) about $H_a=H$, $H_b=0$, '
        '$\\beta=\\beta_{\\rm eq}$, with $\\delta H_a$ removed by the constraint, gives $\\beta\'\'+\\gamma\\beta\'+m_{\\rm eff}^2(\\beta-\\beta_{\\rm eq})=0$ with (for $\\beta_{\\rm eq}\\to0$)')
    d.align([r'\frac{\gamma}{3H}&=\frac{1-16\alpha_2H^2+48\alpha_2^2H^4-48\alpha_3H^4+512\alpha_2^3H^6+192\alpha_2\alpha_3H^6-1024\alpha_2^4H^8}{(1-4\alpha_2H^2+32\alpha_2^2H^4)^2},',
             r'\frac{m_{\rm eff}^2}{m_r^2}&=\frac1{1-4\alpha_2H^2+32\alpha_2^2H^4},\qquad 3m_r^2\beta_{\rm eq}=3H^2-12\alpha_2H^4+O(\beta_{\rm eq}^2).'])
    d.p('For $\\alpha_2=\\alpha_3=0$: $\\gamma=3H$, $m_{\\rm eff}=m_r$. A numerical test in units $\\ell_2=1$ ($\\Lambda=3\\times10^{-4}$, $m_r=0.1$, $\\xi=10^{-6}$, '
        'start 0.05 above equilibrium) relaxes to $\\beta=%s$, $H_a=%s$, the solution of $3H^2=\\Lambda+6m_r^2\\beta^2$, $3m_r^2\\beta=3H^2-12H^4$ '
        '($%s$, $%s$). In the physical far future $\\alpha_2H_\\infty^2\\approx10^{-92}$ and $\\beta_{\\rm eq}\\approx H_\\infty^2/m_r^2$: the auxiliary space '
        'stays frozen. Without the stabilizing potential the far future would instead be eight-dimensional de Sitter, $a=b=e^{Hx_0}$, an exact '
        'vacuum solution when all three equations read' % (g(ds, 0, 10), g(ds, 1, 10), g(ds, 2, 10), g(ds, 3, 10)))
    d.eq(r'21H^2-420\alpha_2H^4+2520\alpha_3H^6=\Lambda ;')
    d.p('the auxiliary space would then inflate and $G_4\\propto b^{-4}$ would fall like $e^{-4Hx_0}$, which experiments exclude '
        '($|\\dot G/G|$ is less than a few percent of $H_0$).')

    # ==============================================================================================
    d.h1('9. The size of every interaction in the three regimes')
    d.p('**Regime I.** The shares of the four terms of the 00 equation at four points of the coupled primordial solution of section 6.3: '
        'each size ($S_1$, $S_2$, $S_3$ of section 6.3, and $\\kappa\\varepsilon$) divided by $S_1+S_2+S_3+\\kappa\\varepsilon$:')
    d.table(['point of regime I', 'Einstein $S_1$', 'Gauss–Bonnet $S_2$', 'cubic $S_3$', 'gas $\\kappa\\varepsilon$'],
            budget_rows(D), widths=[0.32, 0.17, 0.17, 0.17, 0.17], small=True)
    d.p('In the cubic regime the Gauss–Bonnet and cubic terms balance each other; in the Gauss–Bonnet and Einstein regimes one term carries '
        'the 00 equation alone (its monomials cancel among themselves, which is what fixes $u_2$ and $u_1$); near the singularity the gas '
        'joins the Gauss–Bonnet and cubic terms in the ratio $\\alpha_2E^{(2)}:\\alpha_3E^{(3)}:\\kappa\\varepsilon=9/4:-5/4:1$ of section 6.3.')
    # numbers for the table of all interactions (computed here from the printed values)
    MPl, M8, mr2 = 2.435e18, 1.0e4, 100.0
    Hinf = hinf[1] if len(hinf) > 1 else None
    OmL = late4[1] if len(late4) > 1 else None
    H0 = Hinf/math.sqrt(OmL) if (Hinf and OmL) else None
    lam_I = OmL*gbnow[0] if (OmL and gbnow) else None                     # Lambda4/(3H^2) at H = 1/l2 = M8: Omega_L H0^2/M8^2
    fw_I = M8**2/(3*MPl**2)                                                  # kappa f_W/(3 H^2 V_perp) at H = M8, f_W = M8^4
    a2_III = (Hinf/M8)**2 if Hinf else None; a3_III = 1e-6*(Hinf/M8)**4 if Hinf else None
    bet_III = Hinf**2/mr2 if Hinf else None
    gh_II = 2*math.pi*1e-9/H0 if H0 else None
    d.p('**All interactions.** Every term of the action of section 3 with its weight in the field equations of each regime ($M_8=10$ TeV, '
        '$\\alpha_2=1/M_8^2$, $\\alpha_3=10^{-6}/M_8^4$, $m_r=10$ GeV, $f_W=M_8^4$, $M_{\\rm Pl}=2.435\\times10^{18}$ GeV; regime I at $H=1/\\ell_2=M_8$; '
        '"estimate" marks an order-of-magnitude statement):')
    d.table(['interaction', 'regime I (primordial)', 'regime II (today)', 'regime III (far future)'],
            [['$S_0$: $\\Lambda_4/(3H^2)$', '$\\Omega_\\Lambda\\alpha_2H_0^2=%s$ (if $\\Lambda\\approx\\Lambda_4$)' % fmt(lam_I, 2),
              '$\\Omega_\\Lambda=%s$' % g(late4, 1, 4), '$\\to1$'],
             ['$S_1$: Einstein term', 'dominant for $\\tau\\gg\\ell_2$', 'dominant (four-dimensional Friedmann equations)', 'dominant'],
             ['$S_2$: $\\alpha_2H^2$', 'dominant for $\\tau_{32}\\ll\\tau\\ll\\ell_2$', g(gbnow, 0, 3), fmt(a2_III, 3)],
             ['$S_3$: $\\alpha_3H^4$', 'dominant for $\\tau\\ll\\tau_{32}$', g(gbnow, 1, 3), fmt(a3_III, 3)],
             ['$S_\\psi$: dirac16complex (gas; matter)', 'see the table above; dominant near $\\tau=0$', '$\\Omega_m=0.315$ (its cold dark matter, $\\Omega_c=\\Omega_ch^2/h^2=0.120/0.674^2=%s$, if it is the dark matter)' % fmt(0.120/0.674**2, 3),
              '$\\to0$ like $e^{-3H_\\infty t}$'],
             ['$S_m$: $m/H$', '$\\ll1$ (massless gas)', '$\\gg1$ (cold matter)', '$\\gg1$'],
             ['$S_{\\rm self}$: $g_rn^{6/7}\\approx(k_{\\rm phys}/M_8)^6$ (estimate)', '$<1$ below the eight-dimensional Planck scale', 'negligible', 'negligible'],
             ['$S_U$: radion displacement $\\beta$', 'catches $b$ when $H\\approx m_r$, then oscillates', '$\\beta_{\\rm eq}=%s$' % g(radnow, 0, 3),
              '$\\beta_{\\rm eq}\\approx H_\\infty^2/m_r^2=%s$' % fmt(bet_III, 3)],
             ['$S_f$: $\\kappa f_W/(3H^2V_\\perp)$', fmt(fw_I, 3), 'part of $\\Lambda_4$; shifts $\\beta_{\\rm eq}$ by %s' % g(radnow, 1, 3), 'part of $\\Lambda_4$'],
             ['$S_{\\rm SM}$: brane plasma', g(brane, 0, 3) + ' (at $T=M_8$)', '$\\Omega_r=%s$; baryons $\\Omega_b\\approx0.049$' % g(late4, 0, 3), '$\\to0$ like $e^{-4H_\\infty t}$'],
             ['$S_P$: portals (estimate)', 'negligible (on the brane)', 'production of dirac16complex (freeze-in)', 'negligible'],
             ['$\\mathcal L_{\\phi W}$: radion–brane coupling', 'radion not yet caught', 'brane channels alone: $T_{\\rm RH}=%s$ GeV ($m_r=10$ GeV)' % g(decay, 2, 3), '—'],
             ['radion–bulk couplings $g_n$', 'radion not yet caught', 'decay into bulk Kaluza–Klein modes in %s s' % g(nbv(D, 'bulk'), 2, 2), '—'],
             ['pair creation, $2\\pi M/H$ ($M=1$ eV)', 'of order one (Planckian times)', fmt(gh_II, 3), g(gh, 0, 3)]],
            widths=[0.25, 0.25, 0.25, 0.25], small=True)

    # ==============================================================================================
    d.h1('10. What the solutions mean')
    d.ul(['**The primordial regime is under control.** Its vacuum solution is exact (section 6.1), and its three Lovelock sub-regimes have '
          'accurate closed forms. Near $x_0=0$ any dirac16complex gas makes the universe matter-dominated with $a\\propto\\tau^{3/4}$ and a finite '
          'contraction rate $h_*$ of the auxiliary space; pair creation from the singularity supplies such a gas.',
          '**Dynamical compactification is temporary.** In the Gauss–Bonnet and Einstein regimes the auxiliary space contracts ($q_2,q_1<0$) while '
          'ordinary space expands; but once the bulk gas dominates, the auxiliary space re-expands ($b\\propto\\tau^{2/5}$). The Standard Model '
          'brane, at $T\\le M_8$, is too dilute in the bulk to stop it. A stabilizing potential is a necessary part of the model.',
          '**Stabilization costs a moduli problem.** The potential catches $b$ when $H\\approx m_r$, but the radion then carries almost all the energy '
          '(%s). Through brane fields alone it would decay far too late for nucleosynthesis; with its $10^{13}$ channels of gravitational '
          'strength into the bulk Kaluza–Klein modes it decays within %s s, but its energy then sits in the bulk and returns to the brane only '
          'through volume-suppressed couplings. The minimal model needs an additional dilution and reheating of the brane after stabilization.'
          % (g(caprad, 0, 3), g(nbv(D, 'bulk'), 2, 2)),
          '**The modern regime is standard cosmology on the brane.** With $b$ fixed (to $10^{-86}$–$10^{-24}$) gravity is four-dimensional with '
          '$G_4=G_8/V_\\perp$, the Lovelock corrections are $10^{-92}$ and smaller, and the brane Friedmann solution reproduces the age of the '
          'universe; dirac16complex enters only through its energy density (cold matter).',
          '**The far future is de Sitter** with a stable, frozen auxiliary space; without the stabilizer it would be eight-dimensional de Sitter '
          '(excluded). No pairs are created.',
          '**Validity.** The classical solutions need curvatures below the eight-dimensional Planck scale ($\\tau\\gg\\kappa^{1/6}$), which requires '
          '$\\ell_2\\gg\\kappa^{1/6}$ for the cubic and Gauss–Bonnet regimes to be classical; the brane plasma needs $T<M_8$ and $T<T_*$; the massless-gas '
          'description of dirac16complex needs $m\\ll H_a$ and negligible self-couplings.'])

    # ==============================================================================================
    d.h1('11. Verification')
    d.p('Every statement above is a check of 7+1.nb (Wolfram Language 15.0.1) and of the independent Python program tools/ref_coupled.py '
        '(numpy, scipy, sympy, mpmath), which uses its own closed forms (from the independent Lovelock engine tools/lovelock71.py), its own '
        'quadrature, root finding and Runge–Kutta integrator (DOP853). The numbers printed by both programs:')
    trows = []
    for lab, a, b, tol, r in cmp_:
        trows.append([lab, ', '.join(fmt(x, 7) for x in a[:6]) + (' …' if len(a) > 6 else '') if a else '?',
                      ', '.join(fmt(x, 7) for x in b[:6]) + (' …' if len(b) > 6 else '') if b else '?',
                      fmt(r, 2) if r is not None else '—', fmt(tol, 1)])
    d.table(['quantity', '7+1.nb', 'Python', 'max. relative difference', 'tolerance'], trows, widths=[0.2, 0.27, 0.27, 0.13, 0.13], small=True)
    d.p('The largest relative difference is %s. The checks of sections 20–23 of the notebook:' % (fmt(worst, 2) if worst is not None else '?'))
    d.table(['#', 'check', 'result', 'section'],
            [[str(i + 1), '`' + c[0].replace('`', "'") + '`', c[1], (c[3][:4] if c[3] else c[2][:3])] for i, c in enumerate(ck)],
            widths=[0.06, 0.74, 0.08, 0.12], small=True)
    d.p('The Python program’s boolean checks:')
    d.table(['#', 'check', 'result'], [[str(i + 1), '`' + k.replace('`', "'") + '`', v] for i, (k, v) in enumerate(rb)],
            widths=[0.06, 0.82, 0.12], small=True)
    d.h1('12. How to reproduce')
    d.ol(['`python3 tools/make_7p1_nb.py` writes 7+1.nb (unevaluated).',
          'Evaluate 7+1.nb in a fresh kernel of Wolfram 15.0.1 (Evaluation ▸ Quit Kernel, then Evaluation ▸ Evaluate Notebook; or `Get["D:/Developer/github/7+1/tools/run_eval.wl"]`) and save it; the last section summarizes all checks.',
          '`python3 tools/ref_coupled.py` repeats the checks of this document independently (exit code 0 if all pass).',
          '`python3 tools/make_docs.py` regenerates this document, its symbol dictionary and its provenance record.'])
    d.h1('13. References')
    d.ul(['D. Lovelock, J. Math. Phys. 12 (1971) 498 — the Lovelock tensors; F. Müller-Hoissen, Phys. Lett. B 163 (1985) 106 — the cubic density.',
          'E. Kasner, Am. J. Math. 43 (1921) 217 — power-law solutions; N. Deruelle and L. Fariña-Busto, Phys. Rev. D 41 (1990) 3696; '
          'I. V. Kirnos, A. N. Makarenko, S. A. Pavluchenko and A. V. Toporensky, Gen. Rel. Grav. 42 (2010) 2633 — Kasner regimes of '
          'Einstein–Gauss–Bonnet cosmology.',
          'N. Arkani-Hamed, S. Dimopoulos and G. Dvali, Phys. Lett. B 429 (1998) 263 and Phys. Rev. D 59 (1999) 086004 — large extra dimensions, '
          'the normalcy temperature; G. D. Coughlan, W. Fischler, E. W. Kolb, S. Raby and G. G. Ross, Phys. Lett. B 131 (1983) 59 — the '
          'cosmological moduli problem.',
          'G. W. Gibbons and S. W. Hawking, Phys. Rev. D 15 (1977) 2738 — the de Sitter temperature $H/2\\pi$.',
          'Planck Collaboration, Astron. Astrophys. 641 (2020) A6 — $H_0$, $\\Omega_m$, the age; S. Navas et al. (Particle Data Group), Phys. Rev. D '
          '110 (2024) 030001 and its 2025 update — particle masses; HotQCD, Phys. Lett. B 795 (2019) 15 — the QCD crossover temperature.'])
    return d


def budget_rows(D):
    """the shares S_1, S_2, S_3, kappa eps (each divided by their sum) printed by section 23.3 of 7+1.nb"""
    v = nbv(D, 'shares')
    labs = ['cubic ($H_a=3\\times10^5$)', 'Gauss–Bonnet ($H_a=30$)', 'Einstein ($H_a=10^{-2}$)', 'near the singularity ($H_a=10^{15}$)']
    if len(v) != 16:
        return [['(not found)', '?', '?', '?', '?']]
    return [[labs[i]] + [fmt(x, 3) for x in v[4*i:4*i + 4]] for i in range(4)]
