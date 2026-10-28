"""Content of DARK_SECTOR (Markdown + LaTeX): can the field dirac16complex be the dark matter, and can it be the dark energy? An
analysis based on the accurate approximate solutions of the coupled field equations in the three gravity regimes (section 24 of
7+1.nb, which uses sections 20-23).  Every equation used is written here; every number is parsed from the evaluated notebook or
from the output of tools/ref_dark.py, or computed here from them."""
import math
from docgen import Doc
from doc_stage3_common import fmt, nums_in
from doc_common import signature_block

# ------------------------------------------------------------------------------------------------ numbers of both programs
NB = {  # short name -> output key of section 24 of 7+1.nb
    'req': 'Omega_c, rho_c/h^2 (GeV/cm^3), s0 (1/cm^3), m Y (GeV)', 'md': 'M_D (GeV)',
    'pauli': 'Pauli bound for g = 16 internal states (keV)', 'dbeta': 'radion displacement by dark matter today',
    'gas': 'gas share after capture (section 22.1) at tau = 10^5, 3 10^5', 'late': 'without stabilizer (section 21.3): tau Hb at Ha = 10^-14',
    'coup': 'lambda_eff, y_eff, G_r per unit coupling (GeV^-1, 1, GeV^-2)', 'sigsi': 'sigma_SI (cm^2) at c_P = 1, m = 40 GeV',
    'sigm': 'sigma/m (cm^2/g) at c_r = 1, m = 1 GeV; c_r allowed by the Bullet Cluster at m = 1 GeV',
    'life': 'portal mass mu_P (GeV) at c_y = 1; lifetimes chi -> nu gamma, chi -> 3 nu (s) at c_y = 1, m = 1 GeV and 10 keV',
    'mmax': 'largest m (GeV) with lifetime > age for c_y = 1',
    'over': 'gas share after capture f, overabundance O at T_RH = 1 GeV and 4 MeV', 'bulk': 'fraction of radion decays into the bulk; brane share',
    'tstar': 'T* (GeV), S(m_h/T*)',
    'yzero': 'zero mode: Yukawa-portal Y at c_y = 1 (N = 1, S = 1), Higgs-portal Y at c_P = 1 (S = 1)',
    'mstar': 'm* (GeV) for N = 1 and N = 12', 'cy': 'c_y needed and allowed at m = 10 keV, 10 MeV',
    'tower': 'Kaluza-Klein tower below m_h: number, mass sum (GeV), Y_tower and mass per entropy (GeV) at c_y = 1 (N = 1, S = 1); c_y bound at T*',
    'kklife': 'tau(chi_n -> 3 nu) (s) at mu_n = 10 GeV, c_y = 1; largest m (GeV) with tau(chi -> 3 nu) > age for c_y = 10^-5',
    'hp': 'Higgs portal pairs: I_hp, I_hp^E, number of pairs, mass per entropy (GeV) at c_P = 1 (N = 4, S = 1); c_P needed for S = 1 and at T*',
    'hpfrac': 'Higgs portal pairs at c_P = 1, T_RH = T*: fraction of the dark matter', 'vec': 'vector portal: T_V (GeV), (T*/T_V)^11',
    'trh': 'T_RH (GeV) = 5, 8, 10, T*: S(m_h/T_RH), c_y bound of the tower (N = 1), c_P needed by the Higgs-portal pairs (N = 4)',
    'wmin': 'w_min of Kohn-Sham matter and its kF (units m = 1); saturation kF0', 'rhol': 'rho_Lambda (GeV^4)',
    'cas': 'Casimir energy / rho_Lambda: dirac16complex (m = 0), graviton (20 dof)',
    'xl': 'm l_perp and m (GeV) at which the dirac16complex Casimir energy equals rho_Lambda',
    'wtab': 'kF, p, w at kF = 0.1, 0.3, 0.5, 0.6, 0.67, 0.7, 0.75 (units m = 1)',
    'castab': 'dirac16complex Casimir energy / rho_Lambda at m l_perp = 1, 10, 30, 60'}
PY = {  # short name -> output key of tools/ref_dark.py
    'req': 'NUM Omega_c, rho_c/h^2 (GeV/cm^3), s0 (1/cm^3), m Y (GeV)', 'md': 'NUM M_D (GeV) = (2 pi)^(2/3) M8',
    'pauli': 'NUM Pauli bound for g = 16 internal states, 0.13 keV (2/16)^(1/4) (keV)',
    'dbeta': 'NUM radion displacement by dark matter today H0^2 Omega_c/(2 m_r^2) (m_r = 10 GeV)',
    'coupled': 'NUM from tools/ref_coupled.py: gas share at tau = 1e5, 3e5; Gamma_brane, Gamma_bulk (GeV); tau Hb at Ha = 1e-14',
    'coup': 'NUM lambda_eff, y_eff, G_r per unit coupling (GeV^-1, 1, GeV^-2)', 'sigsi': 'NUM sigma_SI (cm^2) at c_P = 1, m = 40 GeV',
    'sigm': 'NUM sigma/m (cm^2/g) at c_r = 1, m = 1 GeV; c_r allowed by the Bullet Cluster at m = 1 GeV',
    'life': 'NUM mu_P (GeV) at c_y = 1; lifetimes chi -> nu gamma, chi -> 3 nu (s) at c_y = 1, m = 1 GeV and 10 keV',
    'mmax': 'NUM largest m (GeV) with lifetime chi -> 3 nu > age for c_y = 1',
    'over': 'NUM gas share after capture f, overabundance O at T_RH = 1 GeV and 4 MeV',
    'bulk': 'NUM fraction of radion decays into the bulk; brane share Gamma_brane/Gamma_bulk', 'tstar': 'NUM T* (GeV), S(m_h/T*)',
    'yzero': 'NUM zero mode: Yukawa-portal Y at c_y = 1 (N = 1, S = 1), Higgs-portal Y at c_P = 1 (S = 1)',
    'mstar': 'NUM m* (GeV) for N = 1 and N = 12', 'cy': 'NUM c_y needed and allowed at m = 10 keV, 10 MeV',
    'lat': 'NUM lattice sums for R = 100: Sum r_4(k)(1 - k/R^2)^2 / (pi^2 R^4/12) - 1, Sum r_4(k) sqrt(k)/R (1 - k/R^2)^2 / (16 pi^2 R^4/315) - 1',
    'tower': 'NUM Kaluza-Klein tower below m_h: number, mass sum (GeV), Y_tower and mass per entropy (GeV) at c_y = 1 (N = 1, S = 1); c_y bound at T*',
    'kklife': 'NUM tau(chi_n -> 3 nu) (s) at mu_n = 10 GeV, c_y = 1; largest m (GeV) with tau(chi -> 3 nu) > age for c_y = 1e-5',
    'hp': 'NUM Higgs portal pairs: I_hp, I_hp^E, number of pairs, mass per entropy (GeV) at c_P = 1 (N = 4, S = 1); c_P needed for S = 1 and at T*',
    'hpfrac': 'NUM Higgs portal pairs at c_P = 1, T_RH = T*: fraction of the dark matter', 'vec': 'NUM vector portal: T_V (GeV), (T*/T_V)^11',
    'trh': 'NUM T_RH (GeV) = 5, 8, 10, T*: S(m_h/T_RH), c_y bound of the tower (N = 1), c_P needed by the Higgs-portal pairs (N = 4)',
    'wmin': 'NUM w_min of Kohn-Sham matter and its kF (units m = 1); saturation kF0', 'rhol': 'NUM rho_Lambda (GeV^4)',
    'cas': 'NUM Casimir energy / rho_Lambda: dirac16complex (m = 0), graviton (20 dof)',
    'xl': 'NUM m l_perp and m (GeV) at which the dirac16complex Casimir energy equals rho_Lambda',
    'wtab': 'NUM kF, p, w at kF = 0.1, 0.3, 0.5, 0.6, 0.67, 0.7, 0.75 (units m = 1)',
    'castab': 'NUM dirac16complex Casimir energy / rho_Lambda at m l_perp = 1, 10, 30, 60'}


def nbv(D, k):
    return nums_in(D['nums'].get('dk', {}).get(NB[k], ''))


def pyv(D, k):
    return nums_in(D['ref'].get('dk', ('', '', {}))[2].get(PY[k], ''))


def _sel(v, idx):
    try:
        return [v[i] for i in idx]
    except (IndexError, TypeError):
        return []


def compare(D):
    """(label, notebook values, Python values, tolerance, max relative difference) for every number printed by both programs"""
    rows = [
        ('$\\Omega_c$, $\\rho_c/h^2$, $s_0$, $mY$', nbv(D, 'req'), pyv(D, 'req'), 1e-12),
        ('$M_D$', nbv(D, 'md'), pyv(D, 'md'), 1e-14),
        ('Pauli bound for $g=16$', nbv(D, 'pauli'), pyv(D, 'pauli'), 1e-14),
        ('radion displacement by dark matter', nbv(D, 'dbeta'), pyv(D, 'dbeta'), 1e-12),
        ('gas share after capture (from the coupled solution)', nbv(D, 'gas'), pyv(D, 'coupled')[:2], 1e-8),
        ('$\\tau H_b$ without stabilizer (gas-dominated)', nbv(D, 'late'), _sel(pyv(D, 'coupled'), (4,)), 1e-5),
        ('$\\lambda_{\\rm eff}$, $y_{\\rm eff}$, $G_r$', nbv(D, 'coup'), pyv(D, 'coup'), 1e-14),
        ('$\\sigma_{\\rm SI}$', nbv(D, 'sigsi'), pyv(D, 'sigsi'), 1e-12),
        ('$\\sigma/m$; $c_r$ allowed', nbv(D, 'sigm'), pyv(D, 'sigm'), 1e-12),
        ('$\\mu_P$; lifetimes', nbv(D, 'life'), pyv(D, 'life'), 1e-12),
        ('largest stable $m$ ($c_y=1$)', nbv(D, 'mmax'), pyv(D, 'mmax'), 1e-9),
        ('gas share $f$, overabundance $O$', nbv(D, 'over'), pyv(D, 'over'), 1e-8),
        ('radion decays: bulk fraction, brane share', nbv(D, 'bulk'), pyv(D, 'bulk'), 1e-9),
        ('$T_*$, $S(m_h/T_*)$', nbv(D, 'tstar'), pyv(D, 'tstar'), 1e-9),
        ('zero mode: $Y$ (Yukawa, Higgs portal)', nbv(D, 'yzero'), pyv(D, 'yzero'), 1e-12),
        ('$m_*$ ($N=1$, 12)', nbv(D, 'mstar'), pyv(D, 'mstar'), 1e-9),
        ('$c_y$ needed and allowed', nbv(D, 'cy'), pyv(D, 'cy'), 1e-9),
        ('Kaluza–Klein tower (Yukawa portal)', nbv(D, 'tower'), pyv(D, 'tower'), 1e-9),
        ('Kaluza–Klein lifetime; largest stable $m$ ($c_y=10^{-5}$)', nbv(D, 'kklife'), pyv(D, 'kklife'), 1e-9),
        ('Higgs-portal pairs', nbv(D, 'hp'), pyv(D, 'hp'), 1e-8),
        ('Higgs-portal pairs: fraction of the dark matter', nbv(D, 'hpfrac'), pyv(D, 'hpfrac'), 1e-8),
        ('vector portal: $T_V$, $(T_*/T_V)^{11}$', nbv(D, 'vec'), pyv(D, 'vec'), 1e-12),
        ('dependence on $T_{\\rm RH}$', nbv(D, 'trh'), pyv(D, 'trh'), 1e-8),
        ('Kohn–Sham matter: $w_{\\min}$, its $k_F$, $k_{F0}$', nbv(D, 'wmin'), pyv(D, 'wmin'), 1e-9),
        ('Kohn–Sham matter: $(k_F,p,w)$ table', nbv(D, 'wtab'), pyv(D, 'wtab'), 1e-9),
        ('$\\rho_\\Lambda$', nbv(D, 'rhol'), pyv(D, 'rhol'), 1e-12),
        ('Casimir energies (massless)', nbv(D, 'cas'), pyv(D, 'cas'), 1e-12),
        ('$m\\ell_\\perp$ and $m$ where the Casimir energy equals $\\rho_\\Lambda$', nbv(D, 'xl'), pyv(D, 'xl'), 1e-9),
        ('Casimir energy at $m\\ell_\\perp=1,10,30,60$', nbv(D, 'castab'), pyv(D, 'castab'), 1e-9)]
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


def V(D):
    """all printed numbers of section 24, by short name"""
    return {k: nbv(D, k) for k in NB}


# ================================================================================================================ document
def build(D, today):
    ck = D['checks'].get('dk', []); ref = D['ref'].get('dk', ('', '', {}))[2]
    nck = len(ck); nok = sum(1 for c in ck if c[1] == 'True')
    rb = [(k, v) for k, v in ref.items() if v in ('True', 'False')]
    cmp_ = compare(D)
    worst = max([r for *_, r in cmp_ if r is not None], default=None)
    v = V(D)
    req, md, pauli, dbeta, gas, late = v['req'], v['md'], v['pauli'], v['dbeta'], v['gas'], v['late']
    coup, sigsi, sigm, life, mmax = v['coup'], v['sigsi'], v['sigm'], v['life'], v['mmax']
    over, bulk, tstar, yzero, mstar, cy = v['over'], v['bulk'], v['tstar'], v['yzero'], v['mstar'], v['cy']
    tower, kklife, hp, hpfrac, vec, trh = v['tower'], v['kklife'], v['hp'], v['hpfrac'], v['vec'], v['trh']
    wmin, rhol, cas, xl, wtab, castab = v['wmin'], v['rhol'], v['cas'], v['xl'], v['wtab'], v['castab']
    pyc = pyv(D, 'coupled')
    gbr, gbulk = (pyc[2], pyc[3]) if len(pyc) > 3 else (None, None)

    d = Doc('dirac16complex as dark matter and dark energy',
            'An analysis based on the accurate approximate solutions of the coupled field equations of the 7+1 model in the three gravity '
            'regimes — section 24 of the notebook 7+1.nb',
            ['Repository: github.com/once-ere/7plus1 — prepared for Nash',
             'Notebook: 7+1.nb (Wolfram Language 15.0.1), evaluated ' + D['log'].get('finished', ''),
             'Document generated by tools/make_docs.py on ' + today],
            'Can the 16-component field dirac16complex of the 7+1 model be the dark matter of the universe, and can it be the dark energy? '
            'This document answers on the basis of the accurate approximate solutions of the coupled field equations in the three gravity '
            'regimes of the model (primordial, modern, far future), the couplings of the model, and present observational data. Every '
            'statement is marked DERIVED (an exact consequence of the field equations, checked), ESTIMATE (an order-of-magnitude statement '
            'with stated assumptions) or OPEN (not decided by the model). The verdict: dirac16complex has every property the solutions and '
            'the observations demand of dark matter — it is cold, dark and collisionless and enters the late-universe solution exactly as '
            'the cold-dark-matter density — provided the auxiliary space is stabilized and its Yukawa portal is very weak; but the minimal '
            'model does not produce the observed abundance with couplings of natural size (OPEN). It cannot be the dark energy: as matter '
            'its equation of state never falls below $w=-0.06$ in a stable state, and its vacuum energy is a constant contribution to the '
            'cosmological constant ($w=-1$ exactly). Section 24 of 7+1.nb contains %d checks (%d True); the independent Python program '
            'tools/ref_dark.py contains %d (%d True); the %d groups of numbers printed by both agree to a largest relative difference of %s.'
            % (nck, nok, len(rb), sum(1 for k, x in rb if x == 'True'), len(cmp_), fmt(worst, 2) if worst is not None else '?'))

    # ==============================================================================================
    d.h1('1. Results at a glance')
    d.table(['question', 'answer', 'status'], [
        ['Is dirac16complex matter cold, and can it be the $\\Omega_c$ of the late-universe solution?',
         'yes: the lightest mode has $p_a/\\varepsilon=k_a^2/(3(m^2a^2+k_a^2))\\to0$, $\\varepsilon\\propto a^{-3}$; it is exactly the term $\\Omega_ca^{-3}$ (section 4)', 'DERIVED'],
        ['Does it disturb the auxiliary space today?', 'no: it displaces the radion by $\\delta\\beta=H_0^2\\Omega_c/(2m_r^2)=%s$ (section 4.3)' % g(dbeta, 0, 2), 'DERIVED'],
        ['Does it need the stabilizing potential?', 'yes: without it a bulk gas re-expands the auxiliary space, $b\\propto\\tau^{2/5}$ (section 4.2)', 'DERIVED'],
        ['Is it dark?', '$\\sigma_{\\rm SI}=%s$ cm$^2$ at $m=40$ GeV for a Higgs portal of natural size (LZ: $2.2\\times10^{-48}$ cm$^2$) (section 5)' % g(sigsi, 0, 2), 'ESTIMATE'],
        ['Is it collisionless?', '$\\sigma/m=%s$ cm$^2$/g at $m=1$ GeV (Bullet Cluster: $<0.7$ cm$^2$/g) (section 5)' % g(sigm, 0, 2), 'ESTIMATE'],
        ['Is it stable?', 'absolutely without the Yukawa portal; with it, longer than the age of the universe only for $m<%s$ MeV$\\times c_y^{-2/3}$ (section 5)' % (fmt(mmax[0]*1e3, 3) if mmax else '?'), 'ESTIMATE'],
        ['Does the model make the observed abundance $mY=%s$ GeV?' % g(req, 3, 3),
         'its primordial gas and radion make far too much (overabundance $%s\\,T_{\\rm RH}$/GeV, brane share $%s$); after the dilution its moduli problem requires, its portals with natural couplings make far too little (the Higgs portal would need $c_P=%s$) (section 6)'
         % (g(over, 1, 2), g(bulk, 1, 2), g(hp, 5, 3)), 'OPEN'],
        ['Can its matter be the dark energy?', 'no: $w=p/\\varepsilon\\ge%s$, and $p<0$ only where uniform matter is unstable or metastable (section 7.1)' % g(wmin, 0, 3), 'DERIVED'],
        ['Can its vacuum energy be the dark energy?', 'only as a constant contribution to $\\Lambda_4$ ($w=-1$); its Casimir energy on the auxiliary torus is $%s\\,\\rho_\\Lambda$ if massless (section 7.3)' % g(cas, 0, 2), 'DERIVED / ESTIMATE'],
        ['Evolving dark energy (DESI DR2)?', 'the model predicts $w=-1$ exactly; a confirmed $w_0>-1$, $w_a<0$ would contradict it (section 7.4)', 'DERIVED'],
    ], widths=[0.27, 0.58, 0.15], small=True)

    # ==============================================================================================
    d.h1('2. The model, the conventions and the objects used')
    d.p('Units $\\hbar=c=1$; coordinates $X=(x_0,x_1,\\dots,x_7)$: $x_0$ is time (the direction of the octonion unit $e_0=1$), $x_1,x_2,x_3$ '
        'ordinary space, $x_4,\\dots,x_7$ the four auxiliary space-like coordinates; $y=(x_4,\\dots,x_7)$. The frame metric is '
        '$\\eta=\\mathrm{diag}(+1,-1,\\dots,-1)$, simply the squares of the octonion basis.')
    signature_block(d, eta=r'\eta', eta71=r'\eta_{71}', gam=r'\Gamma')
    d.p('**The field.** dirac16complex is a 16-component complex spinor field $\\psi(X)$ (anticommuting operators when quantized) with the real '
        '16×16 Dirac matrices $\\Gamma_A$, $\\Gamma_A\\Gamma_B+\\Gamma_B\\Gamma_A=2\\eta_{AB}1_{16}$, $\\Gamma^A=\\eta^{AA}\\Gamma_A$, the Dirac conjugate '
        '$\\bar\\psi=\\psi^\\dagger\\Gamma_0$, the frame matrices $\\Gamma^\\mu=e_A{}^\\mu\\Gamma^A$ ($g_{\\mu\\nu}=\\eta_{AB}e^A{}_\\mu e^B{}_\\nu$) and the spin '
        'connection $\\Omega_\\mu=\\frac14\\eta_{AC}\\omega_\\mu{}^C{}_B\\Gamma^A\\Gamma^B$. **The action of the model** (gravity with the Lovelock terms '
        '$k=0,\\dots,3$, $\\kappa=8\\pi G_8=M_8^{-6}$; dirac16complex with its mass and nine self-couplings; a stabilizing potential of the auxiliary '
        'space; on the world brane $W=\\lbrace y=0\\rbrace$ with induced metric $h$: a tension, the Standard Model and three portals):')
    d.align([r'S&=-\frac1{2\kappa}\int d^8X\sqrt{|g|}\big(2\Lambda+R+\alpha_2L_2+\alpha_3L_3\big)+\int d^8X\sqrt{|g|}\Big[\bar\psi\,i\Gamma^\mu(\partial_\mu+\Omega_\mu)\psi-m\bar\psi\psi',
             r'&\qquad-\frac12\sum_{r=0}^8g_r\sum_{|K|=r}\eta_K(\bar\psi\Theta^K\psi)^2-U(b)\Big]',
             r'&\quad+\int_Wd^4x\sqrt{|h|}\Big[-f_W+\mathcal L_{\rm SM}-\sum_{i,s}y_{is}\bar L_i\tilde\Phi P_R\psi_s+\text{h.c.}-G_P^f(\bar\psi\Gamma^\mu\psi)(\bar f\gamma_\mu f)-\lambda_P\Phi^\dagger\Phi\,\bar\psi\psi\Big],',
             r'L_2&=R^2-4R_{\mu\nu}R^{\mu\nu}+R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma},\qquad L_3=\frac18\delta^{\alpha_1\beta_1\alpha_2\beta_2\alpha_3\beta_3}_{\gamma_1\delta_1\gamma_2\delta_2\gamma_3\delta_3}'
             r'R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}R^{\gamma_2\delta_2}{}_{\alpha_2\beta_2}R^{\gamma_3\delta_3}{}_{\alpha_3\beta_3},',
             r'U(b)&=U_0+\frac{K_U}2\ln^2\frac b{b_0},\qquad b^8=\det(-g_{kl})\ (k,l=4,\dots,7),'])
    d.p('where $\\Theta^K=i^{\\sigma_r}\\Gamma^{A_1}\\cdots\\Gamma^{A_r}$ for $K=(A_1<\\dots<A_r)$, $\\sigma_r=\\frac{r(r-1)}2\\bmod2$, $\\eta_K=\\eta_{A_1A_1}\\cdots\\eta_{A_rA_r}$; '
        '$\\delta^{\\cdots}_{\\cdots}$ is the generalized Kronecker delta (the determinant of Kronecker deltas); $\\mathcal L_{\\rm SM}$ is the Standard '
        'Model Lagrangian in the field of $h$ (Higgs doublet $\\Phi$, $\\tilde\\Phi=i\\sigma_2\\Phi^*$, lepton doublets $L_i$, fermions $f$); $\\psi_s$ '
        'are brane-right-handed parts of $\\psi$ at the brane; $y_{is}$, $G_P^f$, $\\lambda_P$ are the Yukawa (neutrino), vector and Higgs portals.')
    d.p('**The homogeneous metric of the three regimes** and its Hubble rates (a prime is $d/dx_0$; $\\tau$ is the time since the singularity):')
    d.eq(r"ds^2=dx_0^2-a^2(dx_1^2+dx_2^2+dx_3^2)-b^2(dx_4^2+\cdots+dx_7^2),\qquad H_a=\frac{a'}a,\qquad H_b=\frac{b'}b,\qquad \beta=\ln\frac b{b_0}.")
    d.p('Every source is $T^\\mu{}_\\nu=\\mathrm{diag}(\\varepsilon,-p_a,-p_a,-p_a,-p_b,-p_b,-p_b,-p_b)$; for matter whose comoving energy (at fixed '
        'occupation numbers) is $W(a,b)$,')
    d.eq(r"\varepsilon=\frac W{a^3b^4},\qquad p_a=-\frac{\partial_aW}{3a^2b^4},\qquad p_b=-\frac{\partial_bW}{4a^3b^3},\qquad \varepsilon'+3H_a(\varepsilon+p_a)+4H_b(\varepsilon+p_b)=0,")
    d.p('and the field equations $\\sum_k\\alpha_kE^{(k)\\mu}{}_\\nu=\\kappa T^\\mu{}_\\nu$ ($\\alpha_0=2\\Lambda$, $\\alpha_1=1$) become three equations, '
        '$\\sum_k\\alpha_kE^{(k)0}{}_0=\\kappa\\varepsilon$, $\\sum_k\\alpha_kE^{(k)1}{}_1=-\\kappa p_a$, $\\sum_k\\alpha_kE^{(k)4}{}_4=-\\kappa p_b$; the first is the constraint')
    d.eq(r'3(H_a^2+4H_aH_b+2H_b^2)-12\alpha_2H_b(4H_a^3+18H_a^2H_b+12H_aH_b^2+H_b^3)+360\alpha_3H_a^2H_b^3(4H_a+3H_b)-\Lambda=\kappa\varepsilon .')
    d.p('**The radion.** In Einstein gravity the 11 and 44 equations and the constraint give exactly')
    d.eq(r"\beta''+3H_a\beta'+4\beta'^2=\kappa\Big(\frac\varepsilon6-\frac{p_a}2+\frac{p_b}3\Big)+\frac\Lambda3,")
    d.p('so the potential (which contributes $\\varepsilon=U$, $p_a=-U$, $p_b=-U-\\frac b4U\'$) gives the restoring term $-(\\kappa K_U/12)\\beta$: the radion '
        'mass is $m_r^2=\\kappa K_U/12$. With the auxiliary space static the 00 and 11 equations are the Friedmann equations of the brane,')
    d.eq(r"3H_a^2=8\pi G_4\rho_4+\Lambda_4,\qquad H_a'=-\frac\kappa2(\varepsilon+p_a),\qquad 8\pi G_4=\frac\kappa{V_\perp}=\frac1{M_{\rm Pl}^2},\qquad \Lambda_4=\Lambda+\kappa U_0+\frac{\kappa f_W}{V_\perp},")
    d.p('with $V_\\perp=b_0^4L_c^4$ the proper volume of the auxiliary torus ($L_c$ its coordinate period), $\\rho_4=V_\\perp\\varepsilon$ the energy '
        'per brane volume and $M_{\\rm Pl}$ the reduced Planck mass. The canonical radion is $\\phi=\\sqrt{12}M_{\\rm Pl}\\beta$; it couples to brane '
        'fields through $(\\phi/(\\sqrt3M_{\\rm Pl}))T^\\mu{}_\\mu$.')
    d.p('**Kaluza–Klein modes.** On the torus a field is a sum of modes $e^{i\\vec q\\cdot\\vec y}$ with $\\vec q=2\\pi\\vec n/L_c$, $\\vec n\\in\\mathbb Z^4$; '
        'at $b=b_0$ their four-dimensional masses are')
    d.eq(r'M_n=\sqrt{m^2+\mu_n^2},\qquad \mu_n=\frac{2\pi|\vec n|}{\ell_\perp},\qquad \ell_\perp=V_\perp^{1/4},\qquad '
         r'\#\{\vec n:|\vec n|^2=k\}=r_4(k)=8\sigma(k)-32\sigma(k/4),')
    d.p('($\\sigma$ the divisor sum, $\\sigma(k/4)=0$ unless $4\\mid k$: Jacobi’s four-square theorem). Each level is four brane Dirac spinors (the '
        '16 components of $\\psi$ split at the brane into four 4-component spinors). The number of levels with $\\mu_n<\\mu$ and their density are, '
        'in the continuum limit,')
    d.eq(r'N(\mu)=\frac{\pi^2}2\Big(\frac{\mu\ell_\perp}{2\pi}\Big)^4=\frac{V_\perp\mu^4}{32\pi^2},\qquad dN=\frac{V_\perp\mu^3\,d\mu}{8\pi^2}.')
    d.p('The LIGHTEST MODE ($\\vec n=0$) is constant over the torus: $\\psi(x,y)=\\chi(x)/\\sqrt{V_\\perp}$ with canonically normalized four-dimensional '
        'fields $\\chi$ (four Dirac spinors of mass $m$).')
    d.p('**Physical parameters** (as in the solutions of the modern regime; the mass $m$ of dirac16complex is free):')
    d.table(['quantity', 'value', 'meaning'], [
        ['$M_8$', '$10^4$ GeV', 'reduced eight-dimensional Planck mass, $\\kappa=M_8^{-6}$'],
        ['$M_{\\rm Pl}$', '$2.435\\times10^{18}$ GeV', 'reduced four-dimensional Planck mass, $M_{\\rm Pl}^2=V_\\perp/\\kappa=V_\\perp M_8^6$'],
        ['$V_\\perp$, $\\ell_\\perp$', '$M_{\\rm Pl}^2/M_8^6=%s$ GeV$^{-4}$; $%s$ GeV$^{-1}$' % (fmt(2.435e18**2/1e24, 4), fmt((2.435e18**2/1e24)**0.25, 6)), 'volume and size of the auxiliary torus'],
        ['$m_r$', '10 GeV', 'radion mass'],
        ['$\\alpha_2$, $\\alpha_3$, $f_W$', '$1/M_8^2$, $10^{-6}/M_8^4$, $M_8^4$', 'Gauss–Bonnet and cubic couplings, brane tension'],
        ['$m_h$, $v$, $G_F$', '125.20 GeV, $(\\sqrt2G_F)^{-1/2}=246.22$ GeV, $1.1663788\\times10^{-5}$ GeV$^{-2}$', 'Higgs mass, Higgs vacuum value, Fermi constant'],
    ], widths=[0.17, 0.43, 0.40], small=True)

    # ==============================================================================================
    d.h1('3. What the observations require (section 24.1)')
    d.p('Data used, as reported by the collaborations (the confidence levels are theirs):')
    d.table(['observation', 'result', 'source'], [
        ['cosmic microwave background', '$\\Omega_ch^2=0.120\\pm0.001$, $\\Omega_bh^2=0.0224\\pm0.0001$, $H_0=67.4\\pm0.5$ km/s/Mpc, $\\Omega_m=0.315\\pm0.007$', 'Planck 2018, Astron. Astrophys. 641 (2020) A6'],
        ['baryon acoustic oscillations', '$w(a)=w_0+w_a(1-a)$ with $w_0>-1$, $w_a<0$ preferred over $\\Lambda$ at $3.1\\sigma$ (DESI + CMB), $2.8\\sigma$–$4.2\\sigma$ with supernovae', 'DESI DR2, arXiv:2503.14738'],
        ['the same, with ACT', 'DESI + Planck PR4 + ACT DR6 + DES-Y5: $w_0=-0.753\\pm0.056$, $w_a=-0.84^{+0.23}_{-0.20}$ ($4.2\\sigma$)', 'Garcia-Quintero et al., arXiv:2504.18464'],
        ['direct detection', '$\\sigma_{\\rm SI}<2.2\\times10^{-48}$ cm$^2$ at 40 GeV (90% CL)', 'LZ, 4.2 tonne-years, arXiv:2410.17036'],
        ['self-interaction', '$\\sigma/m<1.25$ cm$^2$/g (68% CL), $<0.7$ cm$^2$/g with equal mass-to-light ratios', 'Randall et al., ApJ 679 (2008) 1173 (Bullet Cluster)'],
        ['small-scale structure', 'thermal-relic warm dark matter heavier than 5.7 keV (95% CL)', 'Iršič et al., Phys. Rev. D 109 (2024) 043511 (Lyman-$\\alpha$)'],
        ['the same, for decay-produced dark matter', 'warmer than a thermal relic of equal mass: about 7 keV where the thermal-relic bound is 3 keV', 'Ballesteros, Garcia, Pierre, JCAP 03 (2021) 101'],
        ['Pauli principle', 'fermion with $g=2$ internal states heavier than 0.13 keV (95% CL); the bound scales as $g^{-1/4}$', 'Alvey et al., MNRAS 501 (2021) 1188 (dwarf spheroidals)'],
        ['X-ray lines', 'limits on the mixing of sterile neutrinos of 6–40 keV, from $\\nu_s\\to\\nu\\gamma$', 'NuSTAR, Roach et al., Phys. Rev. D 107 (2023) 023009'],
        ['large extra dimensions', '$M_D>5.9$ to 11.2 TeV for $n=6$ to 2 ($M_{\\rm Pl}^2=R^nM_D^{n+2}$, $V=(2\\pi R)^n$)', 'Particle Data Group 2024 (ATLAS, 139 fb$^{-1}$)'],
    ], widths=[0.2, 0.47, 0.33], small=True)
    d.p('**The required abundance.** With the critical density per $h^2$, the entropy density today and its degrees of freedom,')
    d.align([r'mY&=\frac{\Omega_ch^2\,(\rho_c/h^2)}{s_0},\qquad \frac{\rho_c}{h^2}=3\Big(100\,\frac{\rm km}{\rm s\,Mpc}\Big)^2M_{\rm Pl}^2,\qquad '
             r's_0=\frac{2\pi^2}{45}g_{*s0}T_{\rm CMB}^3,\qquad g_{*s0}=2+\frac{21}4\cdot\frac4{11},'])
    d.p('($Y=n/s$ the comoving number density, $T_{\\rm CMB}=2.7255$ K) one finds $\\Omega_c=0.120/0.674^2=%s$, $\\rho_c/h^2=%s$ GeV/cm$^3$, '
        '$s_0=%s$ cm$^{-3}$ and' % (g(req, 0, 4), g(req, 1, 5), g(req, 2, 5)))
    d.eq(r'mY=%s\ \text{GeV}.' % g(req, 3, 4))
    d.p('**Collider bound.** In the convention of the Particle Data Group, $V_\\perp=(2\\pi R)^4$ and $M_{\\rm Pl}^2=R^4M_D^6$; with $M_{\\rm Pl}^2=V_\\perp M_8^6$ '
        'this is $M_D=(2\\pi)^{2/3}M_8=%s$ TeV, above the strongest bound (11.2 TeV). **Pauli bound.** The lightest mode has $g=16$ internal states '
        '(four brane Dirac spinors, each with two spin states, particle and antiparticle): $m>0.13\\,\\text{keV}\\times(2/16)^{1/4}=%s$ keV.'
        % (fmt(md[0]/1e3, 4) if md else '?', g(pauli, 0, 3)))

    # ==============================================================================================
    d.h1('4. dirac16complex as cold matter in the solutions of the three regimes (section 24.2, DERIVED)')
    d.h2('4.1 The stress tensor of a massive particle')
    d.p('A particle of dirac16complex with comoving momenta $k_a$ (ordinary space) and $k_b$ (auxiliary space) has the comoving energy '
        '$\\omega=W(a,b)$:')
    d.align([r'\omega&=\sqrt{m^2+\frac{k_a^2}{a^2}+\frac{k_b^2}{b^2}},\qquad \varepsilon\,a^3b^4=\omega,\qquad p_a\,a^3b^4=\frac{k_a^2}{3a^2\omega},\qquad p_b\,a^3b^4=\frac{k_b^2}{4b^2\omega},',
             r'\frac{p_a}\varepsilon\Big|_{k_b=0}&=\frac{k_a^2}{3(m^2a^2+k_a^2)}\;\xrightarrow{\,k_a/a\ll m\,}\;0,\qquad p_b\big|_{k_b=0}=0,\qquad '
             r'\varepsilon\big|_{k_a=k_b=0}=\frac m{a^3b^4}\propto a^{-3}\ \ (b\ \text{fixed}).'])
    d.p('The lightest mode ($k_b=0$) is therefore cold, pressureless matter once $k_a/a\\ll m$. Every other Kaluza–Klein level ($k_b/b_0=\\mu_n$) '
        'is matter of mass $M_n=(m^2+\\mu_n^2)^{1/2}$ with a static pressure along the auxiliary space, $p_b/\\varepsilon=\\mu_n^2/(4M_n^2)\\le\\frac14$.')
    d.h2('4.2 Regime I (primordial): the gas wins, and needs the stabilizer')
    d.p('Near $x_0=0$ any fluid with $p_a=w_a\\varepsilon$ ($w_a<1$) dominates the Lovelock terms. With $b$ finite and $H_a=p/\\tau$, the leading '
        'orders of the three equations are')
    d.align([r'\tau^3\textstyle\sum\alpha_kE^{(k)0}{}_0&\to48hp^3(30\alpha_3h^2-\alpha_2),\qquad \tau^3\textstyle\sum\alpha_kE^{(k)1}{}_1\to32hp^2(p-1)(18\alpha_3h^2-\alpha_2),',
             r'\tau^4\textstyle\sum\alpha_kE^{(k)4}{}_4&\to p^3(p-1)(648\alpha_3h^2-12\alpha_2),\qquad (H_b\to h),'])
    d.p('so $h=h_*=-\\sqrt{\\alpha_2/(54\\alpha_3)}$, $p=1/(1+w_a)$ ($a\\propto\\tau^{3/4}$ for the gas, $w_a\\to\\frac13$) and '
        '$\\kappa\\varepsilon\\tau^3\\to\\frac{64}3p^3\\alpha_2|h_*|$. Later, in the Einstein regime, the gas fraction grows and the gas dominates; Einstein gravity '
        'with $p_a=w_a\\varepsilon$, $p_b=w_b\\varepsilon$ has the power laws $a\\propto\\tau^P$, $b\\propto\\tau^Q$ with')
    d.align([r'3(P^2+4PQ+2Q^2)&=\mathcal C,\qquad 3P^2+8PQ-2P+10Q^2-4Q=-w_a\mathcal C,',
             r'3(2P^2+3PQ-P+2Q^2-Q)&=-w_b\mathcal C,\qquad 3P(1+w_a)+4Q(1+w_b)=2;',
             r'(w_a,w_b)&=(0,\tfrac14)\ \Rightarrow\ (P,Q)=(0,\tfrac25).'])
    d.p('The gas of the coupled primordial solution has $p_b/\\varepsilon\\to\\frac14$ (its auxiliary momenta are blue-shifted) and reaches '
        '$\\tau H_b=%s$ (exactly $\\frac25$ in the limit): ordinary space stops expanding and the AUXILIARY SPACE RE-EXPANDS. dirac16complex can '
        'therefore be dark matter only if the stabilizing potential $U(b)$ holds the auxiliary space.' % g(late, 0, 7))
    d.h2('4.3 Regime II (modern): the dark matter of the late-universe solution')
    d.p('With $U(b)$ the auxiliary space is caught when $H\\approx m_r$; the radion then oscillates, $\\beta\\approx A\\sin(m_r\\tau+\\vartheta)/(m_r(\\tau-\\tau_E))$ '
        'with $H_a+2H_b\\approx2/(3(\\tau-\\tau_E))$, and the dirac16complex gas — its auxiliary momenta frozen into Kaluza–Klein masses — is '
        'four-dimensional matter: its share of the energy is constant, %s and %s at $\\tau=10^5$ and $3\\times10^5$ (units $\\ell_2=\\sqrt{\\alpha_2}=1$; '
        'for the normalization of the primordial gas). In the late universe the brane obeys'
        % (g(gas, 0, 6), g(gas, 1, 6)))
    d.align([r'H_a&=H_0\sqrt{\Omega_ra^{-4}+\Omega_ma^{-3}+\Omega_\Lambda},\qquad \Omega_m=\Omega_b+\Omega_c,\qquad \Omega_\Lambda=1-\Omega_m-\Omega_r,',
             r't(a)&\approx\frac{2}{3H_0\Omega_m^2}\Big[(\Omega_ma-2\Omega_r)\sqrt{\Omega_r+\Omega_ma}+2\Omega_r^{3/2}\Big]+\frac{2\,\mathrm{arcsinh}\big(\sqrt{\Omega_\Lambda/\Omega_m}\,a^{3/2}\big)}{3H_0\sqrt{\Omega_\Lambda}}-\frac{2a^{3/2}}{3H_0\sqrt{\Omega_m}},'])
    d.p('the composite solution (accurate to $2.2\\times10^{-4}$; age 13.79 Gyr). Cold dirac16complex matter is exactly the term $\\Omega_ca^{-3}$: '
        'dirac16complex CAN BE the cold dark matter of this solution without changing it. Being uniform over the torus, its energy density per '
        'brane volume is $\\rho_{\\rm DM}=V_\\perp\\varepsilon$ and it gravitates exactly like brane matter ($8\\pi G_4=\\kappa/V_\\perp$). Through the radion '
        'equation (source $\\kappa\\varepsilon/6$ for $p_a=p_b=0$) it shifts the radion to')
    d.eq(r'\delta\beta=\frac{\kappa\varepsilon}{6m_r^2}=\frac{\rho_{\rm DM}}{6M_{\rm Pl}^2m_r^2}=\frac{H_0^2\Omega_c}{2m_r^2}=%s\qquad(\rho_{\rm DM}=3H_0^2M_{\rm Pl}^2\Omega_c,\ m_r=10\ \text{GeV}),'
         % g(dbeta, 0, 4))
    d.p('so $G_4\\propto e^{-4\\beta}$ does not change.')
    d.h2('4.4 Regime III (far future): the dark matter dilutes away')
    d.p('For matter and $\\Lambda$ the brane solution is exact,')
    d.align([r'a(t)&=\Big(\frac{\Omega_m}{\Omega_\Lambda}\Big)^{1/3}\sinh^{2/3}\Big(\frac32H_\infty t\Big),\qquad H_\infty=H_0\Omega_\Lambda^{1/2},',
             r'\frac{\rho_{\rm DM}}{\rho_\Lambda}&=\frac{\Omega_ca^{-3}}{\Omega_\Lambda}=\frac{\Omega_c}{\Omega_m}\,\frac1{\sinh^2(\frac32H_\infty t)}\;\to\;4\,\frac{\Omega_c}{\Omega_m}\,e^{-3H_\infty t},'])
    d.p('the identity $\\Omega_ma^{-3}/\\Omega_\\Lambda=\\mathrm{csch}^2(\\frac32H_\\infty t)$ and the limit $\\mathrm{csch}^2z\\,e^{2z}\\to4$ being checked exactly.')

    # ==============================================================================================
    d.h1('5. The lightest mode: dark, collisionless, and stable or not (section 24.3, ESTIMATE)')
    d.p('Every coupling of the lightest mode is diluted by the volume of the torus. With couplings of natural size at the eight-dimensional '
        'Planck scale — $\\lambda_P=c_P/M_8^5$, $y=c_y/M_8^2$, $G_P=c_V/M_8^6$, $g_r=c_r/M_8^6$ (mass dimensions $-5,-2,-6,-6$; $c_P,c_y,c_V,c_r$ '
        'dimensionless) — and $V_\\perp=M_{\\rm Pl}^2/M_8^6$:')
    d.eq(r'\lambda_{\rm eff}=\frac{\lambda_P}{V_\perp}=\frac{c_PM_8}{M_{\rm Pl}^2},\qquad y_{\rm eff}=\frac{y}{\sqrt{V_\perp}}=\frac{c_yM_8}{M_{\rm Pl}},\qquad '
         r'G_{\rm eff}=\frac{G_P}{V_\perp}=\frac{c_V}{M_{\rm Pl}^2},\qquad G_r=\frac{g_r}{V_\perp}=\frac{c_r}{M_{\rm Pl}^2},')
    d.p('numerically $\\lambda_{\\rm eff}=%s\\,c_P$ GeV$^{-1}$, $y_{\\rm eff}=%s\\,c_y$, $G_r=%s\\,c_r$ GeV$^{-2}$.' % (g(coup, 0, 4), g(coup, 1, 4), g(coup, 2, 4)))
    d.p('**Dark.** The operator $\\lambda_{\\rm eff}\\Phi^\\dagger\\Phi\\bar\\chi\\chi$ contains $\\lambda_{\\rm eff}v\\,h\\bar\\chi\\chi$ ($\\Phi^\\dagger\\Phi=(v+h)^2/2$); the '
        'nucleon couples to the Higgs boson with $f_Nm_N/v$ ($f_N\\approx0.30$, $m_N=0.9383$ GeV). Exchanging $h$ gives the contact interaction '
        '$G_{\\chi N}\\bar\\chi\\chi\\bar NN$ with $G_{\\chi N}=\\lambda_{\\rm eff}f_Nm_N/m_h^2$, and for a Dirac fermion')
    d.eq(r'\sigma_{\rm SI}=\frac{G_{\chi N}^2\mu^2}{\pi}=\frac{\lambda_{\rm eff}^2f_N^2m_N^2\mu^2}{\pi m_h^4},\qquad \mu=\frac{mm_N}{m+m_N};')
    d.p('at $m=40$ GeV and $c_P=1$: $\\sigma_{\\rm SI}=%s$ cm$^2$, some 56 orders of magnitude below the LZ limit.' % g(sigsi, 0, 3))
    d.p('**Collisionless.** A contact self-interaction of strength $G_r$ gives $\\sigma\\approx G_r^2m^2/(4\\pi)$: $\\sigma/m=%s$ cm$^2$/g at $m=1$ GeV '
        '($c_r=1$); the Bullet Cluster allows $c_r$ up to $%s$ at $m=1$ GeV.' % (g(sigm, 0, 3), g(sigm, 1, 3)))
    d.p('**Stable or not.** The Higgs and vector portals are bilinear in $\\chi$ and cannot make it decay; without the Yukawa portal the lightest '
        'mode is absolutely stable. Through the Yukawa portal it mixes with the neutrinos: after symmetry breaking the portal is a Dirac mass '
        '$\\mu_P=y_{\\rm eff}v/\\sqrt2$ between $\\nu_L$ and $\\chi$ ($\\mu_P=%s\\,c_y$ GeV), the mixing angle is $\\theta=\\mu_P/m$, and $\\chi$ decays like a '
        'sterile neutrino (rates of a Majorana fermion; a Dirac fermion decays half as fast, immaterial here):' % g(life, 0, 4))
    d.eq(r'\Gamma(\chi\to\nu\gamma)=\frac{9\,\alpha\,G_F^2\,m^5\sin^22\theta}{1024\,\pi^4},\qquad \Gamma(\chi\to3\nu)=\frac{G_F^2m^5\theta^2}{96\pi^3},\qquad '
         r'\sin^22\theta\approx4\theta^2,\qquad \frac{\Gamma(\chi\to\nu\gamma)}{\Gamma(\chi\to3\nu)}=\frac{27\alpha}{8\pi}.')
    d.p('Since $\\theta^2m^5=\\mu_P^2m^3$, $\\tau(\\chi\\to3\\nu)=K/(c_y^2m^3)$ with $K=%s$ s GeV$^3$. Lifetimes at $c_y=1$:' % g(life, 2, 4))
    d.table(['$m$', '$\\tau(\\chi\\to\\nu\\gamma)$ (s)', '$\\tau(\\chi\\to3\\nu)$ (s)'], [
        ['1 GeV', g(life, 1, 3), g(life, 2, 3)], ['10 keV', g(life, 3, 3), g(life, 4, 3)]], widths=[0.2, 0.35, 0.35])
    d.p('The lifetime exceeds the age of the universe ($t_0=13.79$ Gyr) for $m<m_{\\max}=(K/(c_y^2t_0))^{1/3}=%s$ GeV$\\times c_y^{-2/3}$. This is the '
        'minimum requirement; searches for the X-ray line of $\\chi\\to\\nu\\gamma$ demand far longer lifetimes and only strengthen the conclusions '
        'of section 6.' % g(mmax, 0, 4))

    # ==============================================================================================
    d.h1('6. The abundance (section 24.4, ESTIMATE; the outcome is OPEN)')
    d.h2('6.1 From the solutions of regimes I and II: overproduction')
    d.p('**The primordial gas.** After the capture of the auxiliary space the gas carries the fraction $f=%s$ of the energy and the radion the rest '
        '(for the normalization of the primordial gas; $f$ scales with it). Suppose — most favourably — that the radion energy all became brane '
        'radiation at the temperature $T_{\\rm RH}$. For radiation $\\rho_{\\rm rad}=\\frac{\\pi^2}{30}g_*T^4$ and $s=\\frac{2\\pi^2}{45}g_{*s}T^3$, so '
        '$\\rho_{\\rm rad}/s=\\frac34T$ ($g_*=g_{*s}$), and the dark matter needs $\\rho_{\\rm DM}/\\rho_{\\rm rad}=mYs/(\\frac34Ts)=\\frac43mY/T_{\\rm RH}$. The '
        'gas exceeds it by' % g(over, 0, 4))
    d.eq(r'O=\frac{f}{1-f}\,\frac{T_{\rm RH}}{\frac43mY}=%s\,\frac{T_{\rm RH}}{\rm GeV},\qquad O(4\ \text{MeV})=%s,' % (g(over, 1, 4), g(over, 2, 3)))
    d.p('4 MeV being the lowest reheating temperature compatible with nucleosynthesis.')
    d.p('**Radion decay.** The radion decays into brane fields through $(\\phi/(\\sqrt3M_{\\rm Pl}))T^\\mu{}_\\mu$,')
    d.eq(r'\Gamma(\phi\to f\bar f)=\frac{N_cm_f^2m_r}{24\pi M_{\rm Pl}^2}\Big(1-\frac{4m_f^2}{m_r^2}\Big)^{3/2},\qquad '
         r'\Gamma(\phi\to gg)=\frac{b_{\rm QCD}^2\alpha_s^2m_r^3}{96\pi^3M_{\rm Pl}^2}\quad(b_{\rm QCD}=11-\tfrac23n_f),')
    d.p('and into Kaluza–Klein pairs of dirac16complex. In the Einstein frame a level has the mass $M_E(\\beta)$ below, with $M_n=M_E(0)$, and '
        'hence the coupling $g_n$ to $\\phi=\\sqrt{12}M_{\\rm Pl}\\beta$:')
    d.eq(r'M_E(\beta)=e^{-2\beta}\sqrt{m^2+\mu_n^2e^{-2\beta}},\qquad g_n=\Big|\frac{dM_E}{d\phi}\Big|_{\beta=0}=\frac{2M_n+\mu_n^2/M_n}{\sqrt{12}\,M_{\rm Pl}};')
    d.p('with four brane spinors per level and $m\\ll\\mu_n$ the widths are')
    d.align([r'\Gamma_n&=\frac{3\mu_n^2m_r}{8\pi M_{\rm Pl}^2}\Big(1-\frac{4\mu_n^2}{m_r^2}\Big)^{3/2},\qquad '
             r'\Gamma_{\rm bulk}=\sum_n\Gamma_n\approx\frac{3m_r}{8\pi M_{\rm Pl}^2}\,\frac{V_\perp}{8\pi^2}\Big(\frac{m_r}2\Big)^6\frac{B(3,\frac52)}2,'])
    d.p('$B$ the Euler beta function. With $m_r=10$ GeV: $\\Gamma_{\\rm brane}=%s$ GeV, $\\Gamma_{\\rm bulk}=%s$ GeV — the radion does NOT become brane '
        'radiation: a fraction %s of its decays produces dirac16complex, and the brane receives only $\\Gamma_{\\rm brane}/\\Gamma_{\\rm bulk}=%s$ of its '
        'energy. Without a later dilution the universe after the radion decay is made of dirac16complex — the opposite of the observed ratio.'
        % (fmt(gbr, 3) if gbr else '?', fmt(gbulk, 3) if gbulk else '?', g(bulk, 0, 16), g(bulk, 1, 3)))
    d.p('Both results mean that dirac16complex as dark matter needs the dilution and reheating that the moduli problem of the minimal model '
        'already needs (a later period of inflation on the brane, not part of the model), at a reheating temperature below the normalcy '
        'temperature, above which the brane plasma loses its energy to bulk gravitons ($\\dot\\rho/\\rho\\sim T^7/M_8^6>H\\sim T^2/M_{\\rm Pl}$):')
    d.eq(r'T_{\rm RH}<T_*=\Big(\frac{M_8^6}{M_{\rm Pl}}\Big)^{1/5}=%s\ \text{GeV}.' % g(tstar, 0, 4))
    d.h2('6.2 After the dilution: the portals')
    d.p('**Freeze-in.** A particle $\\chi$ made in decays $B\\to\\chi+\\dots$ of a particle $B$ (internal states $g_B$, mass $m_B$, partial width $\\Gamma$) in '
        'equilibrium obeys $\\dot n_\\chi+3Hn_\\chi=g_B\\Gamma m_B^2TK_1(m_B/T)/(2\\pi^2)$ (as long as $n_\\chi$ is far below equilibrium). With '
        '$H=1.66\\sqrt{g_*}T^2/M_{\\rm Pl}^{\\rm (unreduced)}$, $s=\\frac{2\\pi^2}{45}g_{*s}T^3$ and $x=m_B/T$, $Y=n_\\chi/s$ integrates to (Hall, Jedamzik, '
        'March-Russell and West)')
    d.align([r'Y&=\frac{135\,g_B}{8\pi^3\cdot1.66\,g_*^{3/2}}\,\frac{M_{\rm Pl}^{\rm (unreduced)}\,\Gamma}{m_B^2}\,S(x_{\rm RH}),\qquad '
             r'S(x)=\frac2{3\pi}\int_x^\infty t^3K_1(t)\,dt,\qquad \int_0^\infty t^3K_1(t)\,dt=\frac{3\pi}2,'])
    d.p('where $x_{\\rm RH}=m_B/T_{\\rm RH}$ (instantaneous reheating: production only below $T_{\\rm RH}$), $M_{\\rm Pl}^{\\rm (unreduced)}=1.22091\\times10^{19}$ '
        'GeV, $g_*=g_{*s}=106.75$ (below the electroweak bosons $g_*=86.25$, so for $T_{\\rm RH}<m_W$ this understates $Y$ by at most '
        '$(106.75/86.25)^{3/2}=1.4$). For Higgs decays ($B=h$, $g_B=1$) at $T_{\\rm RH}=T_*$: $S(m_h/T_*)=%s$. Other production processes add to '
        'these yields; the estimates below are lower bounds on what each portal makes — except for the lightest mode alone, where $S=1$ is used '
        'as an upper bound (up to factors of order one).' % g(tstar, 1, 4))
    d.p('**Yukawa portal, the lightest mode.** $\\Gamma(h\\to\\bar\\nu\\chi)+\\Gamma(h\\to\\nu\\bar\\chi)=y_{\\rm eff}^2m_h/(16\\pi)$ per channel; with $N$ channels '
        '(three lepton flavours and the brane spinors, $N\\le12$) $Y=%s\\,c_y^2NS$. With $Y_1=Y(c_y=1)$ the observed abundance needs '
        '$c_y^2=mY/(m\\,Y_1)$, which falls like $1/m$, while the lifetime allows only $c_y^2\\le K/(m^3t_0)$, which falls like $1/m^3$. Both hold '
        'only for' % g(yzero, 0, 4))
    d.eq(r'm<m_*=\Big(\frac{Y_1\,K}{mY\,t_0}\Big)^{1/2}=%s\ \text{keV}\ (N=1),\quad %s\ \text{keV}\ (N=12),'
         % (fmt(mstar[0]*1e6, 3) if mstar else '?', fmt(mstar[1]*1e6, 3) if len(mstar) > 1 else '?'))
    d.p('even with $S=1$ — below the Lyman-$\\alpha$ bound. For example ($N=1$): at $m=10$ keV the abundance needs $c_y=%s$ while stability allows '
        '$c_y<%s$; at $m=10$ MeV it needs $c_y=%s$ and allows $c_y<%s$.' % (g(cy, 0, 3), g(cy, 1, 3), g(cy, 2, 3), g(cy, 3, 3)))
    d.p('**Yukawa portal, the Kaluza–Klein tower.** Every level couples to the brane with the same strength as the lightest mode (its mode function '
        'at $y=0$ is also $1/\\sqrt{V_\\perp}$), so Higgs decays make the whole tower below $m_h$, $\\Gamma(h\\to\\bar\\nu\\chi_n)\\propto(1-\\mu_n^2/m_h^2)^2$. With '
        '$dN=V_\\perp\\mu^3d\\mu/(8\\pi^2)$ and $\\int_0^12\\pi^2u^3(1-u^2)^2du=\\pi^2/12$, $\\int_0^12\\pi^2u^4(1-u^2)^2du=16\\pi^2/315$:')
    d.align([r'\sum_n\Big(1-\frac{\mu_n^2}{m_h^2}\Big)^2&=\frac{V_\perp m_h^4}{192\pi^2}=%s,\qquad \sum_n\mu_n\Big(1-\frac{\mu_n^2}{m_h^2}\Big)^2=\frac{V_\perp m_h^5}{315\pi^2}=%s\ \text{GeV},'
             % (g(tower, 0, 4), g(tower, 1, 4)),
             r'Y_{\rm tower}&=%s\,c_y^2NS,\qquad \frac{\rho_{\rm tower}}s=%s\,c_y^2NS\ \text{GeV}\qquad(\text{mean mass }\tfrac{192}{315}m_h=%s\ \text{GeV}).'
             % (g(tower, 2, 4), g(tower, 3, 4), fmt(192/315*125.2, 4))])
    d.p('(The continuum sums are the limits of the lattice sums with $r_4(k)$: for a radius $R=100$, $\\sum_kr_4(k)(1-k/R^2)^2$ and '
        '$\\sum_kr_4(k)\\sqrt k/R\\,(1-k/R^2)^2$ equal $\\pi^2R^4/12$ and $16\\pi^2R^4/315$ to $10^{-8}$.) At $T_{\\rm RH}=T_*$ the tower carries more energy '
        'than all the dark matter ($\\rho_{\\rm tower}/s>mY$) unless')
    d.eq(r'c_y<%s\qquad(N=1,\ T_{\rm RH}=T_*),' % g(tower, 4, 3))
    d.p('and through the same mixing its states decay into neutrinos and photons long after nucleosynthesis: $\\tau(\\chi_n\\to3\\nu)=%s$ s$/c_y^2$ at '
        '$\\mu_n=10$ GeV (other channels shorten it by about ten). The Yukawa portal must therefore be very weak or absent; the precise bound from '
        'nucleosynthesis and the CMB is not derived here (OPEN). With $c_y<10^{-5}$ the lightest mode outlives the universe for $m<%s$ GeV, and its '
        'freeze-in yield is negligible.' % (g(kklife, 0, 3), g(kklife, 1, 3)))
    d.p('**Higgs portal.** The lightest mode alone gets $Y=%s\\,c_P^2S$ from $\\Gamma(h\\to\\chi\\bar\\chi)=(\\lambda_{\\rm eff}v)^2m_h/(8\\pi)$. But the brane '
        'breaks the conservation of Kaluza–Klein momentum: $\\bar\\psi\\psi|_W=\\sum_{n,n\'}\\bar\\chi_n\\chi_{n\'}/V_\\perp$, so $\\lambda_{\\rm eff}v\\,h\\bar\\chi_n\\chi_{n\'}$ '
        'couples EVERY pair of levels. For a scalar coupling to two fermions of masses $\\mu_1=u_1m_h$, $\\mu_2=u_2m_h$,' % g(yzero, 1, 3))
    d.align([r'\Gamma(h\to\chi_n\bar\chi_{n^\prime})&=\frac{(\lambda_{\rm eff}v)^2m_h}{8\pi}\big(1-(u_1+u_2)^2\big)\lambda^{1/2}(1,u_1^2,u_2^2)',
             r'&=\frac{(\lambda_{\rm eff}v)^2m_h}{8\pi}\big(1-(u_1+u_2)^2\big)^{3/2}\big(1-(u_1-u_2)^2\big)^{1/2},',
             r'\lambda(1,x_1,x_2)&=(1-x_1-x_2)^2-4x_1x_2=\big(1-(\sqrt{x_1}+\sqrt{x_2})^2\big)\big(1-(\sqrt{x_1}-\sqrt{x_2})^2\big).'])
    d.p('Summed over all pairs with $dN\\,dN\'$ the sum is $(V_\\perp m_h^4/(8\\pi^2))^2=%s$ times' % g(hp, 2, 3))
    d.align([r'I_{hp}&=\iint_{u_1+u_2<1}u_1^3u_2^3\big(1-(u_1+u_2)^2\big)^{3/2}\big(1-(u_1-u_2)^2\big)^{1/2}du_1\,du_2=%s,' % g(hp, 0, 5),
             r'I_{hp}^E&=\iint_{u_1+u_2<1}(u_1+u_2)\,u_1^3u_2^3\big(1-(u_1+u_2)^2\big)^{3/2}\big(1-(u_1-u_2)^2\big)^{1/2}du_1\,du_2=%s,' % g(hp, 1, 5)])
    d.p('the second with the mass weight $u_1+u_2$ (mean pair mass $I^E_{hp}/I_{hp}=%s\\,m_h$). With four brane spinors per level the pairs carry '
        '$\\rho/s=%s\\,c_P^2S$ GeV: the observed abundance needs' % (fmt(hp[1]/hp[0], 3) if len(hp) > 1 else '?', g(hp, 3, 3)))
    d.eq(r'c_P=%s\ \ (T_{\rm RH}\gg m_h,\ S\to1),\qquad c_P=%s\ \ (T_{\rm RH}=T_*),' % (g(hp, 4, 3), g(hp, 5, 3)))
    d.p('two orders of magnitude above natural size; for $c_P=1$ the pairs make $%s$ of the dark matter. Such dark matter would be a population of '
        'heavy Kaluza–Klein modes (tens of GeV) whose decays into lighter levels and Standard Model particles through the same portal are not '
        'computed here (OPEN).' % g(hpfrac, 0, 2))
    d.p('**Vector portal** (and any coupling of gravitational strength per pair, $G_{\\rm eff}=c_V/M_{\\rm Pl}^2$). The rate per pair is of order '
        '$G_{\\rm eff}^2T^5$ and the number of pairs lighter than $T$ of order $N(T)^2\\propto V_\\perp^2T^8$, so the production rate divided by '
        '$H\\sim T^2/M_{\\rm Pl}$ is of order (numerical factors omitted; they raise $T_V$)')
    d.eq(r'\frac{\Gamma}{H}\sim c_V^2\Big(\frac{T}{T_V}\Big)^{11},\qquad T_V=\Big(\frac{M_8^{12}}{M_{\rm Pl}}\Big)^{1/11}=%s\ \text{GeV},\qquad \Big(\frac{T_*}{T_V}\Big)^{11}=%s:'
         % (g(vec, 0, 4), g(vec, 1, 2)))
    d.p('negligible. **Dependence on the reheating temperature** (the Boltzmann factor $S$ falls steeply below $T_*$, and with it what the portals make):')
    trows = []
    for i in range(0, len(trh) - 3, 4):
        trows.append([fmt(trh[i], 4), fmt(trh[i + 1], 3), fmt(trh[i + 2], 3), fmt(trh[i + 3], 3)])
    d.table(['$T_{\\rm RH}$ (GeV)', '$S(m_h/T_{\\rm RH})$', 'tower: largest $c_y$ ($N=1$)', 'Higgs-portal pairs: $c_P$ needed'], trows,
            widths=[0.2, 0.25, 0.27, 0.28], small=True)
    d.h2('6.3 Conclusion on the abundance')
    d.p('With couplings of natural size the minimal model does not produce the observed abundance after the dilution its moduli problem requires: '
        'the Yukawa portal must be tiny (else the tower is overproduced), the vector portal is negligible, and the Higgs portal falls short by a '
        'factor $%s$ (it would need $c_P\\approx%s$). Before that dilution the model makes far too much. The origin of the dark-matter abundance is '
        'therefore OPEN: it must come from physics outside the minimal model (for example the reheating after the brane inflation that the moduli '
        'problem needs, or couplings much larger than natural).' % (fmt(1/hpfrac[0], 2) if hpfrac else '?', fmt(hp[5], 3) if len(hp) > 5 else '?'))

    # ==============================================================================================
    d.h1('7. dirac16complex and dark energy (section 24.5)')
    d.h2('7.1 Matter cannot be dark energy (DERIVED)')
    d.p('Free particles have $p_a=-\\partial_aW/(3a^2b^4)\\ge0$, because $W=\\omega$ decreases with $a$. Interacting matter at zero temperature has '
        '$p=-\\partial E/\\partial V$ at fixed particle number $N$; with $E=\\varepsilon V$, $\\rho=N/V$:')
    d.eq(r'p=\rho\frac{d\varepsilon}{d\rho}-\varepsilon=\rho^2\frac{d}{d\rho}\Big(\frac\varepsilon\rho\Big)=\mu\rho-\varepsilon,\qquad \mu=\frac{d\varepsilon}{d\rho}=\varepsilon_F\ \ \text{(Hugenholtz–Van Hove)}.')
    d.p('The pressure is negative only where the energy per particle $\\varepsilon/\\rho$ decreases with density — below the saturation density '
        '$\\rho_0$, where $\\varepsilon/\\rho$ has its minimum. There uniform matter has a higher energy than a mixture of droplets at saturation and '
        'empty space (both of pressure zero: energy $\\rho\\,\\varepsilon_0/\\rho_0<\\varepsilon$), and is metastable or, in the spinodal region '
        '$dp/d\\rho<0$, unstable. The uniform Kohn–Sham matter of dirac16complex (section 8 of 7+1.nb: seven space dimensions, Hartree approximation, '
        'units $m=1$, scalar and vector couplings $g_S=-10^5$, $g_V=10^4$) is')
    d.align([r'\rho&=\frac{k_F^7}{105\pi^4},\qquad c_7=\frac1{15\pi^4},\qquad \rho_S=c_7\int_0^{k_F}\frac{k^6m^*}{\sqrt{k^2+m^{*2}}}dk,\qquad m^*=m+g_S\rho_S,',
             r'\varepsilon&=c_7\int_0^{k_F}k^6\sqrt{k^2+m^{*2}}\,dk-\frac{(m^*-m)^2}{2g_S}+\frac{g_V}2\rho^2,\qquad \varepsilon_F=g_V\rho+\sqrt{k_F^2+m^{*2}},'])
    d.p('saturating at $k_{F0}=%s$. Its equation of state $w=p/\\varepsilon$ (with $p=\\varepsilon_F\\rho-\\varepsilon$):' % g(wmin, 2, 8))
    wrows = []
    for i in range(0, len(wtab) - 2, 3):
        wrows.append([fmt(wtab[i], 3), fmt(wtab[i + 1], 4), fmt(wtab[i + 2], 4)])
    d.table(['$k_F$', '$p$ (units $m^8$)', '$w=p/\\varepsilon$'], wrows, widths=[0.2, 0.35, 0.35], small=True)
    d.p('$p>0$ at low density, $p<0$ between about $k_F=0.5$ and $k_{F0}$ (with a spinodal region), $p=0$ at saturation; the most negative value is '
        '$w_{\\min}=%s$ at $k_F=%s$ — nowhere near $-1$, and only in unstable or metastable states.' % (g(wmin, 0, 4), g(wmin, 1, 3)))
    d.h2('7.2 A condensate is only a cosmological constant (DERIVED)')
    d.p('A Lorentz-invariant condensate of dirac16complex ($\\langle\\bar\\psi\\psi\\rangle\\ne0$, which the self-couplings allow) has '
        '$T_{\\mu\\nu}=-\\rho_{\\rm vac}g_{\\mu\\nu}$, i.e. $p=-\\varepsilon$, $w=-1$ exactly: a contribution to $\\Lambda_4$, not a dynamical dark energy.')
    d.h2('7.3 Vacuum energy on the auxiliary torus (DERIVED, with an ESTIMATE of its size)')
    d.p('The zero-point energy of dirac16complex on the torus (periodic boundary conditions) is a Casimir energy. For a real bosonic field of mass '
        '$m$ in eight dimensions, four of them compact with proper size $\\ell_\\perp$, the finite part of the vacuum energy per unit eight-volume is '
        '$-\\frac2{(2\\pi)^4}\\sum_{\\vec n\\ne0}\\big(\\frac m{\\ell_\\perp|\\vec n|}\\big)^4K_4(m\\ell_\\perp|\\vec n|)$; per unit brane volume ($\\times V_\\perp$):')
    d.align([r'\rho_{\rm Cas}^{(1)}&=-\frac2{\ell_\perp^4}\sum_{\vec n\in\mathbb Z^4\setminus0}\Big(\frac{m\ell_\perp}{2\pi|\vec n|}\Big)^4K_4(m\ell_\perp|\vec n|)'
             r'\;\xrightarrow{\,m\to0\,}\;-\frac{\Gamma(4)\,Z}{\pi^4\ell_\perp^4},\qquad K_\nu(z)\to\frac{\Gamma(\nu)}2\Big(\frac2z\Big)^\nu,',
             r'Z&=\sum_{\vec n\ne0}\frac1{|\vec n|^8}=\sum_{k\ge1}\frac{r_4(k)}{k^4}=8(1-4^{-3})\zeta(4)\zeta(3),'])
    d.p('the last equality from $\\sum_kr_4(k)k^{-s}=8(1-4^{1-s})\\zeta(s)\\zeta(s-1)$ (checked against the shell sum to $k=20000$). A fermion contributes '
        'with the opposite sign; dirac16complex has 16 real degrees of freedom on shell, the eight-dimensional graviton $8\\cdot5/2=20$. With '
        '$\\rho_\\Lambda=3H_0^2M_{\\rm Pl}^2\\Omega_\\Lambda=%s$ GeV$^4$:' % g(rhol, 0, 4))
    d.eq(r'\rho_{\rm Cas}(\text{dirac16complex},m=0)=+16\,\frac{\Gamma(4)Z}{\pi^4\ell_\perp^4}=%s\,\rho_\Lambda,\qquad \rho_{\rm Cas}(\text{graviton})=%s\,\rho_\Lambda .'
         % (g(cas, 0, 4), g(cas, 1, 4)))
    ctab = [['1', g(castab, 0, 4)], ['10', g(castab, 1, 4)], ['30', g(castab, 2, 4)], ['60', g(castab, 3, 4)],
            [g(xl, 0, 6), '1 (by construction)']]
    d.table(['$m\\ell_\\perp$', '$\\rho_{\\rm Cas}(\\text{dirac16complex})/\\rho_\\Lambda$'], ctab, widths=[0.3, 0.5], small=True)
    d.p('It equals $\\rho_\\Lambda$ only at $m\\ell_\\perp=%s$, i.e. $m=%s$ MeV. Every such term is constant while $b$ is frozen ($w=-1$ up to the radion '
        'displacement, $10^{-86}$), and all of them add to $\\Lambda+\\kappa U_0+\\kappa f_W/V_\\perp$: the observed $\\Lambda_4$ needs their sum to cancel to '
        '$10^{-35}$ of each term — the cosmological-constant problem, which dirac16complex does not solve.' % (g(xl, 0, 5), fmt(xl[1]*1e3, 4) if len(xl) > 1 else '?'))
    d.h2('7.4 The prediction $w=-1$, and DESI')
    d.p('With the radion frozen the model predicts $w=-1$ exactly. The DESI preference for $w_0>-1$, $w_a<0$ ($2.8\\sigma$ to $4.2\\sigma$ depending on '
        'the supernova sample; $w_0=-0.753\\pm0.056$, $w_a=-0.84^{+0.23}_{-0.20}$ for DESI + Planck PR4 + ACT DR6 + DES-Y5), if confirmed, would contradict '
        'the model whatever dirac16complex is. The radion, the only scalar of the model, cannot be the dark energy either: with $m_r\\gg H_0$ it is '
        'frozen, and a radion light enough to roll today would change $G_4\\propto e^{-4\\beta}$ by order one (ESTIMATE).')

    # ==============================================================================================
    d.h1('8. Verdict (section 24.6)')
    d.table(['claim', 'status', 'basis'], [
        ['dirac16complex behaves as cold matter ($w=0$) and can be the $\\Omega_c$ of the late-universe solution', 'DERIVED', 'sections 4.1, 4.3'],
        ['its effect on the radion is negligible ($\\delta\\beta=%s$)' % g(dbeta, 0, 2), 'DERIVED', 'section 4.3'],
        ['as dark matter it needs the stabilizing potential (otherwise decompactification)', 'DERIVED', 'section 4.2'],
        ['dark: $\\sigma_{\\rm SI}=%s$ cm$^2$; collisionless: $\\sigma/m=%s$ cm$^2$/g' % (g(sigsi, 0, 2), g(sigm, 0, 2)), 'ESTIMATE', 'section 5'],
        ['stable over the age of the universe if $m<%s$ MeV$\\times c_y^{-2/3}$ (absolutely stable without the Yukawa portal)' % (fmt(mmax[0]*1e3, 3) if mmax else '?'), 'ESTIMATE', 'section 5'],
        ['mass window: $m>5.7$ keV (Lyman-$\\alpha$), $>%s$ keV (Pauli, $g=16$)' % g(pauli, 0, 2), 'data', 'section 3'],
        ['primordial gas and radion decay overproduce dirac16complex', 'DERIVED from the solutions; ESTIMATE for $O$', 'section 6.1'],
        ['Yukawa portal: the Kaluza–Klein tower is overproduced unless $c_y<%s$' % g(tower, 4, 2), 'ESTIMATE', 'section 6.2'],
        ['Higgs portal: the observed abundance needs $c_P=%s$ ($T_{\\rm RH}=T_*$)' % g(hp, 5, 3), 'ESTIMATE', 'section 6.2'],
        ['observed abundance produced by the minimal model with natural couplings', 'OPEN (not produced)', 'section 6.3'],
        ['dirac16complex matter as dark energy', 'excluded (DERIVED: $w\\ge-0.06$; $p<0$ only where uniform matter is unstable or metastable)', 'section 7.1'],
        ['dirac16complex vacuum energy or condensate as dark energy', 'not an explanation (DERIVED: a contribution to $\\Lambda_4$ with $w=-1$)', 'sections 7.2, 7.3'],
        ['evolving dark energy (DESI DR2: $2.8\\sigma$ to $4.2\\sigma$)', 'in tension with the model ($w=-1$ exactly)', 'section 7.4'],
    ], widths=[0.55, 0.28, 0.17], small=True)

    # ==============================================================================================
    d.h1('9. Validity and limitations')
    d.ul(['The DERIVED statements follow from the field equations of the model and the solutions of the three regimes; they hold for any mass $m$ '
          'and any couplings.',
          'The ESTIMATES assume couplings of natural size at $M_8$ ($c_P,c_y,c_V,c_r$ of order one), $M_8=10$ TeV, $m_r=10$ GeV, instantaneous '
          'reheating at $T_{\\rm RH}$, production of the portals by Higgs decays only (other processes add), $g_*=106.75$ in the freeze-in yield, the '
          'continuum limit of the Kaluza–Klein sums (exact to $10^{-8}$ for $R=100$), and sterile-neutrino decay rates with the Majorana normalization; '
          'order-one factors are not tracked.',
          'The normalization of the primordial gas (section 4.2) is a choice; the overabundance $O$ scales with it.',
          'OPEN: the production mechanism of the observed abundance; the precise bound on the Yukawa portal from nucleosynthesis and the CMB; the '
          'decays of heavy Kaluza–Klein modes made by the Higgs portal.',
          'The Casimir energy assumes periodic boundary conditions on the torus and no other bulk fields; it is one term of $\\Lambda_4$ and does not '
          'change the conclusion that the model has $w=-1$.'])

    # ==============================================================================================
    d.h1('10. Verification')
    d.p('Every statement above is a check of section 24 of 7+1.nb (Wolfram Language 15.0.1) and of the independent Python program tools/ref_dark.py '
        '(numpy, scipy, sympy, mpmath), which uses its own quadratures (the Kohn–Sham integrals by numerical quadrature instead of the closed forms '
        'of the notebook), its own root finding and its own lattice sums; the capture and decay results of the solutions are taken from the '
        'independent program tools/ref_coupled.py, which it runs. The numbers printed by both programs:')
    trows = []
    for lab, a, b, tol, r in cmp_:
        trows.append([lab, ', '.join(fmt(x, 7) for x in a[:5]) + (' …' if len(a) > 5 else '') if a else '?',
                      ', '.join(fmt(x, 7) for x in b[:5]) + (' …' if len(b) > 5 else '') if b else '?',
                      fmt(r, 2) if r is not None else '—', fmt(tol, 1)])
    d.table(['quantity', '7+1.nb', 'Python', 'max. relative difference', 'tolerance'], trows, widths=[0.2, 0.27, 0.27, 0.13, 0.13], small=True)
    d.p('The largest relative difference is %s. The checks of section 24 of the notebook:' % (fmt(worst, 2) if worst is not None else '?'))
    d.table(['#', 'check', 'result', 'section'],
            [[str(i + 1), '`' + c[0].replace('`', "'") + '`', c[1], (c[3][:4] if c[3] else c[2][:3])] for i, c in enumerate(ck)],
            widths=[0.06, 0.74, 0.08, 0.12], small=True)
    d.p('The Python program’s boolean checks:')
    d.table(['#', 'check', 'result'], [[str(i + 1), '`' + k.replace('`', "'") + '`', x] for i, (k, x) in enumerate(rb)],
            widths=[0.06, 0.82, 0.12], small=True)
    d.h1('11. How to reproduce')
    d.ol(['`python3 tools/make_7p1_nb.py` writes 7+1.nb (unevaluated).',
          'Evaluate 7+1.nb in a fresh kernel of Wolfram 15.0.1 (Evaluation ▸ Quit Kernel, then Evaluation ▸ Evaluate Notebook; or `Get["D:/Developer/github/7+1/tools/run_eval.wl"]`) and save it; the last section summarizes all checks.',
          '`python3 tools/ref_dark.py` repeats the checks of this document independently (exit code 0 if all pass).',
          '`python3 tools/make_docs.py` regenerates this document, its symbol dictionary and its provenance record.'])
    d.h1('12. References')
    d.ul(['Planck Collaboration, Astron. Astrophys. 641 (2020) A6 — cosmological parameters.',
          'DESI Collaboration, DESI DR2 Results II, arXiv:2503.14738; C. Garcia-Quintero et al. (DESI), arXiv:2504.18464 — evolving dark energy.',
          'LZ Collaboration, arXiv:2410.17036 — direct detection with 4.2 tonne-years.',
          'S. W. Randall, M. Markevitch, D. Clowe, A. H. Gonzalez, M. Bradač, ApJ 679 (2008) 1173 — the Bullet Cluster.',
          'V. Iršič et al., Phys. Rev. D 109 (2024) 043511 — Lyman-$\\alpha$; G. Ballesteros, M. A. G. Garcia, M. Pierre, JCAP 03 (2021) 101 — '
          'Lyman-$\\alpha$ bounds on non-thermal relics.',
          'J. Alvey et al., MNRAS 501 (2021) 1188 — the Pauli principle in dwarf spheroidal galaxies.',
          'B. M. Roach et al., Phys. Rev. D 107 (2023) 023009 — NuSTAR limits on sterile-neutrino decays.',
          'S. Navas et al. (Particle Data Group), Phys. Rev. D 110 (2024) 030001 — large extra dimensions, masses, $G_F$.',
          'L. J. Hall, K. Jedamzik, J. March-Russell, S. M. West, JHEP 03 (2010) 080 — freeze-in.',
          'P. B. Pal and L. Wolfenstein, Phys. Rev. D 25 (1982) 766; V. Barger, R. J. N. Phillips, S. Sarkar, Phys. Lett. B 352 (1995) 365 — decays '
          'of heavy (sterile) neutrinos.',
          'N. Arkani-Hamed, S. Dimopoulos, G. Dvali, Phys. Lett. B 429 (1998) 263 and Phys. Rev. D 59 (1999) 086004 — large extra dimensions, the '
          'normalcy temperature; N. Arkani-Hamed, S. Dimopoulos, G. Dvali, J. March-Russell, Phys. Rev. D 65 (2001) 024032 — bulk fermions coupled '
          'to the lepton doublets.',
          'N. M. Hugenholtz and L. Van Hove, Physica 24 (1958) 363 — $\\mu=\\varepsilon_F$; H. B. G. Casimir, Proc. K. Ned. Akad. Wet. 51 (1948) 793 — '
          'the Casimir energy; C. G. J. Jacobi (1829) — the four-square theorem; G. D. Coughlan et al., Phys. Lett. B 131 (1983) 59 — the moduli problem.'])
    return d
