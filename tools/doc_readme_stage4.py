"""The parts of README, README_SYMBOLS and README_PROVENANCE that describe the fourth stage of 7+1.nb (sections 20-23: the coupled
field equations of the whole model with every interaction in equation form, and their accurate approximate solutions in the
primordial, modern and far-future gravity regimes).  Every equation used is written in the document that uses it; every number is
parsed from the evaluated notebook or from the output of tools/ref_coupled.py (gathered by doc_stage3_common.gather3)."""
from doc_stage3_common import fmt, nums_in
from doc_coupled import nbv, compare


def present(D3):
    return bool(D3 and D3.get('checks', {}).get('cp'))


def counts(D3):
    c = D3['checks'].get('cp', [])
    r = [(k, v) for k, v in D3.get('ref', {}).get('cp', ('', '', {}))[2].items() if v in ('True', 'False')]
    return len(c), sum(1 for x in c if x[1] == 'True'), len(r), sum(1 for k, v in r if v == 'True')


def abstract_text(D3):
    return ('; (vii) in sections 20–23 of the notebook, the coupled field equations of the whole model with every interaction written as an '
            'equation, and accurate approximate solutions of them in the three gravity regimes — primordial (with its cubic, Gauss–Bonnet and '
            'Einstein sub-regimes), modern and far future — each compared with a high-accuracy numerical solution (section 13)')


def results_paragraph(d, D3):
    n = counts(D3)
    pq = nbv(D3, 'pq'); acc = nbv(D3, 'acc'); caprad = nbv(D3, 'caprad'); late4 = nbv(D3, 'late4'); comp = nbv(D3, 'comp')
    d.p('**The coupled field equations and their solutions in the three gravity regimes** (sections 20–23 of the notebook; section 13 below). '
        'Every interaction of the model — the four Lovelock terms, dirac16complex with its mass, nine self-couplings and spin connection, a '
        'stabilizing potential of the auxiliary space, the brane tension, the Standard Model on the brane and the portals — is a term of one '
        'action, whose variation gives coupled field equations; for the homogeneous metric $ds^2=dx_0^2-a^2d\\vec x^2-b^2d\\vec y^2$ they become')
    d.eq(r'\sum_{k=0}^{3}\alpha_kE^{(k)0}{}_0=\kappa\varepsilon,\qquad \sum_{k=0}^{3}\alpha_kE^{(k)1}{}_1=-\kappa p_a,\qquad \sum_{k=0}^{3}\alpha_kE^{(k)4}{}_4=-\kappa p_b .')
    d.p('In the primordial regime the vacuum solution is exact (an algebraic curve and one quadrature), with Kasner limits $(p_1,q_1)=(%s,%s)$ '
        '(Einstein) and $(p_2,q_2)=(%s,%s)$ (Gauss–Bonnet) whose closed-form corrections are accurate to %s and %s; a dirac16complex gas makes '
        'the universe matter-dominated near $x_0=0$ and, without a stabilizing potential, makes the auxiliary space re-expand. With the potential '
        'the radion is caught and carries %s of the energy (the moduli problem of the minimal model); the modern brane cosmology has the '
        'composite analytic solution accurate to %s (age %s Gyr), and the far future is de Sitter with a frozen auxiliary space. Sections 20–23 '
        'contain %d checks of their own.' % (fmt(pq[0], 5) if pq else '?', fmt(pq[1], 5) if pq else '?', fmt(pq[2], 5) if pq else '?',
                                             fmt(pq[3], 5) if pq else '?', fmt(acc[2], 2) if len(acc) > 2 else '?', fmt(acc[1], 2) if len(acc) > 1 else '?',
                                             fmt(caprad[0], 3) if caprad else '?', fmt(comp[0], 2) if comp else '?', fmt(late4[3], 5) if len(late4) > 3 else '?', n[0]))


def notation_rows():
    return [
        ['$\\beta$, $m_r$, $U(b)$', '`betR`, `rIIm`, `uPot`', 'radion $\\ln(b/b_0)$, its mass $m_r^2=\\kappa K_U/12$, the stabilizing potential (sections 20–23)'],
        ['$u$, $X$, $P_E,P_G,P_C$', '`uS`, `xS`, `rIpE`, `rIpG`, `rIpC`', 'variables and polynomials of the exact primordial vacuum curve'],
        ['$p_1,q_1$; $p_2,q_2$', '`rIp1`, `rIq1`; `rIp2`, `rIq2`', 'Einstein- and Gauss–Bonnet–Kasner exponents'],
        ['$S_1,S_2,S_3$, $f_{\\rm gas}$', '`rIstate[h]["S"]`, `["frac"]`', 'sizes of the curvature terms of the 00 equation; gas fraction'],
        ['$g_*$, $g_{*s}$, $t(T)$', '`rIIgstars`, `rIItg`', 'degrees of freedom of the brane plasma, time–temperature relation'],
        ['$H_\\infty$, $\\gamma$, $m_{\\rm eff}$', '`rIIIH`, `rIIIgam`, `rIIImeff`', 'final de Sitter rate; radion friction and mass in de Sitter'],
    ]


def verification_text(D3):
    n = counts(D3); cmp_ = compare(D3)
    worst = max([r for *_, r in cmp_ if r is not None], default=None)
    return ('For sections 20–23 the independent program `tools/ref_coupled.py` (%d checks, %d True; numpy, scipy, sympy, mpmath, with its own '
            'closed forms from `tools/lovelock71.py`, its own quadratures and an 8th-order Runge–Kutta integrator) repeats every derivation and '
            'every numerical solution; the %d groups of numbers printed by both programs agree to a largest relative difference of %s, and '
            '`verify_7p1_nb.py` confirms that every cell and output of sections 0–19 is unchanged.' % (n[2], n[3], len(cmp_), fmt(worst, 2) if worst is not None else '?'))


def notebook_items(n_last):
    return ['**20. The coupled field equations: every interaction** — the action term by term (the cubic Lovelock density in invariant form), '
            'the energy–momentum tensors, the homogeneous equations, the radion and the four-dimensional constants (section 13.1).',
            '**21. Regime I: the primordial Lovelock regime** — the exact vacuum curve, the cubic, Gauss–Bonnet and Einstein solutions and their '
            'accuracy, the matter-dominated cubic regime, decompactification, pair creation (section 13.2).',
            '**22. Regime II: the modern regime** — capture of the radion, its decay and the moduli problem, $g_*(T)$ and $t(T)$, the late '
            'universe and its age (section 13.3).',
            '**23. Regime III: the far future** — de Sitter, the radion in de Sitter, the size of every interaction (section 13.4).',
            '**%d. Summary of all checks** of sections 1–%d.' % (n_last + 1, n_last)]


def file_rows():
    return [
        ['`tools/nb_part_coupled.py`', 'generator part of sections 20–23 of 7+1.nb (executed by `make_7p1_nb.py`)'],
        ['`tools/ref_coupled.py`', 'independent Python checks of sections 20–23'],
    ]


def references():
    return ['D. Lovelock, J. Math. Phys. 12 (1971) 498; F. Müller-Hoissen, Phys. Lett. B 163 (1985) 106 — Lovelock gravity, the cubic density '
            '(section 13); N. Deruelle and L. Fariña-Busto, Phys. Rev. D 41 (1990) 3696 — Lovelock cosmology.',
            'G. D. Coughlan, W. Fischler, E. W. Kolb, S. Raby and G. G. Ross, Phys. Lett. B 131 (1983) 59 — the cosmological moduli problem '
            '(section 13.3); N. Arkani-Hamed, S. Dimopoulos and G. Dvali, Phys. Rev. D 59 (1999) 086004 — the normalcy temperature.']


# =====================================================================================================================
def readme_section(d, D3, secno=13):
    """README section 13: the fourth stage, self-contained"""
    n = counts(D3)
    pq = nbv(D3, 'pq'); u12 = nbv(D3, 'u12'); cub = nbv(D3, 'cub'); acc = nbv(D3, 'acc'); near = nbv(D3, 'near'); late = nbv(D3, 'late')
    capfit = nbv(D3, 'capfit'); caprad = nbv(D3, 'caprad'); decay = nbv(D3, 'decay'); late4 = nbv(D3, 'late4'); comp = nbv(D3, 'comp')
    hinf = nbv(D3, 'hinf'); radnow = nbv(D3, 'radnow'); tT = nbv(D3, 'tT')
    d.h1('%d. Sections 20–23: the coupled field equations and their solutions in the three gravity regimes' % secno)
    d.p('Sections 20–23 of 7+1.nb were appended in a fourth stage (27 October 2026); sections 0–19 are unchanged, and the summary of all checks, '
        'formerly section 20, is now section 24. The coordinates are those of the whole project: $x_0$ time (the octonion unit), $x_1,x_2,x_3$ '
        'ordinary space, $x_4,\\dots,x_7$ auxiliary space; $\\eta=\\mathrm{diag}(+1,-1,\\dots,-1)$ is simply the squares of the octonion basis, and '
        'xAct declares the spacetime with the signature list $\\{7,1,0\\}$ (seven space-like, one time-like dimension). A prime is $d/dx_0$; '
        '$\\kappa=8\\pi G_8=M_8^{-6}$. The sections contain %d checks, of which %d are True.' % (n[0], n[1]))
    d.h2('%d.1 Every interaction in equation form (section 20)' % secno)
    d.p('With the Lovelock densities $L_k=2^{-k}\\delta^{\\alpha_1\\beta_1\\cdots\\alpha_k\\beta_k}_{\\gamma_1\\delta_1\\cdots\\gamma_k\\delta_k}R^{\\gamma_1\\delta_1}{}_{\\alpha_1\\beta_1}\\cdots R^{\\gamma_k\\delta_k}{}_{\\alpha_k\\beta_k}$ '
        '($\\delta^{\\cdots}_{\\cdots}$ the generalized Kronecker delta), the real 16×16 Dirac matrices $\\Gamma^A$, the frame matrices '
        '$\\Gamma^\\mu=e_A{}^\\mu\\Gamma^A$, the spin connection $\\Omega_\\mu=\\frac14\\eta_{AC}\\omega_\\mu{}^C{}_B\\Gamma^A\\Gamma^B$, the Clifford products '
        '$\\Theta^K=i^{\\sigma_r}\\Gamma^{A_1}\\cdots\\Gamma^{A_r}$ ($\\sigma_r=\\frac{r(r-1)}2\\bmod2$, $\\eta_K=\\eta_{A_1A_1}\\cdots\\eta_{A_rA_r}$), the '
        'world brane $W=\\lbrace x_4=\\dots=x_7=0\\rbrace$ with induced metric $h$, and $b^8=\\det(-g_{kl})$ ($k,l=4,\\dots,7$), the action is')
    d.align([r'S&=-\frac1{2\kappa}\int d^8X\sqrt{|g|}\big(2\Lambda+R+\alpha_2L_2+\alpha_3L_3\big)',
             r'&\quad+\int d^8X\sqrt{|g|}\Big[\bar\psi\,i\Gamma^\mu(\partial_\mu+\Omega_\mu)\psi-m\bar\psi\psi',
             r'&\qquad\qquad-\frac12\sum_{r=0}^8g_r\sum_{|K|=r}\eta_K(\bar\psi\Theta^K\psi)^2-U(b)\Big]',
             r'&\quad+\int_Wd^4x\sqrt{|h|}\Big[-f_W+\mathcal L_{\rm SM}-\sum_{i,s}y_{is}\bar L_i\tilde\Phi P_R\psi_s+\text{h.c.}',
             r'&\qquad\qquad-G^f_P(\bar\psi\Gamma^\mu\psi)(\bar f\gamma_\mu f)-\lambda_P\Phi^\dagger\Phi\,\bar\psi\psi\Big],',
             r'L_2&=R^2-4R_{\mu\nu}R^{\mu\nu}+R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma},\qquad U(b)=U_0+\frac{K_U}2\ln^2\frac b{b_0},',
             r'L_3&=R^3-12RR_{\mu\nu}R^{\mu\nu}+16R^\mu{}_\nu R^\nu{}_\rho R^\rho{}_\mu+24R^{\mu\rho}R^{\nu\sigma}R_{\mu\nu\rho\sigma}+3RR_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}',
             r'&\quad-24R^{\mu\nu}R_\mu{}^{\rho\sigma\tau}R_{\nu\rho\sigma\tau}+4R^{\mu\nu}{}_{\rho\sigma}R^{\rho\sigma}{}_{\lambda\tau}R^{\lambda\tau}{}_{\mu\nu}-8R_\mu{}^\rho{}_\nu{}^\sigma R_\rho{}^\lambda{}_\sigma{}^\tau R_\lambda{}^\mu{}_\tau{}^\nu,'])
    d.p('where $\\mathcal L_{\\rm SM}$ is the Standard Model Lagrangian (gauge fields $G^a_\\mu,W^I_\\mu,B_\\mu$, Higgs doublet $\\Phi$, three generations, '
        'Yukawa couplings) in the gravitational field of $h$, $y_{is}$, $G^f_P$, $\\lambda_P$ are the Yukawa, vector and Higgs portal couplings, $L_i$ '
        'the lepton doublets and $\\psi_s$ the brane-right-handed parts of dirac16complex at the brane. Varying it gives '
        '$\\sum_k\\alpha_kE^{(k)\\mu}{}_\\nu=\\kappa T^\\mu{}_\\nu$ ($\\alpha_0=2\\Lambda$, $\\alpha_1=1$, $E^{(k)\\mu}{}_\\nu=-2^{-k-1}\\delta^{\\mu\\alpha_1\\cdots}_{\\nu\\gamma_1\\cdots}R\\cdots R$), '
        'the nonlinear Dirac equation of $\\psi$ with the portal sources at the brane, and the Standard Model equations with portal terms. For '
        'the homogeneous metric every source is $T^\\mu{}_\\nu=\\mathrm{diag}(\\varepsilon,-p_a\\times3,-p_b\\times4)$, and for a comoving energy $W(a,b)$')
    d.align([r"\varepsilon&=\frac W{a^3b^4},\quad p_a=-\frac{\partial_aW}{3a^2b^4},\quad p_b=-\frac{\partial_bW}{4a^3b^3},",
             r"\beta''+3H_a\beta'+4\beta'^2&=\kappa\Big(\frac\varepsilon6-\frac{p_a}2+\frac{p_b}3\Big)+\frac\Lambda3\quad(\beta=\ln\tfrac b{b_0},\ \alpha_2=\alpha_3=0),"])
    d.p('so the radion $\\beta$ has the mass $m_r^2=\\kappa K_U/12$, and with $b$ static the brane obeys $3H_a^2=8\\pi G_4\\rho_4+\\Lambda_4$ with '
        '$8\\pi G_4=\\kappa/V_\\perp$ ($V_\\perp=b^4L_c^4$ the volume of the auxiliary torus) and $\\Lambda_4=\\Lambda+\\kappa U_0+\\kappa f_W/V_\\perp$. The '
        'canonical radion $\\phi=\\sqrt{12}M_{\\rm Pl}\\beta$ couples to the brane through $(\\phi/(\\sqrt3M_{\\rm Pl}))T^\\mu{}_\\mu$.')
    d.h2('%d.2 Regime I: the primordial Lovelock regime (section 21)' % secno)
    d.p('In units $\\ell_2=\\sqrt{\\alpha_2}=1$ with $\\xi=\\alpha_3/\\alpha_2^2=10^{-6}$ and $u=H_b/H_a$, $X=H_a^2$ ($H_a=a\'/a$, $H_b=b\'/b$), the vacuum '
        'constraint is quadratic in $X$:')
    d.align([r'0&=3(2u^2+4u+1)-12u(u^3+12u^2+18u+4)X+360\xi u^3(3u+4)X^2,',
             r'\tau(H_a)&=\int_{H_a}^\infty\frac{dH}{F(H)},\qquad F=-\frac{dH_a}{dx_0}>0,'])
    d.p('an exact solution as an explicit curve plus one quadrature ($\\tau$ the time since the singularity). Its limits: the cubic regime '
        '$H_a=1/\\tau+A_0+A_1\\tau$, $H_b=h_0+h_1\\tau+h_2\\tau^2$ with $h_0=-\\sqrt{\\alpha_2/(30\\alpha_3)}=%s$; the Gauss–Bonnet Kasner regime '
        '$a\\propto\\tau^{p_2}$, $b\\propto\\tau^{q_2}$ ($p_2=3/(3+4u_2)$, $u_2^3+12u_2^2+18u_2+4=0$); the Einstein Kasner regime $p_1=(1+2\\sqrt2)/7=%s$, '
        '$q_1=(2-3\\sqrt2)/14=%s$, with $\\tau=p_1/H_a-e_1H_a+s_1$, $e_1=\\alpha_2(1+9\\sqrt2)/7$. Inside their regimes these formulas are accurate to '
        '%s (cubic), %s (Gauss–Bonnet) and %s (Einstein). With a dirac16complex gas, $p_a=w_a\\varepsilon$, the universe near $\\tau=0$ is matter-'
        'dominated: $a\\propto\\tau^{1/(1+w_a)}$ ($\\tau^{3/4}$ for the gas; measured %s) and $H_b\\to h_*=-\\sqrt{\\alpha_2/(54\\alpha_3)}$. Once the gas '
        'dominates the Einstein regime, the auxiliary space re-expands, $b\\propto\\tau^{2/5}$ (measured $\\tau H_b=%s$): a stabilizing potential is '
        'needed.' % (fmt(cub[0], 7) if cub else '?', fmt(pq[0], 7) if pq else '?', fmt(pq[1], 7) if pq else '?', fmt(acc[0], 2) if acc else '?',
                     fmt(acc[1], 2) if len(acc) > 1 else '?', fmt(acc[2], 2) if len(acc) > 2 else '?', fmt(near[0], 7) if near else '?',
                     fmt(late[1], 6) if len(late) > 1 else '?'))
    d.h2('%d.3 Regime II: the modern regime (section 22)' % secno)
    d.p('The potential catches the auxiliary space when $H\\approx m_r$; the radion then oscillates as')
    d.eq(r'\beta\approx\frac{A\sin(m_r\tau+\vartheta)}{m_r(\tau-\tau_E)},\qquad H_a+2H_b\approx\frac2{3(\tau-\tau_E)},')
    bulk = nbv(D3, 'bulk')
    d.p('(fitted: $A=%s$, largest deviation %s of the amplitude) and carries %s of the energy. Into brane fields a 10 GeV radion ($M_8=10$ TeV) '
        'decays through $\\Gamma(\\phi\\to f\\bar f)=N_cm_f^2m_r(1-4m_f^2/m_r^2)^{3/2}/(24\\pi M_{\\rm Pl}^2)$ and '
        '$\\Gamma(\\phi\\to gg)=b_{\\rm QCD}^2\\alpha_s^2m_r^3/(96\\pi^3M_{\\rm Pl}^2)$, which alone would reheat the brane only to %s GeV; but its '
        '$10^{13}$ decay channels of gravitational strength into Kaluza–Klein pairs of dirac16complex ($\\Gamma_n=3\\mu_n^2m_r(1-4\\mu_n^2/m_r^2)^{3/2}/(8\\pi M_{\\rm Pl}^2)$, '
        '$\\mu_n=2\\pi|n|/V_\\perp^{1/4}$) make it decay into the BULK within %s s. The minimal model therefore has a moduli problem and needs a later '
        'reheating of the brane. After it, the brane plasma ($g_*(T)$ from ideal gases of all Standard Model particles) gives '
        '$t(1\\,\\mathrm{MeV})=%s$ s, and the late universe ($H_a=H_0\\sqrt{\\Omega_ra^{-4}+\\Omega_ma^{-3}+\\Omega_\\Lambda}$, Planck 2018) has the '
        'composite analytic solution'
        % (fmt(capfit[0], 5) if capfit else '?', fmt(capfit[3], 2) if len(capfit) > 3 else '?', fmt(caprad[0], 3) if caprad else '?',
           fmt(decay[2], 3) if len(decay) > 2 else '?', fmt(bulk[2], 2) if len(bulk) > 2 else '?', fmt(tT[3], 5) if len(tT) > 3 else '?'))
    d.eq(r't(a)\approx\frac{2}{3H_0\Omega_m^2}\Big[(\Omega_ma-2\Omega_r)\sqrt{\Omega_r+\Omega_ma}+2\Omega_r^{3/2}\Big]+\frac{2\,\mathrm{arcsinh}\big(\sqrt{\Omega_\Lambda/\Omega_m}\,a^{3/2}\big)}{3H_0\sqrt{\Omega_\Lambda}}-\frac{2a^{3/2}}{3H_0\sqrt{\Omega_m}},')
    d.p('accurate to %s at all times; the age is %s Gyr. The radion sits at $\\beta_{\\rm eq}=%s$ today, so $G_4$ does not change.'
        % (fmt(comp[0], 2) if comp else '?', fmt(late4[3], 5) if len(late4) > 3 else '?', fmt(radnow[0], 2) if radnow else '?'))
    d.h2('%d.4 Regime III: the far future (section 23)' % secno)
    d.p('De Sitter expansion with $H_\\infty=H_0\\sqrt{\\Omega_\\Lambda}=%s$ s$^{-1}$, exactly $H_a=H_\\infty\\coth(\\tfrac32H_\\infty t)$ for matter and '
        '$\\Lambda$; the radion obeys $\\beta\'\'+\\gamma\\beta\'+m_{\\rm eff}^2(\\beta-\\beta_{\\rm eq})=0$ with' % (fmt(hinf[0], 4) if hinf else '?'))
    d.align([r'\frac{\gamma}{3H}&=\frac{1-16\alpha_2H^2+48\alpha_2^2H^4-48\alpha_3H^4+512\alpha_2^3H^6+192\alpha_2\alpha_3H^6-1024\alpha_2^4H^8}{(1-4\alpha_2H^2+32\alpha_2^2H^4)^2},',
             r'\frac{m_{\rm eff}^2}{m_r^2}&=\frac1{1-4\alpha_2H^2+32\alpha_2^2H^4},'])
    d.p('so the auxiliary space stays frozen; without the potential the far future would be eight-dimensional de Sitter, '
        '$21H^2-420\\alpha_2H^4+2520\\alpha_3H^6=\\Lambda$, with an inflating auxiliary space. **Checks.** Sections 20–23 contain %d checks, all True; '
        '`tools/ref_coupled.py` repeats them with %d checks, all True.' % (n[0], n[2]))


# =====================================================================================================================
def symbols_section(d, D3, secno=11):
    d.h1('%d. Sections 20–23 of 7+1.nb: the coupled field equations and their solutions' % secno)
    d.p('The symbols of the fourth stage (sections 20–23 of 7+1.nb). A prime is $d/dx_0$; $\\tau=x_0-x_s$ is the time since the singularity; '
        'regime I uses units $\\ell_2=\\sqrt{\\alpha_2}=1$.')
    d.table(['symbol', 'Mathematica', 'definition', 'meaning'], [
        ['$L_2$, $L_3$', '`feL3coef`', '(%d.1)' % secno, 'Gauss–Bonnet and cubic Lovelock densities'],
        ['$\\alpha_k$, $\\Lambda$, $\\kappa$', '`alf[k]`, `lamb`, `kapS`', '$\\alpha_0=2\\Lambda$, $\\alpha_1=1$, $\\kappa=8\\pi G_8$', 'Lovelock couplings, bulk cosmological constant'],
        ['$\\ell_2$, $\\xi$', '`rIxi`', '$\\ell_2=\\sqrt{\\alpha_2}$, $\\xi=\\alpha_3/\\alpha_2^2$', 'Gauss–Bonnet length, cubic ratio'],
        ['$\\varepsilon$, $p_a$, $p_b$, $W(a,b)$', '`epsM`, `paM`, `pbM`, `wW`', '(%d.2)' % secno, 'energy density, pressures, comoving energy'],
        ['$U(b)$, $U_0$, $K_U$, $b_0$', '`uPot`, `u0B`, `kU`', '$U_0+\\frac{K_U}2\\ln^2(b/b_0)$', 'stabilizing potential'],
        ['$\\beta$, $m_r$, $\\phi$', '`betR`, `rIIm`', '$\\ln(b/b_0)$, $m_r^2=\\kappa K_U/12$, $\\phi=\\sqrt{12}M_{\\rm Pl}\\beta$', 'radion, its mass, canonical radion'],
        ['$f_W$, $V_\\perp$, $G_4$, $\\Lambda_4$', '`fWb`, `vPerpS`', '$V_\\perp=b^4L_c^4$, $8\\pi G_4=\\kappa/V_\\perp$, $\\Lambda_4=\\Lambda+\\kappa U_0+\\kappa f_W/V_\\perp$', 'brane tension and four-dimensional constants'],
        ['$u$, $X$, $P_E,P_G,P_C$', '`uS`, `xS`, `rIpE`, `rIpG`, `rIpC`', '(%d.3)' % secno, 'vacuum curve of regime I'],
        ['$p_1,q_1,p_2,q_2$', '`rIp1`, `rIq1`, `rIp2`, `rIq2`', '(%d.3)' % secno, 'Kasner exponents of the Einstein and Gauss–Bonnet regimes'],
        ['$h_0$, $h_*$', '`rIcub[[1]]`, `rIhStar`', '$-\\sqrt{\\alpha_2/(30\\alpha_3)}$, $-\\sqrt{\\alpha_2/(54\\alpha_3)}$', 'auxiliary contraction rates near $\\tau=0$ (vacuum, matter)'],
        ['$p_{\\rm eff}$', '`rIpeff`', '$-H_a^2/H_a\'$', 'effective Kasner exponent'],
        ['$S_1,S_2,S_3$, $f_{\\rm gas}$', '`rIstate[h]["S"]`, `["frac"]`', '(%d.4)' % secno, 'sizes of the curvature terms; gas fraction'],
        ['$A$, $\\vartheta$, $\\tau_E$', '`rIIaF`, `rIIphF`, `rIItEF`', '$\\beta\\approx A\\sin(m_r\\tau+\\vartheta)/(m_r(\\tau-\\tau_E))$', 'radion oscillation after capture'],
        ['$g_*$, $g_{*s}$, $t(T)$', '`rIIgstars`, `rIItg`', '$\\rho=\\frac{\\pi^2}{30}g_*T^4$, $s=\\frac{2\\pi^2}{45}g_{*s}T^3$', 'brane plasma'],
        ['$H_\\infty$, $\\gamma$, $m_{\\rm eff}$', '`rIIIH`, `rIIIgam`, `rIIImeff`', '$H_0\\Omega_\\Lambda^{1/2}$; (%d.5)' % secno, 'de Sitter rate; radion friction and mass'],
    ], widths=[0.17, 0.22, 0.33, 0.28], small=True)
    d.align([r'L_2&=R^2-4R_{\mu\nu}R^{\mu\nu}+R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma},&&(%d.1)' % secno,
             r'L_3&=\frac18\delta^{\alpha_1\beta_1\alpha_2\beta_2\alpha_3\beta_3}_{\gamma_1\delta_1\gamma_2\delta_2\gamma_3\delta_3}R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}R^{\gamma_2\delta_2}{}_{\alpha_2\beta_2}R^{\gamma_3\delta_3}{}_{\alpha_3\beta_3},',
             r'\varepsilon&=\frac W{a^3b^4},\quad p_a=-\frac{\partial_aW}{3a^2b^4},\quad p_b=-\frac{\partial_bW}{4a^3b^3},&&(%d.2)' % secno,
             r'ds^2&=dx_0^2-a^2d\vec x^2-b^2d\vec y^2,',
             r'P_E&=3(2u^2+4u+1),\quad P_G=-12u(u^3+12u^2+18u+4),\quad P_C=360u^3(3u+4),&&(%d.3)' % secno,
             r'u&=\frac{H_b}{H_a},\quad X=H_a^2,\quad p_1=\frac{1+2\sqrt2}7,\quad q_1=\frac{2-3\sqrt2}{14},\quad p_2=\frac3{3+4u_2},\quad q_2=u_2p_2,',
             r'S_1&=3(H_a^2+4H_a|H_b|+2H_b^2),\qquad S_2=12\alpha_2|H_b|(4H_a^3+18H_a^2|H_b|+12H_aH_b^2+|H_b|^3),&&(%d.4)' % secno,
             r'S_3&=360\alpha_3H_a^2|H_b|^3(4H_a+3|H_b|),\qquad f_{\rm gas}=\frac{\kappa\varepsilon}{S_1+S_2+S_3+\kappa\varepsilon},',
             r'\frac\gamma{3H}&=\frac{1-16\alpha_2H^2+48\alpha_2^2H^4-48\alpha_3H^4+512\alpha_2^3H^6+192\alpha_2\alpha_3H^6-1024\alpha_2^4H^8}{(1-4\alpha_2H^2+32\alpha_2^2H^4)^2},&&(%d.5)' % secno,
             r'\frac{m_{\rm eff}^2}{m_r^2}&=\frac1{1-4\alpha_2H^2+32\alpha_2^2H^4}.'])

# =====================================================================================================================
def provenance_section(d, D, D3, runs, secno=8):
    n = counts(D3); cmp_ = compare(D3)
    worst = max([r for *_, r in cmp_ if r is not None], default=None)
    d.h1('%d. The fourth stage: sections 20–23 of 7+1.nb' % secno)
    d.p('This section records the provenance of the coupled field equations and their solutions in the three gravity regimes (sections 20–23 of '
        '7+1.nb, appended on 27 October 2026; the summary of all checks became section 24).')
    d.h2('%d.1 The task' % secno)
    d.ol(['Compute accurate approximate solutions of the coupled field equations of the model for all three gravity regimes; write each and '
          'every interaction in equation form; define all symbols; explain and analyze the solutions (summarized; the original wording is the '
          'author’s private working material).',
          'Document them in a documentation file with a provenance file and a symbol file, in Markdown and LaTeX, the LaTeX compiled to PDF; no '
          'document refers the reader to another file, and every equation used is written in the file that uses it.',
          'Then analyze whether dirac16complex can be the dark matter and the dark energy, on the basis of these solutions (a separate stage).',
          'Push every completed stage to github.com/once-ere/7plus1 and verify the repository.'])
    d.p('**Interpretation.** The "three gravity regimes" are the three epochs of the model: primordial (near $x_0=0$; inside it the cubic, '
        'Gauss–Bonnet and Einstein terms dominate in turn), modern (a stabilized auxiliary space, four-dimensional gravity on the brane) and far '
        'future (de Sitter).')
    d.h2('%d.2 Inputs, data and decisions' % secno)
    d.p('Inputs from the earlier sections: the closed forms of the Lovelock tensors of the cosmological metric (section 12.6), the pair-set '
        'densities (section 12.5), the Kasner forms (section 14.1). Data: Particle Data Group 2024/2025 (masses), Planck 2018 ($H_0=67.4$ km/s/Mpc, '
        '$\\Omega_m=0.315$), $T_{\\rm CMB}=2.7255$ K, $N_{\\rm eff}=3.044$, $M_{\\rm Pl}=2.435\\times10^{18}$ GeV. Choices: $M_8=10$ TeV, $\\alpha_2=1/M_8^2$, '
        '$\\alpha_3=10^{-6}/M_8^4$, $m_r=10$ GeV, $f_W=M_8^4$; in regime I units $\\ell_2=1$, $\\xi=10^{-6}$, and a massless collisionless dirac16complex '
        'gas normalized to $\\kappa\\varepsilon=10^{-6}\\times3H_a^2$ at $a=b=1$, $H_a=1$; the size of a curvature term is the sum of the absolute values '
        'of its monomials (in a Kasner regime the net value cancels). Before the new sections were appended the evaluated notebook was copied to '
        '`backups/7+1_2026-10-27_before_coupled.nb` (excluded by .gitignore).')
    d.h2('%d.3 Verification of the fourth stage' % secno)
    d.p('Sections 20–23 of 7+1.nb contain %d checks (%d True); they are rows of section 4.2 above. The independent program `tools/ref_coupled.py` '
        '(%d boolean checks, %d True) repeats every derivation and numerical solution with its own code; the %d groups of numbers printed by both '
        'programs agree to a largest relative difference of %s (tolerances from $10^{-14}$ for exact constants to $2\\times10^{-3}$ for quantities '
        'read off long numerical integrations); its complete output is printed in section 4.3.' % (n[0], n[1], n[2], n[3], len(cmp_),
                                                                                                 fmt(worst, 2) if worst is not None else '?'))
