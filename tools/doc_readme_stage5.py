"""The parts of README, README_SYMBOLS and README_PROVENANCE that describe section 24 of 7+1.nb (dirac16complex as dark matter and dark
energy: an analysis based on the accurate approximate solutions of sections 20-23).  Every equation used is written in the document that
uses it; every number is parsed from the evaluated notebook or from the output of tools/ref_dark.py (gathered by
doc_stage3_common.gather3)."""
from doc_stage3_common import fmt
from doc_dark import nbv, compare


def present(D3):
    return bool(D3 and D3.get('checks', {}).get('dk'))


def counts(D3):
    c = D3['checks'].get('dk', [])
    r = [(k, v) for k, v in D3.get('ref', {}).get('dk', ('', '', {}))[2].items() if v in ('True', 'False')]
    return len(c), sum(1 for x in c if x[1] == 'True'), len(r), sum(1 for k, v in r if v == 'True')


def _g(v, i, n=3):
    try:
        return fmt(v[i], n)
    except (IndexError, TypeError):
        return '?'


def abstract_text(D3):
    return ('; (viii) in section 24 of the notebook, an analysis — based on those solutions — of whether dirac16complex can be the dark '
            'matter and the dark energy of the universe (section 14)')


def results_paragraph(d, D3):
    n = counts(D3)
    dbeta = nbv(D3, 'dbeta'); over = nbv(D3, 'over'); bulk = nbv(D3, 'bulk'); tower = nbv(D3, 'tower'); hp = nbv(D3, 'hp')
    wmin = nbv(D3, 'wmin'); cas = nbv(D3, 'cas'); sigsi = nbv(D3, 'sigsi')
    d.p('**dirac16complex as dark matter and dark energy** (section 24 of the notebook; section 14 below). DERIVED from the solutions: the '
        'lightest Kaluza–Klein mode of dirac16complex is cold matter ($p_a/\\varepsilon=k_a^2/(3(m^2a^2+k_a^2))\\to0$) and is exactly the term '
        '$\\Omega_ca^{-3}$ of the late-universe solution, displacing the radion by only $\\delta\\beta=H_0^2\\Omega_c/(2m_r^2)=%s$; as dark matter it '
        'needs the stabilizing potential, without which a bulk gas re-expands the auxiliary space. ESTIMATE with couplings of natural size: it is '
        'dark ($\\sigma_{\\rm SI}=%s$ cm$^2$) and collisionless. The abundance is OPEN: the primordial gas and the radion overproduce it (by '
        '$%s\\,T_{\\rm RH}$/GeV; the radion gives the brane only $%s$ of its energy), and after the dilution required by the moduli problem the '
        'portals underproduce it (the Higgs portal would need $c_P=%s$; the Yukawa portal must satisfy $c_y<%s$ or it overproduces the '
        'Kaluza–Klein tower). It cannot be the dark energy: as matter $w\\ge%s$, and its vacuum energy ($%s\\,\\rho_\\Lambda$ on the torus if '
        'massless) is a constant contribution to $\\Lambda_4$, $w=-1$. Section 24 contains %d checks of its own.'
        % (_g(dbeta, 0, 2), _g(sigsi, 0, 2), _g(over, 1, 2), _g(bulk, 1, 2), _g(hp, 5, 3), _g(tower, 4, 2), _g(wmin, 0, 3), _g(cas, 0, 2), n[0]))


def notation_rows():
    return [
        ['$\\chi$, $\\mu_n$, $r_4(k)$', '`rIIr4`', 'lightest Kaluza–Klein mode; Kaluza–Klein masses $2\\pi|\\vec n|/V_\\perp^{1/4}$; lattice multiplicities (section 24)'],
        ['$c_P,c_y,c_V,c_r$', '`cPS`, `cyS`, `cVS`, `crS`', 'portal and self-couplings in units of $M_8$; natural size 1'],
        ['$mY$, $S(x)$, $T_*$', '`dkYm`, `dkS`, `dkTstar`', 'required abundance; Boltzmann factor of freeze-in; normalcy temperature'],
        ['$w$, $\\rho_{\\rm Cas}$', '`dkWtab`, `dkCas`', 'equation of state of Kohn–Sham matter; Casimir energy on the torus'],
    ]


def verification_text(D3):
    n = counts(D3); cmp_ = compare(D3)
    worst = max([r for *_, r in cmp_ if r is not None], default=None)
    return ('For section 24 the independent program `tools/ref_dark.py` (%d checks, %d True; numpy, scipy, sympy, mpmath, with its own '
            'quadratures, root finding and lattice sums, and the capture and decay results of `tools/ref_coupled.py`) repeats every number; the '
            '%d groups of numbers printed by both programs agree to a largest relative difference of %s, and `verify_7p1_nb.py` confirms that '
            'every cell and output of sections 0–23 is unchanged except one documented correction (section 14).'
            % (n[2], n[3], len(cmp_), fmt(worst, 2) if worst is not None else '?'))


def notebook_items(n_last, with_summary=True):
    return ['**24. dirac16complex as dark matter and dark energy** — the data, cold matter in the three regimes, darkness and stability, the '
            'abundance (overproduction by the solutions, the portals after the dilution, the Kaluza–Klein tower), dark energy (Kohn–Sham matter, '
            'condensate, Casimir energy, $w=-1$), the verdict (section 14).'] + (
           ['**%d. Summary of all checks** of sections 1–%d.' % (n_last + 1, n_last)] if with_summary else [])


def file_rows():
    return [
        ['`tools/nb_part_dark.py`', 'generator part of section 24 of 7+1.nb (executed by `make_7p1_nb.py`)'],
        ['`tools/ref_dark.py`', 'independent Python checks of section 24'],
    ]


def references():
    return ['L. J. Hall, K. Jedamzik, J. March-Russell and S. M. West, JHEP 03 (2010) 080 — freeze-in (section 14); Planck Collaboration, '
            'Astron. Astrophys. 641 (2020) A6; DESI Collaboration, arXiv:2503.14738 — dark-matter density and evolving dark energy (section 14).']


# =====================================================================================================================
def readme_section(d, D3, secno=14):
    """README section 14: the dark-sector analysis, self-contained"""
    n = counts(D3)
    v = {k: nbv(D3, k) for k in ('req', 'dbeta', 'sigsi', 'sigm', 'mmax', 'over', 'bulk', 'tstar', 'yzero', 'mstar', 'tower', 'hp', 'vec',
                                 'wmin', 'cas', 'xl', 'pauli')}
    d.h1('%d. Section 24: dirac16complex as dark matter and dark energy' % secno)
    d.p('Section 24 of 7+1.nb was appended on 27 October 2026 (the summary of all checks became section 25%s). It asks whether the field '
        'dirac16complex can be the dark matter and the dark energy, on the basis of the solutions of sections 20–23, the couplings of the model '
        'and present data, and labels every statement DERIVED (exact consequence of the field equations), ESTIMATE (order of magnitude, '
        'assumptions stated) or OPEN. Coordinates: $x_0$ time (the octonion unit), $x_1,x_2,x_3$ ordinary space, $x_4,\\dots,x_7$ auxiliary space; '
        '$\\eta=\\mathrm{diag}(+1,-1,\\dots,-1)$ is simply the squares of the octonion basis, and xAct declares the spacetime with $\\{7,1,0\\}$. '
        'Physical parameters: $M_8=10$ TeV ($\\kappa=M_8^{-6}$), $M_{\\rm Pl}=2.435\\times10^{18}$ GeV, $V_\\perp=M_{\\rm Pl}^2/M_8^6$, $m_r=10$ GeV. The '
        'section contains %d checks, of which %d are True. In this stage one cell of section 23.3 was corrected: the cold-dark-matter density '
        'in its table read "0.265" and now reads "Omega_c = 0.264" ($0.120/0.674^2=0.26416$).'
        % ('' if D3['n_last'] + 1 == 25 else ', and is now section %d' % (D3['n_last'] + 1), n[0], n[1]))
    d.h2('%d.1 Cold matter in the three regimes (DERIVED)' % secno)
    d.p('With the homogeneous metric $ds^2=dx_0^2-a^2d\\vec x^2-b^2d\\vec y^2$, a source with comoving energy $W(a,b)$ has '
        '$\\varepsilon=W/(a^3b^4)$, $p_a=-\\partial_aW/(3a^2b^4)$, $p_b=-\\partial_bW/(4a^3b^3)$. A massive particle of dirac16complex has '
        '$W=\\omega=\\sqrt{m^2+k_a^2/a^2+k_b^2/b^2}$, so')
    d.eq(r'\frac{p_a}{\varepsilon}\Big|_{k_b=0}=\frac{k_a^2}{3(m^2a^2+k_a^2)}\to0,\qquad \varepsilon\big|_{k_a=k_b=0}=\frac m{a^3b^4}\propto a^{-3}:')
    d.p('the lightest Kaluza–Klein mode (constant over the auxiliary torus) is cold matter. In regime II it is exactly the term $\\Omega_ca^{-3}$ of '
        'the brane Friedmann equation $H_a=H_0\\sqrt{\\Omega_ra^{-4}+\\Omega_ma^{-3}+\\Omega_\\Lambda}$; through the radion equation '
        '$\\beta\'\'+3H_a\\beta\'+4\\beta\'^2=\\kappa(\\varepsilon/6-p_a/2+p_b/3)+\\Lambda/3$ it shifts $\\beta=\\ln(b/b_0)$ by '
        '$\\delta\\beta=\\kappa\\varepsilon/(6m_r^2)=H_0^2\\Omega_c/(2m_r^2)=%s$. In regime I a bulk gas of dirac16complex dominates near $x_0=0$ and, '
        'without the stabilizing potential, re-expands the auxiliary space ($b\\propto\\tau^{2/5}$): dark matter of this kind needs the stabilizer. '
        'In regime III $\\rho_{\\rm DM}/\\rho_\\Lambda=(\\Omega_c/\\Omega_m)\\,\\mathrm{csch}^2(\\frac32H_\\infty t)\\to4(\\Omega_c/\\Omega_m)e^{-3H_\\infty t}$.'
        % _g(v['dbeta'], 0, 3))
    d.h2('%d.2 Dark, collisionless, stable (ESTIMATE)' % secno)
    d.p('A coupling of the lightest mode is diluted by the volume: with natural brane and bulk couplings $\\lambda_P=c_P/M_8^5$, $y=c_y/M_8^2$, '
        '$g_r=c_r/M_8^6$,')
    d.eq(r'\lambda_{\rm eff}=\frac{c_PM_8}{M_{\rm Pl}^2},\qquad y_{\rm eff}=\frac{c_yM_8}{M_{\rm Pl}},\qquad G_r=\frac{c_r}{M_{\rm Pl}^2},\qquad '
         r'\sigma_{\rm SI}=\frac{\lambda_{\rm eff}^2f_N^2m_N^2\mu^2}{\pi m_h^4},\qquad \sigma\approx\frac{G_r^2m^2}{4\pi},')
    d.p('so $\\sigma_{\\rm SI}=%s$ cm$^2$ at $m=40$ GeV (LZ: $2.2\\times10^{-48}$ cm$^2$) and $\\sigma/m=%s$ cm$^2$/g at 1 GeV (Bullet Cluster: '
        '$0.7$ cm$^2$/g). Without the Yukawa portal the mode is stable; with it, it mixes with the neutrinos ($\\theta=y_{\\rm eff}v/(\\sqrt2m)$) and '
        'decays, $\\Gamma(\\chi\\to3\\nu)=G_F^2m^5\\theta^2/(96\\pi^3)$, outliving the universe only for $m<%s$ MeV$\\times c_y^{-2/3}$. Its Pauli bound '
        '(16 internal states) is %s keV; Lyman-$\\alpha$ requires more than 5.7 keV for a thermal relic.'
        % (_g(v['sigsi'], 0, 2), _g(v['sigm'], 0, 2), fmt(v['mmax'][0]*1e3, 3) if v['mmax'] else '?', _g(v['pauli'], 0, 2)))
    d.h2('%d.3 The abundance (OPEN)' % secno)
    d.p('The observed abundance is $mY=\\Omega_ch^2(\\rho_c/h^2)/s_0=%s$ GeV ($Y=n/s$). The solutions of regimes I–II overproduce dirac16complex: '
        'the gas left after the capture of the auxiliary space exceeds it by $O=\\frac f{1-f}T_{\\rm RH}/(\\frac43mY)=%s\\,T_{\\rm RH}/$GeV even if the '
        'radion energy all became brane radiation, and in fact the radion decays into Kaluza–Klein pairs of dirac16complex, the brane receiving only '
        '$\\Gamma_{\\rm brane}/\\Gamma_{\\rm bulk}=%s$. After the dilution that the moduli problem requires, at $T_{\\rm RH}<T_*=(M_8^6/M_{\\rm Pl})^{1/5}=%s$ '
        'GeV, freeze-in by Higgs decays,'
        % (_g(v['req'], 3, 3), _g(v['over'], 1, 3), _g(v['bulk'], 1, 2), _g(v['tstar'], 0, 3)))
    d.eq(r'Y=\frac{135}{8\pi^3\cdot1.66\,g_*^{3/2}}\,\frac{M_{\rm Pl}^{\rm (unreduced)}\Gamma}{m_h^2}\,S\Big(\frac{m_h}{T_{\rm RH}}\Big),\qquad '
         r'S(x)=\frac2{3\pi}\int_x^\infty t^3K_1(t)\,dt,')
    d.p('makes too little through the portals with natural couplings: the lightest mode alone needs a Yukawa coupling that makes it decay '
        '(compatible only for $m<%s$ keV, below the Lyman-$\\alpha$ bound); the whole Kaluza–Klein tower, which couples like the lightest mode, '
        'is overproduced unless $c_y<%s$ ($\\sum_n(1-\\mu_n^2/m_h^2)^2=V_\\perp m_h^4/(192\\pi^2)=%s$); the Higgs portal, which couples every pair of '
        'levels, needs $c_P=%s$ at $T_{\\rm RH}=T_*$; the vector portal is of gravitational strength per pair and negligible. The origin of the '
        'abundance lies outside the minimal model.'
        % (fmt(v['mstar'][0]*1e6, 2) if v['mstar'] else '?', _g(v['tower'], 4, 2), _g(v['tower'], 0, 3), _g(v['hp'], 5, 3)))
    d.h2('%d.4 Dark energy (DERIVED)' % secno)
    d.p('Matter has $p=\\rho^2\\frac d{d\\rho}(\\varepsilon/\\rho)=\\mu\\rho-\\varepsilon$ at zero temperature: negative only below saturation, where uniform '
        'matter is metastable or unstable; for the Kohn–Sham matter of section 8 the most negative value is $w=%s$. A Lorentz-invariant condensate '
        'has $w=-1$. The Casimir energy on the torus, per real bosonic degree of freedom' % _g(v['wmin'], 0, 3))
    d.eq(r'\rho_{\rm Cas}^{(1)}=-\frac2{\ell_\perp^4}\sum_{\vec n\ne0}\Big(\frac{m\ell_\perp}{2\pi|\vec n|}\Big)^4K_4(m\ell_\perp|\vec n|)\;\xrightarrow{m\to0}\;'
         r'-\frac{\Gamma(4)}{\pi^4\ell_\perp^4}\,8(1-4^{-3})\zeta(4)\zeta(3)\qquad(\ell_\perp=V_\perp^{1/4}),')
    d.p('gives for dirac16complex (16 fermionic degrees of freedom) $%s\\,\\rho_\\Lambda$ if massless, $\\rho_\\Lambda$ at $m\\ell_\\perp=%s$; it is constant, '
        '$w=-1$, one more term of $\\Lambda_4$. The model therefore predicts $w=-1$ exactly; the DESI DR2 preference for $w_0>-1$, $w_a<0$ '
        '($2.8\\sigma$–$4.2\\sigma$), if confirmed, would contradict it. **Verdict:** dirac16complex is a viable dark-matter particle in every '
        'property the solutions test, if the auxiliary space is stabilized and the Yukawa portal is very weak, but the minimal model does not '
        'produce its abundance; it is not the dark energy. `tools/ref_dark.py` repeats the section with %d checks, all True.'
        % (_g(v['cas'], 0, 2), _g(v['xl'], 0, 4), n[2]))


# =====================================================================================================================
def symbols_section(d, D3, secno=12):
    d.h1('%d. Section 24 of 7+1.nb: dirac16complex as dark matter and dark energy' % secno)
    d.p('The symbols of section 24 (the Mathematica names begin with `dk`).')
    d.table(['symbol', 'Mathematica', 'definition', 'meaning'], [
        ['$\\chi$', '—', '$\\psi(x,y)=\\chi(x)/\\sqrt{V_\\perp}$', 'lightest Kaluza–Klein mode (four brane Dirac spinors)'],
        ['$\\mu_n$, $r_4(k)$', '`rIIr4`', '$2\\pi|\\vec n|/\\ell_\\perp$; $8\\sigma(k)-32\\sigma(k/4)$', 'Kaluza–Klein mass; number of $\\vec n\\in\\mathbb Z^4$ with $|\\vec n|^2=k$'],
        ['$mY$, $s_0$, $\\rho_c$', '`dkYm`, `dks0`, `dkRhoc`', '(%d.1)' % secno, 'required abundance, entropy density today, critical density'],
        ['$c_P,c_y,c_V,c_r$', '`cPS`, `cyS`, `cVS`, `crS`', '(%d.2)' % secno, 'couplings in units of $M_8$'],
        ['$\\lambda_{\\rm eff}$, $y_{\\rm eff}$, $G_r$', '`dkLam`, `dkYeff`, `dkGr`', '(%d.2)' % secno, 'four-dimensional couplings of the lightest mode'],
        ['$\\sigma_{\\rm SI}$, $\\theta$, $\\mu_P$', '`dkSigSI`, `dkMuP`', '(%d.3)' % secno, 'direct-detection cross section; neutrino mixing; portal mass'],
        ['$S(x)$, $T_*$, $Y$', '`dkS`, `dkTstar`, `dkYfi`', '(%d.4)' % secno, 'Boltzmann factor, normalcy temperature, freeze-in yield'],
        ['$I_{hp}$', '`dkIhp`', '(%d.4)' % secno, 'phase-space integral of Higgs-portal pairs'],
        ['$p$, $w$, $\\rho_{\\rm Cas}$', '`dkPres`, `dkWtab`, `dkCas`', '(%d.5)' % secno, 'pressure, equation of state, Casimir energy'],
    ], widths=[0.17, 0.22, 0.33, 0.28], small=True)
    d.align([r'mY&=\frac{\Omega_ch^2(\rho_c/h^2)}{s_0},\quad \frac{\rho_c}{h^2}=3\Big(100\,\frac{\rm km}{\rm s\,Mpc}\Big)^2M_{\rm Pl}^2,\quad s_0=\frac{2\pi^2}{45}\Big(2+\frac{21}4\cdot\frac4{11}\Big)T_{\rm CMB}^3,&&(%d.1)' % secno,
             r'\lambda_{\rm eff}&=\frac{\lambda_P}{V_\perp}=\frac{c_PM_8}{M_{\rm Pl}^2},\quad y_{\rm eff}=\frac{y}{\sqrt{V_\perp}}=\frac{c_yM_8}{M_{\rm Pl}},\quad G_r=\frac{g_r}{V_\perp}=\frac{c_r}{M_{\rm Pl}^2},&&(%d.2)' % secno,
             r'\sigma_{\rm SI}&=\frac{\lambda_{\rm eff}^2f_N^2m_N^2\mu^2}{\pi m_h^4},\quad \theta=\frac{\mu_P}m,\quad \mu_P=\frac{y_{\rm eff}v}{\sqrt2},&&(%d.3)' % secno,
             r'S(x)&=\frac2{3\pi}\int_x^\infty t^3K_1(t)dt,\quad T_*=\Big(\frac{M_8^6}{M_{\rm Pl}}\Big)^{1/5},&&(%d.4)' % secno,
             r'I_{hp}&=\iint_{u_1+u_2<1}u_1^3u_2^3\big(1-(u_1+u_2)^2\big)^{3/2}\big(1-(u_1-u_2)^2\big)^{1/2}du_1du_2,',
             r'p&=\mu\rho-\varepsilon,\quad w=\frac p\varepsilon,\quad \rho_{\rm Cas}^{(1)}=-\frac2{\ell_\perp^4}\sum_{\vec n\ne0}\Big(\frac{m\ell_\perp}{2\pi|\vec n|}\Big)^4K_4(m\ell_\perp|\vec n|).&&(%d.5)' % secno])


# =====================================================================================================================
def provenance_section(d, D, D3, runs, secno=9):
    n = counts(D3); cmp_ = compare(D3)
    worst = max([r for *_, r in cmp_ if r is not None], default=None)
    d.h1('%d. Section 24 of 7+1.nb: dirac16complex as dark matter and dark energy' % secno)
    d.p('This section records the provenance of section 24 of 7+1.nb (appended on 27 October 2026; the summary of all checks became '
        'section 25%s).' % ('' if D3['n_last'] + 1 == 25 else ', and is now section %d' % (D3['n_last'] + 1)))
    d.h2('%d.1 The task' % secno)
    d.ol(['On the basis of the accurate approximate solutions of the coupled field equations in the three gravity regimes, analyze the validity '
          'of associating dark matter (and dark energy) with dirac16complex (summarized; the original wording is the author’s private working '
          'material).',
          'Document the analysis in a documentation file with a provenance file and a symbol file, in Markdown and LaTeX, the LaTeX compiled to '
          'PDF; no document refers the reader to another file, and every equation used is written in the file that uses it.',
          'Push every completed stage to github.com/once-ere/7plus1 and verify the repository.'])
    d.h2('%d.2 Inputs, data and decisions' % secno)
    d.p('Inputs from the earlier sections: the stress tensor from the comoving energy (section 20.2), the radion equation (20.3), the gas-dominated '
        'regime without stabilizer (21.3), the gas share after capture (22.1), the radion decay rates into brane fields and into the bulk with '
        'Jacobi’s $r_4(k)$ (22.2), the late universe (22.4), the Standard Model constants (17.3), the uniform Kohn–Sham matter (8.7). Data: Planck '
        '2018; DESI DR2 (arXiv:2503.14738) and Garcia-Quintero et al. (arXiv:2504.18464); LZ (arXiv:2410.17036); Randall et al. (Bullet Cluster); '
        'Iršič et al. (Lyman-$\\alpha$); Ballesteros, Garcia and Pierre (decay-produced dark matter); Alvey et al. (Pauli bound, $g=2$, scaling '
        '$g^{-1/4}$); Roach et al. (NuSTAR); Particle Data Group 2024 (large extra dimensions). Decisions: natural couplings at $M_8$; Majorana '
        'normalization of the sterile-neutrino decay rates; instantaneous reheating at $T_{\\rm RH}\\le T_*$; production by Higgs decays as a lower '
        'bound; continuum Kaluza–Klein sums checked against lattice sums; Pauli bound for $g=16$. Before section 24 was appended the evaluated '
        'notebook was copied to `backups/7+1_2026-10-27_before_dark.nb` (excluded by .gitignore). One correction to section 23.3 was made: the '
        'cold-dark-matter density in its table read "0.265" and now reads "Omega_c = 0.264" ($0.120/0.674^2=0.26416$).')
    d.h2('%d.3 Verification' % secno)
    d.p('Section 24 of 7+1.nb contains %d checks (%d True); they are rows of section 4.2 above. The independent program `tools/ref_dark.py` '
        '(%d boolean checks, %d True) recomputes every number with its own code; the %d groups of numbers printed by both programs agree to a '
        'largest relative difference of %s; its complete output is printed in section 4.3. `tools/verify_7p1_nb.py` confirms that, with section 24 '
        'removed, the notebook equals the backup except the documented correction.' % (n[0], n[1], n[2], n[3], len(cmp_),
                                                                                       fmt(worst, 2) if worst is not None else '?'))
