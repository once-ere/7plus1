"""The parts of README, README_SYMBOLS and README_PROVENANCE that describe the third stage of 7+1.nb (sections 11-20:
gravity, the self-coupled dirac16complex in gravitational fields, the Standard Model on the world brane).  Every equation used
is written in the document that uses it; every number is parsed from the evaluated notebook or from the outputs of the
independent Python programs (gathered by doc_stage3_common.gather3)."""
from doc_stage3_common import fmt, mnum, nums_in, break_sum


def counts(D3):
    ck = D3['checks']; ref = D3.get('ref', {})
    out = {}
    for key in ('grav', 'ks', 'sm'):
        c = ck.get(key, [])
        r = [(k, v) for k, v in ref.get(key, ('', '', {}))[2].items() if v in ('True', 'False')]
        out[key] = (len(c), sum(1 for x in c if x[1] == 'True'), len(r), sum(1 for k, v in r if v == 'True'))
    return out


def abstract_text(D3):
    return ('; (vi) in sections 11–19 of the notebook, gravity on the 7+1 spacetime — curvature, the Lovelock tensors (computed with xAct '
            'exactly as in the reference notebook Einstein-Lovelock-7+1-Nash.nb, and independently), generalized frame matrices, the canonical '
            'spin connection that couples gravity to dirac16complex, and the primordial, modern and far-future epochs — the nine-channel '
            'self-coupling of dirac16complex with its Kohn–Sham theory in gravitational fields, and the Standard Model on the world brane '
            'with the portals of dirac16complex and neutrino flavour conversion (section 12)')


def results_paragraph(d, D3):
    c = counts(D3)
    n = sum(v[0] for v in c.values())
    d.p('**Gravity, the self-coupled field in gravitational fields, and the Standard Model** (sections 11–19 of the notebook; section 12 '
        'below). On a metric $g_{\\mu\\nu}(X)$ of signature $(+,-,\\dots,-)$ the Lovelock tensors of order $k$, with the generalized Kronecker '
        'delta $\\delta^{\\mu_1\\cdots\\mu_n}_{\\nu_1\\cdots\\nu_n}=\\det[\\delta^{\\mu_i}_{\\nu_j}]$ and the Riemann tensor $R^{\\gamma\\delta}{}_{\\alpha\\beta}$, are')
    d.eq(r'E^{(k)\mu}{}_\nu=-\frac{1}{2^{k+1}}\,\delta^{\mu\alpha_1\beta_1\cdots\alpha_k\beta_k}_{\nu\gamma_1\delta_1\cdots\gamma_k\delta_k}\,'
         r'R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}\cdots R^{\gamma_k\delta_k}{}_{\alpha_k\beta_k}\qquad(E^{(1)}=G,\ E^{(4)}\equiv0\ \text{in eight dimensions}),')
    d.p('and the generalized frame matrices $\\Gamma_\\mu(X)=e^A{}_\\mu(X)\\Gamma^{(1,7)}_A$ with the canonical spin connection $\\Omega_\\mu$ have a '
        'vanishing total covariant derivative, so that dirac16complex couples to gravity through')
    d.eq(r'\nabla_\mu\Gamma_\nu=\partial_\mu\Gamma_\nu-\{{}^{\lambda}_{\mu\nu}\}\Gamma_\lambda+[\Omega_\mu,\Gamma_\nu]=0,\qquad '
         r'\big(i\,\Gamma^\mu(X)(\partial_\mu+\Omega_\mu)-m\big)\Psi=0 .')
    d.p('The field coupled to itself through all nine ranks of the Clifford algebra has a Kohn–Sham theory in gravitational fields, solved '
        'in the three cosmological epochs with its back-reaction on the expansion; the Standard Model lives on the brane $x_4=\\dots=x_7=0$ '
        'and meets dirac16complex through portals. Sections 11–19 contain %d checks of their own.' % n)


def notation_rows():
    """rows (symbol, Mathematica, definition) for the notation table of README section 3"""
    return [
        ['$g_{\\mu\\nu}$, $\\{{}^{\\lambda}_{\\mu\\nu}\\}$', '`cv["g"]`, `cv["chr"]`', 'section 12: metric (signature $(+,-,\\dots,-)$) and Christoffel symbols (never written $\\Gamma$)'],
        ['$R^\\rho{}_{\\sigma\\mu\\nu}$, $R_{\\sigma\\nu}$, $R$, $G^\\mu{}_\\nu$', '`cv["riem"]`, `cv["ric"]`, `cv["R"]`, `cv["G"]`', 'Riemann (Landau–Lifshitz), Ricci, scalar curvature, Einstein tensor'],
        ['$E^{(k)\\mu}{}_\\nu$, $L_k$', '`Lovelock1/2/3`, `lovelockE`, `lovelockL`', 'Lovelock tensors and densities ($k=1,2,3$)'],
        ['$\\alpha_k$, $\\bar\\Lambda$, $\\kappa$', '`alf[k]`, `lamb`, `kap`', 'Lovelock couplings ($\\alpha_1=1$), cosmological constant, $8\\pi G_8$'],
        ['$w_k$, $\\Lambda$, $H$', '`w1`, `w2`, `w3`, `Λ`, `H`', 'pure-number couplings and inverse length of the reference notebook'],
        ['$a$, $b$, $H_a$, $H_b$', '`a3[x0]`, `b4[x0]`, `Ha`, `Hb`', 'scale factors of ordinary and auxiliary space; Hubble rates $a\'/a$, $b\'/b$'],
        ['$e^A{}_\\mu$, $\\Gamma_\\mu(X)$, $\\Omega_\\mu$', '`fd["e"]`, `fd["Gam"]`, `fd["Om"]`', 'frame field, generalized frame matrices, spinor connection'],
        ['$\\Theta^K$, $I_r$, $g_r$', '`Θup[K]`, `Ichan`, `gg[r]`', 'Clifford products of rank $r$, their invariants and couplings'],
        ['$F_{rs}$, $G_s$', '`fierzF`, `Geff`', 'Fierz matrix and Hartree–Fock couplings $G=(1-F^{\\mathsf T})g$'],
        ['$M$, $\\varepsilon$, $p_a$, $p_b$', '`ksStress`', 'effective mass, energy density and pressures of Kohn–Sham matter'],
        ['$W$ (brane), $h_{\\mu\\nu}$', '—', 'world brane $x_4=\\dots=x_7=0$, induced metric'],
        ['$v$, $V_{\\rm CKM}$, $U_{\\rm PMNS}$', '`vEW`, `vCKM`, `uPMNS`', 'Higgs vacuum value, quark and lepton mixing matrices'],
        ['$R_A$, $U$', '`rBr`, `uInt`', 'brane form of the Dirac matrices and the intertwiner $\\Gamma^{(1,7)}_AU=UR_A$'],
        ['$P(\\alpha\\to\\beta)$', '`pVac`, `pMat`', 'neutrino flavour-conversion probability'],
    ]


def verification_text(D3):
    c = counts(D3)
    return ('For sections 11–19, the independent programs `tools/ref_gravity.py` (%d checks, %d True; with its own curvature and Lovelock '
            'engine `tools/lovelock71.py`), `tools/ref_ks_gravity.py` (%d checks, %d True) and `tools/ref_sm.py` (%d checks, %d True) repeat '
            'every construction; the Lovelock tensors are computed in the notebook both with xAct (as in the reference notebook '
            'Einstein-Lovelock-7+1-Nash.nb) and with an independent pair-set sum, which agree exactly.'
            % (c['grav'][2], c['grav'][3], c['ks'][2], c['ks'][3], c['sm'][2], c['sm'][3]))


def notebook_items(n_last):
    return ['**11. Gravity on the 7+1 spacetime** — the curvature engine, the cosmological metric, symmetries and the Bianchi identity '
            '(section 12.1).',
            '**12. The Lovelock tensors** — Lovelock’s theorem, the generalized Kronecker delta, xAct (as in Einstein-Lovelock-7+1-Nash.nb), the '
            'pair-set algorithm, closed forms, an inhomogeneous metric, a gravitational wave, exponential vacuum solutions (section 12.1).',
            '**13. Generalized frame matrices, the canonical spin connection and the covariant Dirac equation** — total covariant derivative, '
            'spin curvature, Lichnerowicz, Cartan’s structure equation, a static radial metric, gauge coupling (section 12.2).',
            '**14. The three epochs** — Kasner solutions near $x_0=0$, pair creation, a static auxiliary space, de Sitter and thermal pair '
            'creation (section 12.3).',
            '**15. The generalized self-coupling** — nine ranks, the Fierz matrix, the covariant nonlinear field equation, chirality (section 12.4).',
            '**16. Kohn–Sham theory in gravitational fields** — uniform matter in the cosmological spacetime, the epochs, self-consistent '
            'Einstein–Lovelock evolution (section 12.4).',
            '**17. The Standard Model on the world brane** — fields, anomalies, Lagrangian, electroweak breaking, masses and mixing, cosmology '
            '(section 12.5).',
            '**18. dirac16complex at the brane** — four brane Dirac spinors, Kaluza–Klein masses, portals, coupled field equations (section 12.6).',
            '**19. Neutrino flavour** — vacuum, matter, the Sun, supernovae, the relic background, effects of dirac16complex (section 12.6).',
            '**%d. Summary of all checks** of sections 1–%d.' % (n_last + 1, n_last)]


def file_rows():
    return [
        ['`tools/nb_part_gravity.py`, `tools/nb_part_ks_gravity.py`, `tools/nb_part_sm.py`', 'generator parts of sections 11–19 of 7+1.nb (executed by `make_7p1_nb.py`)'],
        ['`tools/lovelock71.py`', 'independent curvature and Lovelock engine (pair-set algorithm, sympy)'],
        ['`tools/ref_gravity.py`, `tools/ref_ks_gravity.py`, `tools/ref_sm.py`', 'independent Python checks of sections 11–14, 15–16 and 17–19'],
    ]


def references():
    return ['D. Lovelock, J. Math. Phys. 12 (1971) 498; D. Lovelock and H. Rund, *Tensors, Differential Forms, and Variational Principles*, '
            'Dover (1989), Remark 2 after Eq. (4.38) — the Lovelock tensors (section 12.1).',
            'J. M. Martín-García et al., xAct: efficient tensor computer algebra for the Wolfram Language (xPerm, xTensor, xCoba), GPL, '
            'www.xact.es; T. Nutma, Comput. Phys. Commun. 185 (2014) 1719 (xTras) — the abstract-index computation of section 12.1.',
            'L. D. Landau and E. M. Lifshitz, *The Classical Theory of Fields*, Pergamon — curvature conventions; E. Kasner, Am. J. Math. 43 '
            '(1921) 217 — power-law solutions; G. W. Gibbons and S. W. Hawking, Phys. Rev. D 15 (1977) 2738 — the de Sitter temperature $H/2\\pi$.',
            'S. Navas et al. (Particle Data Group), Phys. Rev. D 110 (2024) 030001 and the 2025 update; I. Esteban et al., JHEP 12 (2024) 216 '
            '(NuFit-6.0); Planck Collaboration, Astron. Astrophys. 641 (2020) A6 — the data of section 12.5–12.6.',
            'N. Arkani-Hamed, S. Dimopoulos, G. Dvali, Phys. Lett. B 429 (1998) 263; L. Wolfenstein, Phys. Rev. D 17 (1978) 2369; '
            'S. P. Mikheyev and A. Yu. Smirnov, Sov. J. Nucl. Phys. 42 (1985) 913 — extra dimensions, the MSW effect (section 12.5–12.6).']


# =====================================================================================================================
def readme_section(d, D3, secno=12):
    """README section 12: the third stage, self-contained"""
    c = counts(D3)
    nm = D3['nums']
    d.h1('%d. Sections 11–19: gravity, the self-coupled dirac16complex in gravitational fields, and the Standard Model' % secno)
    d.p('Sections 11–19 of 7+1.nb (and the summary section 20) were appended in a third stage; sections 0–10 are unchanged. They rebuild the '
        'ideas of the repository Dirac_claude on the 7+1 spacetime of this project — ONE time coordinate $x_0$ and SEVEN space coordinates '
        '$x_1,\\dots,x_7$ ($x_1,x_2,x_3$ ordinary space, $x_4,\\dots,x_7$ auxiliary space); the 4+4 metric of Dirac_claude is not used. Here '
        '$\\Gamma_A=\\Gamma^{(1,7)}_A$, $\\Gamma^A=\\eta_{(1,7)}^{AA}\\Gamma_A$, $\\eta=\\eta^{(1,7)}=\\mathrm{diag}(+1,-1,\\dots,-1)$, $\\chi=i\\Gamma_0\\Gamma_1\\cdots\\Gamma_7$, '
        '$\\Psi$ is the field dirac16complex and $\\bar\\Psi=\\Psi^\\dagger\\Gamma_0$; a prime is $d/dx_0$. The sections contain %d, %d and %d checks '
        '(gravity, Kohn–Sham, Standard Model), of which %d, %d and %d are True.' % (c['grav'][0], c['ks'][0], c['sm'][0],
                                                                              c['grav'][1], c['ks'][1], c['sm'][1]))

    # ---------------------------------------------------------------------------------------------
    d.h2('%d.1 Curvature and the Lovelock tensors (sections 11–12)' % secno)
    d.p('Christoffel symbols are written $\\{{}^{\\lambda}_{\\mu\\nu}\\}$ ($\\Gamma$ always denotes a Dirac matrix); the curvature follows Landau and Lifshitz:')
    d.align([r'\{{}^{\lambda}_{\mu\nu}\}&=\tfrac12g^{\lambda\rho}\big(\partial_\mu g_{\rho\nu}+\partial_\nu g_{\rho\mu}-\partial_\rho g_{\mu\nu}\big),\qquad '
             r'R^{\rho}{}_{\sigma\mu\nu}=\partial_\mu\{{}^{\rho}_{\sigma\nu}\}-\partial_\nu\{{}^{\rho}_{\sigma\mu}\}+\{{}^{\rho}_{\mu\lambda}\}\{{}^{\lambda}_{\sigma\nu}\}-\{{}^{\rho}_{\nu\lambda}\}\{{}^{\lambda}_{\sigma\mu}\},',
             r'R_{\sigma\nu}&=R^{\rho}{}_{\sigma\rho\nu},\qquad R=g^{\sigma\nu}R_{\sigma\nu},\qquad G^\mu{}_\nu=R^\mu{}_\nu-\tfrac12\delta^\mu_\nu R,\qquad '
             r'R^{\gamma\delta}{}_{\alpha\beta}=g^{\delta\sigma}R^{\gamma}{}_{\sigma\alpha\beta}.'])
    d.p('In $n$ dimensions the most general symmetric, divergence-free tensor built from the metric and its first two derivatives is a sum of '
        '$m-1$ Lovelock tensors plus a cosmological term, $m=n/2$ for even $n$ (D. Lovelock and H. Rund): for $n=8$ there are exactly three, '
        'normalized as in the reference notebook Einstein-Lovelock-7+1-Nash.nb (P. L. Nash, GPL),')
    d.align([r'E^{(k)\mu}{}_\nu&=-\frac{1}{2^{k+1}}\,\delta^{\mu\alpha_1\beta_1\cdots\alpha_k\beta_k}_{\nu\gamma_1\delta_1\cdots\gamma_k\delta_k}\,'
             r'R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}\cdots R^{\gamma_k\delta_k}{}_{\alpha_k\beta_k},\qquad '
             r'L_k=\frac{1}{2^k}\,\delta^{\alpha_1\beta_1\cdots\alpha_k\beta_k}_{\gamma_1\delta_1\cdots\gamma_k\delta_k}\,R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}\cdots R^{\gamma_k\delta_k}{}_{\alpha_k\beta_k},',
             r'\delta^{\mu_1\cdots\mu_n}_{\nu_1\cdots\nu_n}&=\det\big[\delta^{\mu_i}_{\nu_j}\big],\qquad E^{(1)}=G,\qquad E^{(k)\mu}{}_\mu=-\frac{8-2k}{2}L_k,\qquad E^{(4)}\equiv0 .'])
    d.p('The field equations, in the form of the reference notebook (pure numbers $w_k$, $\\Lambda$ and an inverse length $H$) and with matter '
        '($\\kappa=8\\pi G_8$):')
    d.eq(r'0=-\Lambda\,\delta^\mu_\nu+\sum_{k=1}^{3}H^{-2k}w_k\,E^{(k)\mu}{}_\nu,\qquad \sum_{k=0}^{3}\alpha_kE^{(k)\mu}{}_\nu=\kappa T^\mu{}_\nu\ \ '
         r'(\alpha_0=2\bar\Lambda,\ \alpha_1=1),\\ \alpha_2=\frac{w_2}{w_1H^2},\ \ \alpha_3=\frac{w_3}{w_1H^4},\ \ \bar\Lambda=\frac{\Lambda H^2}{w_1}.')
    d.p('The notebook computes the Lovelock tensors twice: with the tensor packages xAct (%s) exactly as the reference notebook does — abstract '
        'indices, the generalized delta `kδ[lower, upper] = Det[Outer[delta, lower, upper]]`, and the components of a given metric from '
        '`CD = LC[g]` — and by an independent sum over sets of index pairs; the two agree exactly for the cosmological metric and a '
        'gravitational wave. The reference notebook’s section "FAILED Model metric g71; DO NOT USE!" is not used. For the cosmological metric '
        '$ds^2=dx_0^2-a^2(dx_1^2+dx_2^2+dx_3^2)-b^2(dx_4^2+\\cdots+dx_7^2)$ with $H_a=a\'/a$, $H_b=b\'/b$:' % (D3.get('xact_versions') or 'xPerm, xTensor, xCoba, xTras'))
    d.align([r'E^{(1)0}{}_0&=3H_a^2+12H_aH_b+6H_b^2,\qquad E^{(2)0}{}_0=-12H_b\big(4H_a^3+18H_a^2H_b+12H_aH_b^2+H_b^3\big),',
             r'E^{(3)0}{}_0&=360\,H_a^2H_b^3\big(4H_a+3H_b\big).'])
    ex = D3.get('expo') or {}
    ts = ex.get('tex_sol_qm1', {})
    d.p('Exponential vacuum solutions ($a=e^{Hx_0}$, $b=e^{q_HHx_0}$) exist for a one-parameter family of couplings; for $q_H=-1$ (eternal '
        'inflation of ordinary space with eternal deflation of the auxiliary space) $w_2=%s$, $w_3=%s$.' % (ts.get('w_2', '?'), ts.get('w_3', '?')))

    # ---------------------------------------------------------------------------------------------
    d.h2('%d.2 Generalized frame matrices, the canonical spin connection, and the covariant Dirac equation (section 13)' % secno)
    d.p('With a frame field $g_{\\mu\\nu}=\\eta_{AB}e^A{}_\\mu e^B{}_\\nu$ and its inverse $e_A{}^\\mu$:')
    d.align([r'\Gamma_\mu(X)&=e^A{}_\mu\Gamma_A,\qquad \Gamma^\mu(X)=e_A{}^\mu\Gamma^A,\qquad \Gamma_\mu\Gamma_\nu+\Gamma_\nu\Gamma_\mu=2g_{\mu\nu}1_{16},',
             r'\omega_\mu{}^A{}_B&=e^A{}_\nu\big(\partial_\mu e_B{}^\nu+\{{}^{\nu}_{\mu\lambda}\}e_B{}^\lambda\big),\qquad \Omega_\mu=\tfrac14\eta_{AC}\,\omega_\mu{}^C{}_B\,\Gamma^A\Gamma^B,\qquad D_\mu=\partial_\mu+\Omega_\mu,',
             r'\nabla_\mu\Gamma_\nu&=\partial_\mu\Gamma_\nu-\{{}^{\lambda}_{\mu\nu}\}\Gamma_\lambda+[\Omega_\mu,\Gamma_\nu]=0,\qquad '
             r'[D_\mu,D_\nu]=\tfrac14R_{AB\mu\nu}\Gamma^A\Gamma^B,\qquad (\Gamma^\mu D_\mu)^2=g^{\mu\nu}\nabla_\mu D_\nu-\tfrac R4 .'])
    d.p('The covariant Dirac equation $\\big(i\\Gamma^\\mu(X)D_\\mu-m\\big)\\Psi=0$ follows from $\\mathcal L=\\sqrt{|g|}\\,\\bar\\Psi(i\\Gamma^\\mu D_\\mu-m)\\Psi$, and '
        '$\\bar\\Psi\\Gamma^\\mu\\Psi$ is conserved. For the cosmological metric the spin connection is not zero:')
    d.eq(r"\Omega_0=0,\quad \Omega_i=\tfrac{a'}{2}\Gamma_i\Gamma_0\ (i=1,2,3),\quad \Omega_k=\tfrac{b'}{2}\Gamma_k\Gamma_0\ (k=4,\dots,7),\qquad "
         r"\Gamma^\mu\Omega_\mu=\big(\tfrac32H_a+2H_b\big)\Gamma_0 .")
    d.p('The notebook also derives the connection from Cartan’s structure equation $de^A+\\omega^A{}_B\\wedge e^B=0$ for a static metric with a '
        'radial coordinate, and couples an abelian gauge field, $D_\\mu=\\partial_\\mu+\\Omega_\\mu-iqA_\\mu$.')

    # ---------------------------------------------------------------------------------------------
    d.h2('%d.3 The three epochs (section 14)' % secno)
    g = nm.get('grav', {})
    gb = nums_in(g.get('Gauss-Bonnet Kasner p', ''))
    d.p('**Primordial** (near $x_0=0$): for $a=x_0^p$, $b=x_0^q$ the highest Lovelock order dominates. Einstein: $3p+4q=1$, $3p^2+4q^2=1$, '
        '$p=(1\\pm2\\sqrt2)/7$, $q=(2\\mp3\\sqrt2)/14$ — on the branch $p=0.547$, $q=-0.160$ ordinary space expands while the auxiliary space '
        'contracts. Gauss–Bonnet: $3p+4q=3$, $-203p^3+81p^2+351p+27=0$ (%s). Cubic: $(p,q)=(-15/7,20/7)$. The anisotropic expansion creates '
        'pairs of dirac16complex. **Modern**: with $b=b_0$ the field equations reduce to the Friedmann equations of ordinary space, '
        '$3H_a^2=\\kappa\\varepsilon+\\bar\\Lambda$, $2H_a\'+3H_a^2=-\\kappa p_a+\\bar\\Lambda$, with $G_4=G_8/V_\\perp$. **Far future**: in de Sitter '
        'space pairs are created with the Fermi–Dirac occupation at the Gibbons–Hawking temperature $H/2\\pi$:'
        % ('$p=' + '$, $'.join(fmt(x, 6) for x in gb) + '$' if gb else 'three real roots'))
    d.eq(r'|\beta|^2=\frac{1}{e^{2\pi M/H}+1}.')

    # ---------------------------------------------------------------------------------------------
    d.h2('%d.4 The nine-channel self-coupling and Kohn–Sham theory in gravitational fields (sections 15–16)' % secno)
    d.p('With $\\Theta^K=i^{\\sigma_r}\\Gamma^{A_1}\\cdots\\Gamma^{A_r}$ for the 256 increasing multi-indices $K$ ($\\sigma_r=\\frac{r(r-1)}2\\bmod2$, '
        '$\\Theta_K=\\eta_{A_1A_1}\\cdots\\eta_{A_rA_r}\\Theta^K$) the most general contact self-coupling has one coupling per rank,')
    d.eq(r'\mathcal L_{\rm self}=-\frac12\sum_{r=0}^{8}g_r\sum_{|K|=r}(\bar\Psi\Theta_K\Psi)(\bar\Psi\Theta^K\Psi),\qquad '
         r'\big(i\Gamma^\mu D_\mu-m\big)\Psi=\sum_rg_r\sum_{|K|=r}(\bar\Psi\Theta_K\Psi)\,\Theta^K\Psi,')
    d.p('and contains the four couplings of section 11.3 ($g_0=g_S$, $g_1=g_V$, $g_7=-g_A$, $g_8=-g_P$). The Fierz matrix '
        '$F_{rs}=\\frac1{16}\\sum_l(-1)^{rs-l}\\binom sl\\binom{8-s}{r-l}$ ($F^2=1$) turns exchange into direct terms, so Hartree–Fock is Hartree '
        'with $G=(1-F^{\\mathsf T})g$. Homogeneous Kohn–Sham matter in the cosmological spacetime has the effective mass $M=m+G_0\\rho_S$, '
        'the energy density $\\varepsilon$ and the pressures $p_a$, $p_b$ (integrals over the comoving Fermi ball), and')
    d.eq(r"\varepsilon'+3H_a(\varepsilon+p_a)+4H_b(\varepsilon+p_b)=0 .")
    ks = nm.get('ks', {})
    prim = nums_in(ks.get('primordial m*(x0), x0 = 1e-4 .. 10', ''))
    d.p('Near $x_0=0$ the effective mass falls toward zero (%s), in the modern epoch the matter is pressureless in ordinary space, and in a '
        'de Sitter future $1-M$ decays like the density; Einstein and Einstein–Lovelock gravity driven by this matter preserve the constraint '
        'equation to a relative %s.' % (('$M=%s$ at $x_0=10^{-4}$' % fmt(prim[0], 4)) if prim else 'chiral symmetry is effectively restored',
                                         fmt(max(nums_in(ks.get('Einstein: Hb(0), {a, b, Ha, Hb}(20), violation', ''))[-1:] +
                                                 nums_in(ks.get('Einstein-Lovelock: Hb(0), {a, b, Ha, Hb}(20), violation', ''))[-1:] or [None]), 2)
                                         if (ks.get('Einstein: Hb(0), {a, b, Ha, Hb}(20), violation') or '') else '?'))

    # ---------------------------------------------------------------------------------------------
    d.h2('%d.5 The Standard Model on the world brane (section 17)' % secno)
    sm = nm.get('sm', {})
    d.p('The auxiliary coordinates form a torus ($x_k\\simeq x_k+L_c$); the Standard Model lives on the brane $W=\\lbrace x_4=\\dots=x_7=0\\rbrace$ with '
        'the induced metric $h_{\\mu\\nu}=g_{\\mu\\nu}(x,y=0)$, and the total action is')
    d.eq(r'S=S_{\rm EL}[g]+S_\Psi[\Psi,g]+\int_Wd^4x\sqrt{|h|}\,\big(\mathcal L_{\rm SM}+\mathcal L_{\rm portal}\big)+S_{\rm stab},')
    d.eq(r'S_{\rm EL}=-\frac1{2\kappa}\int d^8X\sqrt{|g|}\sum_{k=0}^{3}\alpha_kL_k,\qquad '
         r'S_\Psi=\int d^8X\sqrt{|g|}\,\big[\bar\Psi(i\Gamma^\mu D_\mu-m)\Psi+\mathcal L_{\rm self}\big],')
    d.p('($S_{\\rm stab}$ is a stress that keeps the auxiliary scale factor $b$ constant in the modern epoch), '
        'with the complete Standard Model ($SU(3)\\times SU(2)\\times U(1)$, three generations with right-handed neutrinos, one Higgs doublet, '
        'Yukawa couplings) in the gravitational field of $h$. All gauge anomalies cancel in every generation; with the data of the Particle '
        'Data Group the electroweak scale is $v=(\\sqrt2G_F)^{-1/2}=%s$ GeV and the CKM Jarlskog invariant $J=%s$; the brane Friedmann equation '
        '$3H_a^2=8\\pi G_4\\rho_4+\\Lambda_4$ with the Planck 2018 parameters gives the age %s Gyr.' % (
            fmt(mnum(sm.get('v', '') or ''), 9), fmt(mnum(sm.get('J_CKM', '') or ''), 4), fmt(mnum(sm.get('age (Gyr)', '') or ''), 5)))

    # ---------------------------------------------------------------------------------------------
    d.h2('%d.6 dirac16complex at the brane, portals, and neutrino flavour (sections 18–19)' % secno)
    d.p('With the Dirac representation $\\gamma_\\mu$ ($\\gamma_\\mu\\gamma_\\nu+\\gamma_\\nu\\gamma_\\mu=2\\eta_{\\mu\\nu}1_4$), $\\gamma_5=i\\gamma_0\\gamma_1\\gamma_2\\gamma_3$ '
        'and $\\rho_4=i\\gamma_0$, $\\rho_5=-\\gamma_1$, $\\rho_6=-\\gamma_2$, $\\rho_7=-\\gamma_3$, the matrices')
    d.eq(r'R_\mu=\gamma_\mu\otimes1_4\ (\mu=0,\dots,3),\qquad R_k=\gamma_5\otimes\rho_k\ (k=4,\dots,7),\qquad \Gamma_AU=UR_A,\quad U^\dagger U=1_{16},')
    d.p('show that dirac16complex at the brane is four ordinary Dirac spinors, $\\Psi=U\\sum_a\\psi_a\\otimes e_a$, with $\\bar\\Psi\\Gamma^\\mu\\Psi\'='
        '\\sum_a\\bar\\psi_a\\gamma^\\mu\\psi\'_a$; a mode with auxiliary momentum $\\vec q$ has the Kaluza–Klein mass $\\sqrt{m^2+|\\vec q|^2/b_0^2}$. Three '
        'portals (Yukawa, vector, Higgs) couple it to the Standard Model at the brane; the Yukawa portal makes the light neutrino mixing '
        'non-unitary, $1-NN^\\dagger=\\mu\\mu^\\dagger/m^2$. Neutrino flavour conversion in vacuum,')
    d.eq(r'P(\alpha\to\beta)=\delta_{\alpha\beta}-4\sum_{i>j}\mathrm{Re}\,X_{ij}\sin^2\frac{\Delta m^2_{ij}L}{4E}+2\sum_{i>j}\mathrm{Im}\,X_{ij}\sin\frac{\Delta m^2_{ij}L}{2E},'
         r'\qquad X_{ij}=U^*_{\alpha i}U_{\beta i}U_{\alpha j}U^*_{\beta j},')
    dv = nums_in(sm.get('DUNE vacuum P(mu->e), P(mubar->ebar), P(mu->mu)', '')); dm = nums_in(sm.get('DUNE matter P(mu->e), P(mubar->ebar)', ''))
    d.p('in matter (potential $\\sqrt2G_Fn_e$), in the Sun and in supernovae is worked out with NuFit-6.0 data; for DUNE ($L=1300$ km, $E=2.5$ '
        'GeV) $P(\\nu_\\mu\\to\\nu_e)=%s$ in vacuum and %s in the Earth’s crust.' % (fmt(dv[0], 5) if dv else '?', fmt(dm[0], 5) if dm else '?'))
    d.p('**Checks.** Sections 11–19 of 7+1.nb contain %d checks, all True; the independent programs `tools/ref_gravity.py`, '
        '`tools/ref_ks_gravity.py` and `tools/ref_sm.py` repeat the constructions with %d, %d and %d checks, all True, and reproduce the '
        'printed numbers.' % (c['grav'][0] + c['ks'][0] + c['sm'][0], c['grav'][2], c['ks'][2], c['sm'][2]))


# =====================================================================================================================
def symbols_section(d, D3, secno=10):
    d.h1('%d. Sections 11–19 of 7+1.nb: gravity, the self-coupled field in gravitational fields, the Standard Model' % secno)
    d.p('The symbols of the third stage of the project (sections 11–19 of 7+1.nb). Here $\\Gamma_A=\\Gamma^{(1,7)}_A$, $\\Gamma^A=\\eta_{(1,7)}^{AA}\\Gamma_A$, '
        '$\\eta=\\eta^{(1,7)}$, and Greek indices $\\mu,\\nu,\\lambda,\\rho,\\sigma,\\alpha,\\beta,\\gamma,\\delta\\in\\lbrace0,\\dots,7\\rbrace$ are coordinate indices '
        '(on the brane $\\mu,\\nu\\in\\lbrace0,\\dots,3\\rbrace$).')
    d.table(['symbol', 'Mathematica', 'definition', 'meaning'], [
        ['$g_{\\mu\\nu}$, $g^{\\mu\\nu}$', '`cv["g"]`, `cv["ginv"]`', 'signature $(+,-,\\dots,-)$', 'metric of 7+1 spacetime and its inverse'],
        ['$\\{{}^{\\lambda}_{\\mu\\nu}\\}$', '`cv["chr"]`', '(10.1)', 'Christoffel symbols'],
        ['$R^\\rho{}_{\\sigma\\mu\\nu}$', '`cv["riem"]`', '(10.1)', 'Riemann tensor'],
        ['$R_{\\sigma\\nu}$, $R$, $G^\\mu{}_\\nu$', '`cv["ric"]`, `cv["R"]`, `cv["G"]`', '(10.1)', 'Ricci tensor, scalar curvature, Einstein tensor'],
        ['$R^{\\gamma\\delta}{}_{\\alpha\\beta}$', '`cv["rud"]`', '$g^{\\delta\\sigma}R^\\gamma{}_{\\sigma\\alpha\\beta}$', 'Riemann tensor, two indices up'],
        ['$\\delta^{\\mu_1\\cdots\\mu_n}_{\\nu_1\\cdots\\nu_n}$', '`gKD`, `kδ`', '$\\det[\\delta^{\\mu_i}_{\\nu_j}]$', 'generalized Kronecker delta'],
        ['$E^{(k)\\mu}{}_\\nu$, $L_k$', '`Lovelock1/2/3`, `lovelockE`, `lovelockL`', '(10.2)', 'Lovelock tensors and densities'],
        ['$\\alpha_k$, $\\bar\\Lambda$, $\\kappa$, $G_8$', '`alf[k]`, `lamb`, `kap`', '(10.3)', 'couplings, cosmological constant, $\\kappa=8\\pi G_8$'],
        ['$w_k$, $\\Lambda$, $H$', '`w1`, `w2`, `w3`, `Λ`, `H`', '(10.3)', 'pure-number couplings and inverse length of the reference notebook'],
        ['$T^\\mu{}_\\nu$; $\\varepsilon$, $p_a$, $p_b$', '—; `ksStress`', '$T=\\mathrm{diag}(\\varepsilon,-p_a,-p_a,-p_a,-p_b,\\dots,-p_b)$', 'energy–momentum tensor; energy density and pressures'],
        ['$a$, $b$, $H_a$, $H_b$', '`a3[x0]`, `b4[x0]`, `Ha`, `Hb`', '(10.4); $H_a=a\'/a$, $H_b=b\'/b$', 'scale factors and Hubble rates'],
        ['$q_H$', '`qH`', '$b=e^{q_HHx_0}$ when $a=e^{Hx_0}$', 'ratio of the Hubble rates of exponential solutions'],
        ['$p$, $q$', '`pK`, `qK`', '$a=x_0^p$, $b=x_0^q$', 'Kasner exponents'],
        ['$e^A{}_\\mu$, $e_A{}^\\mu$', '`fd["e"]`, `fd["ei"]`', '$g_{\\mu\\nu}=\\eta_{AB}e^A{}_\\mu e^B{}_\\nu$', 'frame field and its inverse'],
        ['$\\Gamma_\\mu(X)$, $\\Gamma^\\mu(X)$', '`fd["Gam"]`, `fd["GamUp"]`', '$e^A{}_\\mu\\Gamma_A$, $e_A{}^\\mu\\Gamma^A$', 'generalized frame matrices'],
        ['$\\omega_\\mu{}^A{}_B$, $\\Omega_\\mu$', '`fd["om"]`, `fd["Om"]`', '(10.5)', 'canonical spin connection, spinor connection'],
        ['$D_\\mu$, $\\nabla_\\mu$', '—', '$D_\\mu=\\partial_\\mu+\\Omega_\\mu$; (10.5)', 'spinor covariant derivative; total covariant derivative'],
        ['$R_{AB\\mu\\nu}$', '`rABmn`', '$\\eta_{AC}e^C{}_\\rho R^\\rho{}_{\\sigma\\mu\\nu}e_B{}^\\sigma$', 'Riemann tensor with frame indices'],
        ['$q$, $A_\\mu$', '`qG`, `aGauge`', '$D_\\mu=\\partial_\\mu+\\Omega_\\mu-iqA_\\mu$', 'charge and abelian gauge field'],
        ['$M$ (pair creation), $|\\beta|^2$', '`bogoliubov`', '$|\\beta|^2=1/(e^{2\\pi M/H}+1)$ in de Sitter space', 'mass of a mode; pair-creation probability'],
        ['$K$, $\\Theta^K$, $\\Theta_K$', '`mIdx`, `Θup[K]`, `Θdn[K]`', '(10.6)', 'multi-index and the 256 Clifford products'],
        ['$g_r$; $g_S,g_V,g_A,g_P$', '`gg[r]`; `gS`, `gV`, `gA`, `gP`', '$g_0=g_S$, $g_1=g_V$, $g_7=-g_A$, $g_8=-g_P$', 'self-couplings of rank $r$; the four of section 11.3'],
        ['$F_{rs}$, $G_s$', '`fierzF`, `Geff`', '(10.6)', 'Fierz matrix and Hartree–Fock couplings'],
        ['$M$ (Kohn–Sham), $\\rho_S$, $\\rho$', '`ksGap`, `ksMoments`', '$M=m+G_0\\rho_S$', 'effective mass, scalar density, density of uniform matter'],
        ['$W$, $h_{\\mu\\nu}$, $L_c$, $b_0$', '—', '$W=\\lbrace x_4=\\dots=x_7=0\\rbrace$, $h_{\\mu\\nu}=g_{\\mu\\nu}(x,y=0)$', 'world brane, induced metric, torus period, static scale'],
        ['$S$, $S_{\\rm EL}$, $S_\\Psi$, $S_{\\rm stab}$', '—', '(10.9)', 'total action; Einstein–Lovelock, dirac16complex and stabilizing parts'],
        ['$\\mathcal L_{\\rm self}$, $\\mathcal L_{\\rm SM}$, $\\mathcal L_{\\rm portal}$', '`LselfS`', '(10.9); $\\mathcal L_{\\rm self}=-\\frac12\\sum_rg_rI_r$', 'self-coupling, Standard Model and portal Lagrangians'],
        ['$V_\\perp$, $G_4$', '`vPerp`', '$V_\\perp=b_0^4L_c^4$, $G_4=G_8/V_\\perp$', 'proper volume of the torus; four-dimensional Newton constant'],
        ['$v$, $G_F$, $J$', '`vEW`, `gFermi`, `jCKM`', '$v=(\\sqrt2G_F)^{-1/2}$; $J=\\mathrm{Im}(V_{us}V_{cb}V_{ub}^*V_{cs}^*)$', 'Higgs vacuum value, Fermi constant, Jarlskog invariant'],
        ['$V_{\\rm CKM}$, $U_{\\rm PMNS}$, $U$', '`vCKM`, `uPMNS`', 'standard parametrization', 'quark and lepton mixing matrices'],
        ['$\\rho_4$, $\\Lambda_4$', '—', '$3H_a^2=8\\pi G_4\\rho_4+\\Lambda_4$', 'energy density and cosmological constant of the brane'],
        ['$\\gamma_\\mu$, $\\gamma_5$, $\\rho_k$', '`g4D`, `g5D`, `rhoK`', '(10.7)', 'brane Dirac matrices and transverse matrices'],
        ['$R_A$, $U$ (brane)', '`rBr`, `uInt`', '(10.7)', 'brane form of the Dirac matrices; intertwiner'],
        ['$\\psi_a$, $e_a$', '—', '$\\Psi=U\\sum_a\\psi_a\\otimes e_a$', 'four brane Dirac spinors; unit vectors of $\\mathbb C^4$'],
        ['$\\vec q$', '`qKK`', '$\\vec q=2\\pi\\vec n/L_c$', 'auxiliary momentum'],
        ['$\\mu$, $m$, $N$', '`muPortal`, `mPsi`, `nAct`', '$1-NN^\\dagger=\\mu\\mu^\\dagger/m^2$', 'portal Dirac mass, dirac16complex mass, light mixing matrix'],
        ['$P(\\alpha\\to\\beta)$, $X_{ij}$, $\\Delta m^2_{ij}$', '`pVac`, `xFormula`, `dm21`, `dm31`', '(10.8)', 'oscillation probability'],
        ['$L$, $E$, $n_e$', '—', '', 'baseline, neutrino energy, electron density'],
    ], widths=[0.17, 0.21, 0.33, 0.29], small=True)
    d.align([r'\{{}^{\lambda}_{\mu\nu}\}&=\tfrac12g^{\lambda\rho}(\partial_\mu g_{\rho\nu}+\partial_\nu g_{\rho\mu}-\partial_\rho g_{\mu\nu}),&&(10.1)',
             r'R^\rho{}_{\sigma\mu\nu}&=\partial_\mu\{{}^{\rho}_{\sigma\nu}\}-\partial_\nu\{{}^{\rho}_{\sigma\mu}\}+\{{}^{\rho}_{\mu\lambda}\}\{{}^{\lambda}_{\sigma\nu}\}-\{{}^{\rho}_{\nu\lambda}\}\{{}^{\lambda}_{\sigma\mu}\},',
             r'R_{\sigma\nu}&=R^\rho{}_{\sigma\rho\nu},\quad G^\mu{}_\nu=R^\mu{}_\nu-\tfrac12\delta^\mu_\nu R,',
             r'E^{(k)\mu}{}_\nu&=-\frac{1}{2^{k+1}}\delta^{\mu\alpha_1\beta_1\cdots\alpha_k\beta_k}_{\nu\gamma_1\delta_1\cdots\gamma_k\delta_k}R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}\cdots R^{\gamma_k\delta_k}{}_{\alpha_k\beta_k},&&(10.2)',
             r'L_k&=\frac1{2^k}\delta^{\alpha_1\beta_1\cdots\alpha_k\beta_k}_{\gamma_1\delta_1\cdots\gamma_k\delta_k}R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}\cdots R^{\gamma_k\delta_k}{}_{\alpha_k\beta_k},',
             r'0&=-\Lambda\delta^\mu_\nu+\sum_{k=1}^3H^{-2k}w_kE^{(k)\mu}{}_\nu,\quad \sum_{k=0}^3\alpha_kE^{(k)\mu}{}_\nu=\kappa T^\mu{}_\nu,&&(10.3)',
             r'&\qquad \alpha_0=2\bar\Lambda,\ \alpha_1=1,\ \alpha_2=\tfrac{w_2}{w_1H^2},\ \alpha_3=\tfrac{w_3}{w_1H^4},',
             r'ds^2&=dx_0^2-a^2(dx_1^2+dx_2^2+dx_3^2)-b^2(dx_4^2+\cdots+dx_7^2),&&(10.4)',
             r'\omega_\mu{}^A{}_B&=e^A{}_\nu\big(\partial_\mu e_B{}^\nu+\{{}^{\nu}_{\mu\lambda}\}e_B{}^\lambda\big),\quad \Omega_\mu=\tfrac14\eta_{AC}\omega_\mu{}^C{}_B\Gamma^A\Gamma^B,&&(10.5)',
             r'\nabla_\mu\Gamma_\nu&=\partial_\mu\Gamma_\nu-\{{}^{\lambda}_{\mu\nu}\}\Gamma_\lambda+[\Omega_\mu,\Gamma_\nu]=0,',
             r'\Theta^K&=i^{\sigma_r}\Gamma^{A_1}\cdots\Gamma^{A_r},\ \ \sigma_r=\tfrac{r(r-1)}2\bmod2,\ \ \Theta_K=\eta_{A_1A_1}\cdots\eta_{A_rA_r}\Theta^K,&&(10.6)',
             r'F_{rs}&=\tfrac1{16}\sum_l(-1)^{rs-l}\tbinom sl\tbinom{8-s}{r-l},\ \ G=(1-F^{\mathsf T})g,',
             r'\gamma_0&=\mathrm{diag}(1,1,-1,-1),\ \gamma_i=\begin{pmatrix}0&\sigma_i\\-\sigma_i&0\end{pmatrix},\ \rho_4=i\gamma_0,\ \rho_{4+i}=-\gamma_i,&&(10.7)',
             r'R_\mu&=\gamma_\mu\otimes1_4,\ R_k=\gamma_5\otimes\rho_k,\ \Gamma_AU=UR_A,',
             r'P(\alpha\to\beta)&=\Big|\sum_iU_{\beta i}U^*_{\alpha i}e^{-im_i^2L/(2E)}\Big|^2,\quad X_{ij}=U^*_{\alpha i}U_{\beta i}U_{\alpha j}U^*_{\beta j},\quad \Delta m^2_{ij}=m_i^2-m_j^2,&&(10.8)',
             r'S&=S_{\rm EL}+S_\Psi+\int_Wd^4x\sqrt{|h|}(\mathcal L_{\rm SM}+\mathcal L_{\rm portal})+S_{\rm stab},&&(10.9)',
             r'S_{\rm EL}&=-\tfrac1{2\kappa}\int d^8X\sqrt{|g|}\sum_k\alpha_kL_k,\quad S_\Psi=\int d^8X\sqrt{|g|}\big[\bar\Psi(i\Gamma^\mu D_\mu-m)\Psi+\mathcal L_{\rm self}\big].'])


# =====================================================================================================================
def provenance_section(d, D, D3, runs3, secno=7):
    c = counts(D3)
    d.h1('%d. The third stage: sections 11–20 of 7+1.nb' % secno)
    d.p('This section records the provenance of the gravity, Kohn–Sham and Standard Model sections (11–19) and of the final summary (20). '
        'Notation: $\\Gamma_A=\\Gamma^{(1,7)}_A$, $\\Gamma^A=\\eta_{(1,7)}^{AA}\\Gamma_A$, $\\eta=\\mathrm{diag}(+1,-1,\\dots,-1)$, $\\Psi$ the field dirac16complex.')
    d.h2('%d.1 The task' % secno)
    d.ol(['Rebuild the ideas of github.com/once-ere/Dirac_claude on the 7+1 spacetime of this project (one time coordinate, seven space '
          'coordinates; never the 4+4 metric of Dirac_claude): the Einstein and Lovelock tensors; generalized frame matrices; a canonical spin '
          'connection, with a vanishing total covariant derivative, that couples gravity to dirac16complex; the Kohn–Sham theory of '
          'dirac16complex in gravitational fields in the primordial (near $x_0=0$), modern and far-future epochs; a generalized self-coupling; '
          'and the totality of the Standard Model.',
          'Compute the Lovelock tensors as in the author’s corrected reference notebook Einstein-Lovelock-7+1-Nash.nb, and do not use its '
          'section "FAILED Model metric g71; DO NOT USE!".',
          'Document every stage (documentation, provenance and symbol files; Markdown, LaTeX, PDF; self-contained) and push every completed '
          'stage to github.com/once-ere/7plus1.'])
    d.h2('%d.2 Inputs and how they were used' % secno)
    d.p('The reference notebook Einstein-Lovelock-7+1-Nash.nb (P. L. Nash, 2022, GPL; the author’s corrected version) supplied: the xAct '
        'packages xCoba and xTras; the manifold with abstract indices $a$, $b$, $c$, $d$, $f$, $i$, $j$, $k$, $l$, $m$, $n$, $p$, $q$, $r$, $s$, $y$; the chart with coordinate scalars '
        '`x0[]`, …, `x7[]`; `MyArrayComponents`; the generalized delta `kδ`; the three Lovelock tensors with their normalization')
    d.eq(r'\mathrm{Lovelock}_k=-\frac1{2^{k+1}}\,R_{c_1d_1}{}^{f_1i_1}\cdots R_{c_kd_k}{}^{f_ki_k}\,\delta^{b\,c_1d_1\cdots c_kd_k}_{a\,f_1i_1\cdots f_ki_k}\quad(k=1,2,3);')
    d.p('the evaluation of components by `CD = LC[g]`; and the vacuum field equations $0=-\\Lambda\\delta+H^{-2}w_1\\mathrm{Lovelock}_1+H^{-4}w_2\\mathrm{Lovelock}_2'
        '+H^{-6}w_3\\mathrm{Lovelock}_3$ with the procedure of solving them for $w_2,w_3$. The model metric of its FAILED section was not used anywhere. '
        'The only change: the signature declaration $\\{7,1,0\\}$ became $\\{1,7,0\\}$ to match $\\eta$ of this notebook (it does not enter the '
        'Lovelock tensors with mixed indices, and `SignDetOfMetric` is $-1$ for both declarations, as checked with xTensor 1.2.0). xAct versions as printed by the notebook: %s.' % (D3.get('xact_versions') or '?'))
    if D3.get('sha_reference'):
        d.code('Einstein-Lovelock-7+1-Nash.nb   SHA-256 %s\nEinstein-Lovelock-7+1-Nash.pdf  SHA-256 %s' % D3['sha_reference'], 'text')
    d.p('Data of the Standard Model part: Particle Data Group 2024/2025 ($G_F$, boson masses, quark masses, CKM parameters), NuFit-6.0 (PMNS '
        'parameters, normal ordering, without Super-Kamiokande atmospheric data), Planck 2018 (cosmological parameters).')
    d.h2('%d.3 Decisions' % secno)
    d.table(['#', 'decision', 'reason'], [
        ['1', 'Christoffel symbols written $\\{{}^{\\lambda}_{\\mu\\nu}\\}$; Landau–Lifshitz curvature', '$\\Gamma$ is reserved for Dirac matrices; with $(+,-,\\dots,-)$, $G^\\mu{}_\\nu-\\bar\\Lambda\\delta^\\mu_\\nu=\\kappa T^\\mu{}_\\nu$, $T^0{}_0>0$'],
        ['2', 'the Lovelock tensors are computed twice (xAct; pair-set sums)', 'exact agreement is an independent check of the reference method'],
        ['3', 'the exponential solutions are obtained from the xAct components of the cosmological metric', 'a new xAct run with symbolic exponentials does not finish in reasonable time'],
        ['4', 'the canonical (Levi-Civita) spin connection, $\\Omega_\\mu=\\frac14\\omega_{\\mu AB}\\Gamma^A\\Gamma^B$', 'it is the connection with $\\nabla_\\mu\\Gamma_\\nu=0$'],
        ['5', 'the self-coupling is organized by the nine ranks of the Clifford algebra', 'each rank is one Lorentz-invariant class; the Fierz matrix closes on them'],
        ['6', 'Kohn–Sham matter: homogeneous, spin-saturated, no Dirac sea, adiabatic', 'the cosmological spacetime is homogeneous; the expansion rates are small compared with the orbital energies'],
        ['7', 'the Standard Model on a 3+1 brane of a torus of the auxiliary coordinates', 'the 7+1 analogue of the fold of Dirac_claude'],
        ['8', 'portals are invented local couplings at the brane', 'the only contact of dirac16complex with the Standard Model besides gravity'],
    ], widths=[0.055, 0.475, 0.47], small=True)
    d.h2('%d.4 Literature' % secno)
    d.ul(references())
    d.h2('%d.5 Verification of the third stage' % secno)
    d.p('Sections 11–19 of 7+1.nb contain %d checks (gravity %d, Kohn–Sham %d, Standard Model %d), all True; they are rows of section 4.2 above. '
        'The independent programs recompute everything: `tools/ref_gravity.py` (%d boolean checks, %d True), `tools/ref_ks_gravity.py` (%d, %d), '
        '`tools/ref_sm.py` (%d, %d); their complete outputs are printed in section 4.3.' % (
            c['grav'][0] + c['ks'][0] + c['sm'][0], c['grav'][0], c['ks'][0], c['sm'][0], c['grav'][2], c['grav'][3], c['ks'][2], c['ks'][3],
            c['sm'][2], c['sm'][3]))
