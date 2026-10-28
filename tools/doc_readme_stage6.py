"""The parts of README, README_SYMBOLS and README_PROVENANCE that describe sections 25 and 26 of 7+1.nb (the auxiliary torus and the
world brane, defined and explained with worked examples; the radion -- the size of the extra space -- on the 4-torus and on the 4-sphere,
with a check of the author's note on the radion).  Every equation used is written in the document that uses it; every number is parsed from
the evaluated notebook or from the outputs of tools/ref_torus_brane.py and tools/ref_radion.py (gathered by doc_stage3_common.gather3)."""
from doc_stage3_common import fmt
import doc_torus_brane_common as TB
import doc_radion_common as RD


def present(D3):
    return bool(D3 and (D3.get('checks', {}).get('tb') or D3.get('checks', {}).get('rd')))


def counts(D3):
    return TB.counts(D3), RD.counts(D3)


def abstract_text(D3):
    return ('; (ix) in sections 25 and 26 of the notebook, the auxiliary torus and the world brane, defined precisely and explained with worked '
            'examples, and the radion — the size of the four extra dimensions — on the 4-torus and on the 4-sphere, with a check of the author’s '
            'note on the radion (sections 15 and 16)')


def results_paragraph(d, D3):
    (a, b, c_, e), (a6, b6, c6, e6) = counts(D3)
    tod = TB.nbv(D3, 'today')
    d.p('**The auxiliary torus, the world brane and the radion** (sections 25–26 of the notebook; sections 15–16 below). The auxiliary space is the '
        'flat torus $T^4=\\mathbb R^4/(L_c\\mathbb Z)^4$ of proper volume $V_\\perp=b^4L_c^4$; today ($M_8=10$ TeV) $\\ell_\\perp=V_\\perp^{1/4}=%s$ m, and its '
        'Kaluza–Klein levels start at $2\\pi/\\ell_\\perp=%s$ MeV. The world brane $W=\\{x_4=\\dots=x_7=0\\}$ carries the Standard Model, sees the induced '
        'metric $h_{\\mu\\nu}=g_{\\mu\\nu}(x,0)$, is totally geodesic for the metrics of the model (they do not depend on $\\vec y$), and couples to every '
        'Kaluza–Klein level with the same strength $1/\\sqrt{V_\\perp}$. '
        'The radion $\\beta=\\ln(b/b_0)$, canonically $\\phi=\\sqrt{12}M_{\\rm Pl}\\beta$, couples to brane matter with $1/(\\sqrt3M_{\\rm Pl})$ on the torus and on '
        'a 4-sphere alike; a sphere needs a stiffer stabilizer ($m_r^2=\\kappa K_U/12-2/r_0^2$) and has no Dirac zero modes. The author’s note on '
        'the radion is right in its ideas; the checks correct a factor 2 in its canonical normalization, a missing Einstein-frame factor, the '
        'sign of a junction condition, its warped normalization and mass, and five-dimensional equations used in eight dimensions. Sections 25 '
        'and 26 contain %d and %d checks of their own (%d and %d True).' % (TB.g(tod, 2, 4), TB.g(tod, 4, 4), a, a6, b, b6))


def notation_rows():
    return [
        ['$T^4$, $L_c$, $b$, $V_\\perp$, $\\ell_\\perp$', '`tbLc`, `b4[x0]`, `tbV`, `tbLp`', 'the auxiliary torus, its period, size, proper volume and size (section 15)'],
        ['$W$, $h_{\\mu\\nu}$', '`tbOnW`, `tbInduced`', 'the world brane and its induced metric (section 15)'],
        ['$\\beta$, $\\phi$, $V(\\beta)$', '`rdBt`, `rdV`', 'the radion, the canonical radion, its potential (section 16)'],
        ['$S^4$, $r$, $d_l$, $Z_l$', '`rdR0`, `rdDl`, `rdZon`', 'the 4-sphere of the comparison, its radius, multiplicities, zonal harmonics (section 16)'],
    ]


def verification_text(D3):
    (a, b, c_, e), (a6, b6, c6, e6) = counts(D3)
    w5, w6 = TB.worst(D3), RD.worst(D3)
    return ('For sections 25 and 26 the independent programs `tools/ref_torus_brane.py` (%d checks, %d True) and `tools/ref_radion.py` (%d checks, '
            '%d True; sympy, numpy, mpmath, with their own lattice sieves, curvature routines and 60-digit arithmetic) repeat every example; the %d '
            'and %d groups of numbers printed by the notebook and by them agree to a largest relative difference of %s and %s, and '
            '`verify_7p1_nb.py` confirms that every cell and output of the earlier sections is unchanged.'
            % (c_, e, c6, e6, len(TB.compare(D3)), len(RD.compare(D3)), fmt(w5, 2) if w5 is not None else '?', fmt(w6, 2) if w6 is not None else '?'))


def notebook_items(n_last):
    return ['**25. The auxiliary torus and the world brane: definitions and worked examples** — the torus as a quotient space, its symmetries, '
            'Fourier modes and Kaluza–Klein masses, Jacobi’s count of levels and Weyl’s law, the size today, images and Kaluza–Klein modes in '
            'Newton’s law; the brane, its induced metric and extrinsic curvature, brane sources, dirac16complex as four brane spinors (section 15).',
            '**26. The radion: the size of the extra space, on the 4-torus and on the 4-sphere** — the Einstein frame and the canonical radion, the '
            'potentials and masses, the coupling to the brane and the fifth force, levels and Dirac spectra on $T^4$ and $S^4$, the warped model of '
            'Randall and Sundrum, the mixing with the Higgs boson; a check of the author’s note (section 16).',
            '**%d. Summary of all checks** of sections 1–%d.' % (n_last + 1, n_last)]


def file_rows():
    return [
        ['`tools/nb_part_torus_brane.py`, `tools/nb_part_radion.py`', 'generator parts of sections 25 and 26 of 7+1.nb (executed by `make_7p1_nb.py`)'],
        ['`tools/ref_torus_brane.py`, `tools/ref_radion.py`', 'independent Python checks of sections 25 and 26'],
        ['`AUXILIARY_TORUS*`, `WORLD_BRANE*`, `RADION*` (`.md`, `.tex`, `.pdf`)', 'the documents on the auxiliary torus, the world brane and the radion, '
         'each with its symbol dictionary and its provenance record'],
    ]


def references():
    return ['Th. Kaluza (1921); O. Klein, Z. Phys. 37 (1926) 895; N. Arkani-Hamed, S. Dimopoulos, G. Dvali, Phys. Lett. B 429 (1998) 263 — compact '
            'and large extra dimensions (sections 15–16); C. G. J. Jacobi (1829); H. Weyl, Math. Ann. 71 (1912) 441 — counting the levels.',
            'L. Randall, R. Sundrum, Phys. Rev. Lett. 83 (1999) 3370; W. D. Goldberger, M. B. Wise, Phys. Rev. Lett. 83 (1999) 4922 and Phys. Lett. '
            'B 475 (2000) 275; G. F. Giudice, R. Rattazzi, J. D. Wells, Nucl. Phys. B 595 (2001) 250 — the warped radion, its stabilization, '
            'Higgs–radion mixing (section 16); R. Camporesi, A. Higuchi, J. Geom. Phys. 20 (1996) 1 — Dirac spectra on spheres; J. G. Lee et '
            'al., Phys. Rev. Lett. 124 (2020) 101101 — the Eöt-Wash bound.']


# =====================================================================================================================
def readme_section(d, D3, secno=15):
    """README sections secno (section 25 of the notebook) and secno + 1 (section 26), self-contained"""
    (a, b, c_, e), (a6, b6, c6, e6) = counts(D3)
    g = TB.g
    tod = TB.nbv(D3, 'today'); nmu = TB.nbv(D3, 'nmu')
    d.h1('%d. Section 25: the auxiliary torus and the world brane' % secno)
    d.p('Section 25 of 7+1.nb was appended on 28 October 2026 (the summary of all checks became section 26, and is now section %d). It defines '
        'the two places of the model — the compact space of the auxiliary coordinates and the submanifold where the Standard Model lives — and '
        'checks a worked example for every statement. Coordinates: $x_0$ time (the octonion unit), $x_1,x_2,x_3$ ordinary space, $\\vec y=(x_4,\\dots,x_7)$ '
        'the auxiliary space; $\\eta=\\mathrm{diag}(+1,-1,\\dots,-1)$ is simply the squares of the octonion basis, and xAct declares the spacetime with '
        '$\\{7,1,0\\}$. The section contains %d checks (%d True).' % (D3['n_last'] + 1, a, b))
    d.h2('%d.1 The auxiliary torus' % secno)
    d.eq(r'x_k\simeq x_k+L_cn_k\ (n_k\in\mathbb Z):\quad T^4=\mathbb R^4/(L_c\mathbb Z)^4,\qquad b^8=\det(-g_{kl}),\qquad V_\perp=b^4L_c^4,\qquad \ell_\perp=V_\perp^{1/4}=bL_c .')
    d.p('A field on the torus is a Fourier series in $\\phi_{\\vec n}=e^{2\\pi i\\vec n\\cdot\\vec y/L_c}/\\sqrt{V_\\perp}$ ($\\vec n\\in\\mathbb Z^4$); seen from ordinary spacetime it '
        'is a Kaluza–Klein tower of masses $\\sqrt{m^2+\\mu_{\\vec n}^2}$,')
    d.eq(r'\mu_{\vec n}=\frac{2\pi|\vec n|}{\ell_\perp},\qquad r_4(k)=8\sigma(k)-32\sigma(k/4)\ \text{states with }|\vec n|^2=k,\qquad N(\mu)\simeq\frac{V_\perp\mu^4}{32\pi^2},')
    d.p('($\\sigma$ the divisor sum; Weyl’s law for the number of levels below $\\mu$). The volume fixes Newton’s constant, $8\\pi G_4=\\kappa/V_\\perp$; today '
        '$V_\\perp=M_{\\rm Pl}^2/M_8^6=%s$ GeV$^{-4}$, $\\ell_\\perp=%s$ m, the first level is at %s MeV and there are $%s$ levels below 5 GeV. By Poisson '
        'summation the images of a mass on the torus equal its Kaluza–Klein modes,' % (g(tod, 0, 4), g(tod, 2, 4), g(tod, 4, 4), g(nmu, 0, 3)))
    d.eq(r'\sum_{\vec n\in\mathbb Z^4}\frac1{(r^2+|\vec n|^2)^{5/2}}=\frac{4\pi^2}{3r}\sum_{\vec m\in\mathbb Z^4}e^{-2\pi r|\vec m|}\quad(\text{lengths in units of }\ell_\perp),')
    d.p('so gravity changes from the eight-dimensional $1/r^5$ law to $\\frac53G_4M/r$ — the graviton plus the massless radion, $1+\\frac23$ — and to '
        'Newton’s law once the stabilizer gives the radion its mass.')
    d.h2('%d.2 The world brane' % secno)
    d.eq(r'W=\{X:\ x_4=x_5=x_6=x_7=0\},\qquad h_{\mu\nu}(x)=g_{\mu\nu}(x,\vec y=0),\qquad \int_Wd^4x\sqrt{|h|}\,\mathcal L=\int d^8X\sqrt{|g|}\,\frac{\delta^4(\vec y)}{b^4}\,\mathcal L .')
    d.p('For the metrics of the model, which do not depend on $\\vec y$, the brane is totally geodesic (its extrinsic curvature $\\frac1{2b}\\partial_kg_{\\mu\\nu}$ '
        'vanishes, so a free particle moving along it stays on it; a $\\vec y$-dependent metric can push particles off it); a tension $f_W$ has '
        '$T^W_{\\mu\\nu}=f_Wh_{\\mu\\nu}$, $w=-1$; smeared over the torus a brane source has $\\varepsilon=\\rho_4/V_\\perp$ and no pressure along the torus. With '
        'the intertwiner $U$ ($\\Gamma_AU=UR_A$, $R_\\mu=\\gamma_\\mu\\otimes1_4$, $R_k=\\gamma_5\\otimes\\rho_k$) dirac16complex at the brane is four brane Dirac spinors, '
        '$\\psi=U\\sum_a\\psi_a\\otimes e_a$; the brane breaks the translations of the torus, so the portals couple the Standard Model to every '
        'Kaluza–Klein level with the same strength $1/\\sqrt{V_\\perp}$.')

    # -----------------------------------------------------------------------------------------------
    hr = RD.nbv(D3, 'hr71')
    d.h1('%d. Section 26: the radion, on the 4-torus and on the 4-sphere' % (secno + 1))
    d.p('Section 26 of 7+1.nb was appended on 28 October 2026 (the summary of all checks became section %d). It defines the radion precisely, '
        'works out every consequence for two extra spaces — the 4-torus of the model and, for comparison, a 4-sphere — and checks a note of the '
        'author on the radion statement by statement. The section contains %d checks (%d True).' % (D3['n_last'] + 1, a6, b6))
    d.h2('%d.1 The definition' % (secno + 1))
    d.eq(r'ds^2=e^{-4\beta(x)}g_{E\mu\nu}dx^\mu dx^\nu-e^{2\beta(x)}\gamma_{ij}dy^idy^j,\qquad \gamma=b_0^2\delta\ (T^4),\quad \gamma=r_0^2\hat\gamma\ (S^4),')
    d.p('with $\\hat\\gamma=\\mathrm{diag}(1,\\sin^2x_4,\\sin^2x_4\\sin^2x_5,\\sin^2x_4\\sin^2x_5\\sin^2x_6)$ the unit sphere. The factor $e^{-4\\beta}$ keeps the coefficient of the '
        'four-dimensional curvature constant (Einstein frame, $2\\alpha+n\\beta=0$). The reduction, checked for arbitrary homogeneous $N$, $a$, $\\beta$, gives')
    d.eq(r'\sqrt{|g_E|}\Big[-\frac{M_{\rm Pl}^2}2R_E+6M_{\rm Pl}^2(\partial\beta)^2-V(\beta)\Big],\qquad M_{\rm Pl}^2=\frac{V_0}\kappa,\qquad \phi=\sqrt{12}M_{\rm Pl}\beta,')
    d.eq(r"V(\beta)=M_{\rm Pl}^2\Big[e^{-4\beta}\big(\Lambda+\kappa U(\beta)\big)-\frac{6\sigma}{r_0^2}e^{-6\beta}\Big]\qquad(\sigma=0:\ T^4,\ \sigma=1:\ S^4),\qquad m_r^2=\frac{V''(0)}{12M_{\rm Pl}^2} .")
    d.h2('%d.2 Torus and sphere' % (secno + 1))
    d.table(['', '$T^4$', '$S^4$'], [
        ['without stabilizer', 'modulus', 'collapses (curvature potential)'],
        ['flat vacuum with $U=U_0+U_1\\beta+\\frac{K_U}2\\beta^2$', '$\\Lambda+\\kappa U_0=0$, $U_1=0$', '$(\\Lambda+\\kappa U_0)r_0^2=%s$, $\\kappa U_1r_0^2=%s$'
         % (RD.ival(D3, 'sphU', 0), RD.ival(D3, 'sphU', 1))],
        ['radion mass', '$\\kappa K_U/12$', '$\\kappa K_U/12%s/r_0^2$' % RD.ival(D3, 'sphU', 2)],
        ['held by $\\Lambda$ alone', '—', '$dS_4\\times S^4$, $\\Lambda r_0^2=%s$, unstable ($m^2=%s H^2$)' % (RD.ival(D3, 'sphL', 0), RD.ival(D3, 'sphL', 2))],
        ['levels', '$2\\pi|\\vec n|/\\ell_\\perp$, $r_4(k)$ states', '$\\sqrt{l(l+3)}/r$, $d_l=(2l+3)(l+1)(l+2)/6$ states'],
        ['seen by a brane', '1 state per level, $\\sqrt{r_4(k)/V}$', '1 zonal state per level, $\\sqrt{d_l/V}$'],
        ['Dirac zero modes', '4', 'none (lowest $|\\lambda|=2/r$)'],
        ['today ($V=V_\\perp$)', '$\\ell_\\perp=%s$ m, first level %s MeV' % (g(tod, 2, 4), RD.val(D3, 'today', 2, 4)),
         '$r=%s$ m, first level %s MeV' % (RD.val(D3, 'today', 1, 4), RD.val(D3, 'today', 3, 4))],
    ], widths=[0.3, 0.32, 0.38], small=True)
    d.h2('%d.3 The radion and the brane' % (secno + 1))
    d.p('On both spaces the brane sees $h=e^{-4\\beta}g_E$, so the radion couples as $\\phi\\,T^\\mu{}_\\mu/(\\sqrt3M_{\\rm Pl})$; a massless radion would add a force of '
        'strength $\\alpha=\\frac23$, and the Eöt-Wash bound requires $m_r\\gtrsim\\hbar c/38.6\\,\\mu{\\rm m}=%s$ eV (the radion of the model has 10 GeV). The '
        'Higgs boson on the brane mixes with the radion through $\\xi R(h)H^\\dagger H$ with $\\gamma=v/(\\sqrt3M_{\\rm Pl})=%s$: for $\\xi=1$ the mixing angle '
        'is $%s$ and the Higgs rates change by $%s$. In the warped model of Randall and Sundrum the radion couples with $1/\\langle\\mu\\rangle$ instead, '
        '$%s$ times more strongly for $\\langle\\mu\\rangle=3$ TeV, and the brane tensions are $\\pm6k/\\kappa$ (UV positive).'
        % (RD.val(D3, 'ew', 0, 4), RD.val(D3, 'ew', 2, 4), RD.val(D3, 'hr71', 1, 3), RD.val(D3, 'hr71', 3, 2), RD.val(D3, 'strength', 0, 3)))
    d.h2('%d.4 The author’s note' % (secno + 1))
    d.p('The note’s ideas — the Kaluza–Klein ansatz, the Einstein frame, the physical volume, the modulus and the fifth force, stabilization, the '
        'warped alternative, Higgs–radion mixing — are confirmed. The checks correct: (1) the reduced action with its $\\beta_n=\\sqrt{2/(n(n+2))}$ is '
        '$\\frac{V_n}{2\\kappa^2}[R-(\\partial\\varphi)^2]$, not $[R-\\frac12(\\partial\\varphi)^2]$ (the kinetic coefficient is $\\frac{n(n+2)}2$; the form with $\\frac12$ needs '
        '$\\beta_n=1/\\sqrt{n(n+2)}$); (2) a bulk potential carries the Einstein-frame factor $e^{-4\\beta}$ (gradient energy $e^{-6\\beta}$); (3) the Israel '
        'junction condition needs $-\\kappa^2S_{\\mu\\nu}$, so the UV tension is $+6k/\\kappa^2$; (4) $f_\\mu=\\sqrt{6/(\\kappa k)}$ and $m_r^2=\\frac43\\kappa k^2\\epsilon^2v_B^2(\\langle\\mu\\rangle/f)^2$ '
        '(the note: $\\sqrt{12/(\\kappa k)}$, $\\frac23$); (5) the five-dimensional normalizations $M_5^3\\int R_5$ and $f=\\sqrt{24M_5^3/k}$ belong to different '
        'conventions; (6) the warped equations offered for eight dimensions are five-dimensional ($G^z{}_z=4A\'\'-10A\'^2$ forbids the solution); '
        '(7) the Higgs–radion kinetic terms lack the $6\\xi\\gamma^2$ that its normalization $Z^2=1+6\\xi\\gamma^2(1-6\\xi)$ contains.')


# =====================================================================================================================
def symbols_section(d, D3, secno=13):
    d.h1('%d. Sections 25–26 of 7+1.nb: the auxiliary torus, the world brane and the radion' % secno)
    d.p('The symbols of sections 25 and 26 (the Mathematica names begin with `tb` and `rd`).')
    d.table(['symbol', 'Mathematica', 'definition', 'meaning'], [
        ['$T^4$, $L_c$', '`tbLc`', '(%d.1)' % secno, 'the auxiliary torus; the coordinate period'],
        ['$b$, $V_\\perp$, $\\ell_\\perp$', '`b4[x0]`, `tbV`, `rIIVp`, `tbLp`', '(%d.1)' % secno, 'size, proper volume and proper size of the torus'],
        ['$\\phi_{\\vec n}$, $\\mu_{\\vec n}$, $r_4(k)$', '`tbPhi`, `rIIr4`', '(%d.2)' % secno, 'Fourier modes, Kaluza–Klein masses, multiplicities'],
        ['$W$, $h_{\\mu\\nu}$', '`tbOnW`, `tbInduced`', '(%d.3)' % secno, 'the world brane; its induced metric'],
        ['$f_W$', '`tbfW`', '(%d.3)' % secno, 'brane tension'],
        ['$U$, $\\psi_a$', '`uInt`, `tbU8`', '(%d.3)' % secno, 'intertwiner; the four brane spinors'],
        ['$\\beta$, $\\phi$, $g_E$', '`rdBt`, `rdBe`', '(%d.4)' % secno, 'the radion, the canonical radion, the Einstein-frame metric'],
        ['$V(\\beta)$, $\\sigma$, $U(\\beta)$', '`rdV`, `rdUq`', '(%d.4)' % secno, 'radion potential; curvature switch; stabilizing energy'],
        ['$S^4$, $r$, $\\hat\\gamma$', '`rdR0`, `rdGam4`', '(%d.4)' % secno, 'the 4-sphere of the comparison'],
        ['$d_l$, $Z_l$', '`rdDl`, `rdZon`', '(%d.5)' % secno, 'multiplicities; zonal harmonics'],
        ['$\\gamma$, $\\xi$, $Z$, $\\theta$', '`rdGam`, `rdXi`, `rdZ`, `rdth`', '(%d.6)' % secno, 'Higgs–radion mixing'],
    ], widths=[0.17, 0.22, 0.2, 0.41], small=True)
    d.align([r'T^4&=\mathbb R^4/(L_c\mathbb Z)^4,\qquad b^8=\det(-g_{kl}),\qquad V_\perp=b^4L_c^4,\qquad \ell_\perp=bL_c,&&(%d.1)' % secno,
             r'\phi_{\vec n}&=\frac{e^{2\pi i\vec n\cdot\vec y/L_c}}{\sqrt{V_\perp}},\qquad \mu_{\vec n}=\frac{2\pi|\vec n|}{\ell_\perp},\qquad r_4(k)=8\sigma(k)-32\sigma(k/4),&&(%d.2)' % secno,
             r'W&=\{x_4=\dots=x_7=0\},\qquad h_{\mu\nu}=g_{\mu\nu}(x,0),\qquad T^W_{\mu\nu}=f_Wh_{\mu\nu},\qquad \psi=U\sum_a\psi_a\otimes e_a,&&(%d.3)' % secno,
             r'ds^2&=e^{-4\beta}g_E-e^{2\beta}\gamma,\quad \phi=\sqrt{12}M_{\rm Pl}\beta,\quad V=M_{\rm Pl}^2\Big[e^{-4\beta}(\Lambda+\kappa U)-\frac{6\sigma}{r_0^2}e^{-6\beta}\Big],&&(%d.4)' % secno,
             r'd_l&=\frac{(2l+3)(l+1)(l+2)}6,\qquad Z_l=\sqrt{\frac{d_l}V}\frac{C_l^{(3/2)}(\cos x_4)}{C_l^{(3/2)}(1)},&&(%d.5)' % secno,
             r'\gamma&=\frac v{\sqrt3M_{\rm Pl}},\qquad Z^2=1+6\xi\gamma^2(1-6\xi),\qquad \tan2\theta=\frac{12\gamma\xi Zm_{h0}^2}{m_{r0}^2-m_{h0}^2(Z^2-36\xi^2\gamma^2)}.&&(%d.6)' % secno])


# =====================================================================================================================
def provenance_section(d, D, D3, runs, secno=10):
    (a, b, c_, e), (a6, b6, c6, e6) = counts(D3)
    w5, w6 = TB.worst(D3), RD.worst(D3)
    d.h1('%d. Sections 25–26 of 7+1.nb: the auxiliary torus, the world brane and the radion' % secno)
    d.p('This section records the provenance of sections 25 and 26 of 7+1.nb (appended on 28 October 2026; the summary of all checks became '
        'section 26 and then section %d).' % (D3['n_last'] + 1))
    d.h2('%d.1 The task' % secno)
    d.ol(['Define and explain the "auxiliary torus" and the "world brane" precisely and correctly, each in its own documentation, symbol and '
          'provenance files (Markdown, LaTeX, PDF), for students, with examples (summarized; the original wording is the author’s private working '
          'material).',
          'Read a note of the author on the radion and use its ideas to verify and clarify the radion and brane documents; create examples in which '
          'the four extra space dimensions are a 4-sphere and a 4-torus.',
          'Push every completed stage to github.com/once-ere/7plus1 and verify the repository.'])
    d.h2('%d.2 Inputs and decisions' % secno)
    d.p('Inputs: the cosmological and inhomogeneous metrics and the curvature routine (section 11), the Standard Model and the brane '
        'representation (17), the Kaluza–Klein operator (18), the stress tensor, the radion and the stabilizer (20), the regimes (21–23), the '
        'portals (24). The note was read from the author’s PDF (not part of the repository). Decisions: the brane at $\\vec y=0$, not dynamical; the '
        'sphere compared at the same volume $V_\\perp$; the warped comparison in five dimensions; all printed numbers computed, never typed; the '
        'Higgs–radion numbers in 60-digit arithmetic. Before each section was appended the evaluated notebook was copied to '
        '`backups/7+1_2026-10-28_before_torus_brane.nb` and `backups/7+1_2026-10-28_before_radion.nb` (excluded by .gitignore).')
    d.h2('%d.3 Verification' % secno)
    d.p('Sections 25 and 26 contain %d and %d checks (%d and %d True); they are rows of section 4.2 above. The independent programs '
        '`tools/ref_torus_brane.py` (%d boolean checks, %d True) and `tools/ref_radion.py` (%d, %d True) recompute every number with their own code; '
        'the %d and %d groups of numbers printed by the notebook and by them agree to a largest relative difference of %s and %s; their complete '
        'outputs are printed in section 4.3. `tools/verify_7p1_nb.py` confirms that, with each new section removed, the notebook equals its backup.'
        % (a, a6, b, b6, c_, e, c6, e6, len(TB.compare(D3)), len(RD.compare(D3)), fmt(w5, 2) if w5 is not None else '?',
           fmt(w6, 2) if w6 is not None else '?'))
