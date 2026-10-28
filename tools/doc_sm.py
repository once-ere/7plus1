"""Content of STANDARD_MODEL (Markdown + LaTeX): the Standard Model on the world brane of 7+1 spacetime, dirac16complex at
the brane (four brane Dirac spinors, Kaluza-Klein masses), the portals of dirac16complex, the coupled field equations, and
neutrino flavour conversion (sections 17-19 of 7+1.nb).  Every equation used is written here; every number is parsed from
the evaluated notebook or from the output of tools/ref_sm.py."""
from docgen import Doc
from doc_stage3_common import fmt, mnum, nums_in
from doc_common import signature_block

PY = {
    'ew': "NUM v = (sqrt2 G_F)^(-1/2) [GeV], g, g', sin^2 theta_W (on shell), e, 1/alpha(tree), lambda = m_H^2/(2v^2)",
    'ckm': 'NUM CKM |V| (rows ud us ub / cd cs cb / td ts tb)',
    'J': 'NUM Jarlskog J = Im(V_us V_cb V_ub* V_cs*)',
    'wolf': 'NUM Wolfenstein lambda = s12, A = s23/s12^2, rhobar + i etabar = -V_ud V_ub*/(V_cd V_cb*)',
    'yt': 'NUM top Yukawa y_t = sqrt2 m_t/v',
    'pmns': 'NUM PMNS |U| and J_PMNS',
    'mnu': 'NUM Dirac neutrino masses (lightest m1 = 0, eV) and their sum',
    'ynu': 'NUM largest neutrino Yukawa sqrt2 m3/v',
    'dvac': 'NUM DUNE vacuum (L = 1300 km, E = 2.5 GeV): P(mu->e), P(mubar->ebar), P(mu->mu)',
    'vcc': 'NUM V_CC = sqrt2 G_F n_e for rho Y_e = 1 g/cm^3 [eV] (textbook 7.63e-14 eV)',
    'dmat': 'NUM DUNE in matter (rho = 2.848 g/cm^3, Y_e = 0.5): P(mu->e), P(mubar->ebar)',
    'sun': 'NUM Sun: P_ee (numerical through the solar density profile, adiabatic formula) for E = 0.3, 1, 5, 10 MeV',
    'sn': 'NUM supernova MSW resonance densities for E = 20 MeV, Y_e = 0.5: H (dm31, theta13), L (dm21, theta12) [g/cm^3]',
    'cnb': 'NUM cosmic neutrino background: T_nu = (4/11)^(1/3) T_gamma [K], number density per flavour [cm^-3]',
    'mix': 'portal mixing: light masses (eV) and heavy masses (GeV) of the 7x7 Dirac mass matrix',
    'dev': 'portal mixing: N N^dag = 1 - mu mu^dag/m^2 to second order (non-unitary light mixing), (1 - N N^dag)_ee',
    'vp': 'NUM vector portal (V_P,e = 0.2 V_CC): P(mu->e) at DUNE with and without the dirac16complex potential',
    'age': 'NUM age of the universe from the brane Friedmann equation (Planck 2018: H0 = 67.4, Omega_m = 0.315) [Gyr]',
    'add': 'NUM ADD example (reduced 8D Planck mass 10 TeV): V_perp [GeV^-4], V_perp^(1/4) [m], lowest KK mass 2 pi/V_perp^(1/4) [MeV]',
}


def compare(D):
    """(label, notebook values, Python values, note) for every number printed by both programs"""
    nm = D['nums']['sm']; ref = D['ref'].get('sm', ('', '', {}))[2]
    g = lambda k: nums_in(nm.get(k, '') or '')
    p = lambda k: nums_in(ref.get(PY[k], '') or '')
    out = []
    out.append(("$v,g,g',\\sin^2\\theta_W,e,1/\\alpha,\\lambda$", sum((g(k) for k in ('v', 'g', "g'", 'sin^2 thetaW', 'e', '1/alpha tree', 'lambda')), []), p('ew'), ''))
    out.append(('$|V_{\\rm CKM}|$', g('|V_CKM|'), p('ckm'), 'Python prints 6 decimals'))
    out.append(('$J_{\\rm CKM}$', g('J_CKM'), p('J'), ''))
    out.append(('$\\bar\\rho,\\bar\\eta$', g('rhobar + i etabar'), p('wolf')[2:4], ''))
    pm = p('pmns')
    out.append(('$|U_{\\rm PMNS}|$', g('|U_PMNS|'), pm[:9], 'Python prints 5 decimals'))
    out.append(('$J_{\\rm PMNS}$', g('J_PMNS'), pm[9:10], ''))
    out.append(('$m_1,m_2,m_3$ (eV)', g('neutrino masses (eV)'), p('mnu')[:3], ''))
    out.append(('$y_t$', g('y_t'), p('yt'), ''))
    out.append(('age (Gyr)', g('age (Gyr)'), p('age'), ''))
    out.append(('$V_\\perp$, $V_\\perp^{1/4}$, lowest KK mass', g('V_perp (GeV^-4)') + g('V_perp^(1/4) (m)') + g('lowest KK mass (MeV)'), p('add'), ''))
    out.append(('light (eV) and heavy (GeV) masses with the portal', g('light masses (eV)') + g('heavy masses (GeV)'), p('mix'), ''))
    out.append(('$(1-NN^\\dagger)_{ee}$', g('(1 - N.N^dag)_ee'), p('dev')[-1:], ''))
    out.append(('DUNE vacuum $P(\\mu\\to e)$, $P(\\bar\\mu\\to\\bar e)$, $P(\\mu\\to\\mu)$', g('DUNE vacuum P(mu->e), P(mubar->ebar), P(mu->mu)'), p('dvac'), ''))
    out.append(('DUNE matter $P(\\mu\\to e)$, $P(\\bar\\mu\\to\\bar e)$', g('DUNE matter P(mu->e), P(mubar->ebar)'), p('dmat'), ''))
    out.append(('vector portal: $P(\\mu\\to e)$ with, without', g('DUNE P(mu->e) with vector portal V_P,e = 0.2 V_CC, without'), p('vp'), ''))
    out.append(('Sun $P_{ee}$ ($E$, numerical, adiabatic)', g('Sun P_ee {E, numerical, adiabatic}'), p('sun'), ''))
    out.append(('supernova resonance densities H, L', g('SN resonance densities H, L (g/cm^3)'), p('sn'), ''))
    out.append(('$T_\\nu$, $n_\\nu$', g('T_nu, n_nu'), p('cnb'), ''))
    return out


def max_rel(a, b):
    if not a or not b or len(a) != len(b):
        return None
    return max(abs(x - y)/max(abs(x), abs(y), 1e-300) for x, y in zip(a, b))


def build(D, today):
    ck = D['checks']['sm']; ref = D['ref'].get('sm', ('', '', {}))[2]; nm = D['nums']['sm']
    nck = len(ck); nok = sum(1 for c in ck if c[1] == 'True')
    rb = [(k, v) for k, v in ref.items() if v in ('True', 'False')]
    g = lambda k: nums_in(nm.get(k, '') or '')
    d = Doc('The Standard Model on the world brane of 7+1 spacetime, the portals of dirac16complex, and neutrino flavour',
            'The complete Standard Model in the gravitational field of 7+1 spacetime, electroweak symmetry breaking, quark and lepton '
            'masses and mixing, dirac16complex at the brane as four Dirac spinors, Kaluza–Klein masses, the portals and the coupled field '
            'equations, and neutrino flavour conversion in vacuum, in matter, in the Sun and in supernovae — sections 17–19 of the '
            'notebook 7+1.nb',
            ['Repository: github.com/once-ere/7plus1 — prepared for Nash',
             'Notebook: 7+1.nb (Wolfram Language 15.0.1), evaluated ' + D['log'].get('finished', ''),
             'Document generated by tools/make_docs.py on ' + today],
            'The 7+1 spacetime of this project has one time coordinate $x_0$, three ordinary space coordinates $x_1,x_2,x_3$ and four '
            'auxiliary space coordinates $x_4,\\dots,x_7$. This document places the complete Standard Model of particle physics on the '
            '3+1-dimensional WORLD BRANE $x_4=x_5=x_6=x_7=0$, in the gravitational field of the whole spacetime, and couples it to the '
            '16-component field dirac16complex, which lives in the whole spacetime. It states the field content and the Lagrangian, '
            'checks anomaly cancellation and electroweak symmetry breaking with current data, decomposes dirac16complex at the brane into '
            'four ordinary Dirac spinors, derives its Kaluza–Klein masses, defines three portals and the coupled field equations, and '
            'works out neutrino flavour conversion with and without dirac16complex. Sections 17–19 of 7+1.nb contain %d checks (%d '
            'True); the independent Python program tools/ref_sm.py contains %d (%d True).' % (nck, nok, len(rb), sum(1 for k, v in rb if v == 'True')))

    # ==============================================================================================
    d.h1('1. Results at a glance')
    v = mnum(nm.get('v', '') or ''); age = mnum(nm.get('age (Gyr)', '') or '')
    J = mnum(nm.get('J_CKM', '') or '')
    dv = g('DUNE vacuum P(mu->e), P(mubar->ebar), P(mu->mu)'); dm = g('DUNE matter P(mu->e), P(mubar->ebar)')
    d.ul(['**Field content.** Gauge group $SU(3)_C\\times SU(2)_L\\times U(1)_Y$, three generations including right-handed neutrinos, '
          'one Higgs doublet; all gauge anomalies cancel in every generation (section 3).',
          '**Electroweak scale.** $v=(\\sqrt2G_F)^{-1/2}=%s$ GeV. At tree level the masses are $m_W=gv/2$, $m_Z=\\sqrt{g^2+g\'^2}\\,v/2$ '
          'and $m_H^2=2\\lambda v^2$ (section 5).' % fmt(v, 10),
          '**Mixing.** CKM matrix with the Jarlskog invariant $J=%s$; PMNS matrix from NuFit-6.0 (section 6).' % fmt(J, 4),
          '**Cosmology on the brane.** The brane Friedmann equation gives the age %s Gyr (section 7).' % fmt(age, 5),
          '**dirac16complex at the brane** is four brane Dirac spinors; a mode with auxiliary momentum $\\vec q$ has the Kaluza–Klein '
          'mass $\\sqrt{m^2+|\\vec q|^2/b_0^2}$ (sections 8–9).',
          '**Portals** (Yukawa, vector, Higgs) couple dirac16complex to the Standard Model at the brane; the Yukawa portal makes the light '
          'neutrino mixing non-unitary (section 10).',
          '**Neutrino flavour.** At DUNE ($L=1300$ km, $E=2.5$ GeV) $P(\\nu_\\mu\\to\\nu_e)=%s$ in vacuum and %s in matter (section 11).'
          % (fmt(dv[0], 5) if dv else '?', fmt(dm[0], 5) if dm else '?')])

    # ==============================================================================================
    d.h1('2. The world brane and the total action')
    d.p('Units $\\hbar=c=1$. The auxiliary coordinates form a torus $T^4$: $x_k\\simeq x_k+L_c$ ($k=4,\\dots,7$; coordinate period $L_c$, proper '
        'size $bL_c$). The Standard Model lives on the world brane')
    d.eq(r'W=\{x_4=x_5=x_6=x_7=0\},\qquad x=(x_0,x_1,x_2,x_3),\quad y=(x_4,\dots,x_7),\\ h_{\mu\nu}(x)=g_{\mu\nu}(x,y=0)\ \ (\mu,\nu=0,\dots,3),')
    d.p('with the induced metric $h_{\\mu\\nu}$; gravity and dirac16complex live in the whole 7+1 spacetime (the bulk).')
    signature_block(d, eta=r'\eta', eta71=r'\eta_{71}', gam=r'\Gamma')
    d.p('The brane has one time-like dimension ($x_0$) and three space-like dimensions ($x_1,x_2,x_3$) — in the counting of xAct the '
        'list $\\{3,1,0\\}$ — and its flat frame metric $\\mathrm{diag}(+1,-1,-1,-1)$ is $\\eta$ restricted to $x_0,\\dots,x_3$: the squares '
        'of the basis $e_0=1,e_1,e_2,e_3$ of the quaternion subalgebra $\\mathbb H$ of the octonions (checked in section 17 of 7+1.nb). '
        'The total action is')
    d.eq(r'S=S_{\rm EL}[g]+S_{\psi}[\psi,g]+\int_Wd^4x\,\sqrt{|h|}\,\big(\mathcal L_{\rm SM}+\mathcal L_{\rm portal}\big)+S_{\rm stab},\qquad '
         r'\int_Wd^4x\,\sqrt{|h|}\,(\dots)=\int d^8X\,\sqrt{|g|}\,\frac{\delta^4(y)}{b^4}\,(\dots),')
    d.p('where the last identity holds for the cosmological metric $ds^2=dx_0^2-a^2\\sum_{i=1}^3dx_i^2-b^2\\sum_{k=4}^7dx_k^2$ '
        '($\\sqrt{|g|}=a^3b^4$, $\\sqrt{|h|}=a^3$). The pieces:')
    d.align([r'S_{\rm EL}&=-\frac1{2\kappa}\int d^8X\sqrt{|g|}\sum_{k=0}^{3}\alpha_kL_k,\qquad L_k=\frac1{2^k}\delta^{\alpha_1\beta_1\cdots\alpha_k\beta_k}_{\gamma_1\delta_1\cdots\gamma_k\delta_k}R^{\gamma_1\delta_1}{}_{\alpha_1\beta_1}\cdots R^{\gamma_k\delta_k}{}_{\alpha_k\beta_k},',
             r'S_\psi&=\int d^8X\sqrt{|g|}\Big[\bar\psi\big(i\Gamma^\mu(X)D_\mu-m\big)\psi-\frac12\sum_{r=0}^{8}g_r\sum_{|K|=r}(\bar\psi\Theta_K\psi)(\bar\psi\Theta^K\psi)\Big],'])
    d.p('($\\kappa=8\\pi G_8$; $\\alpha_0=2\\bar\\Lambda$, $\\alpha_1=1$; $\\Gamma^\\mu(X)$ the generalized frame matrices, $D_\\mu=\\partial_\\mu+\\Omega_\\mu$ '
        'with the canonical spin connection $\\Omega_\\mu$; $\\Theta^K$ the 256 Clifford products with $g_r$ the nine self-couplings of '
        'dirac16complex), and $S_{\\rm stab}$ is a stress that keeps $b$ constant in the modern epoch.')

    # ==============================================================================================
    d.h1('3. Fields, quantum numbers and anomaly cancellation (section 17.1)')
    d.p('Gauge fields $G^a_\\mu$ ($a=1,\\dots,8$), $W^I_\\mu$ ($I=1,2,3$), $B_\\mu$ with couplings $g_s$, $g$, $g\'$. Matter in each of three '
        'generations (Weyl fields; electric charge $Q=T_3+Y$):')
    d.table(['field', '$SU(3)_C$', '$SU(2)_L$', '$Y$', 'chirality $\\chi_f$', 'charges'], [
        ['$Q_L=(u_L,d_L)$', '3', '2', '$1/6$', '$+1$ (left)', '$2/3$, $-1/3$'],
        ['$u_R$', '3', '1', '$2/3$', '$-1$ (right)', '$2/3$'],
        ['$d_R$', '3', '1', '$-1/3$', '$-1$', '$-1/3$'],
        ['$L_L=(\\nu_L,e_L)$', '1', '2', '$-1/2$', '$+1$', '$0$, $-1$'],
        ['$e_R$', '1', '1', '$-1$', '$-1$', '$-1$'],
        ['$\\nu_R$', '1', '1', '$0$', '$-1$', '$0$'],
        ['$\\Phi$ (Higgs)', '1', '2', '$1/2$', '—', '$\\tilde\\Phi=i\\sigma_2\\Phi^*$'],
    ], widths=[0.2, 0.12, 0.12, 0.1, 0.18, 0.28])
    d.p('With $N_f$ the number of colours and $d_f$ the $SU(2)$ dimension of field $f$, the anomaly coefficients of one generation are')
    d.align([r'A_{SU(3)^2U(1)}&=\sum_{f:\,N_f=3}\chi_fd_fY_f,\qquad A_{SU(2)^2U(1)}=\sum_{f:\,d_f=2}\chi_fN_fY_f,\qquad A_{SU(3)^3}=\sum_{f:\,N_f=3}\chi_fd_f,',
             r'A_{U(1)^3}&=\sum_f\chi_fN_fd_fY_f^3,\qquad A_{\rm grav^2U(1)}=\sum_f\chi_fN_fd_fY_f,'])
    d.p('and all five vanish, with or without $\\nu_R$; the number of $SU(2)$ doublets per generation ($3+1=4$) is even (no Witten anomaly). '
        'The relativistic degrees of freedom above the electroweak scale (without $\\nu_R$) are')
    d.eq(r'g_*=\underbrace{16}_{G}+\underbrace{6}_{W}+\underbrace{2}_{B}+\underbrace{4}_{\Phi}+\frac78\cdot3\cdot\big(12+6+6+4+2\big)=28+\frac78\cdot90=106.75 .')

    # ==============================================================================================
    d.h1('4. The Lagrangian on the brane, in the gravitational field (section 17.2)')
    d.p('With the frame (tetrad) of $h_{\\mu\\nu}$, its spin connection and complex 4×4 Dirac matrices $\\gamma^\\mu(x)$ of the brane '
        '($\\gamma^\\mu\\gamma^\\nu+\\gamma^\\nu\\gamma^\\mu=2h^{\\mu\\nu}$), and the gauge-covariant derivative')
    d.eq(r"D_\mu=\nabla_\mu-ig_s\,G^a_\mu T^a-ig\,W^I_\mu\frac{\tau^I}2-ig'\,Y\,B_\mu")
    d.p('($T^a$ the generators of $SU(3)$ in the representation of the field, $\\tau^I$ the Pauli matrices; $\\nabla_\\mu$ includes the spin '
        'connection of $h$ on spinors), the Standard Model Lagrangian is')
    d.align([r'\mathcal L_{\rm SM}&=-\tfrac14G^a_{\mu\nu}G^{a\mu\nu}-\tfrac14W^I_{\mu\nu}W^{I\mu\nu}-\tfrac14B_{\mu\nu}B^{\mu\nu}+(D_\mu\Phi)^\dagger(D^\mu\Phi)-\mu^2\Phi^\dagger\Phi-\lambda(\Phi^\dagger\Phi)^2',
             r'&\quad+\sum_f\bar f\,i\gamma^\mu D_\mu f-\big[\bar Q_LY_u\tilde\Phi u_R+\bar Q_LY_d\Phi d_R+\bar L_LY_e\Phi e_R+\bar L_LY_\nu\tilde\Phi\nu_R+\mathrm{h.c.}\big],',
             r'G^a_{\mu\nu}&=\partial_\mu G^a_\nu-\partial_\nu G^a_\mu+g_sf^{abc}G^b_\mu G^c_\nu,\quad W^I_{\mu\nu}=\partial_\mu W^I_\nu-\partial_\nu W^I_\mu+g\,\epsilon^{IJK}W^J_\mu W^K_\nu,\quad B_{\mu\nu}=\partial_\mu B_\nu-\partial_\nu B_\mu,'])
    d.p('$f$ runs over all Weyl fields of the three generations, and $Y_u,Y_d,Y_e,Y_\\nu$ are complex 3×3 Yukawa matrices. Its field equations '
        'are the ordinary ones with $\\partial\\to\\nabla$ (the Levi-Civita connection of $h$), e.g. $\\nabla_\\mu G^{a\\mu\\nu}+g_sf^{abc}G^b_\\mu G^{c\\mu\\nu}'
        '=-g_s\\sum_q\\bar q\\gamma^\\nu T^aq$.')

    # ==============================================================================================
    d.h1('5. Electroweak symmetry breaking (section 17.3)')
    d.p('With $\\langle\\Phi\\rangle=(0,v/\\sqrt2)$ the term $(D_\\mu\\langle\\Phi\\rangle)^\\dagger(D^\\mu\\langle\\Phi\\rangle)$ is a mass matrix for '
        '$(W^1,W^2,W^3,B)$, its Hessian $M^2_{ij}=\\partial^2\\big|g\\,W^I\\tfrac{\\tau^I}2\\langle\\Phi\\rangle+\\tfrac{g\'}2B\\langle\\Phi\\rangle\\big|^2/\\partial w_i\\partial w_j$ having the eigenvalues')
    d.eq(r"m_W^2=\frac{g^2v^2}4\ \ (\text{twice}),\qquad m_Z^2=\frac{(g^2+g'^2)v^2}4,\qquad m_\gamma^2=0 .")
    d.p('Inputs (Particle Data Group 2025): $G_F=1.1663788\\times10^{-5}$ GeV$^{-2}$, $m_W=80.3692$ GeV, $m_Z=91.1880$ GeV, $m_H=125.20$ GeV. At '
        'tree level (on-shell weak angle)')
    d.eq(r"v=(\sqrt2G_F)^{-1/2},\quad g=\frac{2m_W}{v},\quad \cos\theta_W=\frac{m_W}{m_Z},\quad g'=g\tan\theta_W,\quad e=g\sin\theta_W,\quad \lambda=\frac{m_H^2}{2v^2},")
    d.table(['$v$ (GeV)', '$g$', "$g'$", '$\\sin^2\\theta_W$', '$e$', '$1/\\alpha=4\\pi/e^2$', '$\\lambda$'],
            [[fmt(mnum(nm.get(k, '') or ''), 10) for k in ('v', 'g', "g'", 'sin^2 thetaW', 'e', '1/alpha tree', 'lambda')]],
            widths=[0.15, 0.14, 0.14, 0.15, 0.14, 0.14, 0.14], small=True)
    d.p('($1/\\alpha$ here is the tree-level value at the weak scale, not the low-energy fine-structure constant.)')

    # ==============================================================================================
    d.h1('6. Quark and lepton masses; CKM and PMNS mixing (section 17.4)')
    d.p('Fermion masses are $m=yv/\\sqrt2$ for the singular values $y$ of the Yukawa matrices; the mixing matrices are $V_{\\rm CKM}=U_{uL}^\\dagger '
        'U_{dL}$ and $U_{\\rm PMNS}=U_{eL}^\\dagger U_{\\nu L}$, where $Y=U_L\\,\\mathrm{diag}(y)\\,U_R^\\dagger$. Both are written in the standard '
        'parametrization ($c_{ij}=\\sqrt{1-s_{ij}^2}$):')
    d.eq(r'U(s_{12},s_{13},s_{23},\delta)=\begin{pmatrix}1&0&0\\0&c_{23}&s_{23}\\0&-s_{23}&c_{23}\end{pmatrix}'
         r'\begin{pmatrix}c_{13}&0&s_{13}e^{-i\delta}\\0&1&0\\-s_{13}e^{i\delta}&0&c_{13}\end{pmatrix}'
         r'\begin{pmatrix}c_{12}&s_{12}&0\\-s_{12}&c_{12}&0\\0&0&1\end{pmatrix}.')
    d.p('**Quarks** (PDG 2024 CKM fit: $s_{12}=0.22501$, $s_{13}=0.003732$, $s_{23}=0.04183$, $\\delta=1.147$; PDG 2025 masses: $m_u=2.16$ MeV, '
        '$m_d=4.70$ MeV, $m_s=93.5$ MeV, $m_c=1.2730$ GeV, $m_b=4.183$ GeV, $m_t=172.56$ GeV). With $Y_u=\\sqrt2\\,\\mathrm{diag}(m_u,m_c,m_t)/v$ and '
        '$Y_d=\\sqrt2\\,V\\,\\mathrm{diag}(m_d,m_s,m_b)/v$ the singular values reproduce the masses and $|U_{uL}^\\dagger U_{dL}|=|V|$. The notebook obtains')
    V = g('|V_CKM|')
    if len(V) == 9:
        d.eq(r'|V_{\rm CKM}|=\begin{pmatrix}%s\end{pmatrix},' % r'\\'.join('&'.join(fmt(V[3*i + j], 6) for j in range(3)) for i in range(3)))
    rb_ = g('rhobar + i etabar')
    d.eq(r'J=\mathrm{Im}\big(V_{us}V_{cb}V_{ub}^*V_{cs}^*\big)=%s,\qquad \bar\rho+i\bar\eta=-\frac{V_{ud}V_{ub}^*}{V_{cd}V_{cb}^*}=%s+%s\,i,\\ '
         r'A=\frac{s_{23}}{s_{12}^2}=%s,\qquad y_t=\frac{\sqrt2m_t}{v}=%s,' % (
             fmt(J, 6), fmt(rb_[0], 6) if len(rb_) > 1 else '?', fmt(rb_[1], 6) if len(rb_) > 1 else '?', fmt(0.04183/0.22501**2, 6),
             fmt(mnum(nm.get('y_t', '') or ''), 8)))
    d.p('in agreement with the PDG 2024 values $J=(3.12^{+0.13}_{-0.12})\\times10^{-5}$, $A=0.826$, $\\bar\\rho=0.1591$, $\\bar\\eta=0.3523$.')
    d.p('**Leptons** (NuFit-6.0, normal ordering, without the atmospheric data of Super-Kamiokande: $\\sin^2\\theta_{12}=0.307$, '
        '$\\sin^2\\theta_{23}=0.561$, $\\sin^2\\theta_{13}=0.02195$, $\\delta=177^\\circ$, $\\Delta m^2_{21}=7.49\\times10^{-5}$ eV², '
        '$\\Delta m^2_{31}=2.534\\times10^{-3}$ eV², $\\Delta m^2_{ij}=m_i^2-m_j^2$). With $\\nu_R$ the neutrinos are Dirac particles; with the '
        'lightest mass $m_1=0$: $m_2=\\sqrt{\\Delta m^2_{21}}$, $m_3=\\sqrt{\\Delta m^2_{31}}$, and $Y_\\nu=\\sqrt2\\,U\\,\\mathrm{diag}(m_1,m_2,m_3)/v$ '
        '(electron Yukawa basis diagonal):')
    U = g('|U_PMNS|'); mnu = g('neutrino masses (eV)')
    if len(U) == 9:
        d.eq(r'|U_{\rm PMNS}|=\begin{pmatrix}%s\end{pmatrix},\qquad J_{\rm PMNS}=%s,\\ (m_1,m_2,m_3)=(%s)\ \text{eV}.' % (
            r'\\'.join('&'.join(fmt(U[3*i + j], 5) for j in range(3)) for i in range(3)), fmt(mnum(nm.get('J_PMNS', '') or ''), 5),
            ', '.join(fmt(x, 6) for x in mnu)))

    # ==============================================================================================
    d.h1('7. The modern epoch on the brane (section 17.5)')
    d.p('With $b=b_0$ the 00 field equation is the Friedmann equation of the brane,')
    d.eq(r'3H_a^2=8\pi G_4\,\rho_4+\Lambda_4,\qquad G_4=\frac{G_8}{V_\perp},\qquad V_\perp=b_0^4L_c^4,\qquad \varepsilon_8=\frac{\rho_4}{V_\perp}\ \ \text{(matter on the brane)},')
    d.p('($H_a=a\'/a$, $\\rho_4$ the four-dimensional energy density). With the Planck 2018 parameters $H_0=67.4$ km/s/Mpc, $\\Omega_m=0.315$, '
        '$\\Omega_r=9.14\\times10^{-5}$, $\\Omega_\\Lambda=1-\\Omega_m-\\Omega_r$ the age of the universe is')
    d.eq(r't_0=\frac1{H_0}\int_0^1\frac{da}{a\sqrt{\Omega_r a^{-4}+\Omega_m a^{-3}+\Omega_\Lambda}}=%s\ \text{Gyr}\qquad(\text{Planck 2018: }13.787\pm0.020).' % fmt(age, 6))
    vp = mnum(nm.get('V_perp (GeV^-4)', '') or ''); vq = mnum(nm.get('V_perp^(1/4) (m)', '') or ''); kk = mnum(nm.get('lowest KK mass (MeV)', '') or '')
    d.p('With the reduced eight-dimensional Planck mass $M_8$ ($8\\pi G_8=M_8^{-6}$) the reduced four-dimensional one is $M_{\\rm Pl}^2=M_8^6V_\\perp$ '
        '($M_{\\rm Pl}=2.435\\times10^{18}$ GeV). For $M_8=10$ TeV:')
    d.eq(r'V_\perp=\frac{M_{\rm Pl}^2}{M_8^6}=%s\ \text{GeV}^{-4},\qquad V_\perp^{1/4}=%s\ \text{m},\qquad \frac{2\pi}{V_\perp^{1/4}}=%s\ \text{MeV},'
         % (fmt(vp, 6), fmt(vq, 6), fmt(kk, 6)))
    d.p('the last being the lowest Kaluza–Klein mass of a field on the auxiliary torus (proper period $V_\\perp^{1/4}=b_0L_c$).')

    # ==============================================================================================
    d.h1('8. dirac16complex at the brane: four brane Dirac spinors (section 18)')
    d.p('At the brane the eight Dirac matrices split into four along $W$ and four transverse to it. With the complex 4×4 Dirac matrices of '
        'the Dirac representation and four anti-Hermitian 4×4 matrices $\\rho_k$,')
    d.align([r'\gamma_0&=\mathrm{diag}(1,1,-1,-1),\quad \gamma_i=\begin{pmatrix}0&\sigma_i\\-\sigma_i&0\end{pmatrix},\quad \gamma_5=i\gamma_0\gamma_1\gamma_2\gamma_3,'
             r'\qquad \gamma_\mu\gamma_\nu+\gamma_\nu\gamma_\mu=2\eta_{\mu\nu}1_4,',
             r'\rho_4&=i\gamma_0,\quad \rho_5=-\gamma_1,\quad \rho_6=-\gamma_2,\quad \rho_7=-\gamma_3,\qquad \rho_k\rho_l+\rho_l\rho_k=-2\delta_{kl}1_4,'])
    d.p('($\\sigma_i$ the Pauli matrices, $\\eta=\\mathrm{diag}(1,-1,-1,-1)$, the squares of the quaternion basis $1,e_1,e_2,e_3$) the '
        'sixteen-dimensional matrices')
    d.eq(r'R_\mu=\gamma_\mu\otimes1_4\ \ (\mu=0,\dots,3),\qquad R_k=\gamma_5\otimes\rho_k\ \ (k=4,\dots,7)')
    d.p('satisfy the same Clifford relations as the real Dirac matrices $\\Gamma_A$ of 7+1 spacetime, $R_AR_B+R_BR_A=2\\eta_{AB}1_{16}$ with '
        '$\\eta=\\mathrm{diag}(+1,-1,\\dots,-1)$, the squares of the octonion basis. Both representations are irreducible, so there is a '
        'unitary $U$ with $\\Gamma_AU=UR_A$; it is '
        'obtained as the Schur average over the 256 Clifford products,')
    d.eq(r'U\propto\sum_{K}\Gamma^{(K)}\,s\,\big(R^{(K)}\big)^{-1},\qquad \Gamma^{(K)}=\Gamma_{A_1}\cdots\Gamma_{A_r},\quad R^{(K)}=R_{A_1}\cdots R_{A_r},')
    d.p('with a fixed matrix $s$ (the first matrix unit for which the sum is non-zero), normalized so that $U^\\dagger U=1_{16}$. Writing '
        '$\\psi=U\\zeta$ with $\\zeta\\in\\mathbb C^4\\otimes\\mathbb C^4$, $\\zeta=\\sum_{a=1}^4\\psi_a\\otimes e_a$, dirac16complex at the brane is FOUR brane Dirac '
        'spinors $\\psi_a$ ($a$ the transverse spinor index), and with $\\bar\\psi_a=\\psi_a^\\dagger\\gamma_0$, $\\gamma^\\mu=\\eta^{\\mu\\mu}\\gamma_\\mu$:')
    d.eq(r"\bar\psi\psi'=\sum_{a=1}^4\bar\psi_a\psi'_a,\qquad \bar\psi\,\Gamma^\mu\psi'=\sum_{a=1}^4\bar\psi_a\gamma^\mu\psi'_a\ \ (\mu=0,\dots,3).")
    d.p('The eight-dimensional chirality $\\chi=i\\Gamma_0\\cdots\\Gamma_7$ is')
    d.eq(r'\chi=U\,(\gamma_5\otimes\rho_{4567})\,U^\dagger,\qquad \rho_{4567}=\rho_4\rho_5\rho_6\rho_7\ \ \text{(Hermitian, eigenvalues }+1,+1,-1,-1),')
    d.p('so the left-handed half $\\psi_L$ ($\\chi=-1$) contains brane-right-handed parts with $\\rho_{4567}=-1$ and brane-left-handed parts with '
        '$\\rho_{4567}=+1$.')

    # ==============================================================================================
    d.h1('9. Kaluza–Klein masses on the auxiliary torus (section 18.1)')
    d.p('With $b=b_0$ and $\\psi\\propto e^{i\\vec q\\cdot\\vec y}$, $\\vec q=2\\pi\\vec n/L_c$ ($\\vec n\\in\\mathbb Z^4$), the eight-dimensional Dirac '
        'equation at the brane becomes, in the variables $\\zeta$,')
    d.eq(r'\big(\gamma\cdot p\otimes1_4+X-m\big)\zeta=0,\qquad \gamma\cdot p=p_0\gamma_0-\sum_{i=1}^3p_i\gamma_i,\qquad X=\gamma_5\otimes\frac{\vec\rho\cdot\vec q}{b_0},\qquad '
         r'X^2=-\frac{|\vec q|^2}{b_0^2}1_{16},')
    d.p('and $X$ anticommutes with every $\\gamma_\\mu\\otimes1_4$. Hence each Fourier mode is four brane Dirac fermions of mass')
    d.eq(r'm_{\vec n}=\sqrt{m^2+\frac{|\vec q|^2}{b_0^2}}\qquad(\text{the zero mode }\vec n=0\text{ has mass }m),')
    d.p('checked by the vanishing of $\\det(\\gamma\\cdot p\\otimes1_4+X-m)$ on $p_0^2-|\\vec p|^2=m_{\\vec n}^2$ (and not off shell), in 40-digit '
        'arithmetic, for $m=7/10$, $b_0=13/10$, $\\vec q=2\\pi(2,-1,\\frac12,3)$, $\\vec p=(0.3,-0.2,0.4)$.')

    # ==============================================================================================
    d.h1('10. The portals and the coupled field equations (section 18.2)')
    d.p('The portals of dirac16complex are local couplings at the brane (an invented interaction; separate couplings for the two '
        'eight-dimensional chiral halves), with $\\psi_a(x)$ the brane spinors of $\\psi(x,y=0)$:')
    d.ul(['**Neutrino (Yukawa) portal:** $-\\sum_{i,s}\\big[y^L_{is}\\,\\bar L_i\\tilde\\Phi P_R\\psi^{(L)}_s+y^R_{is}\\,\\bar L_i\\tilde\\Phi P_R\\psi^{(R)}_s\\big]+\\mathrm{h.c.}$, '
          'where $\\psi^{(L)}_s$ ($\\psi^{(R)}_s$), $s=1,2$, are the brane-right-handed parts of the left (right) half of dirac16complex and '
          '$P_R=\\frac12(1+\\gamma_5)$. After symmetry breaking it is a Dirac mass $\\mu_{ia}=y_{ia}v/(\\sqrt2\\sqrt{V_\\perp})$ between $\\nu_{L,i}$ and every '
          'Kaluza–Klein mode ($1/\\sqrt{V_\\perp}$ normalizes the modes on the torus).',
          '**Vector (current) portal:** $-G^f_P\\,(\\bar\\psi\\Gamma^\\mu\\psi)|_W\\,(\\bar f\\gamma_\\mu f)$ for Standard Model fermions $f$. In a background '
          'of dirac16complex with brane density $n_\\psi$ it adds the potential $V^f_P=G^f_Pn_\\psi$ to the forward scattering of $f$ (for '
          'neutrinos: a non-standard matter effect).',
          '**Higgs (scalar) portal:** $-\\lambda_P\\,\\Phi^\\dagger\\Phi\\,(\\bar\\psi\\psi)|_W$, which shifts the dirac16complex mass by $\\lambda_Pv^2/(2V_\\perp)$ '
          'after symmetry breaking.'])
    d.p('The coupled field equations of all fields:')
    d.align([r'\sum_k\alpha_kE^{(k)\mu}{}_\nu&=\kappa\Big(T_\psi+T_{\rm stab}+\frac{\delta^4(y)}{b^4}T_W\Big)^\mu{}_\nu,',
             r'\big(i\Gamma^\mu D_\mu-m\big)\psi-\sum_rg_r\sum_{|K|=r}(\bar\psi\Theta_K\psi)\Theta^K\psi&=\frac{\delta^4(y)}{b^4}\,J_{\rm portal},'])
    d.p('with $E^{(k)}$ the Lovelock tensors, $T_W$ the energy–momentum tensor of the Standard Model fields and portals on the brane, '
        '$J_{\\rm portal}=\\partial\\mathcal L_{\\rm portal}/\\partial\\bar\\psi$ (for the Yukawa portal $U\\big(\\sum_{i}y^*_{is}\\tilde\\Phi^\\dagger L_i\\otimes e_s+\\dots\\big)$), '
        'and the Standard Model field equations of section 4 with the portal terms as additional sources. In the neutrino sector the Yukawa '
        'portal gives, in the basis $(\\nu_L,\\psi_L)\\times(\\nu_R,\\psi_R)$, the Dirac mass matrix of the zero mode')
    d.eq(r'\mathcal M=\begin{pmatrix}m_D&\mu\\0&m\,1_4\end{pmatrix},\qquad m_D=U_{\rm PMNS}\,\mathrm{diag}(m_1,m_2,m_3),')
    d.p('and for $m\\gg\\mu$ the light mixing matrix $N=(1-\\eta)U_{\\rm PMNS}$ is no longer unitary: $1-NN^\\dagger=2\\eta=\\mu\\mu^\\dagger/m^2$ to second '
        'order, while the heavy states are mostly dirac16complex. Example: $m=1$ GeV and $\\mu$ (in eV)')
    d.eq(r'\mu=\begin{pmatrix}3&1&0&0\\1&2&0&0\\-2&1&0&0\end{pmatrix}\times10^{4}:')
    lm = g('light masses (eV)'); hm = g('heavy masses (GeV)'); dev = g('(1 - N.N^dag)_ee')
    d.eq(r'\text{light masses}=(%s)\ \text{eV},\\ \text{heavy masses}=(%s)\ \text{GeV},\\ (1-NN^\dagger)_{ee}=%s .' % (
        ', '.join(fmt(x, 8) for x in lm), ', '.join(fmt(x, 12) for x in hm), fmt(dev[0], 6) if dev else '?'))

    # ==============================================================================================
    d.h1('11. Neutrino flavour and flavour conversion (section 19)')
    d.p('A neutrino of flavour $\\alpha$ produced with energy $E$ is the superposition $\\sum_iU^*_{\\alpha i}|\\nu_i\\rangle$. In vacuum (with '
        '$\\hbar c=1.973269804\\times10^{-16}$ GeV m) the amplitude after a distance $L$ is $A(\\alpha\\to\\beta)=\\sum_iU_{\\beta i}U^*_{\\alpha i}e^{-im_i^2L/(2E)}$, so')
    d.eq(r'P(\alpha\to\beta)=\delta_{\alpha\beta}-4\sum_{i>j}\mathrm{Re}\,X_{ij}\sin^2\Delta_{ij}+2\sum_{i>j}\mathrm{Im}\,X_{ij}\sin2\Delta_{ij},\qquad '
         r'X_{ij}=U^*_{\alpha i}U_{\beta i}U_{\alpha j}U^*_{\beta j},\quad \Delta_{ij}=\frac{\Delta m^2_{ij}L}{4E}.')
    d.p('In matter of electron density $n_e$ the charged-current forward scattering adds $V_{CC}=\\sqrt2G_Fn_e$ to the $\\nu_e$ entry of the flavour '
        'Hamiltonian ($V_{CC}\\to-V_{CC}$ and $U\\to U^*$ for antineutrinos; a potential common to all flavours drops out):')
    d.eq(r'H=\frac1{2E}\,U\,\mathrm{diag}(0,\Delta m^2_{21},\Delta m^2_{31})\,U^\dagger+\mathrm{diag}(V_{CC},0,0),\qquad '
         r'V_{CC}\big(\rho Y_e=1\ \text{g/cm}^3\big)=7.63\times10^{-14}\ \text{eV}.')
    d.p('In two flavours the mixing angle in matter obeys')
    d.eq(r'\tan2\theta_m=\frac{\sin2\theta}{\cos2\theta-2EV/\Delta m^2},')
    d.p('maximal at the MSW resonance $2EV=\\Delta m^2\\cos2\\theta$. **Long baseline** (DUNE: $L=1300$ km, $E=2.5$ GeV, $\\rho=2.848$ g/cm³, $Y_e=0.5$):')
    pp = g('DUNE P(mu->e) with vector portal V_P,e = 0.2 V_CC, without')
    d.table(['', '$P(\\nu_\\mu\\to\\nu_e)$', '$P(\\bar\\nu_\\mu\\to\\bar\\nu_e)$', '$P(\\nu_\\mu\\to\\nu_\\mu)$'], [
        ['vacuum'] + [fmt(x, 8) for x in dv[:3]],
        ['matter'] + [fmt(x, 8) for x in dm[:2]] + ['—'],
        ['matter + vector portal $V_{P,e}=0.2V_{CC}$', fmt(pp[0], 8) if pp else '?', '—', '—'],
    ], widths=[0.34, 0.22, 0.22, 0.22])
    d.p('**The Sun.** With the electron density $n_e(r)=245\\,N_A\\,e^{-10.54\\,r/R_\\odot}$ cm$^{-3}$ ($R_\\odot=6.957\\times10^8$ m, $N_A$ Avogadro’s '
        'number) the two-flavour evolution (effective potential $c_{13}^2V_{CC}$) is integrated from the centre to the surface; the result '
        'agrees with the adiabatic MSW formula')
    d.eq(r'P_{ee}=c_{13}^4\,\frac{1+\cos2\theta_{12}\cos2\theta_m(0)}2+s_{13}^4,\qquad \text{vacuum average: }c_{13}^4\Big(1-\frac{\sin^22\theta_{12}}2\Big)+s_{13}^4,')
    sun = g('Sun P_ee {E, numerical, adiabatic}')
    if len(sun) == 12:
        d.table(['$E$ (MeV)', '$P_{ee}$ (numerical)', '$P_{ee}$ (adiabatic formula)'],
                [[fmt(sun[3*i], 3), fmt(sun[3*i + 1], 8), fmt(sun[3*i + 2], 8)] for i in range(4)], widths=[0.2, 0.4, 0.4])
    d.p('(at 0.3 MeV the residual matter effect at the centre, $2EV/\\Delta m^2_{21}\\approx0.15$, still lowers $P_{ee}$ below the vacuum average by '
        'about 0.025).')
    sn = g('SN resonance densities H, L (g/cm^3)'); cnb = g('T_nu, n_nu')
    d.p('**Supernovae.** The MSW resonance densities $\\rho_{\\rm res}=\\Delta m^2\\cos2\\theta/(2EV_{CC}(1\\,\\text{g/cm}^3))$ for $E=20$ MeV, $Y_e=0.5$ '
        'are %s g/cm³ (H: $\\Delta m^2_{31}$, $\\theta_{13}$) and %s g/cm³ (L: $\\Delta m^2_{21}$, $\\theta_{12}$).' % (
            fmt(sn[0], 6) if sn else '?', fmt(sn[1], 6) if len(sn) > 1 else '?'))
    d.p('**Cosmic neutrino background.** $T_\\nu=(4/11)^{1/3}T_\\gamma$ with $T_\\gamma=2.7255$ K gives $T_\\nu=%s$ K and the number density per '
        'flavour (neutrinos and antineutrinos) $n_\\nu=\\frac{3\\zeta(3)}{2\\pi^2}T_\\nu^3=%s$ cm$^{-3}$.' % (
            fmt(cnb[0], 6) if cnb else '?', fmt(cnb[1], 6) if len(cnb) > 1 else '?'))
    d.p('**dirac16complex and neutrino flavour.** With the Yukawa portal the oscillation probabilities of the light states, computed with $N$ in '
        'place of $U$, have $\\sum_\\beta P(\\alpha\\to\\beta)<1$ (a zero-distance effect, compensated by the heavy states); for $\\mu\\to0$ they become '
        'the three-flavour vacuum probabilities above. The vector portal acts as a non-standard matter potential: $V_{P,e}=0.2\\,V_{CC}$ changes '
        '$P(\\nu_\\mu\\to\\nu_e)$ at DUNE from %s to %s.' % (fmt(pp[1], 6) if len(pp) > 1 else '?', fmt(pp[0], 6) if pp else '?'))

    # ==============================================================================================
    d.h1('12. Verification')
    d.p('Every statement above is a check of 7+1.nb (Wolfram Language 15.0.1) and of the independent Python program tools/ref_sm.py '
        '(numpy, scipy, sympy). The numbers printed by both programs:')
    rows = []
    rels = []
    for lab, a, b, note in compare(D):
        r = max_rel(a, b); rels.append((r, note))
        rows.append([lab, str(len(a)), fmt(r, 2) if r is not None else '—', note])
    d.table(['quantity', 'values', 'max. relative difference', 'note'], rows, widths=[0.42, 0.1, 0.2, 0.28], small=True)
    full = [r for r, note in rels if r is not None and not note]
    d.p('Apart from the matrices that the Python program prints rounded, the largest relative difference is %s.' % (fmt(max(full), 2) if full else '?'))
    d.p('The checks of sections 17–19 of the notebook:')
    d.table(['#', 'check', 'result', 'section'],
            [[str(i + 1), '`' + c[0].replace('`', "'") + '`', c[1], (c[3][:4] if c[3] else c[2][:3])] for i, c in enumerate(ck)],
            widths=[0.06, 0.74, 0.08, 0.12], small=True)
    d.p('The Python program’s boolean checks:')
    d.table(['#', 'check', 'result'], [[str(i + 1), '`' + k.replace('`', "'") + '`', v_] for i, (k, v_) in enumerate(rb)],
            widths=[0.06, 0.82, 0.12], small=True)

    d.h1('13. How to reproduce')
    d.ol(['`python3 tools/make_7p1_nb.py` writes 7+1.nb (unevaluated).',
          'Evaluate 7+1.nb in a fresh kernel of Wolfram 15.0.1 (Evaluation ▸ Quit Kernel, then Evaluation ▸ Evaluate Notebook; or `Get["D:/Developer/github/7+1/tools/run_eval.wl"]`) and save it; the last section summarizes all checks.',
          '`python3 tools/ref_sm.py` repeats the checks of this document independently (exit code 0 if all pass).',
          '`python3 tools/make_docs.py` regenerates this document, its symbol dictionary and its provenance record.'])

    d.h1('14. References')
    d.ul(['S. Navas et al. (Particle Data Group), *Review of Particle Physics*, Phys. Rev. D 110 (2024) 030001, and its 2025 update '
          '(pdg.lbl.gov) — $G_F$, $m_W$, $m_Z$, $m_H$, quark masses, the CKM parameters and $J$.',
          'I. Esteban, M. C. Gonzalez-Garcia, M. Maltoni, I. Martinez-Soler, J. P. Pinheiro, T. Schwetz, *NuFit-6.0: updated global '
          'analysis of three-flavor neutrino oscillations*, JHEP 12 (2024) 216, arXiv:2410.05380 — the PMNS parameters.',
          'Planck Collaboration, N. Aghanim et al., *Planck 2018 results. VI. Cosmological parameters*, Astron. Astrophys. 641 (2020) A6.',
          'N. Arkani-Hamed, S. Dimopoulos, G. Dvali, Phys. Lett. B 429 (1998) 263 — gravity in large extra dimensions, $G_4=G_8/V_\\perp$.',
          'L. Wolfenstein, Phys. Rev. D 17 (1978) 2369; S. P. Mikheyev, A. Yu. Smirnov, Sov. J. Nucl. Phys. 42 (1985) 913 — the MSW effect.',
          'J. N. Bahcall, *Neutrino Astrophysics*, Cambridge University Press (1989) — the exponential solar electron density.',
          'C. Jarlskog, Phys. Rev. Lett. 55 (1985) 1039 — the invariant $J$.',
          'E. Witten, Phys. Lett. B 117 (1982) 324 — the global $SU(2)$ anomaly.'])
    return d
