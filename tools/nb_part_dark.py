# Section 24 of 7+1.nb: can the field dirac16complex be the dark matter, and can it be the dark energy?  The analysis uses the
# accurate approximate solutions of sections 20-23 (the three gravity regimes) and the couplings of sections 8 and 15-19, together
# with present observational data.  Appended on 27 October 2026.  Executed by tools/make_7p1_nb.py, which provides section(),
# subsection(), text(), item(), inp() and B.

section('24. dirac16complex as dark matter and dark energy')
text([('About section 24. ', B),
      'Section 24 was appended on 27 October 2026. The notebook as it was before (sections 0–23, with the summary of all checks '
      'numbered 24; evaluated) is kept in the workspace as backups/7+1_2026-10-27_before_dark.nb (excluded from the repository by '
      '.gitignore).'])
text([('The question. ', B),
      'Can the field dirac16complex be the dark matter of the universe, and can it be the dark energy? This section answers on the '
      'basis of the accurate approximate solutions of sections 21–23 (the three gravity regimes), the couplings of sections 8 and 15–20, '
      'and present observational data. Every statement is marked DERIVED (an exact consequence of the field equations of the model, '
      'checked below), ESTIMATE (an order-of-magnitude statement whose assumptions are stated), or OPEN (not decided by the model). '
      'Physical parameters as in section 22: M8 = 10 TeV (κ = 8πG8 = M8⁻⁶), reduced M_Pl = 2.435 × 10¹⁸ GeV, V⊥ = M_Pl²/M8⁶, '
      'm_r = 10 GeV; the mass m of dirac16complex is free. Coordinates as everywhere: x0 time (the octonion unit), x1, x2, x3 ordinary '
      'space, x4 … x7 auxiliary space; η17 = diag(+1, −1, …, −1) is simply the squares of the octonion basis, and xAct declares the '
      'spacetime with the signature list {7, 1, 0} (seven space-like and one time-like dimension).'])

subsection('24.1 What the observations require')
text('Data used (as reported by the collaborations; the confidence levels are theirs):')
item('Planck 2018 (Astron. Astrophys. 641 (2020) A6): Ω_c h² = 0.120 ± 0.001, Ω_b h² = 0.0224 ± 0.0001, H0 = 67.4 ± 0.5 km/s/Mpc, '
     'Ω_m = 0.315 ± 0.007.')
item('DESI DR2 (DESI Collaboration, arXiv:2503.14738): the evolving dark energy w(a) = w0 + wa (1 − a) with w0 > −1, wa < 0 is preferred '
     'over a cosmological constant at 3.1σ (DESI BAO + CMB) and at 2.8σ to 4.2σ when supernovae are added; with Planck PR4 + ACT DR6 + '
     'DES-Y5 supernovae, w0 = −0.753 ± 0.056, wa = −0.84 (+0.23, −0.20), 4.2σ (Garcia-Quintero et al., arXiv:2504.18464).')
item('LZ, 4.2 tonne-years (LZ Collaboration, arXiv:2410.17036): spin-independent WIMP–nucleon cross section below 2.2 × 10⁻⁴⁸ cm² at '
     '40 GeV (90% CL).')
item('Bullet Cluster (Randall et al., ApJ 679 (2008) 1173): self-interaction cross section per mass σ/m < 1.25 cm²/g (68% CL), '
     '< 0.7 cm²/g with equal mass-to-light ratios.')
item('Lyman-α forest (Iršič et al., Phys. Rev. D 109 (2024) 043511): the mass of a thermal-relic warm dark-matter particle exceeds '
     '5.7 keV (95% CL). Dark matter produced by decays of particles in equilibrium is warmer than a thermal relic of the same mass and '
     'its bound is higher (about 7 keV where the thermal-relic bound is 3 keV: Ballesteros, Garcia and Pierre, JCAP 03 (2021) 101).')
item('Pauli principle in dwarf spheroidal galaxies (Alvey et al., MNRAS 501 (2021) 1188): a fermion with g = 2 internal states is heavier '
     'than 0.13 keV (95% CL); the bound scales as g^(−1/4).')
item('X-ray lines (NuSTAR, Roach et al., Phys. Rev. D 107 (2023) 023009): limits on the mixing of sterile neutrinos of 6–40 keV with the '
     'active ones, from the line of ν_s → νγ.')
item('Large extra dimensions (Particle Data Group 2024): with M_Pl² = R^n M_D^(n+2) and V = (2πR)^n, colliders exclude M_D below 5.9 to '
     '11.2 TeV for n = 6 to 2 (ATLAS, 139 fb⁻¹).')
text('The dark matter must have the comoving abundance Y = n/s with')
text('      m Y = Ω_c h² (ρ_c/h²)/s0,   ρ_c/h² = 3 (100 km/s/Mpc)² M_Pl²,   s0 = (2π²/45) g*s0 T_CMB³,   g*s0 = 2 + (21/4)(4/11),')
text('and in the convention of the Particle Data Group the 7+1 model has M_D = (2π)^(4/6) M8 (V⊥ = (2πR)⁴, M_Pl² = V⊥ M8⁶ = R⁴ M_D⁶). '
     'The lightest Kaluza–Klein level of dirac16complex is four brane Dirac spinors (section 18), i.e. g = 16 internal states (two spin '
     'states, particle and antiparticle, four spinors): its Pauli bound is 0.13 keV × (2/16)^(1/4).')
inp(r'''
dkOch2 = 0.120; dkh = 0.674; dkOc = dkOch2/dkh^2; dkHbarC = 1.973269804*^-14;          (* hbar c in GeV cm *)
dkRhoc = 3 (100./rIImpc rIIhbar)^2 rIIMPl^2/dkHbarC^3;                                   (* rho_c/h^2 in GeV/cm^3 *)
dkGs0 = 2 + 21/4 (4/11); dks0 = 2 Pi^2/45 dkGs0 rIITcmb^3/dkHbarC^3;                    (* entropy density today, 1/cm^3 *)
dkYm = dkOch2 dkRhoc/dks0;                                                               (* m Y in GeV *)
dkMD = (2 Pi)^(2/3) rIIM8; dkPauli = 0.13 (2/16)^(1/4);                                  (* GeV; keV *)
check["required abundance: m Y = Omega_c h^2 (rho_c/h^2)/s0 = 4.37 10^-10 GeV with rho_c/h^2 = 1.053 10^-5 GeV/cm^3 and s0 = 2891 cm^-3 (to 10^-3)",
  Abs[dkRhoc/1.0534*^-5 - 1] < 10^-3 && Abs[dks0/2891.2 - 1] < 10^-3 && Abs[dkYm/4.372*^-10 - 1] < 10^-3]
check["collider bound: with M_Pl^2 = R^n M_D^(n+2) and V_perp = (2 Pi R)^4 = M_Pl^2/M8^6, M_D = (2 Pi)^(2/3) M8 = 34 TeV, above the strongest bound 11.2 TeV",
  Simplify[((2 Pi)^(2/3) m8S)^6 rRS^4 - m8S^6 (2 Pi rRS)^4] === 0 && dkMD > 11.2*^3]
{"Omega_c, rho_c/h^2 (GeV/cm^3), s0 (1/cm^3), m Y (GeV)" -> {dkOc, dkRhoc, dks0, dkYm}, "M_D (GeV)" -> dkMD,
 "Pauli bound for g = 16 internal states (keV)" -> dkPauli}
''')

subsection('24.2 dirac16complex as cold matter in the solutions of sections 21–23 (DERIVED)')
text('A massive dirac16complex particle with comoving momenta k_a (ordinary space) and k_b (auxiliary space) has the comoving energy '
     'ω = (m² + k_a²/a² + k_b²/b²)^(1/2); the stress tensor of section 20.2 (W = ω) gives εV = ω, p_a V = k_a²/(3a²ω), '
     'p_b V = k_b²/(4b²ω), with V = a³b⁴. The lightest mode has k_b = 0; once k_a/a ≪ m it is pressureless matter, ε ∝ a⁻³ at fixed b. '
     'Every other Kaluza–Klein level is matter of mass (m² + k_b²/b0²)^(1/2) with a static pressure p_b along the auxiliary space. '
     'Consequences in the three regimes:')
item('Regime I: a bulk gas of dirac16complex dominates near x0 = 0 (a ∝ τ^(3/4), section 21.3) and, without the stabilizing potential, '
     'decompactifies the auxiliary space (τHb → 2/5 in the gas-dominated Einstein regime). dirac16complex as dark matter therefore '
     'REQUIRES the stabilizer of section 22.1.')
item('Regime II: after the capture of the auxiliary space the frozen auxiliary momenta are Kaluza–Klein masses and the gas is '
     'four-dimensional matter, its share of the energy constant (section 22.1). Matter of this kind is exactly the term Ω_c a⁻³ of the '
     'late-universe solution of section 22.4, so dirac16complex CAN BE the Ω_c of that solution without changing it. Through the radion '
     'equation of section 20.3 (κ(ε/6 − p_a/2 + p_b/3) with p_a = p_b = 0) it displaces the radion by δβ = κε/(6 m_r²) = H0² Ω_c/(2 m_r²) '
     '(κε = ρ_DM/M_Pl² = 3H0²Ω_c): G4 = G8/V⊥ does not change.')
item('Regime III: with a(t) = (Ω_m/Ω_Λ)^(1/3) sinh^(2/3)(3H∞t/2) (section 23.1), ρ_DM/ρ_Λ = Ω_c a⁻³/Ω_Λ = (Ω_c/Ω_m)/sinh²(3H∞t/2) → '
     '4 (Ω_c/Ω_m) e^(−3H∞t): the dark matter, like all matter, dilutes away.')
inp(r'''
dkSt = feStress[Function[{aa, rIIIbb}, Sqrt[mDM^2 + kaP^2/aa^2 + kbP^2/rIIIbb^2]]];
check["a massive dirac16complex particle: eps V == omega, pa V == k_a^2/(3 a^2 omega), pb V == k_b^2/(4 b^2 omega), omega = Sqrt[m^2 + k_a^2/a^2 + k_b^2/b^2]",
  With[{om = Sqrt[mDM^2 + kaP^2/aS^2 + kbP^2/bS^2]}, Simplify[dkSt aS^3 bS^4 - {om, kaP^2/(3 aS^2 om), kbP^2/(4 bS^2 om)}, aS > 0 && bS > 0 && mDM > 0] === {0, 0, 0}]]
check["the lightest mode (k_b = 0) is cold matter: pa/eps == k_a^2/(3 (m^2 a^2 + k_a^2)) -> 0, pb == 0, and eps ~ a^-3 at fixed b for k_a = 0",
  Simplify[(dkSt[[2]]/dkSt[[1]] /. kbP -> 0) - kaP^2/(3 (mDM^2 aS^2 + kaP^2)), aS > 0 && mDM > 0] === 0 && (dkSt[[3]] /. kbP -> 0) === 0 &&
  Simplify[aS D[Log[dkSt[[1]] /. {kaP -> 0, kbP -> 0}], aS], aS > 0 && mDM > 0] === -3]
dkDbeta = rIIH0g^2 dkOc/(2 100.);
check["dark matter displaces the radion today by H0^2 Omega_c/(2 m_r^2) < 10^-80 (m_r = 10 GeV): the four-dimensional Newton constant is constant",
  dkDbeta < 10^-80]
check["regime III: Omega_m a^-3/Omega_Lambda == Csch[3 H t/2]^2 for a = (Omega_m/Omega_Lambda)^(1/3) Sinh[3 H t/2]^(2/3), and -> 4 Exp[-3 H t]",
  Simplify[(omS/olS)/((omS/olS)^(1/3) Sinh[3/2 hS tS]^(2/3))^3 - Csch[3/2 hS tS]^2, omS > 0 && olS > 0 && hS > 0 && tS > 0] === 0 &&
  Limit[Csch[zS]^2 Exp[2 zS], zS -> Infinity] == 4]
{"radion displacement by dark matter today" -> dkDbeta, "gas share after capture (section 22.1) at tau = 10^5, 3 10^5" -> rIIfrVals,
 "without stabilizer (section 21.3): tau Hb at Ha = 10^-14" -> rIlate["tau"] rIlate["Hb"]}
''')

subsection('24.3 The lightest mode: darkness, collisionlessness and stability (ESTIMATE)')
text('The lightest mode is constant over the auxiliary torus: ψ(x, y) = χ(x)/√V⊥ with canonically normalized four-dimensional fields '
     'χ (four brane Dirac spinors, section 18). Every coupling at the brane or in the bulk is therefore diluted by the volume. With '
     'couplings of natural size at the eight-dimensional Planck scale — Higgs portal λ_P = c_P/M8⁵, Yukawa portal y = c_y/M8², vector '
     'portal G_P = c_V/M8⁶, self-couplings g_r = c_r/M8⁶ (mass dimensions −5, −2, −6, −6) — and V⊥ = M_Pl²/M8⁶ the four-dimensional '
     'couplings of χ are')
text('      λ_eff = λ_P/V⊥ = c_P M8/M_Pl²,   y_eff = y/√V⊥ = c_y M8/M_Pl,   G_eff = G_P/V⊥ = c_V/M_Pl²,   G_r = g_r/V⊥ = c_r/M_Pl².')
text('Direct detection (Higgs exchange; the operator λ_eff Φ†Φ χ̄χ contains λ_eff v h χ̄χ, the nucleon couples to h with f_N m_N/v, '
     'f_N ≈ 0.30): σ_SI = λ_eff² f_N² m_N² μ²/(π m_h⁴), μ = m m_N/(m + m_N). Self-interaction (contact): σ = G_r² m²/(4π). Stability: '
     'through the Yukawa portal χ mixes with the neutrinos by θ = μ_P/m, μ_P = y_eff v/√2 (the portal Dirac mass of section 18.2), and '
     'decays like a sterile neutrino:')
text('      Γ(χ → νγ) = 9 α G_F² m⁵ sin²2θ/(1024 π⁴),     Γ(χ → 3ν) = G_F² m⁵ θ²/(96 π³)     (sin²2θ ≈ 4θ²),')
text('the rates of a Majorana fermion (a Dirac fermion decays half as fast; immaterial for these estimates); Γ(χ → νγ)/Γ(χ → 3ν) = '
     '27α/(8π). The Higgs and vector portals are bilinear in χ and cannot make it decay; without the Yukawa portal χ is absolutely '
     'stable. A lifetime longer than the age of the universe is the minimum requirement; searches for the X-ray line of χ → νγ demand '
     'far longer lifetimes and only strengthen the conclusions of section 24.4.')
inp(r'''
check["couplings of the lightest mode: with V_perp = M_Pl^2/M8^6, lambda_P/V_perp == c_P M8/M_Pl^2, y/Sqrt[V_perp] == c_y M8/M_Pl, G_P/V_perp == c_V/M_Pl^2, g_r/V_perp == c_r/M_Pl^2",
  With[{vp = mPlS^2/m8S^6}, Simplify[{(cPS/m8S^5)/vp - cPS m8S/mPlS^2, (cyS/m8S^2)/Sqrt[vp] - cyS m8S/mPlS, (cVS/m8S^6)/vp - cVS/mPlS^2,
     (crS/m8S^6)/vp - crS/mPlS^2}, mPlS > 0 && m8S > 0] === {0, 0, 0, 0}]]
dkLam = rIIM8/rIIMPl^2; dkYeff = rIIM8/rIIMPl; dkGr = 1/rIIMPl^2; dkGeV2g = 1.78266192*^-24;     (* per unit c_P, c_y, c_r; grams per GeV *)
dkSigSI[cpx_, mx_] := With[{mu = mx 0.9383/(mx + 0.9383)}, (cpx dkLam)^2 0.30^2 0.9383^2 mu^2/(Pi mHb^4) dkHbarC^2];    (* cm^2 *)
dkSigM[crx_, mx_] := (crx dkGr)^2 mx/(4 Pi) dkHbarC^2/dkGeV2g;                                                      (* cm^2/g *)
dkCrMax[mx_] := Sqrt[4 Pi 0.7 dkGeV2g/(dkHbarC^2 mx)]/dkGr;
dkMuP[cyx_] := cyx dkYeff vEW/Sqrt[2];
dkTauG[cyx_, mx_] := rIIhbar/(9/137.035999 gFermi^2 mx^5 4 (dkMuP[cyx]/mx)^2/(1024 Pi^4));
dkTau3[cyx_, mx_] := rIIhbar/(gFermi^2 mx^5 (dkMuP[cyx]/mx)^2/(96 Pi^3));
dkAge = rIIage rIIgyr;
dkMmax = (dkTau3[1, 1.]/dkAge)^(1/3);                            (* tau(chi -> 3 nu) = tau(m = 1 GeV)/(c_y^2 m^3) *)
check["natural couplings make the lightest mode dark and collisionless: sigma_SI(c_P = 1, m = 40 GeV) < 10^-50 x (LZ limit 2.2 10^-48 cm^2), sigma/m(c_r = 1, m = 1 GeV) < 10^-70 cm^2/g (Bullet Cluster: 0.7 cm^2/g)",
  dkSigSI[1, 40.] < 10^-50 2.2*^-48 && dkSigM[1, 1.] < 10^-70]
check["stability through the Yukawa portal: for c_y = 1 the lifetime chi -> 3 nu equals the age of the universe at m = 40 MeV; it scales as 1/(c_y^2 m^3), and Gamma(nu gamma)/Gamma(3 nu) == 27 alpha/(8 Pi)",
  0.03 < dkMmax < 0.05 && Abs[dkTau3[1, dkMmax]/dkAge - 1] < 10^-12 && Abs[dkTau3[2, 0.01]/dkTau3[1, 0.01] - 1/4] < 10^-12 &&
  Abs[dkTau3[1, 0.02]/dkTau3[1, 0.01] - 1/8] < 10^-12 &&
  Abs[dkTau3[1, 0.01]/dkTauG[1, 0.01] - 27/(8 Pi 137.035999)] < 10^-12]
{"lambda_eff, y_eff, G_r per unit coupling (GeV^-1, 1, GeV^-2)" -> {dkLam, dkYeff, dkGr}, "sigma_SI (cm^2) at c_P = 1, m = 40 GeV" -> dkSigSI[1, 40.],
 "sigma/m (cm^2/g) at c_r = 1, m = 1 GeV; c_r allowed by the Bullet Cluster at m = 1 GeV" -> {dkSigM[1, 1.], dkCrMax[1.]},
 "portal mass mu_P (GeV) at c_y = 1; lifetimes chi -> nu gamma, chi -> 3 nu (s) at c_y = 1, m = 1 GeV and 10 keV" -> {dkMuP[1], dkTauG[1, 1.], dkTau3[1, 1.], dkTauG[1, 10^-5], dkTau3[1, 10^-5]},
 "largest m (GeV) with lifetime > age for c_y = 1" -> dkMmax}
''')

subsection('24.4 The abundance (ESTIMATE; the outcome is OPEN)')
text([('From the solutions of sections 21–22: overproduction. ', B),
      'The minimal model makes far too much dirac16complex, not too little:'])
item('The primordial gas (sections 21.3 and 22.1). After capture the gas carries the fraction f = 0.0087 of the energy and the radion the '
     'rest (for the normalization of section 21.3; f scales with it). Even if the radion energy all became brane radiation at T_RH, the '
     'gas-to-radiation ratio f/(1 − f) would have to equal ρ_DM/ρ_rad = (4/3) m Y/T_RH (ρ_rad/s = 3T/4); it exceeds it by '
     'O = (f/(1 − f)) T_RH/((4/3) m Y) = 1.5 × 10⁷ (T_RH/GeV), i.e. 6 × 10⁴ even for the lowest T_RH = 4 MeV allowed by nucleosynthesis.')
item('Radion decay (section 22.2). The radion does NOT become brane radiation: it decays 1.3 × 10¹³ times faster into Kaluza–Klein pairs '
     'of dirac16complex (if 2m < m_r) than into brane fields. Without a later dilution the universe after the radion decay is made of '
     'dirac16complex, with a brane share of only Γ_brane/Γ_bulk ≈ 8 × 10⁻¹⁴ — the opposite of the observed ratio.')
text('Both results mean that dirac16complex as dark matter needs the dilution and reheating that the moduli problem already needs '
     '(section 22.2), at T_RH below the normalcy temperature T* = (M8⁶/M_Pl)^(1/5) = 13.3 GeV. After it the abundance is set by the '
     'portals. With instantaneous reheating, production by decays of Higgs bosons in equilibrium below T_RH carries the factor')
text('      S(x) = (2/(3π)) ∫_x^∞ t³ K_1(t) dt,   x = m_h/T_RH   (S(0) = 1 because ∫_0^∞ t³ K_1(t) dt = 3π/2;  S = 0.0078 at T_RH = T*),')
text('and the freeze-in yield of a particle made in h decays with partial width Γ is Y = (135/(8π³·1.66)) (M_Pl^(unreduced)/(g*^(3/2) m_h²)) Γ S '
     '(Hall, Jedamzik, March-Russell and West, JHEP 03 (2010) 080; M_Pl^(unreduced) = 1.221 × 10¹⁹ GeV, g* = 106.75; below the electroweak '
     'bosons g* = 86.25, so for T_RH < m_W this understates Y by at most (106.75/86.25)^(3/2) = 1.4). Production by other processes adds to '
     'this; the estimates below are therefore lower bounds on what each portal makes, except for the zero mode alone (an upper bound up to '
     'factors of order one, S = 1).')
item('Yukawa portal, the zero mode. Γ(h → ν̄χ) + Γ(h → νχ̄) = y_eff² m_h/(16π) per channel (N channels: three lepton flavours and the '
     'brane spinors, N ≤ 12): Y = 9.7 × 10⁻¹⁸ c_y² N S. The abundance grows like c_y² m while the lifetime τ(χ → 3ν) = K/(c_y² m³) falls; '
     'both conditions hold together only for m < m* = (Y(c_y = 1) K/(m Y t0))^(1/2) = 1.2 keV (N = 1), 4.2 keV (N = 12), even with S = 1 '
     '— below the Lyman-α bound.')
item('Yukawa portal, the Kaluza–Klein tower. Every Kaluza–Klein level couples to the brane with the same strength as the zero mode '
     '(section 18.2), so the Higgs decays make the whole tower below m_h: the yield is multiplied by Σ_n (1 − μ_n²/m_h²)² = '
     'V⊥ m_h⁴/(192π²) = 7.7 × 10¹⁷ (continuum limit of the lattice sum; for the lattice sum with r_4(k) see the check below), with mean '
     'mass (192/315) m_h = 76 GeV: Y_tower = 7.5 c_y² N S and a mass per entropy 571 c_y² N S GeV. At T_RH = T* the tower carries more '
     'energy than all the dark matter unless c_y < 9.9 × 10⁻⁶ (N = 1); and through the same mixing its states decay into neutrinos and '
     'photons long after nucleosynthesis (τ(χ_n → 3ν) = 2.8 × 10¹⁰ s/c_y² at μ_n = 10 GeV; other channels shorten it by about ten). The '
     'Yukawa portal must therefore be very weak or absent; the precise bound from nucleosynthesis and the CMB is not derived here (OPEN). '
     'With c_y < 10⁻⁵ the zero mode lives longer than the universe for m < 86 GeV, and its freeze-in yield is negligible.')
item('Higgs portal. The zero mode gets Y = 2.0 × 10⁻⁴⁹ c_P² S (Γ(h → χχ̄) = (λ_eff v)² m_h/(8π)). But the brane breaks the conservation '
     'of Kaluza–Klein momentum, so λ_eff v h χ̄_n χ_m couples EVERY pair of levels: Γ(h → χ_n χ̄_m) = (λ_eff v)² m_h (1 − (u_n + u_m)²)^(3/2) '
     '(1 − (u_n − u_m)²)^(1/2)/(8π), u = μ/m_h. With the density of levels dN = V⊥ μ³ dμ/(8π²) the sum over all pairs is '
     '(V⊥ m_h⁴/(8π²))² = 3.4 × 10³⁸ times I_hp = ∫∫_(u1+u2<1) u1³ u2³ (1 − (u1 + u2)²)^(3/2) (1 − (u1 − u2)²)^(1/2) du1 du2 = 9.5 × 10⁻⁵ '
     '(with the mass weight u1 + u2 under the integral: I_hp^E). For four brane spinors the pairs carry 2.5 × 10⁻¹² c_P² S GeV per '
     'entropy (mean pair mass I_hp^E/I_hp = 0.77 m_h): the observed abundance needs c_P = 13 for T_RH ≫ m_h (S → 1) and c_P = 150 at '
     'T_RH = T* — two orders of magnitude above natural size; for c_P = 1 it is 4.5 × 10⁻⁵ of the dark matter. Such dark matter would be '
     'a population of heavy Kaluza–Klein modes (tens of GeV), whose decays into lighter levels and Standard Model particles through the '
     'same portal are not computed here (OPEN).')
item('Vector portal (and any coupling of gravitational strength per pair, G_eff = c_V/M_Pl²). Summed over the pairs lighter than T the '
     'production rate divided by H is of order c_V² (T/T_V)¹¹ with T_V = (M8¹²/M_Pl)^(1/11) = 490 GeV (numerical factors omitted; they raise '
     'T_V): 5 × 10⁻¹⁸ at T*. Negligible.')
text('Hence, with couplings of natural size, the minimal model does not produce the observed abundance after the dilution its moduli '
     'problem requires: the Yukawa portal must be tiny, the vector portal is negligible, and the Higgs portal falls short by a factor '
     '2 × 10⁴ (it would need c_P ≈ 150). The abundance must come from physics outside the minimal model (OPEN).')
inp(r'''
dkF = rIIfrVals[[2]]; dkOver[trh_] := (dkF/(1 - dkF)) trh/(4/3 dkYm);
dkFracBulk = rIIgBulk/(rIIgBulk + rIIgam);
check["primordial gas: overabundance O = (f/(1 - f)) T_RH/((4/3) m Y) = 1.5 10^7 T_RH/GeV (6 10^4 at T_RH = 4 MeV)",
  Abs[dkOver[1.]/1.51*^7 - 1] < 0.01 && 5*^4 < dkOver[0.004] < 7*^4]
check["the radion puts more than 1 - 10^-12 of its decays into the bulk (Kaluza-Klein pairs of dirac16complex)", 1 - dkFracBulk < 10^-12]
{"gas share after capture f, overabundance O at T_RH = 1 GeV and 4 MeV" -> {dkF, dkOver[1.], dkOver[0.004]}, "fraction of radion decays into the bulk; brane share" -> {dkFracBulk, rIIgam/rIIgBulk}}
''')
inp(r'''
dkTstar = (rIIM8^6/rIIMPl)^(1/5);
dkS[x_?NumericQ] := 2/(3 Pi) NIntegrate[tV^3 BesselK[1, tV], {tV, x, Infinity}, PrecisionGoal -> 12, AccuracyGoal -> 30];
dkSstar = dkS[mHb/dkTstar];
check["Boltzmann factor of production below T_RH: Int_0^Infinity t^3 K_1(t) dt == 3 Pi/2 (exact), S(0) = 1, S(m_h/T*) < 10^-2 with T* = (M8^6/M_Pl)^(1/5) = 13.3 GeV",
  Integrate[tV^3 BesselK[1, tV], {tV, 0, Infinity}] == 3 Pi/2 && Abs[dkS[0] - 1] < 10^-10 && dkSstar < 10^-2 && 13 < dkTstar < 14]
dkC = 135/(8 Pi^3 1.66) 1/106.75^(3/2) 1.22091*^19;            (* Y = dkC Gamma S/m_h^2 *)
dkYfi[cyx_, nch_] := dkC nch ((cyx dkYeff)^2 mHb/(16 Pi))/mHb^2;   (* zero mode, Yukawa portal, S = 1 *)
dkYhp[cpx_] := dkC ((cpx dkLam vEW)^2 mHb/(8 Pi))/mHb^2;           (* zero mode, Higgs portal, S = 1 *)
dkK = dkTau3[1, 1.];                                               (* tau(chi -> 3 nu) = dkK/(c_y^2 m^3), m in GeV *)
dkMstar[nch_] := Sqrt[dkYfi[1, nch] dkK/(dkYm dkAge)];
dkCyNeed[mx_, nch_] := Sqrt[dkYm/(mx dkYfi[1, nch])]; dkCyMax[mx_] := Sqrt[dkK/(mx^3 dkAge)];
check["zero mode through the Yukawa portal: the abundance (~ c_y^2 m) and the stability (lifetime ~ 1/(c_y^2 m^3) > age) are compatible only for m < m*, and m* < 5.7 keV (Lyman-alpha) for N = 1 and N = 12 channels even with S = 1",
  dkMstar[1] < 5.7*^-6 && dkMstar[12] < 5.7*^-6 && Abs[dkCyNeed[dkMstar[1], 1]/dkCyMax[dkMstar[1]] - 1] < 10^-9]
dkTowN = rIIVp mHb^4/(192 Pi^2); dkTowE = rIIVp mHb^5/(315 Pi^2);
dkLat = 1 + Total[rIIr4[[1 ;; 10^4]] (1 - Range[10^4]/10.^4)^2]; dkLatE = Total[rIIr4[[1 ;; 10^4]] Sqrt[Range[10^4]/10.^4] (1 - Range[10^4]/10.^4)^2];
check["Kaluza-Klein tower: Sum_n (1 - mu_n^2/m_h^2)^2 -> V_perp m_h^4/(192 Pi^2) and Sum_n mu_n (1 - mu_n^2/m_h^2)^2 -> V_perp m_h^5/(315 Pi^2) (d^4n = 2 Pi^2 |n|^3 d|n|: exact integrals); the lattice sums with r_4(k) for R = 100 equal Pi^2 R^4/12 and 16 Pi^2 R^4/315 to 10^-8",
  Integrate[2 Pi^2 uV^3 (1 - uV^2)^2, {uV, 0, 1}] == Pi^2/12 && Integrate[2 Pi^2 uV^4 (1 - uV^2)^2, {uV, 0, 1}] == 16 Pi^2/315 &&
  Simplify[(Pi^2/12) (mhS rLS/(2 Pi))^4 - rLS^4 mhS^4/(192 Pi^2)] === 0 && Simplify[(16 Pi^2/315) mhS (mhS rLS/(2 Pi))^4 - rLS^4 mhS^5/(315 Pi^2)] === 0 &&
  Abs[dkLat/(Pi^2 10.^8/12) - 1] < 10^-8 && Abs[dkLatE/(16 Pi^2 10.^8/315) - 1] < 10^-8]
dkCyTow = Sqrt[dkYm/(dkYfi[1, 1] dkTowE dkSstar)];
check["Kaluza-Klein tower through the Yukawa portal: 7.7 10^17 states below m_h with mean mass (192/315) m_h; at T_RH = T* it carries more energy than all the dark matter unless c_y < 2 10^-5 (N = 1)",
  Abs[dkTowN/7.688*^17 - 1] < 10^-3 && Abs[dkTowE/(dkTowN 192/315 mHb) - 1] < 10^-12 && dkCyTow < 2*^-5]
dkIhp = NIntegrate[u1V^3 u2V^3 (1 - (u1V + u2V)^2)^(3/2) Sqrt[1 - (u1V - u2V)^2], {u1V, 0, 1}, {u2V, 0, 1 - u1V}, PrecisionGoal -> 10, AccuracyGoal -> 20];
dkIhpE = NIntegrate[(u1V + u2V) u1V^3 u2V^3 (1 - (u1V + u2V)^2)^(3/2) Sqrt[1 - (u1V - u2V)^2], {u1V, 0, 1}, {u2V, 0, 1 - u1V}, PrecisionGoal -> 10, AccuracyGoal -> 20];
dkPairs = (rIIVp mHb^4/(8 Pi^2))^2;
dkRhp[cpx_, nch_, x_] := dkYhp[cpx] nch dkPairs mHb dkIhpE dkS[x];        (* mass of the Kaluza-Klein pairs per entropy, GeV *)
dkCpNeed[x_] := Sqrt[dkYm/dkRhp[1, 4, x]];
check["Higgs portal: the zero mode gets Y < 10^-40 (c_P = 1); the pairs of Kaluza-Klein levels (lambda(1, u1^2, u2^2) == (1 - (u1 + u2)^2)(1 - (u1 - u2)^2)) give the observed abundance only for c_P between 10 and 20 (T_RH >> m_h, S -> 1) or between 100 and 200 (T_RH = T*)",
  dkYhp[1] < 10^-40 && Expand[(1 - u1V^2 - u2V^2)^2 - 4 u1V^2 u2V^2 - (1 - (u1V + u2V)^2) (1 - (u1V - u2V)^2)] === 0 &&
  10 < dkCpNeed[0] < 20 && 100 < dkCpNeed[mHb/dkTstar] < 200]
dkTV = (rIIM8^12/rIIMPl)^(1/11);
check["vector portal: gravitational strength per pair; the production rate of all pairs lighter than T divided by H ~ (T/T_V)^11 with T_V = (M8^12/M_Pl)^(1/11) = 490 GeV is below 10^-17 at T*",
  480 < dkTV < 500 && (dkTstar/dkTV)^11 < 10^-17]
{"T* (GeV), S(m_h/T*)" -> {dkTstar, dkSstar}, "zero mode: Yukawa-portal Y at c_y = 1 (N = 1, S = 1), Higgs-portal Y at c_P = 1 (S = 1)" -> {dkYfi[1, 1], dkYhp[1]},
 "m* (GeV) for N = 1 and N = 12" -> {dkMstar[1], dkMstar[12]}, "c_y needed and allowed at m = 10 keV, 10 MeV" -> {dkCyNeed[10^-5, 1], dkCyMax[10^-5], dkCyNeed[0.01, 1], dkCyMax[0.01]},
 "Kaluza-Klein tower below m_h: number, mass sum (GeV), Y_tower and mass per entropy (GeV) at c_y = 1 (N = 1, S = 1); c_y bound at T*" -> {dkTowN, dkTowE, dkYfi[1, 1] dkTowN, dkYfi[1, 1] dkTowE, dkCyTow},
 "tau(chi_n -> 3 nu) (s) at mu_n = 10 GeV, c_y = 1; largest m (GeV) with tau(chi -> 3 nu) > age for c_y = 10^-5" -> {dkTau3[1, 10.], (dkK/(10^-10 dkAge))^(1/3)},
 "Higgs portal pairs: I_hp, I_hp^E, number of pairs, mass per entropy (GeV) at c_P = 1 (N = 4, S = 1); c_P needed for S = 1 and at T*" -> {dkIhp, dkIhpE, dkPairs, dkRhp[1, 4, 0], dkCpNeed[0], dkCpNeed[mHb/dkTstar]},
 "Higgs portal pairs at c_P = 1, T_RH = T*: fraction of the dark matter" -> dkRhp[1, 4, mHb/dkTstar]/dkYm,
 "vector portal: T_V (GeV), (T*/T_V)^11" -> {dkTV, (dkTstar/dkTV)^11},
 "T_RH (GeV) = 5, 8, 10, T*: S(m_h/T_RH), c_y bound of the tower (N = 1), c_P needed by the Higgs-portal pairs (N = 4)" ->
   Table[{trh, dkS[mHb/trh], Sqrt[dkYm/(dkYfi[1, 1] dkTowE dkS[mHb/trh])], dkCpNeed[mHb/trh]}, {trh, {5., 8., 10., dkTstar}}]}
''')

subsection('24.5 dirac16complex and dark energy')
text([('Matter cannot be dark energy (DERIVED). ', B),
      'Free particles have p_a = −∂_aW/(3a²b⁴) ≥ 0 because W decreases with a. The interacting Kohn–Sham matter of section 8 '
      '(units m = 1, gS = −10⁵, gV = 10⁴) has, at zero temperature, p = μρ − ε with μ = ε_F (Hugenholtz–Van Hove, section 8.7), so '
      'p = ρ² d(ε/ρ)/dρ: the pressure is negative only where the energy per particle decreases with density — below the saturation '
      'density, where uniform matter has a higher energy than droplets at saturation in empty space (both have p = 0) and is metastable or, '
      'in the spinodal region dp/dρ < 0, unstable. Its most negative equation of state there is w = p/ε ≈ −0.058 (at kF ≈ 0.67): nowhere '
      'near −1.'])
text([('A condensate is only a cosmological constant (DERIVED). ', B),
      'A Lorentz-invariant condensate of dirac16complex (⟨ψbar ψ⟩ ≠ 0, which the self-couplings of section 15 allow) has T_μν ∝ g_μν, '
      'i.e. w = −1 exactly: a contribution to Λ4, not a dynamical dark energy.'])
text([('Vacuum energy is only a contribution to Λ4 (DERIVED, with an ESTIMATE of its size). ', B),
      'The zero-point energy of dirac16complex on the auxiliary torus (periodic boundary conditions) is a Casimir energy. For a field '
      'of mass m in eight dimensions with four compact directions of proper size ℓ⊥ = V⊥^(1/4) the four-dimensional vacuum energy per '
      'real bosonic degree of freedom is −(2/ℓ⊥⁴) Σ_(n ∈ Z⁴, n ≠ 0) (mℓ⊥/(2π|n|))⁴ K_4(mℓ⊥|n|) → −Γ(4) Z/(π⁴ ℓ⊥⁴) for m → 0, with the Epstein sum '
      'Z = Σ_(n ≠ 0) |n|⁻⁸ = 8(1 − 4⁻³) ζ(4) ζ(3) (from r_4(k) of section 22.2); a fermion contributes with the opposite sign, and dirac16complex has '
      '16 real degrees of freedom on shell. Massless it gives +6.8 × 10³⁴ ρ_Λ; it equals ρ_Λ only for mℓ⊥ = 92.2 (m = 59 MeV), and the graviton '
      '(20 degrees of freedom) contributes −8.5 × 10³⁴ ρ_Λ. Every such term is constant while b is frozen (w = −1 up to the radion '
      'displacement, 10⁻⁸⁶), and they add to Λ + κU0 + κf_W/V⊥: the observed Λ4 needs their sum to cancel to 10⁻³⁵ of each — the '
      'cosmological-constant problem, not solved by dirac16complex.'])
text([('The prediction w = −1. ', B),
      'With the radion frozen the model predicts w = −1 exactly; the DESI preference for w0 > −1, wa < 0 (2.8σ to 4.2σ), if confirmed, would '
      'contradict it whatever dirac16complex is. (The radion, the only scalar of the model, cannot be the dark energy either: with '
      'm_r ≫ H0 it is frozen; a radion light enough to roll today would change G4 = G8 e^(−4β)/V⊥ by order one — ESTIMATE.)'])
inp(r'''
dkPE[kf_?NumericQ] := With[{e = eps7[gS0, gV0, kf]}, {fermi7[gS0, gV0, kf] rho7[kf] - e, e}];   (* {p, eps}: p = mu rho - eps with mu = eps_F *)
dkPres[kf_?NumericQ] := dkPE[kf][[1]];
dkWtab = Table[With[{pe = dkPE[kf]}, {kf, pe[[1]], pe[[1]]/pe[[2]]}], {kf, 5/100, 85/100, 1/100}];
dkWmin = Min[dkWtab[[All, 3]]]; dkKmin = dkWtab[[Ordering[dkWtab[[All, 3]], 1][[1]], 1]];
check["Kohn-Sham matter of section 8: p > 0 at low density, p < 0 below saturation, p = 0 at kF0, a spinodal region (dp/dkF < 0) in between, and w = p/eps >= -0.06 everywhere",
  dkWtab[[1, 2]] > 0 && dkPres[N[kF0, 40] - 1/100] < 0 && Abs[dkPres[kF0]] < 10^-20 && AnyTrue[Differences[dkWtab[[All, 2]]], # < 0 &] && -0.06 < dkWmin < -0.04]
dkZ48 = 8 (1 - 4^-3) Zeta[4] Zeta[3];
check["Epstein sum over Z^4 without 0: Sum |n|^-8 == 8 (1 - 4^-3) zeta(4) zeta(3), compared with the shell sum Sum r_4(k)/k^4 to k = 20000 (tail < 3 10^-8)",
  Abs[N[dkZ48, 20] - Total[rIIr4[[1 ;; 20000]]/N[Range[20000]]^4]] < 3*^-8]
dkLp = rIIVp^(1/4); dkRhoL = 3 rIIH0g^2 rIIMPl^2 rIIOL;
dkC0 = Gamma[4]/Pi^4 dkZ48;                                   (* massless, per real degree of freedom, times 1/l^4 *)
dkCas[xx_?NumericQ] := 16 Total[rIIr4[[1 ;; 400]] 2 (xx/(2 Pi Sqrt[N[Range[400]]]))^4 BesselK[4, xx Sqrt[N[Range[400]]]]]/dkLp^4;
dkXL = xxV /. FindRoot[Log[dkCas[xxV]/dkRhoL] == 0, {xxV, 90}];
check["Casimir energy of dirac16complex on the auxiliary torus: massless 16 Gamma(4) Z/(Pi^4 l^4) > 10^34 rho_Lambda, the massive sum -> the massless value for m l -> 0 (to 10^-4), and equals rho_Lambda at m l = 92.2",
  16 dkC0/dkLp^4/dkRhoL > 10^34 && Abs[dkCas[10^-3]/(16 dkC0/dkLp^4) - 1] < 10^-4 && Abs[dkXL - 92.2] < 0.05]
{"w_min of Kohn-Sham matter and its kF (units m = 1); saturation kF0" -> {SetPrecision[dkWmin, MachinePrecision], N[dkKmin], N[kF0]}, "rho_Lambda (GeV^4)" -> dkRhoL,
 "Casimir energy / rho_Lambda: dirac16complex (m = 0), graviton (20 dof)" -> {16 dkC0/dkLp^4/dkRhoL, -20 dkC0/dkLp^4/dkRhoL},
 "m l_perp and m (GeV) at which the dirac16complex Casimir energy equals rho_Lambda" -> {dkXL, dkXL/dkLp},
 "kF, p, w at kF = 0.1, 0.3, 0.5, 0.6, 0.67, 0.7, 0.75 (units m = 1)" -> N[Select[dkWtab, MemberQ[{1/10, 3/10, 1/2, 3/5, 67/100, 7/10, 3/4}, #[[1]]] &]],
 "dirac16complex Casimir energy / rho_Lambda at m l_perp = 1, 10, 30, 60" -> (dkCas[#]/dkRhoL & /@ {1., 10., 30., 60.})}
''')
inp(r'''
ListPlot[{dkWtab[[All, {1, 3}]]}, Joined -> True, Frame -> True, FrameLabel -> {"kF (units m = 1)", "w = p/eps"}, GridLines -> {{N[kF0]}, {0}},
 PlotLabel -> "Kohn-Sham matter of section 8: the equation of state", ImageSize -> 600]
''')

subsection('24.6 Verdict')
text('The table collects the findings. dirac16complex has every property that the solutions of sections 21–23 and the observations '
     'demand of dark matter — it is cold, dark, collisionless, and enters the late-universe solution exactly as Ω_c — provided the '
     'auxiliary space is stabilized and the Yukawa portal is very weak or absent. The minimal model does not, however, produce the '
     'observed abundance: its primordial gas and its radion make far too much (which the dilution required by the moduli problem must '
     'remove), and after that dilution its portals with couplings of natural size make far too little — the Higgs portal would need '
     'c_P ≈ 150. The origin of the abundance is therefore OPEN. As dark energy dirac16complex fails: as matter it never has w < −0.06 in '
     'a stable state, and as vacuum energy or condensate it is only one more contribution to Λ4, with w = −1 exactly.')
inp(r'''
Grid[{{"claim", "status", "basis"},
  {"dirac16complex behaves as cold matter (w = 0) and can be Omega_c of the late-universe solution", "DERIVED", "sections 20.2, 22.1, 22.4, 24.2"},
  {Row[{"its effect on the radion is negligible (delta beta = ", ScientificForm[dkDbeta, 2], ")"}], "DERIVED", "sections 20.3, 24.2"},
  {"as dark matter it needs the stabilizing potential (otherwise decompactification)", "DERIVED", "section 21.3"},
  {Row[{"dark: sigma_SI = ", ScientificForm[dkSigSI[1, 40.], 2], " cm^2 << LZ; collisionless: sigma/m = ", ScientificForm[dkSigM[1, 1.], 2], " cm^2/g"}], "ESTIMATE", "volume dilution, section 24.3"},
  {Row[{"stable over the age of the universe if m < ", NumberForm[1000 dkMmax, 2], " MeV c_y^(-2/3) (absolutely stable without the Yukawa portal)"}], "ESTIMATE", "section 24.3"},
  {Row[{"mass window: m > 5.7 keV (Lyman-alpha), > ", NumberForm[dkPauli, 2], " keV (Pauli, g = 16)"}], "data", "section 24.1"},
  {Row[{"primordial gas and radion decay overproduce dirac16complex (O = ", ScientificForm[dkOver[1.], 2], " T_RH/GeV; brane share ", ScientificForm[rIIgam/rIIgBulk, 2], ")"}], "DERIVED from the solutions; ESTIMATE for O", "sections 21.3, 22.1, 22.2, 24.4"},
  {Row[{"Yukawa portal: the Kaluza-Klein tower is overproduced unless c_y < ", ScientificForm[dkCyTow, 2]}], "ESTIMATE", "section 24.4"},
  {Row[{"Higgs portal: the observed abundance needs c_P = ", NumberForm[dkCpNeed[mHb/dkTstar], 3], " (T_RH = T*)"}], "ESTIMATE", "section 24.4"},
  {"observed abundance m Y = 4.37 10^-10 GeV produced by the minimal model with natural couplings", "OPEN (not produced)", "section 24.4"},
  {"dirac16complex matter as dark energy", "excluded (DERIVED: w >= -0.06; p < 0 only where uniform matter is unstable or metastable)", "sections 8.7, 24.5"},
  {Row[{"dirac16complex vacuum energy or condensate as dark energy", " (Casimir: ", ScientificForm[16 dkC0/dkLp^4/dkRhoL, 2], " rho_Lambda)"}], "not an explanation (DERIVED: a contribution to Lambda4 with w = -1)", "section 24.5"},
  {"evolving dark energy (DESI DR2: 2.8 to 4.2 sigma)", "in tension with the model (w = -1 exactly)", "section 24.5"}}, Frame -> All, Alignment -> Left]
''')
