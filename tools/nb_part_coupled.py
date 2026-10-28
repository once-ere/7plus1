# Sections 20-23 of 7+1.nb (the coupled field equations with every interaction, and their accurate approximate solutions in
# the three gravity regimes: primordial, modern, far future); appended on 27 October 2026.
# Executed by tools/make_7p1_nb.py, which provides section(), subsection(), text(), item(), inp() and B.

section('20. The coupled field equations: every interaction')
text([('About sections 20–23. ', B),
      'Sections 20–23 were appended on 27 October 2026. The notebook as it was before (sections 0–20, evaluated; section 20 was '
      'then the summary of all checks) is kept in the workspace as backups/7+1_2026-10-27_before_coupled.nb (excluded from the '
      'repository by .gitignore). Section 20 writes the complete action of the model, interaction by interaction, derives the '
      'coupled field equations and the energy–momentum tensor of every sector, and reduces them to the homogeneous cosmological '
      'equations. Sections 21–23 solve the coupled equations in the three gravity regimes of the model: the PRIMORDIAL regime '
      '(section 21: near x0 = 0 the higher Lovelock terms dominate; inside it the cubic, Gauss–Bonnet and Einstein terms take over '
      'one after another), the MODERN regime (section 22: the auxiliary space x4 … x7 is held fixed by its stabilizing potential and '
      'gravity on the brane is four-dimensional Einstein gravity), and the FAR FUTURE (section 23: de Sitter expansion of ordinary '
      'space). In each regime the solution is given in closed form (exact where possible, asymptotic otherwise), and every closed '
      'form is compared with a high-accuracy numerical solution of the full equations; the size of every interaction in every regime '
      'is listed in section 23.3.'])
text([('Conventions (as in sections 11–19). ', B),
      'ħ = c = 1; coordinates X = {x0, …, x7}: x0 is time (the unit e0 = 1 of the octonions), x1, x2, x3 ordinary space, x4 … x7 '
      'the four auxiliary space-like coordinates. η17 = diag(+1, −1, …, −1) is simply the squares of the octonion basis (e0² = +1, '
      'e_i² = −1); in xAct the spacetime is declared with the signature list {7, 1, 0}: seven space-like and one time-like '
      'dimension. Curvature conventions, the Lovelock tensors E^(k) and densities L_k are those of sections 11–12; κ = 8πG8.'])

subsection('20.1 The action, term by term')
text('The total action is the sum of the following terms (√|g| d⁸X is the volume element of the bulk, √|h| d⁴x that of the '
     'world brane W = {x4 = … = x7 = 0}, h_μν = g_μν(x, y = 0) its induced metric; ψ is dirac16complex, ψbar = ψ† Γ17[0]):')
item('Gravity, k = 0 (cosmological term):  S_0 = −(1/(2κ)) ∫ d⁸X √|g| α_0,  α_0 = 2Λ (Λ the bulk cosmological constant).')
item('Gravity, k = 1 (Einstein–Hilbert):  S_1 = −(1/(2κ)) ∫ d⁸X √|g| R.')
item('Gravity, k = 2 (Gauss–Bonnet):  S_2 = −(α_2/(2κ)) ∫ d⁸X √|g| (R² − 4 R_μν R^μν + R_μνρσ R^μνρσ).')
item('Gravity, k = 3 (cubic Lovelock):  S_3 = −(α_3/(2κ)) ∫ d⁸X √|g| L_3,  L_3 = R³ − 12 R R_μν R^μν + 16 R^μ_ν R^ν_ρ R^ρ_μ '
     '+ 24 R^μρ R^νσ R_μνρσ + 3 R R_μνρσ R^μνρσ − 24 R^μν R_μ^ρστ R_νρστ + 4 R^μν_ρσ R^ρσ_λτ R^λτ_μν − 8 R_μ^ρ_ν^σ R_ρ^λ_σ^τ R_λ^μ_τ^ν '
     '(= (1/8) δ^{α1β1α2β2α3β3}_{γ1δ1γ2δ2γ3δ3} R^{γ1δ1}_{α1β1} R^{γ2δ2}_{α2β2} R^{γ3δ3}_{α3β3}, section 12; the invariant form is checked '
     'below). The quartic term would be the Euler density, a total derivative in eight dimensions.')
item('dirac16complex and its coupling to gravity:  S_ψ = ∫ d⁸X √|g| ψbar i Γ^μ(X) (∂_μ + Ω_μ) ψ, with the frame matrices '
     'Γ^μ = e_A^μ Γ17up[A] and the spin connection Ω_μ = (1/4) ω_μAB Γ17up[A].Γ17up[B] of section 13. Expanded about flat space '
     '(g_μν = η17_μν + h_μν) this is the graviton vertex −(1/2) h_μν T_ψ^μν + O(h²).')
item('The mass of dirac16complex:  S_m = −∫ d⁸X √|g| m ψbar ψ.')
item('The nine self-couplings of dirac16complex (section 15):  S_self = −(1/2) Σ_{r=0}^{8} g_r ∫ d⁸X √|g| Σ_{|K|=r} η_K (ψbar Θ^K ψ)².')
item('The stabilizing potential of the auxiliary space:  S_U = −∫ d⁸X √|g| U(b), a function of the proper size b of the auxiliary '
     'torus (b⁸ = det(−g_kl), k, l = 4 … 7). Sections 21–23 use U(b) = U_0 + (K_U/2) (ln(b/b0))², whose minimum b0 fixes the size of '
     'the auxiliary space today.')
item('The brane tension:  S_f = −∫_W d⁴x √|h| f_W.')
item('The Standard Model on the brane (section 17.2):  S_SM = ∫_W d⁴x √|h| L_SM, with the gauge interactions (vertices '
     'gs q̄ γ^μ T^a q G^a_μ, g f̄_L γ^μ (τ^I/2) f_L W^I_μ, g′ Y f̄ γ^μ f B_μ, the cubic and quartic self-couplings of G and W in '
     '−(1/4) F², and the Higgs–gauge couplings in (D_μΦ)†(D^μΦ)), the Higgs potential −μ² Φ†Φ − λ (Φ†Φ)², and the Yukawa couplings '
     '−[Q̄_L Y_u Φ̃ u_R + Q̄_L Y_d Φ d_R + L̄_L Y_e Φ e_R + L̄_L Y_ν Φ̃ ν_R + h.c.]. Its fields feel gravity through h_μν, the brane '
     'tetrad and the brane spin connection.')
item('The three portals between the Standard Model and dirac16complex (section 18.2):  S_P = ∫_W d⁴x √|h| [−Σ y_is L̄_i Φ̃ P_R ψ_s + h.c. '
     '− G_P^f (ψbar Γ^μ ψ)|_W (f̄ γ_μ f) − λ_P Φ†Φ (ψbar ψ)|_W].')
text('In eight-dimensional form a brane integral is ∫_W d⁴x √|h| (…) = ∫ d⁸X √|g| (δ⁴(y)/b⁴) (…). Two checks below make the '
     'gravitational part explicit: the invariant form of L_3 and of the Gauss–Bonnet density L_2 is compared EXACTLY with the '
     'generalized-delta (pair-set) densities lovelockL of section 12.5, on random integer tensors with all the algebraic symmetries '
     'of the Riemann tensor (an algebraic identity needs no metric). Such tensors are sums of products '
     'R_abcd = S_ac S_bd − S_ad S_bc of random symmetric integer matrices S.')
inp(r'''
feRandCurv[seed_] := Module[{S, R4 = ConstantArray[0, {8, 8, 8, 8}]}, SeedRandom[seed];
   Do[S = RandomInteger[{-3, 3}, {8, 8}]; S = S + Transpose[S];
    R4 += Table[S[[a, c]] S[[b, d]] - S[[a, d]] S[[b, c]], {a, 8}, {b, 8}, {c, 8}, {d, 8}], {3}]; R4];
feEtaD = Diagonal[\[Eta]17];
feInvariants[RL_] := Module[{rud4, ric, rs, ricU, ricM, rU, tC, raU, m64},
   rud4 = Table[feEtaD[[a]] feEtaD[[b]] RL[[a, b, c, d]], {a, 8}, {b, 8}, {c, 8}, {d, 8}];       (* R^{ab}_{cd} *)
   ric = Table[Sum[feEtaD[[a]] RL[[a, c, a, d]], {a, 8}], {c, 8}, {d, 8}];                       (* R_cd = g^ab R_acbd *)
   rs = Sum[feEtaD[[c]] ric[[c, c]], {c, 8}];
   ricU = \[Eta]17.ric.\[Eta]17; ricM = \[Eta]17.ric;
   rU = Table[feEtaD[[a]] feEtaD[[b]] feEtaD[[c]] feEtaD[[d]] RL[[a, b, c, d]], {a, 8}, {b, 8}, {c, 8}, {d, 8}];
   tC = Table[feEtaD[[c]] feEtaD[[d]] RL[[a, c, b, d]], {a, 8}, {c, 8}, {b, 8}, {d, 8}];       (* R_a^c_b^d *)
   raU = Table[feEtaD[[c]] feEtaD[[d]] feEtaD[[e]] RL[[b, c, d, e]], {b, 8}, {c, 8}, {d, 8}, {e, 8}];
   m64[t_] := ArrayReshape[t, {64, 64}];
   <|"R" -> rs, "Ric.Ric" -> Total[ric ricU, 2], "Riem.Riem" -> Total[rU RL, 4], "R^3" -> rs^3, "R Ric.Ric" -> rs Total[ric ricU, 2],
     "Ric^3" -> Tr[ricM.ricM.ricM], "Ric Ric Riem" -> Sum[ricU[[a, c]] ricU[[b, d]] RL[[a, b, c, d]], {a, 8}, {b, 8}, {c, 8}, {d, 8}],
     "R Riem.Riem" -> rs Total[rU RL, 4], "Ric Riem Riem" -> Total[ricU (ArrayReshape[RL, {8, 512}].Transpose[ArrayReshape[raU, {8, 512}]]), 2],
     "Riem^3 pairs" -> Tr[m64[rud4].m64[rud4].m64[rud4]],
     "Riem^3 crossed" -> Tr[#.#.#] &[ArrayReshape[Transpose[tC, {1, 3, 2, 4}], {64, 64}]]|>];
feL3coef = <|"R^3" -> 1, "R Ric.Ric" -> -12, "Ric^3" -> 16, "Ric Ric Riem" -> 24, "R Riem.Riem" -> 3, "Ric Riem Riem" -> -24, "Riem^3 pairs" -> 4, "Riem^3 crossed" -> -8|>;
feL3test = Table[With[{RL = feRandCurv[seed]}, Module[{inv = feInvariants[RL], ra},
     ra = rudAssoc[Table[feEtaD[[a]] feEtaD[[b]] RL[[a, b, c, d]], {a, 8}, {b, 8}, {c, 8}, {d, 8}]];
     {lovelockL[ra, 3], Total[KeyValueMap[#2 inv[#1] &, feL3coef]], lovelockL[ra, 2], inv["R"]^2 - 4 inv["Ric.Ric"] + inv["Riem.Riem"]}]], {seed, {3108, 2710, 1977, 88}}];
check["cubic Lovelock density in invariant form: L_3 == R^3 - 12 R Ric.Ric + 16 Ric^3 + 24 R^ac R^bd R_abcd + 3 R Riem.Riem - 24 R^ab R_a^cde R_bcde + 4 R^ab_cd R^cd_ef R^ef_ab - 8 R_a^c_b^d R_c^e_d^f R_e^a_f^b (exact, four random integer algebraic curvature tensors)",
  AllTrue[feL3test, #[[1]] == #[[2]] &]]
check["Gauss-Bonnet density in invariant form: L_2 == R^2 - 4 Ric.Ric + Riem.Riem (the same tensors, exact)", AllTrue[feL3test, #[[3]] == #[[4]] &]]
{"L_3 (pair-set) and the invariant form for the four tensors" -> feL3test[[All, 1 ;; 2]]}
''')

subsection('20.2 The coupled field equations and the energy–momentum tensors')
text('Varying the action gives the coupled field equations of all fields:')
text('      Σ_{k=0}^{3} α_k E^(k)μ_ν = κ [T_ψ + T_U + (δ⁴(y)/b⁴) (T_f + T_SM + T_P)]^μ_ν                (gravity, E^(0) = −δ/2, α_1 = 1),')
text('      (i Γ^μ D_μ − m) ψ − Σ_r g_r Σ_{|K|=r} (ψbar Θ_K ψ) Θ^K ψ = (δ⁴(y)/b⁴) [Σ y*_is (Φ̃† L_i) ⊗ (portal spinor s) + G_P^f (f̄γ_μ f) Γ^μ ψ + λ_P Φ†Φ ψ]   (dirac16complex),')
text('      ∇_μ F^aμν = (gauge currents) + (portal currents), D_μD^μΦ + μ²Φ + 2λ(Φ†Φ)Φ + (Yukawa terms) + λ_P (ψbar ψ)|_W Φ = 0, '
     'i γ^μ D_μ f = (Yukawa masses) + G_P^f (ψbar Γ_μ ψ)|_W γ^μ f   (Standard Model on W).')
text('Every energy–momentum tensor is T_μν = (2/√|g|) δS/δg^μν (sections 12 and 16: δS_matter = (1/2) ∫ √|g| T_μν δg^μν, T^0_0 = ε). '
     'For the homogeneous cosmological metric ds² = n² dx0² − a² dx⃗² − b² dy⃗² every source is diag(ε, −pa, −pa, −pa, −pb, …, −pb). '
     'For matter whose comoving energy depends only on the scale factors, W(a, b) = a³ b⁴ ε (at fixed occupation numbers), the '
     'action is −∫ n W dx0 and')
text('      ε = W/(a³ b⁴),    pa = −(∂W/∂a)/(3 a² b⁴),    pb = −(∂W/∂b)/(4 a³ b³),')
text('which gives: the stabilizing potential (W = a³b⁴ U): ε = U, pa = −U, pb = −U − (b/4) U′(b); the brane tension smeared over the '
     'torus (δ⁴(y) → 1/Lc⁴, V⊥ = b⁴ Lc⁴): ε = f_W/V⊥, pa = −f_W/V⊥, pb = 0; a brane fluid with ρ_4 ∝ a^(−3(1+w)): ε = ρ_4/V⊥, pa = w ε, '
     'pb = 0; a massless dirac16complex particle of comoving momentum (k_a, k_b): εV = ω = √(k_a²/a² + k_b²/b²), paV = k_a²/(3a²ω), '
     'pbV = k_b²/(4b²ω). Each such source is conserved, ε′ + 3Ha(ε + pa) + 4Hb(ε + pb) = 0, and the minisuperspace variation of '
     '−(1/(2κ)) ∫ n a³b⁴ Σ α_k L_k − ∫ n W (with the identities of section 12.6) gives exactly the three cosmological field equations')
text('      Σ_k α_k E^(k)0_0 = κ ε,    Σ_k α_k E^(k)1_1 = −κ pa,    Σ_k α_k E^(k)4_4 = −κ pb.')
inp(r'''
feW = Function[{aa, rIIIbb}, wW[aa, rIIIbb]];
feStress[W_] := Simplify[{W[aS, bS]/(aS^3 bS^4), -D[W[aS, bS], aS]/(3 aS^2 bS^4), -D[W[aS, bS], bS]/(4 aS^3 bS^3)}];
check["stabilizing potential: W = a^3 b^4 U(b) gives eps = U, pa = -U, pb = -U - (b/4) U'(b)",
  Simplify[feStress[Function[{aa, rIIIbb}, aa^3 rIIIbb^4 uPot[rIIIbb]]] - {uPot[bS], -uPot[bS], -uPot[bS] - bS uPot'[bS]/4}] === {0, 0, 0}]
check["brane tension smeared over the torus (V_perp = b^4 Lc^4): eps = f_W/V_perp, pa = -f_W/V_perp, pb = 0",
  Simplify[feStress[Function[{aa, rIIIbb}, aa^3 fWb/lcT^4]] - {fWb/(bS^4 lcT^4), -fWb/(bS^4 lcT^4), 0}] === {0, 0, 0}]
check["brane fluid rho_4 = N a^(-3(1+w)), smeared: eps = rho_4/V_perp, pa = w eps, pb = 0",
  Simplify[feStress[Function[{aa, rIIIbb}, aa^3 nW aa^(-3 (1 + wW0))/lcT^4]] - {nW aS^(-3 (1 + wW0))/(bS^4 lcT^4), wW0 nW aS^(-3 (1 + wW0))/(bS^4 lcT^4), 0}] === {0, 0, 0}]
check["massless dirac16complex particle (comoving momentum (k_a, k_b)): eps V = omega, pa V = k_a^2/(3 a^2 omega), pb V = k_b^2/(4 b^2 omega)",
  With[{st = feStress[Function[{aa, rIIIbb}, Sqrt[kaP^2/aa^2 + kbP^2/rIIIbb^2]]], om = Sqrt[kaP^2/aS^2 + kbP^2/bS^2]},
   Simplify[st aS^3 bS^4 - {om, kaP^2/(3 aS^2 om), kbP^2/(4 bS^2 om)}, aS > 0 && bS > 0 && kaP > 0 && kbP > 0] === {0, 0, 0}]]
consW = With[{e = wW[a3[x0], b4[x0]]/(a3[x0]^3 b4[x0]^4), pa = -(Derivative[1, 0][wW][a3[x0], b4[x0]])/(3 a3[x0]^2 b4[x0]^4),
     pb = -(Derivative[0, 1][wW][a3[x0], b4[x0]])/(4 a3[x0]^3 b4[x0]^3)},
   D[e, x0] + 3 a3'[x0]/a3[x0] (e + pa) + 4 b4'[x0]/b4[x0] (e + pb)];
check["every such source is conserved: eps' + 3Ha(eps + pa) + 4Hb(eps + pb) == 0 for arbitrary W(a, b), a(x0), b(x0)", Simplify[consW] === 0]
feELn = -(1/(2 kapS)) (-2 aS^3 bS^4 s00) - wW[aS, bS]; feELa = -(1/(2 kapS)) (-6 aS^2 bS^4 s11) - D[wW[aS, bS], aS];
feELb = -(1/(2 kapS)) (-8 aS^3 bS^3 s44) - D[wW[aS, bS], bS];
check["minisuperspace: EL_n == (a^3 b^4/kappa)(S00 - kappa eps), EL_a == (3 a^2 b^4/kappa)(S11 + kappa pa), EL_b == (4 a^3 b^3/kappa)(S44 + kappa pb) with S = Sum alpha_k E^(k) (identities of section 12.6)",
  With[{st = feStress[Function[{aa, rIIIbb}, wW[aa, rIIIbb]]]},
   Simplify[{feELn - aS^3 bS^4/kapS (s00 - kapS st[[1]]), feELa - 3 aS^2 bS^4/kapS (s11 + kapS st[[2]]), feELb - 4 aS^3 bS^3/kapS (s44 + kapS st[[3]])}] === {0, 0, 0}]]
''')

subsection('20.3 The homogeneous equations, the radion and the four-dimensional constants')
text('With the closed forms of section 12.6 the cosmological field equations are (Hubble rates Ha = a′/a, Hb = b′/b, a prime is '
     'd/dx0; the code writes eTot[c] = Σ_k α_k E^(k)c_c − Λ with alf[k] = α_k, lamb = Λ):')
text('      00:  3(Ha² + 4HaHb + 2Hb²) − 12α_2 Hb(4Ha³ + 18Ha²Hb + 12HaHb² + Hb³) + 360α_3 Ha²Hb³(4Ha + 3Hb) − Λ = κ ε,')
text('      11:  Σ_k α_k E^(k)1_1 − Λ = −κ pa,     44:  Σ_k α_k E^(k)4_4 − Λ = −κ pb,')
text('the 11 and 44 equations being linear in Ha′ and Hb′ and the 00 equation a constraint. The RADION is the logarithm of the size '
     'of the auxiliary space, β = ln(b/b0), Hb = β′. In Einstein gravity the 11 and 44 equations combined with the constraint give the '
     'exact radion equation')
text('      β″ + 3 Ha β′ + 4 β′² = κ (ε/6 − pa/2 + pb/3) + Λ/3.')
text('With U = U_0 + (K_U/2) β² the right-hand side contains −(κK_U/12) β: the radion mass is m_r² = κ K_U/12. With the auxiliary space '
     'static (Hb = Hb′ = 0) the 00 and 11 equations are the Friedmann equations of the brane,')
text('      3Ha² = 8πG4 ρ_4 + Λ4,   Ha′ = −(κ/2)(ε + pa),   8πG4 = κ/V⊥,   Λ4 = Λ + κ U_0 + κ f_W/V⊥,')
text('and the 44 equation fixes where the radion sits (the Gauss–Bonnet term enters; the cubic term vanishes for Hb = 0):')
text('      (κ/4) b U′(b) = Λ + κ U + 2κ f_W/V⊥ + (κ/2)(ε − 3pa) + κ pb − 12 α_2 Ha² (Ha′ + Ha²).')
text('Finally the 00 equation in four-dimensional (Einstein-frame) form: with H_E = Ha + 2Hb (the expansion rate of a_E = a b²) it '
     'reads 3H_E² = 6β′² + κε exactly, i.e. 3 M_Pl² H_E² = (1/2) φ′² + ρ with the canonical radion φ = √12 M_Pl β. Brane fields see '
     'the metric h = (b0/b)⁴ g_E, so the radion couples to them through L = (φ/(√3 M_Pl)) T^μ_μ (the trace of their '
     'energy–momentum tensor); this is the coupling through which the radion decays (section 22.2).')
inp(r'''
feTot[c_] := Expand[eTot[c] /. {alf[1] -> 1}];
feSrc = {epsM + uB, paM - uB, pbM - uB - bUp/4};               (* matter + stabilizing potential: eps, pa, pb; uB = U, bUp = b U'(b) *)
feSolE = First[Solve[{(feTot[2] /. {alf[2] -> 0, alf[3] -> 0}) == -kapS feSrc[[2]], (feTot[5] /. {alf[2] -> 0, alf[3] -> 0}) == -kapS feSrc[[3]]}, {dHa, dHb}]];
feRad = Expand[(dHb /. feSolE) /. Ha^2 -> (kapS feSrc[[1]] + lamb - 12 Ha Hb - 6 Hb^2)/3];
check["radion equation of Einstein gravity: Hb' + 3 Ha Hb + 4 Hb^2 == kappa (eps/6 - pa/2 + pb/3) + Lambda/3 (constraint used)",
  Simplify[feRad - (-3 Ha Hb - 4 Hb^2 + kapS (feSrc[[1]]/6 - feSrc[[2]]/2 + feSrc[[3]]/3) + lamb/3)] === 0]
check["radion mass: U = U0 + (K_U/2) beta^2, b U' = K_U beta gives the restoring term -(kappa K_U/12) beta, m_r^2 = kappa K_U/12",
  Simplify[D[(-3 Ha Hb - 4 Hb^2 + kapS (feSrc[[1]]/6 - feSrc[[2]]/2 + feSrc[[3]]/3) + lamb/3) /. {uB -> u0B + kU/2 betR^2, bUp -> kU betR, Hb -> 0}, betR] /. betR -> 0] === -kapS kU/12]
feSt = {epsM + uB + fV, paM - uB - fV, pbM - uB - bUp/4};      (* with the smeared brane tension fV = f_W/V_perp *)
feStatic = Solve[{(feTot[1] /. {Hb -> 0}) == kapS feSt[[1]], (feTot[2] /. {Hb -> 0, dHb -> 0}) == -kapS feSt[[2]]}, {lamb, dHa}][[1]];
check["static auxiliary space: (kappa/4) b U' == Lambda + kappa U + 2 kappa f_W/V_perp + (kappa/2)(eps - 3pa) + kappa pb - 12 alpha2 Ha^2 (Ha' + Ha^2)",
  Simplify[(bUp /. First[Solve[((feTot[5] /. {Hb -> 0, dHb -> 0}) == -kapS feSt[[3]]) /. feStatic, bUp]]) -
    (4/kapS (lamb + kapS uB + 2 kapS fV + kapS (epsM - 3 paM)/2 + kapS pbM - 12 alf[2] Ha^2 (dHa + Ha^2)) /. feStatic)] === 0]
check["static auxiliary space: Ha' == -(kappa/2)(eps + pa) and E^(3) does not contribute (Hb = 0)",
  Simplify[(dHa /. feStatic) + kapS (epsM + paM)/2] === 0 && FreeQ[Simplify[{feTot[1], feTot[2], feTot[5]} /. {Hb -> 0, dHb -> 0}], alf[3]]]
check["brane Friedmann equation: 3 Ha^2 == (kappa/V_perp) rho_4 + Lambda4 with Lambda4 = Lambda + kappa U0 + kappa f_W/V_perp, i.e. 8 Pi G4 = kappa/V_perp",
  Simplify[(feTot[1] /. Hb -> 0) - kapS (rho4/vPerpS + u0B + fV) - (3 Ha^2 - (kapS/vPerpS) rho4 - (lamb + kapS u0B + kapS fV))] === 0]
check["Einstein-frame form of the constraint: 3 (Ha + 2 Hb)^2 == 6 Hb^2 + (E^(1)0_0 of section 12.6), i.e. 3 M_Pl^2 H_E^2 == phi'^2/2 + rho with phi = Sqrt[12] M_Pl beta",
  Expand[3 (Ha + 2 Hb)^2 - 6 Hb^2 - toH[ECos[[2, 1, 1]]]] === 0]
''')

# =====================================================================================
section('21. Regime I: the primordial Lovelock regime (near x0 = 0)')
text('Units: the Gauss–Bonnet length ℓ2 = √α_2 = 1 (time and length in units of ℓ2, Hubble rates in units of 1/ℓ2), so that the '
     'only parameter of the vacuum equations is ξ = α_3/α_2². This section uses ξ = 10⁻⁶, which separates the three Lovelock regimes '
     'clearly; the crossover times scale as stated below. τ = x0 − x_s is the time since the singularity x_s. Λ is negligible in '
     'this regime (Λ ℓ2² ≈ (H0/M8)² ≈ 10⁻⁹², section 23.3) and set to 0. If ℓ2 = 1/M8 with M8 = 10 TeV, the unit of time is '
     'ℓ2 = 6.6 × 10⁻²⁹ s.')

subsection('21.1 The exact vacuum solution: an algebraic curve and one quadrature')
text('Write the vacuum constraint (the 00 equation with ε = 0) in the variables u = Hb/Ha and X = Ha². Dividing by Ha² it becomes')
text('      P_E(u) + X P_G(u) + ξ X² P_C(u) = 0,   P_E = 3(2u² + 4u + 1),   P_G = −12u(u³ + 12u² + 18u + 4),   P_C = 360 u³ (3u + 4),')
text('the Einstein, Gauss–Bonnet and cubic parts. For every u in (u1, 0), u1 = −1 + 1/√2, it is a quadratic equation for X with '
     'exactly one positive root (P_C < 0 and P_E > 0 there), so the vacuum solution is the explicit curve')
text('      Ha(u) = √X(u),   Hb(u) = u √X(u),   X(u) = [P_G + √(P_G² − 4ξ P_C P_E)]/(2ξ |P_C|).')
text('The evolution moves along this curve (the constraint is preserved, section 16.3); Ha′ follows from the 11 and 44 equations, '
     'and the time is one quadrature, τ(Ha) = ∫_Ha^∞ dH/F(H) with F = −Ha′ > 0. The effective Kasner exponent p_eff = −Ha²/Ha′ (= p for '
     'a = τ^p) runs from p1 = (1 + 2√2)/7 (u → u1: Einstein regime, Ha → 0) through about p2 (u near u2, the root of '
     'u³ + 12u² + 18u + 4 in (−1/3, 0): Gauss–Bonnet regime) to 1 (u → 0: cubic regime, Ha → ∞, Hb → h0 = −√(α_2/(30α_3))). The two '
     'Kasner solutions are p1 = 1/(3 + 4u1), q1 = u1 p1 (3p + 4q = 1) and p2 = 3/(3 + 4u2), q2 = u2 p2 (3p + 4q = 3), the same '
     'exponents as in section 14.1.')
inp(r'''
rIxi = 10^-6;
rIcon = Expand[Expand[(feTot[1] /. {alf[2] -> 1, alf[3] -> xiS, lamb -> 0, Hb -> uS Ha})/Ha^2] /. Ha -> Sqrt[xS]];
{rIpE, rIpG, rIpC} = {Factor[Coefficient[rIcon, xS, 0]], Factor[Coefficient[rIcon, xS, 1]], Factor[Coefficient[rIcon, xS, 2]/xiS]};
check["vacuum constraint in u = Hb/Ha, X = Ha^2 (alpha2 = 1, alpha3 = xi): P_E + X P_G + xi X^2 P_C == 0 with P_E = 3(2u^2 + 4u + 1), P_G = -12u(u^3 + 12u^2 + 18u + 4), P_C = 360 u^3 (3u + 4)",
  Expand[rIpE - 3 (2 uS^2 + 4 uS + 1)] === 0 && Expand[rIpG + 12 uS (uS^3 + 12 uS^2 + 18 uS + 4)] === 0 && Expand[rIpC - 360 uS^3 (3 uS + 4)] === 0 &&
  Expand[rIcon - (rIpE + xS rIpG + xiS xS^2 rIpC)] === 0]
rIpEf[u_] := Evaluate[rIpE /. uS -> u]; rIpGf[u_] := Evaluate[rIpG /. uS -> u]; rIpCf[u_] := Evaluate[rIpC /. uS -> u];
rIX[u_, xi_] := (rIpGf[u] + Sqrt[rIpGf[u]^2 - 4 xi rIpCf[u] rIpEf[u]])/(2 xi Abs[rIpCf[u]]);
rIu1 = -1 + 1/Sqrt[2]; rIu2 = Root[#^3 + 12 #^2 + 18 # + 4 &, 3];
{rIp1, rIq1, rIp2, rIq2} = {1/(3 + 4 rIu1), rIu1/(3 + 4 rIu1), 3/(3 + 4 rIu2), 3 rIu2/(3 + 4 rIu2)};
check["Kasner constants: P_E(u1) == 0, p1 == (1 + 2 Sqrt[2])/7, q1 == (2 - 3 Sqrt[2])/14; u2 in (-1/3, 0) with P_G(u2) == 0, 3 p2 + 4 q2 == 3",
  Simplify[rIpEf[rIu1]] == 0 && Simplify[rIp1 - (1 + 2 Sqrt[2])/7] == 0 && Simplify[rIq1 - (2 - 3 Sqrt[2])/14] == 0 &&
  -1/3 < N[rIu2] < 0 && RootReduce[rIpGf[rIu2]] == 0 && RootReduce[3 rIp2 + 4 rIq2] == 3]
rIdHsol = First[Solve[{feTot[2] == 0, feTot[5] == 0} /. {alf[2] -> 1, alf[3] -> xiS, lamb -> 0}, {dHa, dHb}]];
rIdHa[ha_, hb_, xi_] := Evaluate[(dHa /. rIdHsol) /. {Ha -> ha, Hb -> hb, xiS -> xi}];
rIdHb[ha_, hb_, xi_] := Evaluate[(dHb /. rIdHsol) /. {Ha -> ha, Hb -> hb, xiS -> xi}];
rIpeff[u_?NumericQ, xi_] := With[{h = Sqrt[rIX[u, xi]]}, -h^2/rIdHa[h, u h, xi]];
check["the explicit curve solves the vacuum constraint: residual < 10^-25 (relative) at u = u1 + 10^-6, -0.28, u2, -0.1, -10^-4 (xi = 10^-6, 40 digits)",
  AllTrue[{N[rIu1, 40] + 10^-6, -28/100, N[rIu2, 40], -1/10, -1/10^4}, With[{x = rIX[SetPrecision[#, 40], rIxi]},
     Abs[rIpEf[SetPrecision[#, 40]] + x rIpGf[SetPrecision[#, 40]] + rIxi x^2 rIpCf[SetPrecision[#, 40]]] < 10^-25 Max[1, x^2]] &]]
rIpeffVals = {rIpeff[N[rIu1, 40] + 10^-10, rIxi], rIpeff[N[rIu2, 40], rIxi], rIpeff[-10^-6, rIxi]};
check["the effective exponent along the curve: p_eff -> p1 (u -> u1), about p2 at u = u2 (within 10^-2), -> 1 (u -> 0)",
  Abs[rIpeffVals[[1]] - rIp1] < 10^-8 && Abs[rIpeffVals[[2]] - rIp2] < 10^-2 && Abs[rIpeffVals[[3]] - 1] < 10^-5]
{"u1, u2" -> N[{rIu1, rIu2}, 16], "p1, q1, p2, q2" -> N[{rIp1, rIq1, rIp2, rIq2}, 16], "p_eff at u -> u1, u2, u -> 0" -> N[rIpeffVals, 12]}
''')

subsection('21.2 Accurate approximate solutions in the three Lovelock regimes')
text([('Cubic regime (τ → 0). ', B),
      'The series Ha = 1/τ + A0 + A1 τ + …, Hb = h0 + h1 τ + h2 τ² + … solves the three equations order by order with'])
text('      h0 = −√(α_2/(30α_3)),   h1 = α_2/(16α_3) − 1/(32α_2),   A0 = √30 (58α_2² − 45α_3)/(1920 α_2^(3/2) α_3^(1/2)),')
text('      A1 = −3647α_2/(92160α_3) − 223/(6144α_2) + 645α_3/(8192α_2³),   h2 = √30 (−1212α_2⁴ + 1180α_2²α_3 + 225α_3²)/(307200 α_2^(5/2) α_3^(3/2)).')
text('So a ∝ τ (ordinary space expands like a Milne universe) while b stays finite: the auxiliary space contracts at the finite rate '
     '|h0| — here the cubic and Gauss–Bonnet terms balance (α_3 E^(3) ≈ −α_2 E^(2)) and the Einstein term is negligible.')
text([('Gauss–Bonnet regime. ', B),
      'Expanding the curve about u2 gives p_eff = p2 + f1/Ha² + f2 ξ Ha² + … (f1 and f2 are computed below), so that'])
text('      τ = p2/Ha + f1/(3Ha³) − f2 ξ Ha + s2,   a ∝ (τ − s2)^p2,   b ∝ (τ − s2)^q2   (p2 = 1.5611…, q2 = −0.42086…),')
text('with the time shift s2 fixed at the middle of the regime, Ha_mid = (f1/(f2 ξ))^(1/4).')
text([('Einstein regime (Ha → 0). ', B),
      'Expanding about u1: u = u1 + d1 Ha² + …, d1 = α_2 (16 − 11√2)/4, and p_eff = p1 + e1 Ha² + …, e1 = α_2 (1 + 9√2)/7, so'])
text('      τ = p1/Ha − e1 Ha + s1,   Ha = p1/(τ − s1) − e1 p1²/(τ − s1)³ + …,   a ∝ (τ − s1)^p1,   b ∝ (τ − s1)^q1,')
text('with s1 = τ(u2) − p1/Ha(u2) + ∫_0^Ha(u2) (1/F − p1/H²) dH from the exact quadrature along the curve. The crossovers are at '
     'τ_32 ≈ |q2|/|h0| = |q2| √(30ξ) ℓ2 (cubic → Gauss–Bonnet) and τ_21 ≈ ℓ2 (Gauss–Bonnet → Einstein).')
text('The checks below derive the coefficients, integrate the vacuum equations from cubic-series data at τ = 10⁻⁸ (in the variable '
     'σ = −ln Ha, with x0 as a dependent variable), and compare: each regime formula is accurate to better than 10⁻³ inside its '
     'regime; across the two crossover windows the exact curve-plus-quadrature solution takes over.')
inp(r'''
rIc = Expand[{feTot[1], feTot[2], feTot[5]} /. lamb -> 0];
rIser = Expand[Expand[rIc /. {Ha -> 1/tS + cA0 + cA1 tS + cA2 tS^2, Hb -> -Sqrt[alf[2]/(30 alf[3])] + ch1 tS + ch2 tS^2 + ch3 tS^3,
        dHa -> -1/tS^2 + cA1 + 2 cA2 tS, dHb -> ch1 + 2 ch2 tS + 3 ch3 tS^2}] tS^4];
rIco = Table[Simplify[Coefficient[rIser[[c]], tS, n], alf[2] > 0 && alf[3] > 0], {c, 3}, {n, 0, 3}];
rIs1 = First[Solve[{rIco[[1, 3]] == 0, rIco[[3, 2]] == 0}, {ch1, cA0}]];
rIs2 = First[Solve[{(rIco[[1, 4]] /. rIs1) == 0, (rIco[[3, 3]] /. rIs1) == 0}, {ch2, cA1}]];
rIsame[e1_, e2_] := AllTrue[{{1, 10^-6}, {7/3, 2/5}, {1/2, 3}}, RootReduce[Together[(e1 - e2) /. {alf[2] -> #[[1]], alf[3] -> #[[2]]}]] === 0 &];      (* exact *)
check["cubic regime: Ha = 1/tau + A0 + A1 tau, Hb = h0 + h1 tau + h2 tau^2 solves the 00, 11, 44 equations through order tau^-1 (11: checked, not used)",
  AllTrue[rIco[[All, 1]], # === 0 &] && Simplify[rIco[[1, 2]]] === 0 && Simplify[rIco[[2, 2]]] === 0 &&
  rIsame[rIco[[2, 3]] /. rIs1, 0] && rIsame[(rIco[[2, 4]] /. rIs1) /. rIs2, 0]]
check["cubic regime coefficients: h1 == alpha2/(16 alpha3) - 1/(32 alpha2), A0 == Sqrt[30](58 alpha2^2 - 45 alpha3)/(1920 alpha2^(3/2) alpha3^(1/2)), A1 and h2 as in the text",
  rIsame[ch1 /. rIs1, alf[2]/(16 alf[3]) - 1/(32 alf[2])] && rIsame[cA0 /. rIs1, Sqrt[30] (58 alf[2]^2 - 45 alf[3])/(1920 alf[2]^(3/2) Sqrt[alf[3]])] &&
  rIsame[cA1 /. rIs2, -3647 alf[2]/(92160 alf[3]) - 223/(6144 alf[2]) + 645 alf[3]/(8192 alf[2]^3)] &&
  rIsame[ch2 /. rIs2, Sqrt[30] (-1212 alf[2]^4 + 1180 alf[2]^2 alf[3] + 225 alf[3]^2)/(307200 alf[2]^(5/2) alf[3]^(3/2))]]
rIcub = N[{-Sqrt[alf[2]/(30 alf[3])], ch1 /. rIs1, ch2 /. rIs2, cA0 /. rIs1, cA1 /. rIs2} /. {alf[2] -> 1, alf[3] -> rIxi}, 20];   (* h0, h1, h2, A0, A1 *)
rIeSer = Normal[Series[(rIc[[1]] /. {Hb -> (rIu1 + cd1 eS^2 + cd2 eS^4) Ha} /. Ha -> eS), {eS, 0, 6}]];
rId1 = cd1 /. First[Solve[Coefficient[rIeSer, eS, 4] == 0, cd1]];
rId2 = cd2 /. First[Solve[(Coefficient[rIeSer, eS, 6] /. cd1 -> rId1) == 0, cd2]];
rIdHsolG = First[Solve[{rIc[[2]] == 0, rIc[[3]] == 0}, {dHa, dHb}]];
rIpeE = Normal[Series[-eS^2/((dHa /. rIdHsolG) /. {Hb -> (rIu1 + rId1 eS^2 + rId2 eS^4) Ha} /. Ha -> eS), {eS, 0, 2}]];
rIe1 = Simplify[Coefficient[rIpeE, eS, 2]];
check["Einstein regime: u = u1 + d1 Ha^2, d1 == alpha2 (16 - 11 Sqrt[2])/4; p_eff = p1 + e1 Ha^2, e1 == alpha2 (1 + 9 Sqrt[2])/7",
  Simplify[rIeSer /. {eS -> 0}] === 0 && Simplify[Coefficient[rIeSer, eS, 2]] === 0 && Simplify[rId1 - alf[2] (16 - 11 Sqrt[2])/4] === 0 &&
  Simplify[rIe1 - alf[2] (1 + 9 Sqrt[2])/7] === 0 && Simplify[Coefficient[rIpeE, eS, 0] - (1 + 2 Sqrt[2])/7] === 0]
rIe1v = N[rIe1 /. alf[2] -> 1, 20];
rIuOn[h_, xi_] := u /. FindRoot[rIpEf[u] + h^2 rIpGf[u] + xi h^4 rIpCf[u] == 0, {u, N[rIu2, 40]}, WorkingPrecision -> 40];
rIpeH[h_, xi_] := With[{uu = rIuOn[h, xi]}, -h^2/rIdHa[h, uu h, xi]];
rIf1 = (rIpeH[10^5, 0] - N[rIp2, 40]) 10^10; rIf2 = (rIpeH[10^8, 10^-24] - N[rIp2, 40]) 10^8;
check["Gauss-Bonnet regime: p_eff = p2 + f1/Ha^2 + f2 xi Ha^2 with f1 = -0.69887346..., f2 = -0.7853407... (differences give the same values at Ha = 10^4 and 10^-6 xi Ha^2 to 10^-5)",
  Abs[(rIpeH[10^4, 0] - N[rIp2, 40]) 10^8 - rIf1] < 10^-5 && Abs[(rIpeH[10^8, 10^-22] - N[rIp2, 40]) 10^6 - rIf2] < 10^-5]
{"cubic: h0, h1, h2, A0, A1 (xi = 10^-6)" -> rIcub, "Einstein: d1, e1 (alpha2 = 1)" -> N[{rId1, rIe1} /. alf[2] -> 1, 16], "Gauss-Bonnet: f1, f2" -> N[{rIf1, rIf2}, 12]}
''')
inp(r'''
rIrhs[sg_?NumericQ, yy_?VectorQ] := With[{ha = Exp[-sg]}, With[{dha = rIdHa[ha, yy[[1]], rIxi], dhb = rIdHb[ha, yy[[1]], rIxi]}, {dhb (-ha/dha), -ha/dha}]];
rItau0 = 10^-8; rIHa0 = 1/rItau0 + rIcub[[4]] + rIcub[[5]] rItau0; rIHb0 = rIcub[[1]] + rIcub[[2]] rItau0 + rIcub[[3]] rItau0^2;
check["cubic-series data at tau = 10^-8 satisfy the vacuum constraint to < 10^-12 (relative to 48 Ha^3 |Hb|)",
  Abs[(rIc[[1]] /. {alf[2] -> 1, alf[3] -> rIxi, Ha -> rIHa0, Hb -> rIHb0})]/(48 rIHa0^3 Abs[rIHb0]) < 10^-12]
rIvacY = NDSolveValue[{rIyV'[sV] == rIrhs[sV, rIyV[sV]], rIyV[-Log[rIHa0]] == {rIHb0, rItau0}}, rIyV, {sV, -Log[rIHa0], -Log[10^-4]},
    Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 8}, PrecisionGoal -> 13, AccuracyGoal -> 40, MaxSteps -> 10^6, InterpolationOrder -> All];
rIvacAt[h_] := rIvacY[-Log[h]];          (* {Hb, tau} on the vacuum solution at Ha = h *)
rIdHdu[u_, xi_] := With[{x = rIX[u, xi]}, -(rIpEf'[u] + x rIpGf'[u] + xi x^2 rIpCf'[u])/(rIpGf[u] + 2 xi x rIpCf[u])/(2 Sqrt[x])];
rIF[u_, xi_] := With[{h = Sqrt[rIX[u, xi]]}, -rIdHa[h, u h, xi]];
rItauU2 = NIntegrate[rIdHdu[-v^2, rIxi]/rIF[-v^2, rIxi] 2 v, {v, 0, 10^-3, 3/100, 3/10, Sqrt[-rIu2]}, WorkingPrecision -> 30, PrecisionGoal -> 16];
rIHu2 = Sqrt[rIX[N[rIu2, 40], rIxi]];
rIJ = NIntegrate[(1/rIF[rIu1 + w^2, rIxi] - rIp1/rIX[rIu1 + w^2, rIxi]) rIdHdu[rIu1 + w^2, rIxi] 2 w,
    {w, 0, Sqrt[N[rIu2 - rIu1, 40]] 10^-3, Sqrt[N[rIu2 - rIu1, 40]]/10, Sqrt[N[rIu2 - rIu1, 40]]}, WorkingPrecision -> 30, PrecisionGoal -> 16];
rIs1v = rItauU2 - rIp1/rIHu2 + rIJ;
rIHmid = (rIf1/(rIf2 rIxi))^(1/4);
rIs2v = rIvacAt[rIHmid][[2]] - rIp2/rIHmid - rIf1/(3 rIHmid^3) + rIf2 rIxi rIHmid;
check["the ODE solution agrees with the exact quadrature: tau at Ha(u2) to 10^-10, and tau - p1/Ha + e1 Ha at Ha = 3 10^-3 equals s1 to 10^-6",
  Abs[rIvacAt[rIHu2][[2]] - rItauU2] < 10^-10 && Abs[rIvacAt[3/1000][[2]] - rIp1/(3/1000) + rIe1v 3/1000 - rIs1v] < 10^-6]
rItauCub[h_] := Nest[1/(h - rIcub[[4]] - rIcub[[5]] #) &, 1/h, 60];
rIacc = Table[With[{h = 10.^lh}, With[{tv = rIvacAt[h][[2]]},
     {lh, tv, Abs[rItauCub[h]/tv - 1], Abs[(rIp2/h + rIf1/(3 h^3) - rIf2 rIxi h + rIs2v)/tv - 1], Abs[(rIp1/h - rIe1v h + rIs1v)/tv - 1]}]], {lh, 8, -3.5, -0.25}];
rIbest[lo_, hi_, col_] := Max[Select[rIacc, lo <= #[[1]] <= hi &][[All, col]]];
check["regime approximations: cubic error < 10^-5 for Ha >= 10^4, Gauss-Bonnet < 10^-3 for 3 <= Ha <= 300, Einstein < 10^-3 for Ha <= 0.1",
  rIbest[4, 8, 3] < 10^-5 && rIbest[0.5, 2.5, 4] < 10^-3 && rIbest[-3.5, -1, 5] < 10^-3]
{"tau(u2), Ha(u2), s1, s2, Ha_mid" -> N[{rItauU2, rIHu2, rIs1v, rIs2v, rIHmid}, 13],
 "max errors: cubic (Ha >= 10^4), GB (3 .. 300), Einstein (<= 0.1); worst of the best approximation in the crossovers" -> {rIbest[4, 8, 3], rIbest[0.5, 2.5, 4], rIbest[-3.5, -1, 5], Max[Min /@ rIacc[[All, 3 ;; 5]]]},
 "tau at Ha = 10^3 (cubic -> GB) and Ha = 1 (GB -> Einstein)" -> {rIvacAt[10^3][[2]], rIvacAt[1][[2]]}}
''')
inp(r'''
ListLogLogPlot[{rIacc[[All, {2, 3}]], rIacc[[All, {2, 4}]], rIacc[[All, {2, 5}]]}, Joined -> True, PlotRange -> {10^-12, 10},
 PlotLegends -> {"cubic series", "Gauss-Bonnet Kasner + corrections", "Einstein Kasner + corrections"}, Frame -> True,
 FrameLabel -> {"tau (units l2)", "relative error of tau(Ha)"}, ImageSize -> 600]
''')

subsection('21.3 Coupled to dirac16complex: the matter-dominated cubic regime and decompactification')
text('Now couple the gravity of section 21.2 to a gas of dirac16complex particles. In the primordial regime the mass and the '
     'self-couplings are negligible (m ≪ Ha and g_r n^(6/7) ≪ 1 with n the number density of the gas, i.e. (k_phys/M8)^6 ≪ 1 for g_r ≈ M8⁻⁶; section 23.3), so the Kohn–Sham matter of section 16 reduces to a '
     'massless collisionless gas with an isotropic comoving momentum distribution (normalized so that a = b = 1 at Ha = 1, where its '
     'energy density is κε = 10⁻⁶ × 3Ha²): with the weight cos²θ sin³θ of the seven-dimensional momentum sphere,')
text('      κε = C (a³b⁴)⁻¹ ⟨√(cos²θ/a² + sin²θ/b²)⟩,   κ pa = C (a³b⁴)⁻¹ ⟨cos²θ/(3a² ω)⟩,   κ pb = C (a³b⁴)⁻¹ ⟨sin²θ/(4b² ω)⟩.')
text([('Matter-dominated cubic regime. ', B),
      'Near τ = 0 any fluid with pa = w_a ε (w_a < 1) wins over the vacuum terms. With b finite and Ha = p/τ the leading 44 '
      'equation is p³(p − 1)(648α_3 h² − 12α_2) = 0, so Hb → h* = −√(α_2/(54α_3)) (not h0), and the 11 and 00 equations give'])
text('      p = 1/(1 + w_a),   κ ε τ³ → (64/3) p³ α_2 |h*|,   α_2E^(2) : α_3E^(3) : κε → 9/4 : −5/4 : 1.')
text('For the gas (and for brane radiation) w_a → 1/3 as a/b → 0: a ∝ τ^(3/4). The GAS FRACTION is κε divided by the sum of the sizes '
     'of all terms of the 00 equation, S_1 + S_2 + S_3 + κε, where the size of a curvature term is the sum of the absolute values of '
     'its monomials in Ha and Hb: S_1 = 3(Ha² + 4Ha|Hb| + 2Hb²), S_2 = 12α_2|Hb|(4Ha³ + 18Ha²|Hb| + 12Ha Hb² + |Hb|³), '
     'S_3 = 360α_3 Ha²|Hb|³(4Ha + 3|Hb|) (in a Kasner regime the monomials of the dominant term cancel each other, so its net value '
     'would hide it). Moving forward in time the gas fraction decreases '
     'in the cubic regime (ε ∝ τ⁻⁴ against τ⁻³) and in the Gauss–Bonnet regime (ε ∝ τ^(−4(p2+q2)) against τ⁻⁴), so the vacuum '
     'regimes of section 21.2 follow; in the Einstein regime it grows like τ^(2 − 3p1 − 5q1) = τ^1.160 (the auxiliary momenta are '
     'blue-shifted, pb → ε/4, pa → 0).')
text([('Decompactification. ', B),
      'Einstein gravity with a fluid pa = w_a ε, pb = w_b ε has the power-law solutions a ∝ τ^P, b ∝ τ^Q with 3P(1 + w_a) + 4Q(1 + w_b) = 2: '
      '(w_a, w_b) = (0, 1/4) → (P, Q) = (0, 2/5); (1/3, 0) → (1/2, 0); (1/7, 1/7) → (1/4, 1/4). Once the gas dominates, ordinary space '
      'stops expanding and the AUXILIARY SPACE EXPANDS, b ∝ τ^(2/5): without the stabilizing potential U(b) the 7+1 universe '
      'decompactifies. Section 22.1 shows how U(b) prevents this.'])
inp(r'''
Needs["NumericalDifferentialEquationAnalysis`"];
rIgl = GaussianQuadratureWeights[128, 0, Pi/2, 30]; rIth = N[rIgl[[All, 1]]]; rIwt = N[rIgl[[All, 2]] Cos[rIgl[[All, 1]]]^2 Sin[rIgl[[All, 1]]]^3];
rIwt = rIwt/Total[rIwt]; rIc2 = Cos[rIth]^2; rIs2 = Sin[rIth]^2;
rIgas[la_?NumericQ, lb_?NumericQ] := Module[{ia2 = Exp[-2 la], ib2 = Exp[-2 lb], root, pref = Exp[-3 la - 4 lb]},
   root = Sqrt[rIc2 ia2 + rIs2 ib2]; pref {rIwt.root, rIwt.(rIc2 ia2/(3 root)), rIwt.(rIs2 ib2/(4 root))}];
rIC00[ha_, hb_] := rIc[[1]] /. {alf[2] -> 1, alf[3] -> rIxi, Ha -> ha, Hb -> hb};
rIlinM[ha_, hb_] := {{Coefficient[rIc[[2]], dHa], Coefficient[rIc[[2]], dHb]}, {Coefficient[rIc[[3]], dHa], Coefficient[rIc[[3]], dHb]}} /. {alf[2] -> 1, alf[3] -> rIxi, Ha -> ha, Hb -> hb};
rIlinB[ha_, hb_] := ({rIc[[2]], rIc[[3]]} /. {dHa -> 0, dHb -> 0}) /. {alf[2] -> 1, alf[3] -> rIxi, Ha -> ha, Hb -> hb};
rIlinMf = Function[{ha, hb}, Evaluate[rIlinM[ha, hb]]]; rIlinBf = Function[{ha, hb}, Evaluate[rIlinB[ha, hb]]];
rIsolve2[m_, v_] := {v[[1]] m[[2, 2]] - v[[2]] m[[1, 2]], m[[1, 1]] v[[2]] - m[[2, 1]] v[[1]]}/(m[[1, 1]] m[[2, 2]] - m[[1, 2]] m[[2, 1]]);   (* 2x2 Cramer rule *)
rIgasD[ha_?NumericQ, hb_?NumericQ, la_?NumericQ, lb_?NumericQ, cg_] := Module[{g = rIgas[la, lb]},
   rIsolve2[rIlinMf[ha, hb], -rIlinBf[ha, hb] - cg {g[[2]], g[[3]]}]];
rIgasRhs[sg_?NumericQ, yy_?VectorQ, cg_] := Module[{ha = Exp[-sg], d, dx},
   d = rIgasD[ha, yy[[1]], yy[[2]], yy[[3]], cg]; dx = -ha/d[[1]]; {d[[2]] dx, ha dx, yy[[1]] dx, dx}];
rICg = 10^-6 3; {rIhbVac1, rIx1} = rIvacAt[1];
rIhb1 = Hb /. FindRoot[rIC00[1, Hb] == rICg, {Hb, rIhbVac1}, WorkingPrecision -> 20];
rIcpl[smax_] := NDSolveValue[{rIyC'[sV] == rIgasRhs[sV, rIyC[sV], rICg], rIyC[0] == N[{rIhb1, 0, 0, rIx1}]}, rIyC, {sV, 0, smax},
    Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 8}, PrecisionGoal -> 12, AccuracyGoal -> 30, MaxSteps -> 10^6, InterpolationOrder -> All];
rIcplB = rIcpl[-Log[10^16]]; rIcplF = rIcpl[-Log[10^-14]];
rIxs = rIcplB[-Log[10^16]][[4]] - 3/4 10^-16;
rIstate[h_] := Module[{sol = If[h >= 1, rIcplB, rIcplF], hb, la, lb, x, g, d, e2, e3, sz, hm},
   {hb, la, lb, x} = sol[-Log[h]]; g = rIgas[la, lb]; d = rIgasD[h, hb, la, lb, rICg]; hm = Abs[hb];
   e2 = -12 hb (4 h^3 + 18 h^2 hb + 12 h hb^2 + hb^3); e3 = 360 rIxi h^2 hb^3 (4 h + 3 hb);
   sz = {3 (h^2 + 4 h hm + 2 hm^2), 12 hm (4 h^3 + 18 h^2 hm + 12 h hm^2 + hm^3), 360 rIxi h^2 hm^3 (4 h + 3 hm), rICg g[[1]]};   (* sizes S1, S2, S3, kappa eps *)
   <|"Hb" -> hb, "tau" -> x - rIxs, "peff" -> -h^2/d[[1]], "keps" -> rICg g[[1]], "E2" -> e2, "E3" -> e3, "S" -> sz, "frac" -> sz[[4]]/Total[sz],
     "pa/e" -> g[[2]]/g[[1]], "pb/e" -> g[[3]]/g[[1]]|>];
rIhStar = -1/Sqrt[54 rIxi];
rInear = rIstate[10^15];
check["matter-dominated cubic regime: p_eff -> 3/4, Hb -> -Sqrt[alpha2/(54 alpha3)], kappa eps tau^3 -> 9 alpha2 |h*|, alpha2 E2 : alpha3 E3 : kappa eps -> 9/4 : -5/4 : 1 (Ha = 10^15)",
  Abs[rInear["peff"] - 3/4] < 10^-3 && Abs[rInear["Hb"]/rIhStar - 1] < 10^-3 && Abs[rInear["keps"] rInear["tau"]^3/(9 Abs[rIhStar]) - 1] < 10^-2 &&
  Abs[rInear["E2"]/rInear["keps"] - 9/4] < 10^-2 && Abs[rInear["E3"]/rInear["keps"] + 5/4] < 10^-2]
rIplat = rIstate /@ {3 10^5, 30, 10^-2};
check["coupled: where the gas is a small fraction the effective exponents are those of the vacuum (1, p2, p1 within 10^-2)",
  Abs[rIplat[[1]]["peff"] - 1] < 10^-2 && Abs[rIplat[[2]]["peff"] - rIp2] < 10^-2 && Abs[rIplat[[3]]["peff"] - rIp1] < 10^-2]
rIlate = rIstate[10^-14];
check["without the stabilizing potential the gas-dominated Einstein regime decompactifies: tau Ha -> 0, tau Hb -> 2/5, pb/eps -> 1/4 (Ha = 10^-14)",
  rIlate["tau"] 10^-14 < 10^-3 && Abs[rIlate["tau"] rIlate["Hb"] - 2/5] < 2 10^-3 && Abs[rIlate["pb/e"] - 1/4] < 10^-3]
rIgrow = With[{gA = rIstate[10^-2], gB = rIstate[10^-4]}, Log[(gB["keps"]/(3 10^-8))/(gA["keps"]/(3 10^-4))]/Log[gB["tau"]/gA["tau"]]];
check["the gas fraction grows like tau^(2 - 3 p1 - 5 q1) = tau^1.160 in the Einstein regime (to 2 10^-2)", Abs[rIgrow - (2 - 3 rIp1 - 5 rIq1)] < 2 10^-2]
rIfl[wa_, wb_] := Solve[{(3 (pP^2 + 4 pP qQ + 2 qQ^2)) == cC, (3 pP^2 + 8 pP qQ - 2 pP + 10 qQ^2 - 4 qQ) == -wa cC, 3 (2 pP^2 + 3 pP qQ - pP + 2 qQ^2 - qQ) == -wb cC,
     3 pP (1 + wa) + 4 qQ (1 + wb) == 2}, {pP, qQ, cC}];
check["Einstein + fluid power laws a ~ tau^P, b ~ tau^Q: (wa, wb) = (0, 1/4) -> (0, 2/5); (1/3, 0) -> (1/2, 0); (1/7, 1/7) -> (1/4, 1/4) (tau^2 E^(1) from EKas of section 14.1)",
  Expand[(EKas[[1]] /. {pK -> pP, qK -> qQ}) - {3 (pP^2 + 4 pP qQ + 2 qQ^2), 3 pP^2 + 8 pP qQ - 2 pP + 10 qQ^2 - 4 qQ, 3 (2 pP^2 + 3 pP qQ - pP + 2 qQ^2 - qQ)}] === {0, 0, 0} &&
  rIfl[0, 1/4] == {{pP -> 0, qQ -> 2/5, cC -> 24/25}} && rIfl[1/3, 0] == {{pP -> 1/2, qQ -> 0, cC -> 3/4}} && rIfl[1/7, 1/7] == {{pP -> 1/4, qQ -> 1/4, cC -> 21/16}}]
rIlead = {Coefficient[Expand[(rIc[[1]] /. {Ha -> pw/tS, Hb -> rIIIhh, dHa -> -pw/tS^2, dHb -> 0}) tS^3], tS, 0],
   Coefficient[Expand[(rIc[[2]] /. {Ha -> pw/tS, Hb -> rIIIhh, dHa -> -pw/tS^2, dHb -> 0}) tS^3], tS, 0],
   Coefficient[Expand[(rIc[[3]] /. {Ha -> pw/tS, Hb -> rIIIhh, dHa -> -pw/tS^2, dHb -> 0}) tS^4], tS, 0]};
check["matter-dominated cubic regime (any fluid, b finite): leading 44 term p^3 (p - 1)(648 alpha3 h^2 - 12 alpha2); with h^2 = alpha2/(54 alpha3): p = 1/(1 + w_a) and kappa eps tau^3 = (64/3) p^3 alpha2 |h|",
  Expand[rIlead[[3]] - pw^3 (pw - 1) (648 alf[3] rIIIhh^2 - 12 alf[2])] === 0 &&
  (Select[pw /. Solve[(rIlead[[2]] + wA rIlead[[1]] /. rIIIhh -> -Sqrt[alf[2]/(54 alf[3])]) == 0, pw], # =!= 0 &] // Simplify) == {1/(1 + wA)} &&
  Simplify[(rIlead[[1]] /. rIIIhh -> -Sqrt[alf[2]/(54 alf[3])]) - 64/3 pw^3 alf[2] Sqrt[alf[2]/(54 alf[3])], alf[2] > 0 && alf[3] > 0] === 0]
{"near the singularity (Ha = 10^15): p_eff, Hb, kappa eps tau^3, alpha2 E2/(kappa eps), alpha3 E3/(kappa eps)" -> {rInear["peff"], rInear["Hb"], rInear["keps"] rInear["tau"]^3, rInear["E2"]/rInear["keps"], rInear["E3"]/rInear["keps"]},
 "(p_eff, gas fraction) at Ha = 3 10^5, 30, 10^-2" -> ({#["peff"], #["frac"]} & /@ rIplat), "late (Ha = 10^-14): tau Ha, tau Hb, pa/eps, pb/eps" -> {rIlate["tau"] 10^-14, rIlate["tau"] rIlate["Hb"], rIlate["pa/e"], rIlate["pb/e"]},
 "Einstein-regime growth exponent: numerical, 2 - 3 p1 - 5 q1" -> {rIgrow, N[2 - 3 rIp1 - 5 rIq1]}}
''')

subsection('21.4 Pair creation in the primordial regime and the brane')
text('Section 14.3 computed the Bogoliubov coefficient |β|² of a massless mode on the Einstein–Kasner background as a function '
     'F(K) of the scaling variable K. For large K, F does NOT fall exponentially, as it would for a smooth expansion history, but as '
     'a power,')
text('      F(K) ∝ K^(−2(p1 − q1)/(1 − p1)) = K^(−(3√2 + 2)/2) = K^(−3.121…),')
text('because near x0 = 0 the phase ∫ω dτ ∝ K τ^(1−p1) while the direction of the mode Hamiltonian turns like τ^(p1−q1): the '
     'non-adiabaticity comes from the singularity itself. The comoving number of pairs created in each logarithmic interval of '
     'crossing time therefore grows toward τ = 0 (∝ τ⁻⁶ in the Kasner phase): the primordial production of dirac16complex is '
     'dominated by the earliest, Planckian times and gives a gas fraction of order one there — the gas used in section 21.3.')
text('The Standard Model plasma on the brane can only matter in the bulk when its smeared density ρ_4/V⊥ is comparable to the '
     'curvature terms. At T = M8 and H = M8 its fraction is (π² g*/30) M8⁴/(3 M_Pl² M8²) ≈ 2 × 10⁻²⁸ (M8 = 10 TeV, '
     'M_Pl = 2.435 × 10¹⁸ GeV, g* = 106.75): the brane is gravitationally negligible in the whole primordial regime.')
inp(r'''
rIsg1 = {{0, 1}, {1, 0}}; rIsg2 = {{0, -I}, {I, 0}};
rIbog[k_, p_, q_] := Module[{c, sol, t0 = 10^-15, t1 = 300, v0, v1, h0m, h1m},
   h0m = N[k t0^-p rIsg1 + k t0^-q rIsg2, 30]; h1m = N[k t1^-p rIsg1 + k t1^-q rIsg2, 30];
   v0 = Normalize[Eigenvectors[h0m][[Ordering[Eigenvalues[h0m]]]][[2]]]; v1 = Normalize[Eigenvectors[h1m][[Ordering[Eigenvalues[h1m]]]][[1]]];
   sol = NDSolveValue[{c'[sB] == -I Exp[sB] (k Exp[-p sB] rIsg1 + k Exp[-q sB] rIsg2).c[sB], c[Log[t0]] == N[v0]}, c[Log[t1]], {sB, Log[t0], Log[t1]},
     Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 8}, PrecisionGoal -> 11, AccuracyGoal -> 13, MaxSteps -> Infinity];
   Abs[Conjugate[N[v1]].sol]^2];
rIFk = {rIbog[1, N[rIp1], N[rIq1]], rIbog[8, N[rIp1], N[rIq1]], rIbog[16, N[rIp1], N[rIq1]]};
rItail = Log[rIFk[[3]]/rIFk[[2]]]/Log[2];
check["Einstein-Kasner pair creation: F(1) == 0.0107493 (section 14.3) and the tail F ~ K^(-(3 Sqrt[2] + 2)/2) (exponent between K = 8 and 16 to 2 10^-2)",
  Abs[rIFk[[1]] - 0.0107493] < 10^-6 && Abs[rItail + (3 Sqrt[2] + 2)/2] < 2 10^-2 && Simplify[2 (rIp1 - rIq1)/(1 - rIp1) - (3 Sqrt[2] + 2)/2] === 0]
rIbrane = (Pi^2 106.75/30) (10.^4)^4/(3 (2.435 10^18)^2 (10.^4)^2);
check["brane plasma at T = M8, H = M8: fraction (Pi^2 g*/30) M8^4/(3 M_Pl^2 M8^2) < 10^-27", rIbrane < 10^-27]
{"F(1), F(8), F(16)" -> rIFk, "tail exponent (K = 8 .. 16), -(3 Sqrt[2] + 2)/2" -> {rItail, N[-(3 Sqrt[2] + 2)/2]}, "brane fraction at T = H = M8" -> rIbrane}
''')

# =====================================================================================
section('22. Regime II: the modern regime (a stabilized auxiliary space, four-dimensional gravity on the brane)')
subsection('22.1 The stabilizing potential catches the auxiliary space: radion oscillations')
text('Section 21.3 showed that the auxiliary space must be held. Add the stabilizing potential U(b) = (K_U/2) (ln(b/b0))² '
     '(U_0 = 0, κK_U = 12 m_r², m_r = 10⁻³/ℓ2) to the coupled solution of section 21.3 at Ha = 10⁻² (Einstein–Kasner regime, τ ≈ 57), '
     'with b0 where the Kasner contraction would arrive at τ = 1/m_r. The equations are integrated in x0 (Ha now oscillates). '
     'When the expansion rate falls below m_r the radion β = ln(b/b0) is caught and oscillates about 0; its energy then behaves '
     'like four-dimensional matter, and the exact solution of β″ + (2/τ)β′ + m_r²β = 0 (the radion equation of section 20.3 in a '
     'matter-dominated universe, H_E = 2/(3τ)) is the accurate approximate solution')
text('      β ≈ A sin(m_r τ + ϑ)/(m_r (τ − τ_E)),   H_E = Ha + 2Hb ≈ 2/(3(τ − τ_E)),   corrections O(1/(m_r τ))   (amplitude A, phase ϑ).')
text('The dirac16complex gas (frozen auxiliary momenta = Kaluza–Klein masses) is then also four-dimensional matter: its share of the '
     'energy stays constant, and the radion carries the rest (99%): this is the cosmological moduli problem of the minimal model.')
inp(r'''
rIIm = 1/1000;
{rIIhb0, rIIla0, rIIlb0, rIIx0} = rIcplF[-Log[1/100]]; rIItau0 = rIIx0 - rIxs;
rIIlnb0 = rIIlb0 + rIq1 Log[1/(rIIm rIItau0)];
rIIrhs[tt_?NumericQ, yy_?VectorQ] := Module[{ha, hb, la, lb, g, bet, kU, d}, {ha, hb, la, lb} = yy; g = rIgas[la, lb]; bet = lb - rIIlnb0; kU = 6 rIIm^2 bet^2;
   d = rIsolve2[rIlinMf[ha, hb], -rIlinBf[ha, hb] - {rICg g[[2]] - kU, rICg g[[3]] - kU - 3 rIIm^2 bet}]; {d[[1]], d[[2]], ha, hb}];
rIIha0 = 1/100;
rIIhbC = Hb /. FindRoot[rIC00[rIIha0, Hb] == rICg rIgas[rIIla0, rIIlb0][[1]] + 6 rIIm^2 (rIIlb0 - rIIlnb0)^2, {Hb, rIIhb0}];
rIIcap = NDSolveValue[{rIIyK'[tt] == rIIrhs[tt, rIIyK[tt]], rIIyK[rIItau0] == N[{rIIha0, rIIhbC, rIIla0, rIIlb0}]}, rIIyK, {tt, rIItau0, 3 10^5},
    Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 8, "StiffnessTest" -> False}, PrecisionGoal -> 11, AccuracyGoal -> 16, MaxSteps -> 10^7, MaxStepSize -> 50];
rIIHEavg[t1_] := {Mean[Table[rIIcap[tt][[1]] + 2 rIIcap[tt][[2]], {tt, t1, t1 + 2 Pi/rIIm, 2 Pi/rIIm/2000}]], t1 + Pi/rIIm};
rIItE = With[{h = rIIHEavg[2.9 10^5]}, h[[2]] - 2/(3 h[[1]])];          (* H_E = 2/(3 (tau - tau_E)) over the last period *)
rIIdata = Table[{tt, rIIcap[tt][[4]] - rIIlnb0}, {tt, 10^5, 3 10^5, 50}];
rIIcs = LeastSquares[Table[{Sin[rIIm d[[1]]], Cos[rIIm d[[1]]]}, {d, rIIdata}], Table[d[[2]] rIIm (d[[1]] - rIItE), {d, rIIdata}]];
rIIfit = {rIIaF -> Norm[rIIcs], rIIphF -> ArcTan[rIIcs[[1]], rIIcs[[2]]], rIItEF -> rIItE};
rIIres = Max[Abs[rIIdata[[All, 2]] - ((rIIaF Sin[rIIm # + rIIphF]/(rIIm (# - rIItEF))) & /@ rIIdata[[All, 1]] /. rIIfit)]]/Max[Abs[rIIdata[[All, 2]]]];
rIIHE[t1_] := With[{h = rIIHEavg[t1]}, h[[1]] (h[[2]] - rIItE)];
rIIfr[tt_] := With[{y = rIIcap[tt]}, rICg rIgas[y[[3]], y[[4]]][[1]]/(3 (y[[1]] + 2 y[[2]])^2)];
rIIheVals = {rIIHE[10^5], rIIHE[2 10^5]}; rIIfrVals = {rIIfr[10^5], rIIfr[3 10^5]};
check["capture: the radion oscillates like a four-dimensional massive field in a matter-dominated universe: with tau_E from <H_E> = 2/(3 (tau - tau_E)) over the last period, beta = A sin(m_r tau + vartheta)/(m_r (tau - tau_E)) to 2 10^-2, <H_E>(tau - tau_E) = 2/3 to 10^-2 at tau = 10^5 and 2 10^5, the gas share constant to 2 10^-2",
  rIIres < 2 10^-2 && AllTrue[rIIheVals, Abs[# - 2/3] < 10^-2 &] && Abs[rIIfrVals[[2]]/rIIfrVals[[1]] - 1] < 2 10^-2]
{"start tau, beta(start)" -> {rIItau0, rIIlb0 - rIIlnb0}, "fit A, vartheta, tau_E; max residual/amplitude" -> {rIIaF, rIIphF, rIItEF, rIIres} /. rIIfit,
 "<H_E>(tau - tau_E) at tau = 10^5, 2 10^5" -> rIIheVals, "gas share at tau = 10^5, 3 10^5" -> rIIfrVals, "radion share 1 - gas share" -> 1 - rIIfrVals[[2]]}
''')
inp(r'''
Plot[{rIIcap[tt][[4]] - rIIlnb0, (rIIaF Sin[rIIm tt + rIIphF]/(rIIm (tt - rIItEF))) /. rIIfit}, {tt, 2 10^4, 10^5}, PlotPoints -> 400, Frame -> True,
 FrameLabel -> {"tau (units l2)", "beta = ln(b/b0)"}, PlotLegends -> {"coupled field equations (numerical)", "A sin(m_r tau + vartheta)/(m_r (tau - tau_E))"}, ImageSize -> 600]
''')

subsection('22.2 Radion decay, the moduli problem and the Kaluza–Klein channels')
text('In physical units (M8 = 10 TeV, reduced M_Pl = 2.435 × 10¹⁸ GeV, V⊥ = M_Pl²/M8⁶) the radion couples to the brane fields through '
     'L = (φ/(√3 M_Pl)) T^μ_μ (section 20.3). For a fermion of mass m_f (N_c colours) and, through the trace anomaly '
     'T^μ_μ ⊃ −(b_QCD α_s/(8π)) G^a_μν G^aμν (b_QCD = 11 − 2n_f/3, the first coefficient of the QCD beta function), for two gluons:')
text('      Γ(φ → f f̄) = N_c m_f² m_r (1 − 4m_f²/m_r²)^(3/2)/(24π M_Pl²),   Γ(φ → gg) = b_QCD² α_s² m_r³/(96π³ M_Pl²).')
text('For m_r = 10 GeV (b, c, τ, μ pairs and gluons, α_s = 0.18, n_f = 5) the brane channels give Γ_brane ≈ 4.6 × 10⁻³⁷ GeV, a lifetime '
     'of 1.4 × 10¹² s, and would reheat the brane only to T_RH = [Γ_brane M_Pl/√(π² g*/90)]^(1/2) ≈ 10⁻⁹ GeV (g* = 10.75), far below the '
     '4 MeV that nucleosynthesis needs.')
text([('Decay into the bulk. ', B),
      'The radion also couples to every Kaluza–Klein mode of the bulk fields. A mode of dirac16complex with Kaluza–Klein momentum '
      'μ_n = 2π|n|/ℓ⊥ (n in Z⁴, ℓ⊥ = V⊥^(1/4) the size of the auxiliary torus) has, in the four-dimensional Einstein frame '
      '(h = e^(−4β) g_E), the mass M_E = e^(−2β) (m² + μ_n² e^(−2β))^(1/2); it therefore couples to φ = √12 M_Pl β with '
      'g_n = (2M_n + μ_n²/M_n)/(√12 M_Pl), M_n = (m² + μ_n²)^(1/2). For m ≪ μ_n each Kaluza–Klein level (four brane Dirac spinors) gives'])
text('      Γ_n = 4 g_n² m_r (1 − 4M_n²/m_r²)^(3/2)/(8π) = 3 μ_n² m_r (1 − 4μ_n²/m_r²)^(3/2)/(8π M_Pl²),')
text('and the sum over the 1.17 × 10¹³ momenta with 2μ_n < m_r, done shell by shell with the number r_4(k) = 8σ(k) − 32σ(k/4) of '
     'lattice vectors with |n|² = k (σ the divisor sum, σ(k/4) = 0 unless 4 divides k), agrees with its continuum limit')
text('      Γ_bulk = Σ_n Γ_n ≈ (3 m_r/(8π M_Pl²)) (V⊥/(8π²)) (m_r/2)⁶ B(3, 5/2)/2 ≈ 6.0 × 10⁻²⁴ GeV   (lifetime 0.11 s):')
text('thirteen orders of magnitude faster than the brane channels (the Kaluza–Klein gravitons add channels of the same strength; for '
     'm > m_r/2 the dirac16complex channels close). The radion decays into the BULK: its energy goes to Kaluza–Klein modes, which return '
     'it to the brane only through couplings suppressed by the volume (a Kaluza–Klein graviton of mass M decays into brane fields at a '
     'rate of order M³/M_Pl²). This is the cosmological moduli problem of the minimal model in the form it takes with large extra '
     'dimensions. The minimal model therefore needs the radion and bulk energy left at the end of regime I to be removed before '
     'nucleosynthesis (for example by a later period of inflation on the brane, not part of the model); sections 22.3–22.4 start after '
     'such a reheating, below the normalcy temperature T* ≈ (M8⁶/M_Pl)^(1/5) ≈ 13 GeV.')
inp(r'''
rIIMPl = 2.435 10^18; rIIM8 = 10.^4; rIIVp = rIIMPl^2/rIIM8^6; rIIhbar = 6.582119569 10^-25;
rIIgff[mphi_, mf_, nc_] := If[4 mf^2 < mphi^2, nc mf^2 mphi (1 - 4 mf^2/mphi^2)^(3/2)/(24 Pi rIIMPl^2), 0];
rIIggg[mphi_, as_, nf_] := (11 - 2 nf/3)^2 as^2 mphi^3/(96 Pi^3 rIIMPl^2);
rIIgam = rIIgff[10., 4.183, 3] + rIIgff[10., 1.273, 3] + rIIgff[10., 1.77693, 1] + rIIgff[10., 0.1056583755, 1] + rIIggg[10., 0.18, 5];
rIITRH = Sqrt[rIIgam rIIMPl/Sqrt[Pi^2 10.75/90]];
rIINKK = rIIVp (Pi^2/2) (10./2)^4/(2 Pi)^4;
check["through brane fields alone a 10 GeV radion would decay far too late: T_RH < 10^-6 GeV (nucleosynthesis needs 4 10^-3 GeV); 10^13 Kaluza-Klein channels below m_r/2",
  rIITRH < 10^-6 && 10^13 < rIINKK < 2 10^13]
rIIkmax = Floor[(10. rIIVp^(1/4)/(4 Pi))^2];                        (* |n|^2 = k <= kmax  <=>  2 mu_n < m_r *)
rIIsig = DivisorSigma[1, Range[rIIkmax]];
rIIr4 = 8 rIIsig; rIIr4[[4 ;; ;; 4]] -= 32 rIIsig[[1 ;; Floor[rIIkmax/4]]];       (* r_4(k) = 8 sigma(k) - 32 sigma(k/4) *)
rIImuK = 2 Pi Sqrt[N[Range[rIIkmax]]]/rIIVp^(1/4);
rIIgBulk = 3 10./(8 Pi rIIMPl^2) Total[rIIr4 rIImuK^2 (1 - 4 rIImuK^2/100.)^(3/2)];
rIIgBulkC = 3 10./(8 Pi rIIMPl^2) rIIVp/(8 Pi^2) 5.^6 Beta[3, 5/2]/2;
check["the number of lattice vectors with |n|^2 = k is r_4(k) = 8 sigma(k) - 32 sigma(k/4) (direct count for k <= 60)",
  With[{nn = Total[Tuples[Range[-8, 8], 4]^2, {2}]}, Table[Count[nn, kk], {kk, 60}]] == rIIr4[[1 ;; 60]]]
check["the radion decays into the bulk: the lattice sum over the Kaluza-Klein modes of dirac16complex (m << m_r/2) equals its continuum limit to 10^-6, the number of momenta to 10^-5, and Gamma_bulk/Gamma_brane > 10^12",
  Abs[rIIgBulk/rIIgBulkC - 1] < 10^-6 && Abs[Total[rIIr4]/rIINKK - 1] < 10^-5 && rIIgBulk/rIIgam > 10^12]
{"Gamma (GeV), lifetime (s), T_RH (GeV)" -> {rIIgam, rIIhbar/rIIgam, rIITRH}, "Kaluza-Klein momenta below m_r/2" -> rIINKK,
 "radion into bulk dirac16complex pairs (m << m_r/2): Gamma_bulk lattice, continuum (GeV), lifetime (s), Gamma_bulk/Gamma_brane" -> {rIIgBulk, rIIgBulkC, rIIhbar/rIIgBulk, rIIgBulk/rIIgam},
 "V_perp (GeV^-4), l_perp (m), 2 Pi/l_perp (GeV)" -> {rIIVp, rIIVp^(1/4) 1.973269804 10^-16, 2 Pi/rIIVp^(1/4)}, "normalcy temperature (M8^6/M_Pl)^(1/5) (GeV)" -> (rIIM8^6/rIIMPl)^(1/5)}
''')

subsection('22.3 The Standard Model plasma on the brane: g*(T), g*s(T) and the time–temperature relation')
text('With b static the brane obeys the Friedmann equations of section 20.3, 3 M_Pl² Ha² = ρ_4 + M_Pl² Λ4 and '
     'Ha′ = −(ρ_4 + p_4)/(2M_Pl²). In the radiation era ρ_4 = (π²/30) g*(T) T⁴ and the entropy s = (2π²/45) g*s(T) T³ is conserved, '
     's a³ = constant, so that')
text('      Ha = (π² g*/90)^(1/2) T²/M_Pl,     ln a = −ln T − (1/3) ln g*s(T) + const,     t(T) = t(T1) + ∫ d(ln a)/Ha '
     '= t(T1) + ∫_T^T1 (dT′/T′) [1 + (1/3) d ln g*s/d ln T′]/Ha(T′),')
text('evaluated on 481 points equally spaced in ln T from T1 = 10 GeV (where t(T1) = 1/(2Ha)) to 10 keV; on each step 1/Ha is taken '
     'exponential in ln T (exact when g* is constant), so a step adds Δ(ln a) (1/H_(i+1) − 1/H_i)/ln(H_i/H_(i+1)); halving the grid '
     'changes t by less than 10⁻⁵.')
text('g* and g*s are computed here for ideal gases of all Standard Model particles with their masses (PDG values of section 17), '
     'with the series for each species i (x = m_i/T; +: bosons, −: fermions with alternating signs):')
text('      ρ_i = (g_i T⁴/(2π²)) Σ_{n≥1} (±1)^(n+1) [3(nx)² K_2(nx) + (nx)³ K_1(nx)]/n⁴,   p_i = (g_i T⁴/(2π²)) Σ_{n≥1} (±1)^(n+1) (nx)² K_2(nx)/n⁴,')
text('s_i = (ρ_i + p_i)/T. The QCD crossover is modelled by the switch S(T) = (1 + tanh((T − 0.160 GeV)/0.030 GeV))/2 between quarks '
     'and gluons and the lightest hadrons (π, K, η, ρ, ω, K*, η′, nucleons); this model differs from the lattice equation of state by '
     'up to about 10% near the crossover. The neutrinos decouple while e± are relativistic: T_ν = T (g_γe(T)/(11/2))^(1/3), and '
     'their density is scaled by N_eff/3 = 3.044/3. (The Python check tools/ref_coupled.py integrates the momentum integrals instead '
     'of summing the series.)')
inp(r'''
rIIsp = {{2, 0., -1, "always"}, {6, 80.3692, -1, "always"}, {3, 91.1880, -1, "always"}, {1, 125.20, -1, "always"}, {4, 0.51099895 10^-3, 1, "always"},
   {4, 0.1056583755, 1, "always"}, {4, 1.77693, 1, "always"}, {6, 0., 1, "nu"}, {16, 0., -1, "qgp"}, {12, 2.16 10^-3, 1, "qgp"}, {12, 4.70 10^-3, 1, "qgp"},
   {12, 93.5 10^-3, 1, "qgp"}, {12, 1.2730, 1, "qgp"}, {12, 4.183, 1, "qgp"}, {12, 172.56, 1, "qgp"}, {2, 0.13957039, -1, "hadron"}, {1, 0.1349768, -1, "hadron"},
   {2, 0.493677, -1, "hadron"}, {2, 0.497611, -1, "hadron"}, {1, 0.547862, -1, "hadron"}, {9, 0.77526, -1, "hadron"}, {3, 0.78266, -1, "hadron"},
   {12, 0.8955, -1, "hadron"}, {1, 0.95778, -1, "hadron"}, {4, 0.93827209, 1, "hadron"}, {4, 0.93956542, 1, "hadron"}};
rIInn = N[Range[400]];
rIIrs[t_, g_, m_, st_] := Module[{x = m/t, c, k2, k1},
   If[x > 60, Return[{0., 0.}]];
   c = If[st == 1, (-1.)^(rIInn + 1), 1. + 0 rIInn];
   If[x == 0, Return[With[{r = g t^4 6 Total[c/rIInn^4]/(2 Pi^2)}, {r, 4 r/(3 t)}]]];
   k2 = BesselK[2, rIInn x]; k1 = BesselK[1, rIInn x];
   With[{r = g t^4 Total[c (3 (rIInn x)^2 k2 + (rIInn x)^3 k1)/rIInn^4]/(2 Pi^2), p = g t^4 Total[c (rIInn x)^2 k2/rIInn^4]/(2 Pi^2)}, {r, (r + p)/t}]];
rIIS[t_] := (1 + Tanh[(t - 0.160)/0.030])/2;
rIIgge[t_] := (rIIrs[t, 2, 0., -1][[2]] + rIIrs[t, 4, 0.51099895 10^-3, 1][[2]]) 45/(2 Pi^2 t^3);
rIIgstars[t_?NumericQ] := Module[{tot = {0., 0.}, tn, w},
   Do[Which[sp[[4]] == "nu", tn = t Min[1, (rIIgge[t]/5.5)^(1/3)]; tot += 3.044/3 rIIrs[tn, sp[[1]], 0., 1],
      True, w = Switch[sp[[4]], "qgp", rIIS[t], "hadron", 1 - rIIS[t], _, 1]; tot += w rIIrs[t, sp[[1]], sp[[2]], sp[[3]]]], {sp, rIIsp}];
   {tot[[1]] 30/(Pi^2 t^4), tot[[2]] 45/(2 Pi^2 t^3)}];
rIIgtab = rIIgstars /@ {1000., 10., 1., 0.3, 0.1, 0.01, 0.001, 10.^-4, 10.^-6};
check["g* limits: 106.75 above the electroweak scale (to 0.1), 10.75 near 3 MeV (to 0.1), 2 + (21/4)(4/11)^(4/3)(3.044/3) = 3.3855 after e+e- annihilation (to 10^-3)",
  Abs[rIIgtab[[1, 1]] - 106.75] < 0.1 && Abs[rIIgstars[0.003][[1]] - 10.75] < 0.1 && Abs[rIIgtab[[9, 1]] - (2 + 21/4 (4/11)^(4/3) 3.044/3)] < 10^-3]
rIIlnT = Table[Log[10.] - i (Log[10.] - Log[10.^-5])/480, {i, 0, 480}];
rIIgg = rIIgstars[Exp[#]] & /@ rIIlnT;
rIIHg = MapThread[Sqrt[Pi^2 #2[[1]]/90] Exp[2 #1]/rIIMPl &, {rIIlnT, rIIgg}];
rIIlna = -rIIlnT - Log[rIIgg[[All, 2]]]/3;
rIItg = Join[{0.}, Accumulate[Differences[rIIlna] (1/Rest[rIIHg] - 1/Most[rIIHg])/Log[Most[rIIHg]/Rest[rIIHg]]]] + 1/(2 rIIHg[[1]]);   (* 1/Ha exponential in ln T on each step *)
rIIidx[t_] := Round[(Log[10.] - Log[t])/((Log[10.] - Log[10.^-5])/480)] + 1;
rIItvals = rIItg[[rIIidx[#]]] & /@ {1., 0.1, 0.01, 0.001, 10.^-4, 10.^-5};
rIItgHalf = With[{la = rIIlna[[1 ;; -1 ;; 2]], hg = rIIHg[[1 ;; -1 ;; 2]]}, Join[{0.}, Accumulate[Differences[la] (1/Rest[hg] - 1/Most[hg])/Log[Most[hg]/Rest[hg]]]] + 1/(2 hg[[1]])];
check["the time-temperature quadrature is converged: halving the grid changes t(10 keV) by less than 10^-5", Abs[Last[rIItgHalf]/Last[rIItg] - 1] < 10^-5]
check["radiation era: t(1 MeV) = 0.74 s within 5% (1/(2 Ha) with g* = 10.75)", Abs[rIItvals[[4]] rIIhbar/0.74 - 1] < 0.05]
{"g*, g*s at 1 TeV, 10 GeV, 1 GeV, 300 MeV, 100 MeV, 10 MeV, 1 MeV, 100 keV, 1 keV" -> rIIgtab, "t(T) in s at 1 GeV, 100 MeV, 10 MeV, 1 MeV, 100 keV, 10 keV" -> rIItvals rIIhbar,
 "1/(2 Ha(T)) / t(T)" -> (1/(2 rIIHg[[rIIidx[#]]]) & /@ {1., 0.1, 0.01, 0.001, 10.^-4, 10.^-5})/rIItvals}
''')
inp(r'''
ListLogLinearPlot[Transpose[{Exp[rIIlnT], #}] & /@ Transpose[rIIgg], Joined -> True, PlotLegends -> {"g*(T)", "g*s(T)"}, Frame -> True,
 FrameLabel -> {"T (GeV)", "relativistic degrees of freedom"}, ImageSize -> 600]
''')

subsection('22.4 The late universe: an accurate analytic solution, the age, and the radion today')
text('Below about 1 MeV the brane contains photons and three neutrino species (Ω_r from T_CMB = 2.7255 K and N_eff = 3.044), matter '
     '(Ω_m = 0.315, of which the cold dark matter would be dirac16complex if it is the dark matter) and Λ4 (Ω_Λ = 1 − Ω_m − Ω_r), with '
     'H0 = 67.4 km/s/Mpc (Planck 2018). The Friedmann equation Ha = H0 (Ω_r a⁻⁴ + Ω_m a⁻³ + Ω_Λ)^(1/2) has two exact solutions for two '
     'components:')
text('      t_rm(a) = (2/(3H0 Ω_m²)) [(Ω_m a − 2Ω_r) (Ω_r + Ω_m a)^(1/2) + 2Ω_r^(3/2)]   (radiation + matter),')
text('      t_mΛ(a) = (2/(3H0 Ω_Λ^(1/2))) arcsinh((Ω_Λ/Ω_m)^(1/2) a^(3/2))   (matter + Λ),')
text('and the composite t(a) = t_rm(a) + t_mΛ(a) − t_m(a), t_m = (2/(3H0 Ω_m^(1/2))) a^(3/2), is accurate to 2.2 × 10⁻⁴ at every '
     'time (the error is the product of the small radiation and Λ corrections). The age is 13.79 Gyr. During this regime the radion '
     'sits at β_eq = [Λ4 + κ f_W/V⊥ + (κ/2)(ε − 3pa) + κ pb − 12α_2 Ha²(Ha′ + Ha²)]/(3 m_r²) (section 20.3) — of order H0²/m_r² ≈ 10⁻⁸⁶ '
     'without brane tension, f_W/(3 M_Pl² m_r²) ≈ 6 × 10⁻²⁴ with f_W = M8⁴ (m_r = 10 GeV) — so G4 = G8/V⊥ does not change; the '
     'Gauss–Bonnet and cubic corrections are α_2 H0² ≈ 2 × 10⁻⁹² and α_3 H0⁴ ≈ 4 × 10⁻¹⁹⁰ (α_2 = 1/M8², α_3 = 10⁻⁶/M8⁴).')
inp(r'''
rIImpc = 3.0856775814913673 10^19; rIIgyr = 3.15576 10^16; rIIH0s = 67.4/rIImpc; rIIH0g = rIIH0s rIIhbar;
rIIrhoc = 3 rIIH0g^2 rIIMPl^2; rIITcmb = 2.7255 8.617333262 10^-14;
rIIOg = (Pi^2/15) rIITcmb^4/rIIrhoc; rIIOr = rIIOg (1 + 7/8 (4/11)^(4/3) 3.044); rIIOm = 0.315; rIIOL = 1 - rIIOm - rIIOr;
rIIH[a_] := rIIH0s Sqrt[rIIOr/a^4 + rIIOm/a^3 + rIIOL];
rIIage = NIntegrate[1/(a rIIH[a]), {a, 0, 1}, PrecisionGoal -> 12, MaxRecursion -> 30]/rIIgyr;
check["age of the universe from the brane Friedmann equation: 13.79 Gyr (Planck 2018: 13.787 +- 0.020)", Abs[rIIage - 13.79] < 0.02]
rIIt[a_] := NIntegrate[1/(xa rIIH[xa]), {xa, 0, a}, PrecisionGoal -> 12, MaxRecursion -> 30];
rIItrm[a_] := 2/(3 rIIH0s rIIOm^2) ((rIIOm a - 2 rIIOr) Sqrt[rIIOr + rIIOm a] + 2 rIIOr^(3/2));
rIItmL[a_] := 2/(3 rIIH0s Sqrt[rIIOL]) ArcSinh[Sqrt[rIIOL/rIIOm] a^(3/2)]; rIItm[a_] := 2/(3 rIIH0s Sqrt[rIIOm]) a^(3/2);
rIIcomp = Table[Abs[(rIItrm[a] + rIItmL[a] - rIItm[a])/rIIt[a] - 1], {a, 10^Range[-8, 0.5, 0.25]}];
check["composite analytic solution t(a) = t_rm(a) + t_mLambda(a) - t_m(a) accurate to < 3 10^-4 for 10^-8 <= a <= 10^0.5", Max[rIIcomp] < 3 10^-4]
{"Omega_r, Omega_Lambda, z_eq, age (Gyr)" -> {rIIOr, rIIOL, rIIOm/rIIOr - 1, rIIage}, "max error of the composite solution" -> Max[rIIcomp],
 "radion displacement today (m_r = 10 GeV): Lambda4 + matter, tension f_W = M8^4" -> {(3 rIIH0g^2 rIIOL + 1.5 rIIH0g^2 rIIOm)/(3 100.), (rIIM8^4/rIIMPl^2)/(3 100.)},
 "alpha2 H0^2, alpha3 H0^4" -> {rIIH0g^2/rIIM8^2, 10^-6 rIIH0g^4/rIIM8^4}}
''')

# =====================================================================================
section('23. Regime III: the far future')
subsection('23.1 The de Sitter attractor')
text('Once matter has diluted, Λ4 dominates and ordinary space approaches de Sitter expansion with H∞ = H0 Ω_Λ^(1/2) '
     '(1.81 × 10⁻¹⁸ s⁻¹, Gibbons–Hawking temperature H∞/(2π) ≈ 2.2 × 10⁻³⁰ K). For matter + Λ the solution is exact,')
text('      a(t) = (Ω_m/Ω_Λ)^(1/3) sinh^(2/3)(3H∞t/2),   Ha(t) = H∞ coth(3H∞t/2) → H∞ + 2H∞ e^(−3H∞t),')
text('and it agrees with the full solution (radiation included) to 10⁻⁴ already at a = 3. Pairs of dirac16complex are created with the '
     'Fermi–Dirac occupation 1/(e^(2πM/H∞) + 1) of section 14.5; 2πM/H∞ ≈ 5 × 10³³ for M = 1 eV, so none are.')
inp(r'''
rIIIH = rIIH0s Sqrt[rIIOL];
check["far future: Ha(t) = H_inf Coth[3 H_inf t/2] agrees with the full brane solution at a = 3 and a = 10 to 10^-4",
  AllTrue[{3, 10}, Abs[rIIIH Coth[3/2 rIIIH rIIt[#]]/rIIH[#] - 1] < 10^-4 &]]
{"H_inf (1/s, GeV), T_GH (K)" -> {rIIIH, rIIIH rIIhbar, rIIIH rIIhbar/(2 Pi)/(8.617333262 10^-14)}, "2 Pi M/H_inf for M = 1 eV, 1 GeV" -> {2 Pi 10^-9/(rIIIH rIIhbar), 2 Pi/(rIIIH rIIhbar)}}
''')

subsection('23.2 The radion in de Sitter: Gauss–Bonnet corrections and relaxation')
text('Linearize the full Einstein–Lovelock equations (with U = (K_U/2)β², κK_U = 12m_r²) around Ha = H, Hb = 0, β = β_eq, eliminating '
     'δHa with the constraint. The radion obeys β″ + γ β′ + m_eff² (β − β_eq) = 0 with (for β_eq → 0)')
text('      γ/(3H) = (1 − 16α_2H² + 48α_2²H⁴ − 48α_3H⁴ + 512α_2³H⁶ + 192α_2α_3H⁶ − 1024α_2⁴H⁸)/(1 − 4α_2H² + 32α_2²H⁴)²,')
text('      m_eff²/m_r² = 1/(1 − 4α_2H² + 32α_2²H⁴),     3 m_r² β_eq = 3H² − 12α_2H⁴ (+ O(β_eq²)).')
text('For α_2 = α_3 = 0: γ = 3H and m_eff = m_r. The radion is stable and relaxes to its shifted equilibrium. A numerical test in '
     'units ℓ2 = 1 (Λ = 3 × 10⁻⁴, m_r = 0.1, ξ = 10⁻⁶, start 0.05 above equilibrium) reaches the equilibrium (β_eq, H_eq) of the two '
     'algebraic equations 3H² = Λ + 6m_r²β² and 3m_r²β = 3H² − 12H⁴ to 10⁻⁶. In the physical far future α_2H∞² ≈ 10⁻⁹² and '
     'β_eq ≈ H∞²/m_r²: the auxiliary space is frozen forever.')
text('Without the stabilizing potential the far future would instead be eight-dimensional de Sitter, a = b = e^(Hx0), an exact vacuum '
     'solution when 21H² − 420α_2H⁴ + 2520α_3H⁶ = Λ: the auxiliary space would inflate and G4 ∝ b⁻⁴ would fall like e^(−4Hx0), '
     'which experiments exclude (|Ġ/G| is less than a few percent of H0).')
inp(r'''
rIIIeqC = (feTot[1] /. {Ha -> haV, Hb -> hbV}) - 6 m2V beV^2;
rIIIsl = First[Solve[{(feTot[2] /. {Ha -> haV, Hb -> hbV}) == 6 m2V beV^2, (feTot[5] /. {Ha -> haV, Hb -> hbV}) == 6 m2V beV^2 + 3 m2V beV}, {dHa, dHb}]];
rIIIFb = dHb /. rIIIsl;
rIIIbg = {haV -> hS, hbV -> 0, beV -> bqV, lamb -> (lamb /. First[Solve[(rIIIeqC /. {haV -> hS, hbV -> 0, beV -> bqV}) == 0, lamb]])};
rIIIcd = D[rIIIeqC, #] /. rIIIbg & /@ {haV, hbV, beV}; rIIIfd = D[rIIIFb, #] /. rIIIbg & /@ {haV, hbV, beV};
rIIIgam = Simplify[-(rIIIfd[[2]] - rIIIfd[[1]] rIIIcd[[2]]/rIIIcd[[1]])]; rIIImeff = Simplify[-(rIIIfd[[3]] - rIIIfd[[1]] rIIIcd[[3]]/rIIIcd[[1]])];
check["radion in de Sitter: gamma/(3H) and m_eff^2/m_r^2 as in the text (beta_eq -> 0); gamma = 3H and m_eff^2 = m_r^2 for alpha2 = alpha3 = 0",
  Simplify[(rIIIgam /. bqV -> 0)/(3 hS) - (1 - 16 alf[2] hS^2 + 48 alf[2]^2 hS^4 - 48 alf[3] hS^4 + 512 alf[2]^3 hS^6 + 192 alf[2] alf[3] hS^6 - 1024 alf[2]^4 hS^8)/(1 - 4 alf[2] hS^2 + 32 alf[2]^2 hS^4)^2] === 0 &&
  Simplify[(rIIImeff /. bqV -> 0)/m2V - 1/(1 - 4 alf[2] hS^2 + 32 alf[2]^2 hS^4)] === 0 &&
  Simplify[rIIIgam /. {bqV -> 0, alf[2] -> 0, alf[3] -> 0}] === 3 hS && Simplify[rIIImeff /. {bqV -> 0, alf[2] -> 0, alf[3] -> 0}] === m2V]
check["radion equilibrium in de Sitter: 3 m_r^2 beta_eq == 3 H^2 - 12 alpha2 H^4 + O(beta_eq^2)",
  Simplify[(bqV /. Solve[Normal[Series[rIIIFb /. rIIIbg, {bqV, 0, 1}]] == 0, bqV][[1]]) - (3 hS^2 - 12 alf[2] hS^4)/(3 m2V)] === 0]
rIIIlam = 3 10^-4; rIIImr = 1/10;
rIIIrhs[tt_?NumericQ, yy_?VectorQ] := Module[{ha, hb, la, lb, kU, d}, {ha, hb, la, lb} = yy; kU = 6 rIIImr^2 lb^2;
   d = rIsolve2[rIlinMf[ha, hb], -rIlinBf[ha, hb] + {rIIIlam + kU, rIIIlam + kU + 3 rIIImr^2 lb}]; {d[[1]], d[[2]], ha, hb}];
rIIIeqs = FindRoot[{3 rIIIhh^2 == rIIIlam + 6 rIIImr^2 rIIIbb^2, 3 rIIImr^2 rIIIbb == 3 rIIIhh^2 - 12 rIIIhh^4}, {{rIIIbb, 1/100}, {rIIIhh, 1/100}}, WorkingPrecision -> 30];
rIIIb0 = (rIIIbb /. rIIIeqs) + 1/20;
rIIIds = NDSolveValue[{rIIIyD'[tt] == rIIIrhs[tt, rIIIyD[tt]], rIIIyD[0] == N[{Sqrt[(rIIIlam + 6 rIIImr^2 rIIIb0^2)/3], 0, 0, rIIIb0}]}, rIIIyD, {tt, 0, 3000},
    Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 8, "StiffnessTest" -> False}, PrecisionGoal -> 11, AccuracyGoal -> 14, MaxSteps -> 10^7, MaxStepSize -> 1];
check["the radion relaxes to the equilibrium of the two algebraic equations: beta(3000) to 10^-6, Ha(3000) to 10^-8 (relative)",
  Abs[rIIIds[3000][[4]]/(rIIIbb /. rIIIeqs) - 1] < 10^-6 && Abs[rIIIds[3000][[1]]/(rIIIhh /. rIIIeqs) - 1] < 10^-8]
rIIIds8 = Expand[feTot[1] /. {Ha -> h8, Hb -> h8}];
check["eight-dimensional de Sitter (no stabilizer): all three equations read 21 H^2 - 420 alpha2 H^4 + 2520 alpha3 H^6 == Lambda",
  Expand[rIIIds8 - (21 h8^2 - 420 alf[2] h8^4 + 2520 alf[3] h8^6 - lamb)] === 0 &&
  AllTrue[{feTot[2], feTot[5]}, Expand[(# /. {Ha -> h8, Hb -> h8, dHa -> 0, dHb -> 0}) - rIIIds8] === 0 &]]
{"gamma/(3H), m_eff^2/m_r^2 (beta_eq -> 0)" -> {Factor[(rIIIgam /. bqV -> 0)/(3 hS)], Factor[(rIIImeff /. bqV -> 0)/m2V]},
 "de Sitter test: beta(3000), Ha(3000), beta_eq, H_eq" -> {rIIIds[3000][[4]], rIIIds[3000][[1]], rIIIbb /. rIIIeqs, rIIIhh /. rIIIeqs}}
''')

subsection('23.3 The size of every interaction in the three regimes')
text('The first table gives, at four points of the coupled primordial solution of section 21.3, the share of each term of the 00 '
     'equation: its size (S_1, S_2, S_3 of section 21.3 for the Einstein, Gauss–Bonnet and cubic terms, κε for the dirac16complex gas) '
     'divided by the sum of the four sizes. The second table collects, for each term of the action of section 20.1, its weight in the '
     'field equations of each regime, computed in sections 21–23 (M8 = 10 TeV, α_2 = 1/M8², α_3 = 10⁻⁶/M8⁴, m_r = 10 GeV, f_W = M8⁴; '
     '"estimate" marks an order-of-magnitude statement).')
inp(r'''
rIIIshare[h_] := With[{sz = rIstate[h]["S"]}, sz/Total[sz]];
rIIIbudget = {{"regime I: cubic (Ha = 3 10^5)", Sequence @@ rIIIshare[3 10^5]}, {"regime I: Gauss-Bonnet (Ha = 30)", Sequence @@ rIIIshare[30]},
   {"regime I: Einstein (Ha = 10^-2)", Sequence @@ rIIIshare[10^-2]}, {"regime I: near the singularity (Ha = 10^15)", Sequence @@ rIIIshare[10^15]}};
check["the shares of the four terms of the 00 equation: Einstein > 0.99 at Ha = 10^-2, Gauss-Bonnet > 0.9 at Ha = 30, Gauss-Bonnet and cubic about 1/2 each at Ha = 3 10^5, gas > 0.1 at Ha = 10^15",
  rIIIbudget[[3, 2]] > 0.99 && rIIIbudget[[2, 3]] > 0.9 && Abs[rIIIbudget[[1, 3]] - 1/2] < 0.01 && Abs[rIIIbudget[[1, 4]] - 1/2] < 0.01 && rIIIbudget[[4, 5]] > 0.1]
Grid[Prepend[Map[If[NumberQ[#], NumberForm[#, 4], #] &, rIIIbudget, {2}], {"", "Einstein S_1", "Gauss-Bonnet S_2", "cubic S_3", "dirac16complex gas"}], Frame -> All]
{"shares (Einstein, Gauss-Bonnet, cubic, gas) at Ha = 3 10^5, 30, 10^-2, 10^15" -> rIIIbudget[[All, 2 ;; 5]]}
''')
inp(r'''
Grid[{{"term of the action (section 20.1)", "regime I (primordial)", "regime II (today)", "regime III (far future)"},
  {"S_0: Lambda4/(3 H^2)", "~ 10^-92 (Lambda l2^2)", N[rIIOL, 4], "-> 1"},
  {"S_1: Einstein term", "dominant for tau >> l2", "dominant (four-dimensional Friedmann equations)", "dominant"},
  {"S_2: alpha2 H^2", "dominant for 10^-3 l2 < tau < l2", N[rIIH0g^2/rIIM8^2, 3], N[(rIIIH rIIhbar)^2/rIIM8^2, 3]},
  {"S_3: alpha3 H^4", "dominant for tau < 10^-3 l2", N[10^-6 rIIH0g^4/rIIM8^4, 3], N[10^-6 (rIIIH rIIhbar)^4/rIIM8^4, 3]},
  {"S_psi: dirac16complex (gas; matter)", "gas: first table; dominant near tau = 0", "Omega_m = 0.315 (0.265 if it is the dark matter)", "-> 0 like e^(-3 H t)"},
  {"S_m: m/H", "<< 1 (massless gas)", ">> 1 (cold matter)", ">> 1"},
  {"S_self: g_r n^(6/7) ~ (k_phys/M8)^6 (estimate)", "< 1 below the eight-dimensional Planck scale", "negligible", "negligible"},
  {"S_U: radion displacement beta", "catches b when H ~ m_r, then oscillates", N[(3 rIIH0g^2 rIIOL + 1.5 rIIH0g^2 rIIOm)/(3 100.), 3], N[(rIIIH rIIhbar)^2/100., 3]},
  {"S_f: brane tension kappa f_W/(3 H^2 V_perp)", N[rIIM8^2/(3 rIIMPl^2), 3], "part of Lambda4", "part of Lambda4"},
  {"S_SM: brane plasma", N[rIbrane, 3], N[rIIOr, 3], "-> 0 like e^(-4 H t)"},
  {"S_P: portals (estimate)", "negligible (on the brane)", "production of dirac16complex (freeze-in)", "negligible"},
  {"L_phiW: radion-brane coupling", "radion not yet caught", Row[{"reheating, T_RH = ", NumberForm[rIITRH, 3], " GeV (m_r = 10 GeV)"}], "-"},
  {"pair creation 2 Pi M/H (M = 1 eV)", "order 1 (Planckian times)", N[2 Pi 10^-9/rIIH0g, 3], N[2 Pi 10^-9/(rIIIH rIIhbar), 3]}}, Frame -> All]
''')
