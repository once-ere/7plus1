# Sections 17-19 of 7+1.nb (the Standard Model on the world brane of 7+1 spacetime, the portals of dirac16complex, and
# neutrino flavour conversion); appended on 26 October 2026.  Executed by tools/make_7p1_nb.py.

section('17. The Standard Model on the world brane of 7+1 spacetime')
text([('About sections 17–19. ', B),
      'Dirac_claude places the Standard Model on a seven-dimensional fold of its 4+4 spacetime. Here the same idea is rebuilt for 7+1: '
      'the four auxiliary coordinates x4 … x7 form a torus T⁴ (x_k ≃ x_k + Lc, coordinate period Lc, proper size b·Lc), and the Standard '
      'Model lives on the 3+1-dimensional WORLD BRANE W = {x4 = x5 = x6 = x7 = 0}, spanned by the time x0 and ordinary space x1, x2, x3. '
      'Gravity (sections 11–14) and dirac16complex (sections 7–10, 13, 15–16) live in the whole 7+1 spacetime (the bulk); every Standard '
      'Model field lives on W and feels gravity through the induced metric h_μν = g_μν(x, y = 0) (μ, ν = 0 … 3). The two sectors meet '
      'at W through portals. The total action is'])
text('      S = S_EL[g] + S_dirac16complex[ψ, g] + ∫_W d⁴x √|h| (L_SM + L_portal) + S_stab,')
text('with S_EL the Einstein–Lovelock action of section 12, the self-coupled covariant dirac16complex action of sections 13 and 15, '
     'L_SM the complete Standard Model Lagrangian below (gauge fields, Higgs doublet, three generations of quarks and leptons including '
     'right-handed neutrinos, Yukawa couplings), L_portal the couplings of section 18, and S_stab the stabilizing stress of section 14.4 '
     'that keeps b constant in the modern epoch. In eight-dimensional form ∫_W d⁴x √|h| (…) = ∫ d⁸X √|g| δ⁴(y)/b⁴ (…).')

subsection('17.1 The fields, their quantum numbers, and anomaly cancellation')
text('Gauge group SU(3)_C × SU(2)_L × U(1)_Y with gauge fields G^a_μ (a = 1 … 8), W^I_μ (I = 1, 2, 3), B_μ and couplings gs, g, g′. '
     'Matter in each of three generations (Weyl fields; Q = T3 + Y):')
item('Q_L = (u_L, d_L): (3, 2, 1/6);  u_R: (3, 1, 2/3);  d_R: (3, 1, −1/3);  L_L = (ν_L, e_L): (1, 2, −1/2);  e_R: (1, 1, −1);  ν_R: (1, 1, 0).')
item('Higgs doublet Φ: (1, 2, 1/2), Φ̃ = iσ2 Φ*.')
text('The gauge anomalies (with left-handed fields counted +, right-handed −) cancel generation by generation, with or without ν_R, and '
     'the number of SU(2) doublets per generation is even (no Witten anomaly). The relativistic degrees of freedom of the fields above, '
     'without ν_R, give g* = 28 + (7/8)·90 = 106.75.')
inp(r'''
smF = {{"Q", 3, 2, 1/6, 1}, {"uR", 3, 1, 2/3, -1}, {"dR", 3, 1, -1/3, -1}, {"L", 1, 2, -1/2, 1}, {"eR", 1, 1, -1, -1}, {"nuR", 1, 1, 0, -1}};
anom = <|"SU(3)^2 U(1)" -> Total[(#[[5]] #[[3]] #[[4]]) & /@ Select[smF, #[[2]] == 3 &]],
   "SU(2)^2 U(1)" -> Total[(#[[5]] #[[2]] #[[4]]) & /@ Select[smF, #[[3]] == 2 &]],
   "U(1)^3" -> Total[(#[[5]] #[[2]] #[[3]] #[[4]]^3) & /@ smF], "grav^2 U(1)" -> Total[(#[[5]] #[[2]] #[[3]] #[[4]]) & /@ smF],
   "SU(3)^3" -> Total[(#[[5]] #[[3]]) & /@ Select[smF, #[[2]] == 3 &]]|>;
check["anomaly cancellation per generation (with nu_R): SU(3)^2U(1), SU(2)^2U(1), U(1)^3, grav^2U(1), SU(3)^3 all vanish; doublets per generation even",
  Values[anom] == {0, 0, 0, 0, 0} && EvenQ[Total[Select[smF, #[[3]] == 2 &][[All, 2]]]]]
check["electric charges Q = T3 + Y: (u, d) = (2/3, -1/3), (nu, e) = (0, -1), u_R 2/3, d_R -1/3, e_R -1, nu_R 0",
  ({1/6 + 1/2, 1/6 - 1/2, -1/2 + 1/2, -1/2 - 1/2} == {2/3, -1/3, 0, -1})]
gStar = (8 2 + 3 2 + 2 + 4) + 7/8 3 Total[(2 #[[2]] #[[3]]) & /@ Select[smF, #[[1]] != "nuR" &]];
check["g* = 28 + (7/8) 90 = 106.75 for the Standard Model fields (without nu_R)", gStar == 427/4]
''')

subsection('17.2 The Lagrangian on the brane, in the gravitational field')
text('With the induced metric h_μν, its frame (tetrad) and spin connection, complex 4×4 Dirac matrices γ^μ(x) of the brane '
     '({γ^μ, γ^ν} = 2 h^μν) and the gauge-covariant derivatives D_μ = ∇_μ − i gs G^a_μ T^a − i g W^I_μ τ^I/2 − i g′ Y B_μ:')
text('      L_SM = −(1/4) G^a_μν G^aμν − (1/4) W^I_μν W^Iμν − (1/4) B_μν B^μν + (D_μΦ)†(D^μΦ) − μ² Φ†Φ − λ (Φ†Φ)²')
text('            + Σ_f fbar i γ^μ D_μ f − [Q̄_L Y_u Φ̃ u_R + Q̄_L Y_d Φ d_R + L̄_L Y_e Φ e_R + L̄_L Y_ν Φ̃ ν_R + h.c.],')
text('with G^a_μν = ∂_μG^a_ν − ∂_νG^a_μ + gs f^abc G^b_μ G^c_ν, W^I_μν = ∂_μW^I_ν − ∂_νW^I_μ + g ε^IJK W^J_μ W^K_ν, B_μν = ∂_μB_ν − ∂_νB_μ, '
     'f summing over all Weyl fields of the three generations, and 3×3 Yukawa matrices Y_u, Y_d, Y_e, Y_ν. Its field equations '
     '(∇_μ G^aμν = gs Σ_q q̄γ^ν T^a q + gauge self-coupling, and so on) are the ordinary ones with ∂ → ∇ (Levi-Civita connection of h).')

subsection('17.3 Electroweak symmetry breaking')
text('With ⟨Φ⟩ = (0, v/√2) the term |D_μ⟨Φ⟩|² is a mass matrix for (W¹, W², W³, B) with eigenvalues g²v²/4 (W±, twice), '
     '(g² + g′²)v²/4 (Z) and 0 (photon). The inputs G_F = 1.1663788 × 10⁻⁵ GeV⁻², m_W = 80.3692 GeV, m_Z = 91.1880 GeV, '
     'm_H = 125.20 GeV (PDG 2025) fix v = (√2 G_F)^(−1/2), g = 2m_W/v, cosθ_W = m_W/m_Z (on shell), g′ = g tanθ_W, e = g sinθ_W and '
     'λ = m_H²/(2v²) at tree level.')
inp(r'''
wF = {wW1, wW2, wW3, bB};       (* the components W^1, W^2, W^3 and B of the electroweak gauge fields *) tauP = {PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]};
dPhi = -I (gW Sum[wF[[i]] tauP[[i]]/2, {i, 3}] + gY/2 bB IdentityMatrix[2]).{0, vH/Sqrt[2]};
massTerm = Expand[ComplexExpand[Conjugate[dPhi].dPhi]];
mGB = Table[D[massTerm, wF[[i]], wF[[j]]], {i, 4}, {j, 4}];
check["gauge-boson mass matrix: eigenvalues g^2 v^2/4 (twice), (g^2 + g'^2) v^2/4 and 0 (photon)",
  Sort[Expand[Eigenvalues[mGB]]] == Sort[Expand[{gW^2 vH^2/4, gW^2 vH^2/4, (gW^2 + gY^2) vH^2/4, 0}]]]
gFermi = 1.1663788*^-5; mWb = 80.3692; mZb = 91.1880; mHb = 125.20;
vEW = (Sqrt[2] gFermi)^(-1/2); gEW = 2 mWb/vEW; cW = mWb/mZb; sW2 = 1 - cW^2; gpEW = gEW Sqrt[sW2]/cW; eEW = gEW Sqrt[sW2]; lamH = mHb^2/(2 vEW^2);
check["tree level: m_W == g v/2, m_Z == Sqrt[g^2 + g'^2] v/2, rho == m_W^2/(m_Z^2 cos^2 theta_W) == 1, m_H^2 == 2 lambda v^2",
  Abs[gEW vEW/2 - mWb] < 10^-12 && Abs[Sqrt[gEW^2 + gpEW^2] vEW/2 - mZb] < 10^-10 && Abs[mWb^2/(mZb^2 cW^2) - 1] < 10^-14 && Abs[2 lamH vEW^2 - mHb^2] < 10^-9]
{"v" -> vEW, "g" -> gEW, "g'" -> gpEW, "sin^2 thetaW" -> sW2, "e" -> eEW, "1/alpha tree" -> 4 Pi/eEW^2, "lambda" -> lamH}
''')

subsection('17.4 Quark and lepton masses; CKM and PMNS mixing')
text('Fermion masses are m = y v/√2 for the singular values y of the Yukawa matrices; the mixing matrices are V_CKM = U_uL† U_dL and '
     'U_PMNS = U_eL† U_νL. Both are written in the standard parametrization U(s12, s13, s23, δ) = R23 · U13(δ) · R12. Inputs: PDG 2024 '
     '(s12 = 0.22501, s13 = 0.003732, s23 = 0.04183, δ = 1.147); PDG 2025 quark masses (m_u = 2.16 MeV, m_d = 4.70 MeV, m_s = 93.5 MeV, '
     'm_c = 1.2730 GeV, m_b = 4.183 GeV, m_t = 172.56 GeV); NuFit-6.0, normal ordering, without Super-Kamiokande atmospheric data '
     '(sin²θ12 = 0.307, sin²θ23 = 0.561, sin²θ13 = 0.02195, δ = 177°, Δm²21 = 7.49 × 10⁻⁵ eV², Δm²31 = 2.534 × 10⁻³ eV²). With ν_R '
     'the neutrinos are Dirac particles; with the lightest mass m1 = 0: m2 = √Δm²21, m3 = √Δm²31.')
inp(r'''
mixU[s12_, s13_, s23_, d_] := Module[{c12 = Sqrt[1 - s12^2], c13 = Sqrt[1 - s13^2], c23 = Sqrt[1 - s23^2]},
   {{1, 0, 0}, {0, c23, s23}, {0, -s23, c23}}.{{c13, 0, s13 Exp[-I d]}, {0, 1, 0}, {-s13 Exp[I d], 0, c13}}.{{c12, s12, 0}, {-s12, c12, 0}, {0, 0, 1}}];
vCKM = mixU[0.22501, 0.003732, 0.04183, 1.147];
jCKM = Im[vCKM[[1, 2]] vCKM[[2, 3]] Conjugate[vCKM[[1, 3]]] Conjugate[vCKM[[2, 2]]]];
rbCKM = -vCKM[[1, 1]] Conjugate[vCKM[[1, 3]]]/(vCKM[[2, 1]] Conjugate[vCKM[[2, 3]]]);
check["CKM: unitary; Jarlskog invariant J = 3.12 x 10^-5 (PDG 2024: 3.12 +0.13 -0.12 x 10^-5); A, rhobar, etabar within the PDG uncertainties",
  Norm[vCKM.ConjugateTranspose[vCKM] - IdentityMatrix[3]] < 10^-15 && Abs[jCKM - 3.12*^-5] < 0.01*^-5 &&
  Abs[0.04183/0.22501^2 - 0.826] < 0.016 && Abs[Re[rbCKM] - 0.1591] < 0.0094 && Abs[Im[rbCKM] - 0.3523] < 0.0073]
mUp = {2.16*^-3, 1.2730, 172.56}; mDown = {4.70*^-3, 93.5*^-3, 4.183};
yU = Sqrt[2] DiagonalMatrix[mUp]/vEW; yD = Sqrt[2] vCKM.DiagonalMatrix[mDown]/vEW;
{uuL, suY, uuR} = SingularValueDecomposition[yU]; {udL, sdY, udR} = SingularValueDecomposition[yD];
check["Yukawa matrices: singular values x v/Sqrt[2] are the quark masses and |U_uL^dag.U_dL| == |V_CKM|",
  Norm[Sort[Diagonal[suY]] vEW/Sqrt[2] - Sort[mUp]] < 10^-12 && Norm[Sort[Diagonal[sdY]] vEW/Sqrt[2] - Sort[mDown]] < 10^-12 &&
  Norm[Abs[Reverse[Reverse /@ (ConjugateTranspose[uuL].udL)]] - Abs[vCKM]] < 10^-9]
uPMNS = mixU[Sqrt[0.307], Sqrt[0.02195], Sqrt[0.561], 177 Degree]; dm21 = 7.49*^-5; dm31 = 2.534*^-3;
mNu = {0, Sqrt[dm21], Sqrt[dm31]};
jPMNS = Im[uPMNS[[1, 2]] uPMNS[[2, 3]] Conjugate[uPMNS[[1, 3]]] Conjugate[uPMNS[[2, 2]]]];
check["PMNS unitary; Dirac neutrino Yukawa Sqrt[2] U.diag(m)/v reproduces the masses (singular values)",
  Norm[uPMNS.ConjugateTranspose[uPMNS] - IdentityMatrix[3]] < 10^-15 &&
  Norm[Sort[Diagonal[SingularValueDecomposition[Sqrt[2] uPMNS.DiagonalMatrix[mNu 10^-9]/vEW][[2]]]] vEW/Sqrt[2] 10^9 - Sort[mNu]] < 10^-15]
{"|V_CKM|" -> Abs[vCKM], "J_CKM" -> jCKM, "rhobar + i etabar" -> rbCKM, "|U_PMNS|" -> Abs[uPMNS], "J_PMNS" -> jPMNS, "neutrino masses (eV)" -> N[mNu], "y_t" -> Sqrt[2] 172.56/vEW}
''')

subsection('17.5 The modern epoch on the brane')
text('With b = b0 the 00 equation of section 14.4 is the Friedmann equation of the brane, 3Ha² = 8π G4 ρ4 + Λ4 with G4 = G8/V⊥, '
     'V⊥ = b0⁴ Lc⁴ (for matter on the brane, ε8 = ρ4/V⊥). Integrating it with the Planck 2018 parameters (H0 = 67.4 km/s/Mpc, '
     'Ω_m = 0.315, Ω_r = 9.14 × 10⁻⁵, Ω_Λ = 1 − Ω_m − Ω_r) gives the age of the universe. With a reduced eight-dimensional Planck '
     'mass M8 (8πG8 = M8⁻⁶) the reduced four-dimensional one is MPl² = M8⁶ V⊥ (2.435 × 10¹⁸ GeV); for M8 = 10 TeV this fixes V⊥ and the '
     'lowest Kaluza–Klein mass 2π/V⊥^(1/4) of a field on the auxiliary torus.')
inp(r'''
mpcKm = 3.0856775814913673*^19; gyrS = 3.15576*^16; h0P = 67.4; omM = 0.315; omR = 9.14*^-5; omL = 1 - omM - omR;
ageU = NIntegrate[1/(aa Sqrt[omR/aa^4 + omM/aa^3 + omL]), {aa, 0, 1}, PrecisionGoal -> 12]/(h0P/mpcKm)/gyrS;
check["age of the universe from the brane Friedmann equation: 13.79 Gyr (Planck 2018: 13.787 +- 0.020 Gyr)", Abs[ageU - 13.79] < 0.03]
check["ADD relation: kappa8 rho4/V_perp == 8 Pi G4 rho4 with G4 = G8/V_perp (kappa8 = 8 Pi G8)", Simplify[k8 r4/vP - 8 Pi (k8/(8 Pi)/vP) r4] === 0]
vPerp = (2.435*^18)^2/(1.0*^4)^6; hbarcGeVm = 1.973269804*^-16;
{"age (Gyr)" -> ageU, "V_perp (GeV^-4)" -> vPerp, "V_perp^(1/4) (m)" -> vPerp^(1/4) hbarcGeVm, "lowest KK mass (MeV)" -> 2 Pi/vPerp^(1/4) 10^3}
''')

# =====================================================================================
section('18. dirac16complex at the world brane: the 16 components as four brane Dirac spinors, and the portals')
text('At the brane the eight Dirac matrices split into four along W and four transverse to it. With complex 4×4 Dirac matrices γ^(4)_μ '
     '(Dirac representation, γ5 = iγ0γ1γ2γ3) and four anti-Hermitian 4×4 matrices ρ_k (k = 4 … 7), {ρ_k, ρ_l} = −2δkl '
     '(ρ4 = iγ0, ρ5 = −γ1, ρ6 = −γ2, ρ7 = −γ3), the sixteen matrices')
text('      R_μ = γ^(4)_μ ⊗ 1_4  (μ = 0 … 3),      R_k = γ5 ⊗ ρ_k  (k = 4 … 7)')
text('satisfy the same Clifford relations as Γ17[A]; both representations are irreducible, so there is a unitary U (Schur average over the '
     '256 Clifford products) with Γ17[A].U = U.R_A. Writing ψ = U.ζ with ζ ∈ ℂ⁴ ⊗ ℂ⁴, ζ = Σ_a ψ_a ⊗ e_a, dirac16complex at the brane is '
     'FOUR brane Dirac spinors ψ_a (a = 1 … 4, the transverse spinor index): ψbar ψ′ = Σ_a ψ̄_a ψ′_a and ψbar Γ^μ ψ′ = Σ_a ψ̄_a γ^μ ψ′_a. '
     'The eight-dimensional chirality is χ17 = U.(γ5 ⊗ ρ4567).U† (ρ4567 = ρ4ρ5ρ6ρ7, eigenvalues +1, +1, −1, −1): the left-handed half '
     'ψ_L contains brane-right-handed parts with ρ4567 = −1 and brane-left-handed parts with ρ4567 = +1.')
inp(r'''
g4D = {DiagonalMatrix[{1, 1, -1, -1}]};
Do[AppendTo[g4D, ArrayFlatten[{{0, PauliMatrix[i]}, {-PauliMatrix[i], 0}}]], {i, 3}];
g5D = I g4D[[1]].g4D[[2]].g4D[[3]].g4D[[4]];
rhoK = {I g4D[[1]], -g4D[[2]], -g4D[[3]], -g4D[[4]]};
rBr = Join[Table[KroneckerProduct[g4D[[mu + 1]], IdentityMatrix[4]], {mu, 0, 3}], Table[KroneckerProduct[g5D, rhoK[[k]]], {k, 4}]];
check["brane representation: R_A R_B + R_B R_A == 2 eta17_AB ID16 for R = (gamma_mu (x) 1_4, gamma5 (x) rho_k)",
  And @@ Flatten[Table[rBr[[A]].rBr[[B]] + rBr[[B]].rBr[[A]] == 2 \[Eta]17[[A, B]] ID16, {A, 8}, {B, 8}]]]
clifProd[gs_, K_] := Dot @@ Join[{ID16}, gs[[K + 1]]];
uSchur[seed_] := Sum[clifProd[Table[\[CapitalGamma]17[A], {A, 0, 7}], K].seed.Inverse[clifProd[rBr, K]], {K, mIdx}];
uInt = Module[{s, u}, Do[s = ConstantArray[0, {16, 16}]; s[[Quotient[j, 16] + 1, Mod[j, 16] + 1]] = 1; u = uSchur[s]; If[u =!= ConstantArray[0, {16, 16}], Break[]], {j, 0, 255}];
   u/Sqrt[Tr[ConjugateTranspose[u].u]/16]];
check["intertwiner: Gamma17[A].U == U.R_A for A = 0..7, and U is unitary", And @@ Table[\[CapitalGamma]17[A].uInt == uInt.rBr[[A + 1]], {A, 0, 7}] &&
  Simplify[ConjugateTranspose[uInt].uInt] == ID16]
rho4567 = rhoK[[1]].rhoK[[2]].rhoK[[3]].rhoK[[4]];
check["chi17 == U.(gamma5 (x) rho4567).U^dag; rho4567 Hermitian with eigenvalues {-1, -1, 1, 1}",
  Simplify[\[Chi]17 - uInt.KroneckerProduct[g5D, rho4567].ConjugateTranspose[uInt]] == ConstantArray[0, {16, 16}] &&
  ConjugateTranspose[rho4567] == rho4567 && Sort[Eigenvalues[rho4567]] == {-1, -1, 1, 1}]
{z1B, z2B} = {rndSpinor[], rndSpinor[]};
bar4[z_, a_] := Conjugate[Transpose[Partition[z, 4]][[a]]].g4D[[1]];
check["bilinears at the brane: psibar psi' == Sum_a psibar_a psi'_a and psibar Gamma^mu psi' == Sum_a psibar_a gamma^mu psi'_a (mu = 0..3)",
  Simplify[Conjugate[uInt.z1B].\[CapitalGamma]17[0].(uInt.z2B) - Sum[bar4[z1B, a].Transpose[Partition[z2B, 4]][[a]], {a, 4}]] === 0 &&
  And @@ Table[Simplify[Conjugate[uInt.z1B].\[CapitalGamma]17[0].\[CapitalGamma]17up[mu].(uInt.z2B) -
       Sum[bar4[z1B, a].(\[Eta]17[[mu + 1, mu + 1]] g4D[[mu + 1]]).Transpose[Partition[z2B, 4]][[a]], {a, 4}]] === 0, {mu, 0, 3}]]
''')

subsection('18.1 Kaluza–Klein masses on the auxiliary torus')
text('With b = b0 and ψ ∝ e^{i q⃗·y⃗} (q⃗ = 2π n⃗/Lc, n⃗ ∈ ℤ⁴) the eight-dimensional Dirac equation at the brane becomes, in the variables ζ, '
     '(γ·p ⊗ 1_4 + X − m) ζ = 0 with X = γ5 ⊗ (ρ·q⃗)/b0. X anticommutes with every γ^μ ⊗ 1_4 and X² = −|q⃗|²/b0², so each Fourier mode is '
     'four brane Dirac fermions of mass m_n = √(m² + |q⃗|²/b0²) (a Kaluza–Klein tower; the zero mode n⃗ = 0 has mass m).')
inp(r'''
mPs = 7/10; b0K = 13/10; qKK = 2 Pi {2, -1, 1/2, 3};
xKK = KroneckerProduct[g5D, qKK.rhoK]/b0K;
check["KK operator: X.X == -|q|^2/b0^2 ID16 and X anticommutes with gamma_mu (x) 1_4", Simplify[xKK.xKK + qKK.qKK/b0K^2 ID16] == ConstantArray[0, {16, 16}] &&
  And @@ Table[Simplify[xKK.rBr[[mu]] + rBr[[mu]].xKK] == ConstantArray[0, {16, 16}], {mu, 1, 4}]]
detKK[en_] := Det[N[KroneckerProduct[en g4D[[1]] - (3/10 g4D[[2]] - 2/10 g4D[[3]] + 4/10 g4D[[4]]), IdentityMatrix[4]] + xKK - mPs ID16, 40]];
eOn = Sqrt[mPs^2 + qKK.qKK/b0K^2 + (3/10)^2 + (2/10)^2 + (4/10)^2];
check["KK masses: Det[gamma.p (x) 1 + X - m] vanishes on p0^2 - |p|^2 == m^2 + |q|^2/b0^2 (40-digit arithmetic) and not off shell",
  Abs[detKK[eOn]] < 10^-25 Abs[detKK[11/10 eOn]] && Abs[detKK[11/10 eOn]] > 10^-3]
''')

subsection('18.2 The portals')
text('The portals of dirac16complex (an INVENTED interaction, as in Dirac_claude; separate couplings for the two eight-dimensional chiral '
     'halves) are local couplings at the brane, with ψ_a(x) the brane spinors of ψ(x, y = 0):')
item('Neutrino (Yukawa) portal: −Σ_{i,s} [y^L_is L̄_i Φ̃ P_R ψ^(L)_s + y^R_is L̄_i Φ̃ P_R ψ^(R)_s] + h.c., where ψ^(L)_s (ψ^(R)_s) are the '
     'brane-right-handed parts of the left (right) half of dirac16complex (two each, s = 1, 2). After symmetry breaking it is a Dirac mass '
     'μ_ia = y_ia v/(√2 √V⊥) between ν_L and every Kaluza–Klein mode (the factor 1/√V⊥ normalizes the modes on the torus).')
item('Vector (current) portal: −G_P^f (ψbar Γ^μ ψ)|_W (f̄ γ_μ f) for Standard Model fermions f; in a background of dirac16complex '
     'with brane density n_ψ it adds the potential V_P^f = G_P^f n_ψ to the forward scattering of f (for neutrinos: non-standard matter effects).')
item('Higgs (scalar) portal: −λ_P Φ†Φ (ψbar ψ)|_W, which shifts the dirac16complex mass by λ_P v²/(2 V⊥) after symmetry breaking.')
text('Coupled field equations (all fields): Σ_k α_k E^(k)μ_ν = κ (T_ψ + T_stab + δ⁴(y)/b⁴ T_W)^μ_ν; '
     '(i Γ^μ D_μ − m)ψ − Σ_r g_r Σ_K (ψbar Θ_K ψ) Θ^K ψ = δ⁴(y)/b⁴ × (portal sources: U.(Σ y* Φ̃† L ⊗ e_s + …)); the Standard Model '
     'field equations on W with the portal terms as additional sources. In the neutrino sector the Yukawa portal gives, in the basis '
     '(ν_L, ψ_L) × (ν_R, ψ_R), the Dirac mass matrix {{m_D, μ}, {0, m 1_4}} for the zero mode; with m ≫ μ the light mixing matrix '
     'N = (1 − η) U_PMNS is no longer unitary, η ≈ μ μ†/(2m²), and the heavy states are mostly dirac16complex.')
inp(r'''
mDnu = uPMNS.DiagonalMatrix[N[mNu]]; mPsi = 1.0*^9;
muPortal = Transpose[{{3*^4, 1*^4, -2*^4}, {1*^4, 2*^4, 1*^4}, {0, 0, 0}, {0, 0, 0}}];
mFull = ArrayFlatten[{{mDnu, muPortal}, {ConstantArray[0, {4, 3}], mPsi IdentityMatrix[4]}}];
{ulP, svP, vrP} = SingularValueDecomposition[mFull];
ordP = Ordering[Diagonal[svP]]; ulP = ulP[[All, ordP]]; svL = Diagonal[svP][[ordP]];
nAct = ulP[[1 ;; 3, 1 ;; 3]];
devNN = IdentityMatrix[3] - nAct.ConjugateTranspose[nAct];
check["portal mixing: 1 - N.N^dag == mu.mu^dag/m^2 to second order (non-unitary light mixing)", Norm[devNN - muPortal.ConjugateTranspose[muPortal]/mPsi^2] < 10^-14]
{"light masses (eV)" -> svL[[1 ;; 3]], "heavy masses (GeV)" -> svL[[4 ;; 7]] 10^-9, "(1 - N.N^dag)_ee" -> Re[devNN[[1, 1]]]}
''')

# =====================================================================================
section('19. Neutrino flavour and flavour conversion')
text('A neutrino of flavour α produced with energy E is the superposition Σ_i U*_αi |ν_i⟩. In vacuum, with ħc = 1.973269804 × 10⁻¹⁶ GeV m, '
     'the amplitude after a distance L is A(α → β) = Σ_i U_βi U*_αi exp(−i m_i² L/(2E)), so')
text('      P(α → β) = δαβ − 4 Σ_{i>j} Re X_ij sin²Δ_ij + 2 Σ_{i>j} Im X_ij sin 2Δ_ij,   X_ij = U*_αi U_βi U_αj U*_βj,   Δ_ij = Δm²_ij L/(4E),   Δm²_ij = m_i² − m_j².')
text('In matter of electron density n_e the forward scattering on electrons adds V_CC = √2 G_F n_e to the ν_e row of the flavour '
     'Hamiltonian H = U diag(0, Δm²21, Δm²31) U†/(2E) + diag(V_CC, 0, 0) (V_CC → −V_CC and U → U* for antineutrinos; a common potential '
     'drops out). In two flavours the mixing angle in matter obeys tan 2θ_m = sin 2θ/(cos 2θ − 2EV/Δm²): it is maximal at the MSW '
     'resonance 2EV = Δm² cos 2θ. In the Sun (n_e = 245 N_A e^{−10.54 r/R⊙} cm⁻³) the 8B neutrinos convert adiabatically, '
     'P_ee = c13⁴ (1 + cos 2θ12 cos 2θ_m(0))/2 + s13⁴; low-energy neutrinos approach the vacuum average c13⁴(1 − sin²2θ12/2) + s13⁴ '
     '(at 0.3 MeV the residual matter effect at the centre, 2EV/Δm²21 ≈ 0.15, still lowers P_ee by about 0.025).')
inp(r'''
hbarcGeVmN = 1.973269804*^-16; NAv = 6.02214076*^23; hbarcEVcm = 1.973269804*^-5;
phaseF[L_, en_] := 10^-18 L 10^3/(2 en hbarcGeVmN);
pVac[u_, L_, en_] := Abs[Transpose[u.DiagonalMatrix[Exp[-I {0, dm21, dm31} phaseF[L, en]]].ConjugateTranspose[u]]]^2;
pD = pVac[uPMNS, 1300., 2.5]; pDbar = pVac[Conjugate[uPMNS], 1300., 2.5];
xFormula[u_, a_, b_, L_, en_] := Boole[a == b] + Sum[With[{xx = Conjugate[u[[a, i]]] u[[b, i]] u[[a, j]] Conjugate[u[[b, j]]], dd = {0, dm21, dm31}[[i]] - {0, dm21, dm31}[[j]]},
     -4 Re[xx] Sin[dd phaseF[L, en]/2]^2 + 2 Im[xx] Sin[dd phaseF[L, en]]], {i, 3}, {j, i - 1}];
check["vacuum: P(a->b) doubly stochastic; equals delta - 4 Sum Re X sin^2 Delta + 2 Sum Im X sin 2Delta; CPT P(a->b) == P(bbar->abar)",
  Norm[Total[pD, {2}] - 1] < 10^-13 && Norm[Total[pD] - 1] < 10^-13 && Max[Abs[Table[xFormula[uPMNS, a, b, 1300., 2.5] - pD[[a, b]], {a, 3}, {b, 3}]]] < 10^-13 &&
  Norm[pD - Transpose[pDbar]] < 10^-13]
vCC[rho_, ye_] := Sqrt[2] gFermi 10^-18 rho ye NAv hbarcEVcm^3;
pMat[u_, L_, en_, rho_, ye_, anti_: False, extra_: {0, 0, 0}] := Module[{uu = If[anti, Conjugate[u], u], h, ev, vec},
   h = uu.DiagonalMatrix[{0, dm21, dm31}].ConjugateTranspose[uu]/(2 en 10^9) + If[anti, -1, 1] (DiagonalMatrix[{vCC[rho, ye], 0, 0}] + DiagonalMatrix[extra]);
   Abs[Transpose[MatrixExp[-I h L 10^3/(hbarcGeVmN 10^9)]]]^2];
pM = pMat[uPMNS, 1300., 2.5, 2.848, 0.5];
check["V_CC(rho Y_e = 1 g/cm^3) == 7.63 x 10^-14 eV; matter evolution doubly stochastic; a flavour-universal potential has no effect",
  Abs[vCC[1., 1.] - 7.63*^-14] < 0.01*^-14 && Norm[Total[pM, {2}] - 1] < 10^-12 && Norm[pMat[uPMNS, 1300., 2.5, 2.848, 0.5, False, {1*^-13, 1*^-13, 1*^-13}] - pM] < 10^-11]
th12 = ArcSin[Sqrt[0.307]];
thM[en_, V_] := ArcTan[Cos[2 th12] - 2 en 10^6 V/dm21, Sin[2 th12]]/2;
sinThNum[en_, V_] := Module[{h2 = dm21/(4 en 10^6) {{-Cos[2 th12], Sin[2 th12]}, {Sin[2 th12], Cos[2 th12]}} + DiagonalMatrix[{V, 0}], vv},
   vv = Eigensystem[h2][[2]][[Ordering[Eigenvalues[h2]]]]; Abs[Normalize[vv[[2]]][[1]]]^2];
eRes = dm21 Cos[2 th12]/(2 5*^-12)/10^6;
check["two-flavour MSW: sin^2 theta_m from diagonalization == formula, and sin^2 2theta_m == 1 at 2 E V == dm^2 cos 2theta",
  Max[Table[Abs[sinThNum[e, 5*^-12] - Sin[thM[e, 5*^-12]]^2], {e, {0.5, 2., 10.}}]] < 10^-12 && Abs[Sin[2 thM[eRes, 5*^-12]]^2 - 1] < 10^-12]
{"DUNE vacuum P(mu->e), P(mubar->ebar), P(mu->mu)" -> {pD[[2, 1]], pDbar[[2, 1]], pD[[2, 2]]},
 "DUNE matter P(mu->e), P(mubar->ebar)" -> {pM[[2, 1]], pMat[uPMNS, 1300., 2.5, 2.848, 0.5, True][[2, 1]]}}
''')
inp(r'''
rSun = 6.957*^8; s13sq = 0.02195; c13sq = 1 - s13sq;
vSun[r_] := Sqrt[2] gFermi 10^-18 245 NAv Exp[-10.54 r/rSun] hbarcEVcm^3;
peeSunNum[en_] := Module[{k = dm21/(4 en 10^6), sol, y, ps},
   sol = NDSolveValue[{y'[r] == -I/(hbarcGeVmN 10^9) ({{-k Cos[2 th12] + vSun[r] c13sq/2, k Sin[2 th12]}, {k Sin[2 th12], k Cos[2 th12] - vSun[r] c13sq/2}}.y[r]),
      y[0] == {1., 0.}}, y[rSun], {r, 0, rSun}, MaxSteps -> Infinity, PrecisionGoal -> 9, AccuracyGoal -> 11, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 8}];
   ps = Abs[{{Cos[th12], -Sin[th12]}, {Sin[th12], Cos[th12]}}.sol]^2;
   c13sq^2 (ps[[1]] Cos[th12]^2 + ps[[2]] Sin[th12]^2) + s13sq^2];
peeSunAd[en_] := c13sq^2 (1 + Cos[2 th12] Cos[2 ArcTan[Cos[2 th12] - 2 en 10^6 vSun[0] c13sq/dm21, Sin[2 th12]]/2])/2 + s13sq^2;
sunTab = Table[{e, peeSunNum[e], peeSunAd[e]}, {e, {0.3, 1., 5., 10.}}];
check["Sun: numerical P_ee == adiabatic MSW formula to 2 x 10^-3 at 0.3, 1, 5 and 10 MeV; P_ee(0.3 MeV) within 0.03 of the vacuum average; P_ee falls with energy",
  Max[Abs[sunTab[[All, 2]] - sunTab[[All, 3]]]] < 2*^-3 && Abs[sunTab[[1, 2]] - (c13sq^2 (1 - Sin[2 th12]^2/2) + s13sq^2)] < 0.03 &&
  OrderedQ[Reverse[sunTab[[All, 2]]]]]
rhoRes[dm_, th_, en_] := dm Cos[2 th]/(2 en 10^6 vCC[1., 0.5]);
tNu = (4/11)^(1/3) 2.7255; nNu = 3 Zeta[3]/(4 Pi^2) 2 (8.617333262*^-5 tNu/hbarcEVcm)^3;
check["cosmic neutrino background: T_nu == 1.945 K and n == 112 cm^-3 per flavour", Abs[tNu - 1.945] < 0.001 && Abs[nNu - 112] < 0.5]
pNonU[n_, ms_, L_, en_] := Abs[Transpose[n.DiagonalMatrix[Exp[-I ms^2 phaseF[L, en]]].ConjugateTranspose[n]]]^2;
pNU = pNonU[nAct, svL[[1 ;; 3]], 1300., 2.5];
check["dirac16complex portal: non-unitary oscillations, Sum_b P(a->b) < 1 (zero-distance effect); decoupled limit == three-flavour vacuum result",
  And @@ Table[Total[pNU[[a]]] < 1, {a, 3}] && Norm[pNonU[uPMNS, N[mNu], 1300., 2.5] - pD] < 10^-13]
pMP = pMat[uPMNS, 1300., 2.5, 2.848, 0.5, False, {0.2 vCC[2.848, 0.5], 0, 0}];
{"Sun P_ee {E, numerical, adiabatic}" -> sunTab, "SN resonance densities H, L (g/cm^3)" -> {rhoRes[dm31, ArcSin[Sqrt[s13sq]], 20.], rhoRes[dm21, th12, 20.]},
 "T_nu, n_nu" -> {tNu, nNu}, "DUNE P(mu->e) with vector portal V_P,e = 0.2 V_CC, without" -> {pMP[[2, 1]], pM[[2, 1]]}}
''')
