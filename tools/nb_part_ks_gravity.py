# Sections 15-16 of 7+1.nb (the generalized self-coupling of dirac16complex and its Kohn-Sham theory in gravitational
# fields); appended on 26 October 2026.  Executed by tools/make_7p1_nb.py (section(), subsection(), text(), item(), inp(), B).

section('15. The generalized self-coupling of dirac16complex: all nine ranks')
text([('About sections 15–16. ', B),
      'Section 8 coupled dirac16complex to itself through four contact interactions (scalar, pseudoscalar, vector, axial). The 256 '
      'products of section 8.6 contain nine Lorentz-invariant classes, one for each rank r = 0, 1, …, 8 (the number of different Dirac '
      'matrices in the product), so the most general contact self-coupling without derivatives has nine couplings g0, …, g8. Section 15 '
      'builds it, shows that the four couplings of section 8 are four of the nine, and computes the 9×9 Fierz matrix that turns every '
      'exchange (Fock) term into direct terms. Section 16 writes the Kohn–Sham theory of the self-coupled field in a gravitational field '
      '(using the spin connection of section 13) and solves it in the three epochs of section 14, including the back-reaction of the '
      'Kohn–Sham matter on the expansion through the Einstein–Lovelock equations.'])

subsection('15.1 The nine channels')
text('For a strictly increasing multi-index K = (A1 < A2 < … < Ar) of rank r let Γ^K = Γ17up[A1].….Γ17up[Ar] (Γ^() = ID16) and')
text('      Θ^K = i^σr Γ^K,   σr = (r(r−1)/2) mod 2,      Θ_K = ηK Θ^K,   ηK = η17_A1A1 ⋯ η17_ArAr.')
text('The phase makes every bilinear ψbar Θ^K ψ real (Γ17[0].Θ^K is Hermitian), and Θ_K.Θ^K = ID16, tr(Θ_K.Θ^J) = 16 δKJ. '
     'The invariant of rank r and the general self-coupling are')
text('      I_r = Σ_{|K| = r} (ψbar Θ_K ψ)(ψbar Θ^K ψ) = Σ_{|K| = r} ηK (ψbar Θ^K ψ)²,      L_self = −(1/2) Σ_{r=0}^{8} g_r I_r.')
text('The ranks have 1, 8, 28, 56, 70, 56, 28, 8, 1 members. The couplings of section 8 are g0 = gS, g1 = gV, g7 = −gA, g8 = −gP '
     '(I_8 = −(ψbar iχ17 ψ)², I_7 = −(axial invariant)); g2 … g6 are new (tensor channels of ranks 2 to 6). Each I_r is invariant under '
     'every local Lorentz transformation ψ → S ψ, ψbar → ψbar S⁻¹, so in a gravitational field it is written with the frame matrices '
     'exactly as in flat spacetime; Θ^K ω17 is, up to a phase, Θ of the complementary multi-index (Hodge duality).')
inp(r'''
mIdx = Flatten[Table[Subsets[Range[0, 7], {r}], {r, 0, 8}], 1];
rk[K_] := Length[K];
\[CapitalTheta]up[K_List] := \[CapitalTheta]up[K] = I^Mod[rk[K] (rk[K] - 1)/2, 2] (Dot @@ Join[{ID16}, Table[\[CapitalGamma]17up[A], {A, K}]]);
\[Eta]K[K_List] := Product[\[Eta]17[[A + 1, A + 1]], {A, K}];
\[CapitalTheta]dn[K_List] := \[CapitalTheta]dn[K] = \[Eta]K[K] \[CapitalTheta]up[K];
check["256 multi-indices in 9 ranks with 1, 8, 28, 56, 70, 56, 28, 8, 1 members", Length[mIdx] == 256 && Table[Count[mIdx, K_ /; rk[K] == r], {r, 0, 8}] == {1, 8, 28, 56, 70, 56, 28, 8, 1}]
check["every bilinear psibar Theta^K psi is real: Gamma17[0].Theta^K Hermitian for all 256 K",
  AllTrue[mIdx, ConjugateTranspose[\[CapitalGamma]17[0].\[CapitalTheta]up[#]] == \[CapitalGamma]17[0].\[CapitalTheta]up[#] &]]
check["Theta_K.Theta^K == ID16 and Tr[Theta_K.Theta^J] == 16 delta_KJ",
  AllTrue[mIdx, \[CapitalTheta]dn[#].\[CapitalTheta]up[#] == ID16 &] && Table[Tr[\[CapitalTheta]dn[K].\[CapitalTheta]up[J]], {K, mIdx}, {J, mIdx}] == 16 IdentityMatrix[256]]
check["Hodge duality: Theta^K.omega17 == (phase) Theta^(complement of K), |phase| = 1",
  AllTrue[mIdx, Abs[Tr[ConjugateTranspose[\[CapitalTheta]up[Complement[Range[0, 7], #]]].\[CapitalTheta]up[#].\[Omega]17]]/16 == 1 &]]
SeedRandom[1508];
rndSpinor[] := RandomInteger[{-4, 4}, 16] + I RandomInteger[{-4, 4}, 16];
{psR, pbR} = {rndSpinor[], rndSpinor[]};
Ichan[r_, pb_, ps_] := Sum[\[Eta]K[K] (pb.\[CapitalTheta]up[K].ps)^2, {K, Select[mIdx, rk[#] == r &]}];
bilR[o_] := pbR.o.psR;
check["the couplings of section 8 are four of the nine: I_0 == (psibar psi)^2, I_1 == (psibar Gamma^A psi)(psibar Gamma_A psi), I_8 == -(psibar i chi17 psi)^2, I_7 == -(axial invariant)",
  Ichan[0, pbR, psR] == bilR[ID16]^2 && Ichan[1, pbR, psR] == Sum[bilR[\[CapitalGamma]17up[A]] bilR[\[CapitalGamma]17[A]], {A, 0, 7}] &&
  Ichan[8, pbR, psR] == -bilR[I \[Chi]17]^2 && Ichan[7, pbR, psR] == -Sum[bilR[\[CapitalGamma]17up[A].\[Chi]17] bilR[\[CapitalGamma]17[A].\[Chi]17], {A, 0, 7}]]
wLL = RandomReal[{-1, 1}, {8, 8}]; wLL = wLL - Transpose[wLL];
sLL = MatrixExp[0.3 Sum[wLL[[A + 1, B + 1]] (\[CapitalGamma]17up[A].\[CapitalGamma]17up[B] - \[CapitalGamma]17up[B].\[CapitalGamma]17up[A])/4, {A, 0, 7}, {B, 0, 7}]/2];
sLLi = Inverse[sLL];
check["local Lorentz invariance: I_r(S psi, psibar S^-1) == I_r(psi, psibar) for r = 0..8 (random finite Lorentz transformation), and S^-1 == Gamma17[0].S^dagger.Gamma17[0]",
  And @@ Table[Abs[Ichan[r, pbR.sLLi, sLL.psR] - Ichan[r, pbR, psR]] < 10^-8 Abs[Ichan[r, pbR, psR]] + 10^-8, {r, 0, 8}] &&
  Norm[sLLi - \[CapitalGamma]17[0].ConjugateTranspose[sLL].\[CapitalGamma]17[0]] < 10^-10]
''')

subsection('15.2 The 9×9 Fierz matrix and Hartree–Fock with nine channels')
text('The completeness of the 256 matrices, Σ_K (Θ_K)_cd (Θ^K)_ab = 16 δad δbc, rearranges every product of two bilinears. For '
     'commuting spinors a, b, c, d (for anticommuting fields an extra minus sign appears)')
text('      Σ_{K∈r} (a Θ_K b)(c Θ^K d) = Σ_s F_rs Σ_{J∈s} (a Θ_J d)(c Θ^J b),     F_rs = (1/16) Σ_l (−1)^(rs−l) C(s, l) C(8−s, r−l),')
text('with the binomial coefficients C. F_rs = (1/256) Σ_{K∈r} tr(Θ_K Θ_J Θ^K Θ^J) for any J of rank s, F.F = 1 (five eigenvalues +1, '
     'four −1). With a density matrix R the direct (Hartree) energy is Σ_s (g_s/2) D_s, D_s = Σ_{J∈s} tr(Θ_J R) tr(Θ^J R), and the '
     'exchange (Fock) energy is −Σ_r (g_r/2) Σ_{K∈r} tr(Θ_K R Θ^K R) = −Σ_{r,s} (g_r/2) F_rs D_s. Hartree–Fock is therefore Hartree with '
     'the Fierz-transformed couplings')
text('      G = (1 − Transpose[F]) g,     G_s = g_s − Σ_r F_rs g_r.')
text('For uniform spin-saturated matter R = (ρS ID16 + ρ Γ17[0])/16 only D_0 = ρS² and D_1 = ρ² survive, so only G_0 and G_1 matter; with '
     'the couplings of section 8 (gS = −100000, gV = 10000, all others 0) G_0 = −98750 and G_1 = 20000, the Hartree–Fock couplings of '
     'section 8.7.')
inp(r'''
fierzF = Table[Sum[(-1)^(r s - l) Binomial[s, l] Binomial[8 - s, r - l], {l, 0, Min[r, s]}]/16, {r, 0, 8}, {s, 0, 8}];
fierzFromTraces[r_, J_] := Sum[Tr[\[CapitalTheta]dn[K].\[CapitalTheta]dn[J].\[CapitalTheta]up[K].\[CapitalTheta]up[J]], {K, Select[mIdx, rk[#] == r &]}]/256;
check["Fierz matrix: closed form == (1/256) Sum_{K in r} Tr[Theta_K.Theta_J.Theta^K.Theta^J] for the first and last J of every rank s",
  And @@ Flatten[Table[With[{Js = Select[mIdx, rk[#] == s &]}, fierzFromTraces[r, First[Js]] == fierzF[[r + 1, s + 1]] && fierzFromTraces[r, Last[Js]] == fierzF[[r + 1, s + 1]]], {r, 0, 8}, {s, 0, 8}]]]
check["F.F == IdentityMatrix[9]; eigenvalues: five +1 and four -1", fierzF.fierzF == IdentityMatrix[9] && Sort[Eigenvalues[fierzF]] == {-1, -1, -1, -1, 1, 1, 1, 1, 1}]
{aR, bR, cR, dR} = Table[rndSpinor[], {4}];
check["Fierz identity for commuting spinors (random): Sum_{K in r} (a Theta_K b)(c Theta^K d) == Sum_s F_rs Sum_{J in s} (a Theta_J d)(c Theta^J b), r = 0..8",
  And @@ Table[Sum[(aR.\[CapitalTheta]dn[K].bR) (cR.\[CapitalTheta]up[K].dR), {K, Select[mIdx, rk[#] == r &]}] ==
     Sum[fierzF[[r + 1, s + 1]] Sum[(aR.\[CapitalTheta]dn[J].dR) (cR.\[CapitalTheta]up[J].bR), {J, Select[mIdx, rk[#] == s &]}], {s, 0, 8}], {r, 0, 8}]]
RrF = RandomInteger[{-3, 3}, {16, 16}] + I RandomInteger[{-3, 3}, {16, 16}]; grF = RandomInteger[{-5, 5}, 9];
dS9[R_] := Table[Sum[Tr[\[CapitalTheta]dn[J].R] Tr[\[CapitalTheta]up[J].R], {J, Select[mIdx, rk[#] == s &]}], {s, 0, 8}];
check["Hartree-Fock energy E_H + E_x == Sum_s (G_s/2) D_s with G = (1 - Transpose[F]) g (random density matrix and couplings)",
  grF.dS9[RrF]/2 - Sum[grF[[r + 1]]/2 Sum[Tr[\[CapitalTheta]dn[K].RrF.\[CapitalTheta]up[K].RrF], {K, Select[mIdx, rk[#] == r &]}], {r, 0, 8}] ==
   ((IdentityMatrix[9] - Transpose[fierzF]).grF).dS9[RrF]/2]
gSymb = Array[gg, 9, 0]; Geff = Expand[(IdentityMatrix[9] - Transpose[fierzF]).gSymb];
check["uniform density matrix (rhoS ID16 + rho Gamma17[0])/16: D_0 == rhoS^2, D_1 == rho^2, D_2 .. D_8 == 0",
  Simplify[dS9[(rhoSs ID16 + rhos \[CapitalGamma]17[0])/16]] == Join[{rhoSs^2, rhos^2}, ConstantArray[0, 7]]]
check["section 8 recovered: with g = (gS, gV, 0, 0, 0, 0, 0, -gA, -gP) the exchange energy of uniform matter is that of section 8.7",
  Expand[Sum[(Geff[[s + 1]] - gSymb[[s + 1]])/2 {rhoSs^2, rhos^2, 0, 0, 0, 0, 0, 0, 0}[[s + 1]], {s, 0, 8}] /. {gg[0] -> gS, gg[1] -> gV, gg[7] -> -gA, gg[8] -> -gP, gg[_] -> 0}] ==
   Expand[-(gS/32) (rhoSs^2 + rhos^2) + (gP/32) (rhoSs^2 - rhos^2) - (gV/32) (8 rhoSs^2 - 6 rhos^2) + (gA/32) (8 rhoSs^2 + 6 rhos^2)]]
check["the Hartree-Fock couplings of section 8.7: gS = -100000, gV = 10000 give G_0 = -98750, G_1 = 20000",
  (Geff[[{1, 2}]] /. {gg[0] -> -100000, gg[1] -> 10000, gg[_] -> 0}) == {-98750, 20000}]
{"16 F" -> 16 fierzF, "G_0" -> Geff[[1]], "G_1" -> Geff[[2]], "example g = (-100000, 10000, 3000, 0, ...): {G_0, G_1}" -> (Geff[[{1, 2}]] /. {gg[0] -> -100000, gg[1] -> 10000, gg[2] -> 3000, gg[_] -> 0})}
''')

subsection('15.3 The covariant nonlinear field equation')
text('In a gravitational field the self-coupled field obeys, with the covariant derivative D_μ = ∂_μ + Ω_μ of section 13 and '
     'L = √|g| [ψbar (i Γ^μ D_μ − m) ψ + L_self],')
text('      (i Γ^μ(X) D_μ − m) ψ = Σ_{r=0}^{8} g_r Σ_{|K|=r} (ψbar Θ_K ψ) Θ^K ψ.')
text('The bilinears are Grassmann-even, so the left derivative with respect to ψbar_a acts on them as on commuting variables; the check '
     'differentiates L_self with symbolic ψbar, a random Gaussian-integer ψ and random integer couplings.')
inp(r'''
pbS = Array[pbar, 16]; psN = rndSpinor[]; gN9 = RandomInteger[{-5, 5}, 9];
LselfS = -Sum[gN9[[rk[K] + 1]]/2 \[Eta]K[K] (pbS.\[CapitalTheta]up[K].psN)^2, {K, mIdx}];
check["nonlinear field equation: dL_self/d psibar_a == -Sum_r g_r Sum_K (psibar Theta_K psi) (Theta^K psi)_a",
  Expand[D[LselfS, {pbS}] + Sum[gN9[[rk[K] + 1]] (pbS.\[CapitalTheta]dn[K].psN) \[CapitalTheta]up[K].psN, {K, mIdx}]] === ConstantArray[0, 16]]
''')

subsection('15.4 Chirality in a gravitational field')
text('χ17 commutes with every generator (1/4)[Γ17up[A], Γ17up[B]], hence with every spin connection Ω_μ: gravity does not mix the chiral '
     'halves. In the cosmological spacetime, with ψ = E_L φ_L + E_R φ_R (section 7.10) and (τO·∂)_(a,b) = (1/a) Σ_{i=1..3} τO[i] ∂_i + '
     '(1/b) Σ_{k=4..7} τO[k] ∂_k, the two coupled 8-component equations of section 7.11 become')
text('      (i ∂_0 + i (3Ha/2 + 2Hb) + (τO·∂)_(a,b)) φ_L = m φ_R,      (i ∂_0 + i (3Ha/2 + 2Hb) − (τO·∂)_(a,b)) φ_R = m φ_L.')
inp(r'''
check["chi17 commutes with every generator (Gamma17up[A].Gamma17up[B] - Gamma17up[B].Gamma17up[A])/4, hence with every Omega_mu",
  And @@ Flatten[Table[\[Chi]17.(\[CapitalGamma]17up[A].\[CapitalGamma]17up[B] - \[CapitalGamma]17up[B].\[CapitalGamma]17up[A]) ==
       (\[CapitalGamma]17up[A].\[CapitalGamma]17up[B] - \[CapitalGamma]17up[B].\[CapitalGamma]17up[A]).\[Chi]17, {A, 0, 7}, {B, A + 1, 7}]]]
phiLc = Table[lc[j][x0, x1, x5], {j, 8}]; phiRc = Table[rc[j][x0, x1, x5], {j, 8}];
psiLR = EL.phiLc + ER.phiRc;
dirLR = I Sum[fdCos["GamUp"][[m]].(D[psiLR, X[[m]]] + fdCos["Om"][[m]].psiLR), {m, 8}] - mD psiLR;
hubC = 3 a3'[x0]/(2 a3[x0]) + 2 b4'[x0]/b4[x0];
tauDab[f_] := Sum[\[Tau]O[i].D[f, X[[i + 1]]], {i, 1, 3}]/a3[x0] + Sum[\[Tau]O[k].D[f, X[[k + 1]]], {k, 4, 7}]/b4[x0];
check["curved chiral equations: E_L^dag.Gamma17[0].(i Gamma^mu D_mu - m) psi == (i d0 + i(3Ha/2 + 2Hb) + tauO.d) phi_L - m phi_R, and the same for E_R with -tauO.d",
  Simplify[ConjugateTranspose[EL].\[CapitalGamma]17[0].dirLR - (I D[phiLc, x0] + I hubC phiLc + tauDab[phiLc] - mD phiRc)] === ConstantArray[0, 8] &&
  Simplify[ConjugateTranspose[ER].\[CapitalGamma]17[0].dirLR - (I D[phiRc, x0] + I hubC phiRc - tauDab[phiRc] - mD phiLc)] === ConstantArray[0, 8]]
''')

# =====================================================================================
section('16. The Kohn–Sham theory of dirac16complex in gravitational fields')
text('Hohenberg–Kohn and Kohn–Sham in a gravitational field. For a STATIC metric the energy of the self-coupled field is a functional '
     'of its densities, which are defined with the frame matrices and integrated with √|g| d⁷x; the ground state minimizes it, and the '
     'Kohn–Sham orbitals φ_n solve the covariant Dirac equation with the mean fields,')
text('      (i Γ^μ(X) D_μ − m − Σ(X)) φ_n = 0,      Σ = Σ_{r} g_r Σ_{|K|=r} [tr(Θ_K R) Θ^K − Θ^K R Θ_K],      R = Σ_n occupied φ_n φbar_n,')
text('(direct minus exchange; with the Fierz matrix the exchange part is again of the direct form with G = (1 − Transpose[F]) g). '
     'For a TIME-DEPENDENT metric (cosmology) the same equations with the time-dependent density matrix R(X) are the time-dependent '
     'Kohn–Sham equations; their adiabatic limit — every orbital follows the instantaneous positive-energy eigenstate — holds when the '
     'expansion rates are small compared with the orbital energies (section 14.3 measures the departures). For homogeneous, spin-saturated '
     'matter in the cosmological spacetime the orbitals are the modes of section 13.2 with the comoving momenta k = (p⃗, q⃗) inside a '
     'comoving Fermi ball |k| < kF, and the mean fields reduce to an effective mass M(x0) = m + G_S ρS and a vector potential V0 = G_V ρ:')
text('      i dc/dx0 = (αD·(p⃗/a, q⃗/b) + βD M(x0) + V0(x0)) c.')
text('The physical densities are integrals over the comoving ball, with |k_a| = k cosθ (ordinary space) and |k_b| = k sinθ (auxiliary '
     'space), the solid-angle weight cos²θ sin³θ, and ω = √(k² cos²θ/a² + k² sin²θ/b² + M²):')
text('      ⟨f⟩ = (1/(2π⁴ a³ b⁴)) ∫_0^kF dk k⁶ ∫_0^(π/2) dθ cos²θ sin³θ f,    ρ = ⟨1⟩ = kF⁷/(105 π⁴ a³ b⁴),   ρS = ⟨M/ω⟩,')
text('      ε = ⟨ω⟩ − (M − m)²/(2 G_S) + (G_V/2) ρ²,   pa = ⟨k² cos²θ/(3 a² ω)⟩ + Π,   pb = ⟨k² sin²θ/(4 b² ω)⟩ + Π,   Π = (M − m)²/(2 G_S) + (G_V/2) ρ²,')
text('with the gap equation M = m + G_S ρS(M). T^μ_ν = diag(ε, −pa, −pa, −pa, −pb, −pb, −pb, −pb) is conserved, '
     'ε′ + 3Ha (ε + pa) + 4Hb (ε + pb) = 0, for every expansion history (a, b): the stationarity of ε in M removes the term in M′. '
     'The integrals are evaluated with Gauss–Legendre quadrature (96 nodes in k, 64 in θ). Units: m = 1; couplings of section 8.7 '
     '(G_S = −100000, G_V = 10000, Hartree); comoving kF = kFsat = 0.756636… (the saturation density of section 8.7 at a = b = 1).')

subsection('16.1 Homogeneous Kohn–Sham matter in the cosmological spacetime')
inp(r'''
Needs["NumericalDifferentialEquationAnalysis`"];
glK = GaussianQuadratureWeights[96, 0, 1, 30]; glT = GaussianQuadratureWeights[64, 0, Pi/2, 30];
uK = N[glK[[All, 1]]]; wuK = N[glK[[All, 2]]]; thG = N[glT[[All, 1]]]; wthG = N[glT[[All, 2]]];
kFsat = 0.75663615070149289469`20; gSks = -100000.; gVks = 10000.; mks = 1.;
wGrid[kF_] := wGrid[kF] = kF^7 Outer[Times, wuK uK^6, wthG Cos[thG]^2 Sin[thG]^3];
ksMoments[M_?NumericQ, a_?NumericQ, b_?NumericQ, kF_?NumericQ] := Module[{w = wGrid[kF], pa2, pb2, om, pref},
   pa2 = Outer[Times, (kF uK)^2, Cos[thG]^2]/a^2; pb2 = Outer[Times, (kF uK)^2, Sin[thG]^2]/b^2;
   om = Sqrt[pa2 + pb2 + M^2]; pref = 1/(2 Pi^4 a^3 b^4);
   <|"rho" -> pref Total[w, 2], "rhoS" -> pref Total[w M/om, 2], "ekin" -> pref Total[w om, 2], "pa" -> pref Total[w pa2/om, 2]/3, "pb" -> pref Total[w pb2/om, 2]/4|>];
ksGap[a_?NumericQ, b_?NumericQ, kF_?NumericQ] := Mx /. FindRoot[Mx - mks - gSks ksMoments[Mx, a, b, kF]["rhoS"], {Mx, 10^-14, mks}, Method -> "Brent", AccuracyGoal -> 15, PrecisionGoal -> 15];
ksStress[a_?NumericQ, b_?NumericQ, kF_?NumericQ] := Module[{M = ksGap[a, b, kF], mo, pint},
   mo = ksMoments[M, a, b, kF]; pint = (M - mks)^2/(2 gSks) + gVks/2 mo["rho"]^2;
   <|"M" -> M, "rho" -> mo["rho"], "rhoS" -> mo["rhoS"], "eps" -> mo["ekin"] - (M - mks)^2/(2 gSks) + gVks/2 mo["rho"]^2, "pa" -> mo["pa"] + pint, "pb" -> mo["pb"] + pint|>];
st0ks = ksStress[1., 1., kFsat];
check["flat limit a = b = 1, kF = kFsat: m* == 0.3480973401... and E/A - m == -0.0283214279... (section 8.7), to 10^-11",
  Abs[st0ks["M"] - 0.34809734014624326514] < 10^-11 && Abs[st0ks["eps"]/st0ks["rho"] - 1 + 0.028321427964738575854] < 10^-11]
check["flat limit: pa == pb == rho eps_F - eps (Gibbs-Duhem), which is 0 at the saturation density",
  Abs[st0ks["pa"] - st0ks["pb"]] < 10^-15 && Abs[st0ks["pa"] - (st0ks["rho"] (gVks st0ks["rho"] + Sqrt[kFsat^2 + st0ks["M"]^2]) - st0ks["eps"])] < 10^-15]
aHist[x_] := 1 + 0.3 x + 0.05 x^2; bHist[x_] := 1/(1 + 0.2 x);
epsHist[x_?NumericQ] := ksStress[aHist[x], bHist[x], kFsat];
hst = 10.^-3; xcs = 0.7;
deps = (epsHist[xcs - 2 hst]["eps"] - 8 epsHist[xcs - hst]["eps"] + 8 epsHist[xcs + hst]["eps"] - epsHist[xcs + 2 hst]["eps"])/(12 hst);
consKS = deps + 3 aHist'[xcs]/aHist[xcs] (epsHist[xcs]["eps"] + epsHist[xcs]["pa"]) + 4 bHist'[xcs]/bHist[xcs] (epsHist[xcs]["eps"] + epsHist[xcs]["pb"]);
check["energy-momentum conservation of the Kohn-Sham matter: eps' + 3Ha(eps + pa) + 4Hb(eps + pb) == 0 (a = 1 + 0.3x + 0.05x^2, b = 1/(1 + 0.2x), x0 = 0.7; relative 10^-8)",
  Abs[consKS/deps] < 10^-8]
{"m*" -> st0ks["M"], "rho" -> st0ks["rho"], "eps" -> st0ks["eps"], "p" -> st0ks["pa"], "conservation residual / eps'" -> consKS/deps}
''')

subsection('16.2 The three epochs')
text('Primordial: the Einstein–Kasner background a = x0^0.547, b = x0^(−0.160) of section 14.1, normalized so that a = b = 1 at x0 = 1. '
     'The density grows like 1/x0 (a³b⁴ = x0), but the ordinary momenta are blue-shifted (p/a → ∞), so the scalar density per particle '
     'M/ω falls and the attractive scalar coupling drives the effective mass toward zero: chiral symmetry is effectively restored near '
     'x0 = 0. Modern: a static auxiliary space (b = 1) and growing a: the ordinary momenta redshift, M → m, ε a³ → constant and '
     'pa/ε → 0 — dirac16complex behaves as cold (pressureless) matter in ordinary space, while a stress pb remains along the auxiliary '
     'space, which the stabilizing stress of section 14.4 must balance. Far future: de Sitter a = e^{0.1 x0}, b = 1: ρ = ρ0 e^{−0.3 x0} '
     'and 1 − M decays like ρ.')
inp(r'''
pKs = N[(1 + 2 Sqrt[2])/7, 20]; qKs = N[(2 - 3 Sqrt[2])/14, 20];
primKS = Table[x -> ksStress[x^pKs, x^qKs, kFsat], {x, {10.^-4, 10.^-3, 10.^-2, 10.^-1, 1., 10.}}];
check["primordial: rho x0 == kFsat^7/(105 Pi^4) at every x0, and m* decreases monotonically to below 0.05 at x0 = 10^-4",
  And @@ Table[Abs[(r[[2]]["rho"]) r[[1]] - kFsat^7/(105 Pi^4)] < 10^-12, {r, primKS}] && OrderedQ[primKS[[All, 2]][[All, Key["M"]]]] && primKS[[1, 2]]["M"] < 0.05]
modKS = Table[aa -> ksStress[aa, 1., kFsat], {aa, {1., 2., 4., 10., 100.}}];
check["modern (b = 1): eps a^3 tends to a constant and pa/eps -> 0 (pressureless in ordinary space)",
  Abs[(modKS[[5, 2]]["eps"] 100.^3)/(modKS[[4, 2]]["eps"] 10.^3) - 1] < 0.05 && modKS[[5, 2]]["pa"]/modKS[[5, 2]]["eps"] < 10^-3]
futKS = Table[x -> ksStress[Exp[0.1 x], 1., kFsat], {x, {0., 10., 20., 40., 80.}}];
check["far future (de Sitter, H = 0.1): rho == rho0 e^{-3 H x0} and 1 - m* decays like rho",
  Abs[futKS[[4, 2]]["rho"]/futKS[[1, 2]]["rho"] - Exp[-12.]] < 10^-12 && Abs[(1 - futKS[[5, 2]]["M"])/(1 - futKS[[4, 2]]["M"])/Exp[-12.] - 1] < 0.05]
{"primordial m*(x0), x0 = 1e-4 .. 10" -> primKS[[All, 2]][[All, Key["M"]]], "primordial pa/pb at x0 = 1e-3, 1, 10" -> (#["pa"]/#["pb"] & /@ primKS[[{2, 5, 6}, 2]]),
 "modern m*(a), a = 1, 2, 4, 10, 100" -> modKS[[All, 2]][[All, Key["M"]]], "modern pb/eps at a = 100" -> modKS[[5, 2]]["pb"]/modKS[[5, 2]]["eps"],
 "far future m*(x0), x0 = 0, 10, 20, 40, 80" -> futKS[[All, 2]][[All, Key["M"]]]}
''')

subsection('16.3 Self-consistent evolution: Einstein–Lovelock gravity driven by the Kohn–Sham matter')
text('With T^μ_ν of the Kohn–Sham matter the field equations of section 14 are a closed system for a, b. The 11 and 44 equations are '
     'LINEAR in Ha′ and Hb′ (Lovelock equations are quasi-linear), so they are solved for (Ha′, Hb′) at every step; the 00 equation '
     'contains no second derivatives and is a CONSTRAINT, imposed on the initial data (Hb(0) from Ha(0)) and then preserved by the '
     'evolution because both sides are conserved (Bianchi identity of section 12, conservation law of section 16.1). The runs below '
     '(κ = 1000, Λ = 0, Ha(0) = 0.05, a(0) = b(0) = 1, x0 from 0 to 20) check that preservation for Einstein gravity and for '
     'Einstein–Lovelock gravity with α_2 = 20, α_3 = 500.')
inp(r'''
eTot[c_] := Sum[alf[k] toH[ECos[[k + 1, c, c]]], {k, 1, 3}] - lamb;
coefLin[c_] := {Coefficient[eTot[c], dHa], Coefficient[eTot[c], dHb], eTot[c] /. {dHa -> 0, dHb -> 0}};
check["the 11 and 44 Einstein-Lovelock equations are linear in (Ha', Hb'); the 00 equation contains neither (constraint)",
  Exponent[eTot[2], dHa] <= 1 && Exponent[eTot[2], dHb] <= 1 && Exponent[eTot[5], dHa] <= 1 && Exponent[eTot[5], dHb] <= 1 && FreeQ[eTot[1], dHa | dHb]]
cl11 = coefLin[2]; cl44 = coefLin[5]; con00 = eTot[1];
dHsol[al2_?NumericQ, al3_?NumericQ, kap_?NumericQ, a_?NumericQ, b_?NumericQ, ha_?NumericQ, hb_?NumericQ] := dHsol[al2, al3, kap, a, b, ha, hb] =
   Module[{st = ksStress[a, b, kFsat], rule = {alf[1] -> 1, alf[2] -> al2, alf[3] -> al3, lamb -> 0, Ha -> ha, Hb -> hb}},
    LinearSolve[{cl11[[1 ;; 2]], cl44[[1 ;; 2]]} /. rule, {-kap st["pa"] - (cl11[[3]] /. rule), -kap st["pb"] - (cl44[[3]] /. rule)}]];
dHaF[al2_?NumericQ, al3_?NumericQ, kap_?NumericQ, a_?NumericQ, b_?NumericQ, ha_?NumericQ, hb_?NumericQ] := dHsol[al2, al3, kap, a, b, ha, hb][[1]];
dHbF[al2_?NumericQ, al3_?NumericQ, kap_?NumericQ, a_?NumericQ, b_?NumericQ, ha_?NumericQ, hb_?NumericQ] := dHsol[al2, al3, kap, a, b, ha, hb][[2]];
evolveEL[al2_, al3_, kap_, ha0_] := Module[{st = ksStress[1., 1., kFsat], cpol, hb0, sol, rule, viol},
   cpol = con00 /. {alf[1] -> 1, alf[2] -> al2, alf[3] -> al3, lamb -> 0, Ha -> ha0};
   hb0 = Hb /. First[SortBy[NSolve[cpol == kap st["eps"], Hb, Reals], Abs[Hb /. #] &]];
   sol = NDSolve[{aE'[x] == hE[x] aE[x], bE'[x] == kE[x] bE[x],
       hE'[x] == dHaF[al2, al3, kap, aE[x], bE[x], hE[x], kE[x]], kE'[x] == dHbF[al2, al3, kap, aE[x], bE[x], hE[x], kE[x]],
       aE[0] == 1, bE[0] == 1, hE[0] == ha0, kE[0] == N[hb0]}, {aE, bE, hE, kE}, {x, 0, 20},
      Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 8}, PrecisionGoal -> 10, AccuracyGoal -> 13, MaxSteps -> 10^6][[1]];
   viol = Max[Table[With[{aa = aE[x] /. sol, bb = bE[x] /. sol, ha = hE[x] /. sol, hb = kE[x] /. sol},
       Abs[((con00 /. {alf[1] -> 1, alf[2] -> al2, alf[3] -> al3, lamb -> 0, Ha -> ha, Hb -> hb}) - kap ksStress[aa, bb, kFsat]["eps"])/(kap ksStress[aa, bb, kFsat]["eps"])]], {x, 0, 20, 4}]];
   {N[hb0], {aE[20], bE[20], hE[20], kE[20]} /. sol, viol}];
runE = evolveEL[0., 0., 1000., 0.05];
runL = evolveEL[20., 500., 1000., 0.05];
check["self-consistent evolution: the constraint Sum alpha_k E^(k)0_0 - Lambda == kappa eps is preserved to 10^-7 (relative) for Einstein and for Einstein-Lovelock (alpha2 = 20, alpha3 = 500)",
  runE[[3]] < 10^-7 && runL[[3]] < 10^-7]
{"Einstein: Hb(0), {a, b, Ha, Hb}(20), violation" -> runE, "Einstein-Lovelock: Hb(0), {a, b, Ha, Hb}(20), violation" -> runL}
''')
