# Sections 11-14 of 7+1.nb (gravity on the 7+1 spacetime and dirac16complex in it); appended on 26 October 2026.
# This file is executed by tools/make_7p1_nb.py, which provides section(), subsection(), text(), item(), inp() and B.

section('11. Gravity on the 7+1 spacetime: metric, curvature and the Einstein tensor')
text([('About sections 11–14. ', B),
      'Sections 11–14 were appended on 26 October 2026. The notebook as it was before (sections 0–10, evaluated) is kept in '
      'the workspace as backups/7+1_2026-10-26_before_gravity.nb (excluded from the repository by .gitignore). They take the '
      'IDEAS of the repository github.com/once-ere/Dirac_claude — a 16-component field dirac16complex coupled to gravity through '
      'the canonical spin connection, Einstein–Lovelock gravity, a primordial pair-creation epoch, a Kohn–Sham treatment, the '
      'Standard Model at a lower-dimensional world surface with portals to dirac16complex — and rebuild them on THIS project’s '
      'spacetime. Dirac_claude works on a 4+4 spacetime (four time-like coordinates, a fold, a warped metric); none of that '
      'metric is used here. Here there is ONE time coordinate x0 and SEVEN space coordinates x1 … x7, the frame metric is '
      'η17 = diag(+1, −1, …, −1), and the Dirac matrices are the real Γ17[A] of section 4.2.'])
text([('Task of sections 11–14. ', B),
      '(i) Curvature of a general 7+1 metric: Christoffel symbols, Riemann, Ricci and Einstein tensors; (ii) the Lovelock '
      'tensors of 7+1 spacetime (Einstein, Gauss–Bonnet, cubic; the quartic one vanishes identically in eight dimensions) and the '
      'Einstein–Lovelock field equations; (iii) the generalized frame matrices Γ_μ(X) = e^A_μ(X) Γ17[A], the canonical spin '
      'connection, and the total covariant derivative, under which the frame matrices are constant (∇_μ Γ_ν = 0); (iv) the '
      'covariant Dirac equation of dirac16complex, with which the spin connection couples gravity to the field; (v) the three '
      'cosmological epochs — primordial (near x0 = 0), modern, far future — with the pair creation of dirac16complex in each.'])
text([('Conventions. ', B),
      'ħ = c = 1. Coordinates X = {x0, …, x7}; coordinate (Greek) indices μ, ν, λ, ρ, σ = 0 … 7; frame (Latin capital) indices '
      'A, B, C = 0 … 7, raised and lowered with η17. The metric g_μν has signature (+, −, −, −, −, −, −, −): x0 is time, '
      'x1, x2, x3 are ordinary space and x4 … x7 the four auxiliary space-like coordinates. The Christoffel symbols are written '
      '{λ μν} (never Γ, which is reserved for Dirac matrices); in the code chr[[λ+1, μ+1, ν+1]]. Curvature (Landau–Lifshitz): '
      'R^ρ_σμν = ∂_μ{ρ σν} − ∂_ν{ρ σμ} + {ρ μλ}{λ σν} − {ρ νλ}{λ σμ},  R_σν = R^ρ_σρν,  R = g^σν R_σν,  '
      'G^μ_ν = R^μ_ν − δ^μ_ν R/2. With these conventions the field equations are G^μ_ν − Λ δ^μ_ν = κ T^μ_ν with κ = 8πG8 '
      '(G8 the Newton constant of eight dimensions) and T^0_0 > 0 the energy density.'])

subsection('11.1 The curvature engine')
text('curvature[g, simp] returns, for a metric g given as an 8×8 matrix in the coordinates X, the inverse metric, the Christoffel '
     'symbols, the Riemann tensor riem[[ρ+1, σ+1, μ+1, ν+1]] = R^ρ_σμν, the tensor rud[[c+1, d+1, a+1, b+1]] = R^cd_ab = '
     'g^dσ R^c_σab (two indices up, two down, the form used by the Lovelock tensors), the Ricci tensor, the scalar curvature '
     'and the mixed Einstein tensor G^μ_ν. rudAssoc lists the non-zero components R^cd_ab with c < d and a < b.')
inp(r'''
curvature[g_, simp_: Simplify] := Module[{gi, chr, riem, rud, ric, rs},
  gi = simp[Inverse[g]];
  chr = Table[simp[Sum[gi[[l, r]] (D[g[[r, n]], X[[m]]] + D[g[[r, m]], X[[n]]] - D[g[[m, n]], X[[r]]]), {r, 8}]/2], {l, 8}, {m, 8}, {n, 8}];
  riem = Table[simp[D[chr[[r, s, n]], X[[m]]] - D[chr[[r, s, m]], X[[n]]] + Sum[chr[[r, m, l]] chr[[l, s, n]] - chr[[r, n, l]] chr[[l, s, m]], {l, 8}]],
    {r, 8}, {s, 8}, {m, 8}, {n, 8}];
  rud = Table[simp[Sum[gi[[d, s]] riem[[c, s, a, b]], {s, 8}]], {c, 8}, {d, 8}, {a, 8}, {b, 8}];
  ric = Table[simp[Sum[riem[[r, s, r, n]], {r, 8}]], {s, 8}, {n, 8}];
  rs = simp[Sum[gi[[s, n]] ric[[s, n]], {s, 8}, {n, 8}]];
  <|"g" -> g, "ginv" -> gi, "chr" -> chr, "riem" -> riem, "rud" -> rud, "ric" -> ric, "R" -> rs,
    "G" -> simp[gi.ric - rs IdentityMatrix[8]/2]|>];
rudAssoc[rud_] := Association[Select[Flatten[Table[{c, d, a, b} -> rud[[c + 1, d + 1, a + 1, b + 1]],
     {c, 0, 7}, {d, c + 1, 7}, {a, 0, 7}, {b, a + 1, 7}], 3], #[[2]] =!= 0 &]];
divergenceMixed[cv_, T_] := Table[Sum[D[T[[m, n]], X[[m]]] + Sum[cv["chr"][[m, m, l]] T[[l, n]] - cv["chr"][[l, m, n]] T[[m, l]], {l, 8}], {m, 8}], {n, 8}];
''')

subsection('11.2 The cosmological metric of 7+1 spacetime')
text('The spacetime of the three epochs (section 14) is homogeneous, isotropic within ordinary 3-space and isotropic within the '
     'four-dimensional auxiliary space, but the two spaces may expand at different rates:')
text('      ds² = dx0² − a(x0)² (dx1² + dx2² + dx3²) − b(x0)² (dx4² + dx5² + dx6² + dx7²),')
text('with the scale factor a(x0) (code a3[x0]) of ordinary space and b(x0) (code b4[x0]) of the auxiliary space. The Hubble rates are '
     'Ha = a′/a and Hb = b′/b (a prime is d/dx0; the code writes dHa, dHb for Ha′, Hb′). toH rewrites any expression in Ha, Hb, Ha′, Hb′. '
     'In an orthonormal frame the Riemann tensor of this metric is diagonal on pairs of indices: its only components are the '
     'sectional curvatures K_AB = R^AB_AB (no sum).')
inp(r'''
gCos = DiagonalMatrix[{1, -a3[x0]^2, -a3[x0]^2, -a3[x0]^2, -b4[x0]^2, -b4[x0]^2, -b4[x0]^2, -b4[x0]^2}];
cvCos = curvature[gCos];
toH[e_] := Expand[Simplify[e /. {a3''[x0] -> (dHa + Ha^2) a3[x0], b4''[x0] -> (dHb + Hb^2) b4[x0]} /. {a3'[x0] -> Ha a3[x0], b4'[x0] -> Hb b4[x0]}]];
secK[cv_, i_, j_] := cv["rud"][[i + 1, j + 1, i + 1, j + 1]];
check["cosmological metric: sectional curvatures K_0i = -(Ha'+Ha^2), K_0k = -(Hb'+Hb^2), K_ij = -Ha^2, K_kl = -Hb^2, K_ik = -Ha Hb",
  toH[secK[cvCos, 0, 1]] == -(dHa + Ha^2) && toH[secK[cvCos, 0, 4]] == -(dHb + Hb^2) && toH[secK[cvCos, 1, 2]] == -Ha^2 &&
  toH[secK[cvCos, 4, 5]] == -Hb^2 && toH[secK[cvCos, 1, 4]] == -Ha Hb]
check["cosmological metric: R^cd_ab is non-zero only for {c,d} = {a,b}",
  AllTrue[Keys[rudAssoc[cvCos["rud"]]], Sort[#[[1 ;; 2]]] == Sort[#[[3 ;; 4]]] &]]
check["Einstein tensor: G^0_0 = 3Ha^2 + 12 Ha Hb + 6 Hb^2, G^1_1 = 2Ha' + 3Ha^2 + 4Hb' + 10Hb^2 + 8HaHb, G^4_4 = 3Ha' + 6Ha^2 + 3Hb' + 6Hb^2 + 9HaHb",
  toH[cvCos["G"][[1, 1]]] == 3 Ha^2 + 12 Ha Hb + 6 Hb^2 && toH[cvCos["G"][[2, 2]]] == Expand[2 dHa + 3 Ha^2 + 4 dHb + 10 Hb^2 + 8 Ha Hb] &&
  toH[cvCos["G"][[5, 5]]] == Expand[3 dHa + 6 Ha^2 + 3 dHb + 6 Hb^2 + 9 Ha Hb]]
check["scalar curvature R = -2 (3Ha' + 6Ha^2 + 4Hb' + 10Hb^2 + 12 Ha Hb)", toH[cvCos["R"]] == Expand[-2 (3 dHa + 6 Ha^2 + 4 dHb + 10 Hb^2 + 12 Ha Hb)]]
{"G^0_0" -> toH[cvCos["G"][[1, 1]]], "G^1_1" -> toH[cvCos["G"][[2, 2]]], "G^4_4" -> toH[cvCos["G"][[5, 5]]], "R" -> toH[cvCos["R"]]}
''')

subsection('11.3 Symmetries of the Riemann tensor and the contracted Bianchi identity')
text('An inhomogeneous test metric with three arbitrary positive functions of x0 and x4 (nI, aI, bI) checks the engine away from '
     'homogeneity: ds² = nI² dx0² − aI² (dx1² + dx2² + dx3²) − bI² (dx4² + … + dx7²). With R_ρσμν = g_ρλ R^λ_σμν the Riemann tensor '
     'is antisymmetric in each pair, symmetric under exchange of the pairs, satisfies the first Bianchi identity, and the Einstein '
     'tensor is divergence-free: ∇_μ G^μ_ν = ∂_μ G^μ_ν + {μ μλ} G^λ_ν − {λ μν} G^μ_λ = 0 (contracted Bianchi identity).')
inp(r'''
gInh = DiagonalMatrix[{nI[x0, x4]^2, -aI[x0, x4]^2, -aI[x0, x4]^2, -aI[x0, x4]^2, -bI[x0, x4]^2, -bI[x0, x4]^2, -bI[x0, x4]^2, -bI[x0, x4]^2}];
cvInh = curvature[gInh, Together];
riemLow[cv_] := Table[Sum[cv["g"][[r, l]] cv["riem"][[l, s, m, n]], {l, 8}], {r, 8}, {s, 8}, {m, 8}, {n, 8}];
rlInh = riemLow[cvInh];
check["inhomogeneous metric: R_rsmn = -R_srmn = -R_rsnm = R_mnrs and R^r_[smn] = 0 (first Bianchi identity)",
  AllTrue[Flatten[Table[Together[{rlInh[[r, s, m, n]] + rlInh[[s, r, m, n]], rlInh[[r, s, m, n]] + rlInh[[r, s, n, m]], rlInh[[r, s, m, n]] - rlInh[[m, n, r, s]],
      cvInh["riem"][[r, s, m, n]] + cvInh["riem"][[r, m, n, s]] + cvInh["riem"][[r, n, s, m]]]], {r, 8}, {s, 8}, {m, 8}, {n, 8}]], # === 0 &]]
check["inhomogeneous metric: contracted Bianchi identity nabla_mu G^mu_nu = 0", AllTrue[Together[divergenceMixed[cvInh, cvInh["G"]]], # === 0 &]]
check["cosmological metric: contracted Bianchi identity nabla_mu G^mu_nu = 0", AllTrue[Simplify[divergenceMixed[cvCos, cvCos["G"]]], # === 0 &]]
''')

# =====================================================================================
section('12. The Lovelock tensors of 7+1 spacetime')
text('The generalized Kronecker delta δ^{μ1…μn}_{ν1…νn} is the determinant of the n×n matrix of ordinary Kronecker deltas '
     'δ^{μi}_{νj}: it is +1 (−1) if the lower indices are an even (odd) permutation of the upper ones, all different, and 0 '
     'otherwise. In eight dimensions it vanishes for n > 8. The Lovelock density of order k and the Lovelock tensor of order k are')
text('      L_k = (1/2^k) δ^{α1β1…αkβk}_{γ1δ1…γkδk} R^{γ1δ1}_{α1β1} ⋯ R^{γkδk}_{αkβk},')
text('      E^(k)μ_ν = −(1/2^(k+1)) δ^{μ α1β1…αkβk}_{ν γ1δ1…γkδk} R^{γ1δ1}_{α1β1} ⋯ R^{γkδk}_{αkβk}.')
text('L_0 = 1 and E^(0)μ_ν = −δ^μ_ν/2 (the cosmological term); L_1 = R and E^(1) = G (Einstein); L_2 = R² − 4R_μνR^μν + R_μνρσR^μνρσ '
     '(Gauss–Bonnet) and E^(2) is the Lanczos tensor; L_3 and E^(3) are cubic. E^(k) contains 2k+1 different upper indices, so '
     'E^(4) = 0 identically in eight dimensions; L_4, the Euler density, is not zero, but √|g| L_4 is a total derivative. The trace '
     'is E^(k)μ_μ = −((8 − 2k)/2) L_k. Each E^(k) is symmetric (with both indices down) and divergence-free for every metric, and '
     'δ(√|g| L_k) = √|g| E^(k)_μν δg^μν. The Einstein–Lovelock action and field equations of 7+1 spacetime are therefore')
text('      S = −(1/(2κ)) ∫ d⁸X √|g| Σ_{k=0}^{3} α_k L_k + S_matter,    Σ_{k=0}^{3} α_k E^(k)μ_ν = κ T^μ_ν,')
text('with α_0 = 2Λ (Λ the cosmological constant), α_1 = 1 (Einstein), α_2 and α_3 the Gauss–Bonnet and cubic couplings '
     '(dimensions length² and length⁴), δS_matter = (1/2) ∫ √|g| T_μν δg^μν. For α_2 = α_3 = 0 this is G^μ_ν − Λ δ^μ_ν = κ T^μ_ν.')

subsection('12.1 The generalized Kronecker delta')
inp(r'''
gKD[up_List, dn_List] := Det[Outer[KroneckerDelta, up, dn]];
SeedRandom[7108];
check["generalized Kronecker delta: contraction identity delta^{m1..mp r}_{n1..np r} = (8-p) delta^{m1..mp}_{n1..np} (p = 1..4, 60 random index sets each)",
  And @@ Flatten[Table[With[{u = RandomInteger[{0, 7}, pp], d = RandomInteger[{0, 7}, pp]},
      Sum[gKD[Append[u, r], Append[d, r]], {r, 0, 7}] == (8 - pp) gKD[u, d]], {pp, 1, 4}, {60}]]]
check["generalized Kronecker delta with 9 upper indices vanishes in 8 dimensions (100 random index sets)",
  And @@ Table[gKD[RandomInteger[{0, 7}, 9], RandomInteger[{0, 7}, 9]] == 0, {100}]]
''')

subsection('12.2 Computing the Lovelock tensors: the pair-set algorithm')
text('Because the generalized delta is antisymmetric in its upper and in its lower indices and every curvature factor is antisymmetric '
     'in each of its index pairs, the sum over all index values can be organized by sets of index pairs. Let U run over the sets of k '
     'disjoint pairs (α<β) not containing μ, and L over the ORDERED sequences of k disjoint pairs (γ<δ) of the remaining indices of '
     'U ∪ {μ} without ν. Then')
text('      E^(k)μ_ν = −2^(k−1) k! Σ_U Σ_L sgn(μ U → ν L) Π_i R^{L_i}_{U_i},      L_k = 2^k k! Σ_U Σ_L sgn(U → L) Π_i R^{L_i}_{U_i},')
text('where sgn is the sign of the permutation that takes the upper index sequence to the lower one (Signature[up] Signature[dn]). '
     'The check compares E^(1) with the brute-force δ-contraction −(1/4) δ^{μαβ}_{νγδ} R^γδ_αβ.')
inp(r'''
pairSets[idx_List, 0] := {{}};
pairSets[idx_List, k_Integer?Positive] := pairSets[idx, k] = If[Length[idx] < 2 k, {},
   Join[Flatten[Table[Prepend[#, {First[idx], idx[[j]]}] & /@ pairSets[Delete[Rest[idx], j - 1], k - 1], {j, 2, Length[idx]}], 1],
    pairSets[Rest[idx], k]]];
orderedSeqs[idx_List, k_Integer] := orderedSeqs[idx, k] = Flatten[Permutations /@ pairSets[Sort[idx], k], 1];
lovelockComp[ra_Association, k_Integer, mu_Integer, nu_Integer] := Module[{tot = 0, up, t},
  Do[up = Join[{mu}, Flatten[U]];
   If[MemberQ[up, nu],
    Do[t = Signature[up] Signature[Join[{nu}, Flatten[L]]];
     Do[t = t Lookup[ra, Key[Join[L[[i]], U[[i]]]], 0]; If[t === 0, Break[]], {i, k}];
     tot += t, {L, orderedSeqs[DeleteCases[up, nu], k]}]], {U, pairSets[DeleteCases[Range[0, 7], mu], k]}];
  -2^(k - 1) k! tot];
lovelockE[ra_Association, 0] := -IdentityMatrix[8]/2;
lovelockE[ra_Association, k_Integer?Positive] := Table[lovelockComp[ra, k, mu, nu], {mu, 0, 7}, {nu, 0, 7}];
lovelockL[ra_Association, 0] := 1;
lovelockL[ra_Association, k_Integer?Positive] := Module[{tot = 0, up, t},
  Do[up = Flatten[U];
   Do[t = Signature[up] Signature[Flatten[L]];
    Do[t = t Lookup[ra, Key[Join[L[[i]], U[[i]]]], 0]; If[t === 0, Break[]], {i, k}];
    tot += t, {L, orderedSeqs[up, k]}], {U, pairSets[Range[0, 7], k]}];
  2^k k! tot];
{"number of upper pair-sets (k = 1, 2, 3, 4; mu fixed)" -> Table[Length[pairSets[Range[1, 7], k]], {k, 1, 4}],
 "number of ordered lower sequences of 2k indices (k = 1, 2, 3, 4)" -> Table[Length[orderedSeqs[Range[1, 2 k], k]], {k, 1, 4}]}
''')
inp(r'''
rudCos = rudAssoc[cvCos["rud"]];
E1brute = Table[-1/4 Sum[gKD[{mu, a, b}, {nu, c, d}] cvCos["rud"][[c + 1, d + 1, a + 1, b + 1]], {a, 0, 7}, {b, 0, 7}, {c, 0, 7}, {d, 0, 7}], {mu, 0, 7}, {nu, 0, 7}];
check["pair-set algorithm == brute-force generalized-delta contraction for E^(1) (cosmological metric)", Simplify[E1brute - lovelockE[rudCos, 1]] == ConstantArray[0, {8, 8}]]
''')

subsection('12.3 The Lovelock tensors of the cosmological metric')
text('The checks: E^(1) = G and L_1 = R; the trace identity for k = 0 … 4; E^(4) = 0 while L_4 ≠ 0; every E^(k) is diagonal with '
     'equal entries inside ordinary space and inside the auxiliary space; E^(2) equals the Lanczos tensor '
     'H_μν = 2(R R_μν − 2R_μα R^α_ν − 2R^αβ R_μανβ + R_μ^αβγ R_ναβγ) − (1/2) g_μν L_2; and the closed forms below. '
     'Divergence-freeness becomes, for this metric, d/dx0 E^0_0 + 3Ha (E^0_0 − E^1_1) + 4Hb (E^0_0 − E^4_4) = 0.')
inp(r'''
ECos = Table[Simplify[lovelockE[rudCos, k]], {k, 0, 4}];      (* ECos[[k+1]] = E^(k) *)
LCos = Table[Simplify[lovelockL[rudCos, k]], {k, 0, 4}];
check["E^(1) == Einstein tensor G and L_1 == R (cosmological metric)", Simplify[ECos[[2]] - cvCos["G"]] == ConstantArray[0, {8, 8}] && Simplify[LCos[[2]] - cvCos["R"]] === 0]
check["trace identity E^(k)mu_mu == -((8-2k)/2) L_k for k = 0..4", And @@ Table[Simplify[Tr[ECos[[k + 1]]] + (8 - 2 k)/2 LCos[[k + 1]]] === 0, {k, 0, 4}]]
check["E^(4) == 0 identically in 8 dimensions, while the Euler density L_4 != 0", ECos[[5]] == ConstantArray[0, {8, 8}] && LCos[[5]] =!= 0]
check["every E^(k) is diagonal, E^(k)1_1 = E^(k)2_2 = E^(k)3_3 and E^(k)4_4 = ... = E^(k)7_7",
  And @@ Table[DiagonalMatrix[Diagonal[ECos[[k + 1]]]] == ECos[[k + 1]] && Simplify[ECos[[k + 1, 2, 2]] - ECos[[k + 1, 3, 3]]] === 0 &&
     Simplify[ECos[[k + 1, 2, 2]] - ECos[[k + 1, 4, 4]]] === 0 && And @@ Table[Simplify[ECos[[k + 1, 5, 5]] - ECos[[k + 1, j, j]]] === 0, {j, 6, 8}], {k, 0, 4}]]
lanczosDiag[cv_] := Module[{g = cv["g"], gi = cv["ginv"], ric = cv["ric"], rs = cv["R"], rl = riemLow[cv], ricU, rlU, l2},
  ricU = gi.ric.gi;
  rlU = Table[gi[[a, a]] gi[[b, b]] gi[[c, c]] rl[[m, a, b, c]], {m, 8}, {a, 8}, {b, 8}, {c, 8}];        (* diagonal metric *)
  l2 = rs^2 - 4 Total[ric ricU, 2] + Total[Table[gi[[m, m]] rlU[[m]] rl[[m]], {m, 8}], 4];
  {l2, gi.Table[2 (rs ric[[m, n]] - 2 Sum[ric[[m, a]] gi[[a, a]] ric[[a, n]], {a, 8}] - 2 Sum[ricU[[a, b]] rl[[m, a, n, b]], {a, 8}, {b, 8}]
        + Total[rlU[[m]] rl[[n]], 3]) - g[[m, n]] l2/2, {m, 8}, {n, 8}]}];
lzCos = lanczosDiag[cvCos];
check["E^(2) == Lanczos tensor and L_2 == R^2 - 4 Ric.Ric + Riem.Riem (cosmological metric)", Simplify[lzCos[[1]] - LCos[[3]]] === 0 && Simplify[lzCos[[2]] - ECos[[3]]] == ConstantArray[0, {8, 8}]]
E2closed = -12 Hb (4 Ha^3 + 18 Ha^2 Hb + 12 Ha Hb^2 + Hb^3); E3closed = 360 Ha^2 Hb^3 (4 Ha + 3 Hb);
check["E^(2)0_0 = -12 Hb (4Ha^3 + 18Ha^2 Hb + 12 Ha Hb^2 + Hb^3) and E^(3)0_0 = 360 Ha^2 Hb^3 (4Ha + 3Hb)",
  toH[ECos[[3, 1, 1]]] == Expand[E2closed] && toH[ECos[[4, 1, 1]]] == Expand[E3closed]]
consCos[k_] := Simplify[D[ECos[[k + 1, 1, 1]], x0] + 3 a3'[x0]/a3[x0] (ECos[[k + 1, 1, 1]] - ECos[[k + 1, 2, 2]]) + 4 b4'[x0]/b4[x0] (ECos[[k + 1, 1, 1]] - ECos[[k + 1, 5, 5]])];
check["divergence-free: d0 E^0_0 + 3Ha (E^0_0 - E^1_1) + 4Hb (E^0_0 - E^4_4) == 0 for k = 1, 2, 3 (arbitrary a, b)", And @@ Table[consCos[k] === 0, {k, 1, 3}]]
Grid[Prepend[Table[{k, toH[ECos[[k + 1, 1, 1]]], toH[ECos[[k + 1, 2, 2]]], toH[ECos[[k + 1, 5, 5]]]}, {k, 1, 3}], {"k", "E^(k)0_0", "E^(k)1_1", "E^(k)4_4"}], Frame -> All]
''')
text('The variational origin of E^(k): with a lapse function nL(x0) (g_00 = nL²) the Lagrangian of the cosmological metric is '
     'nL a³ b⁴ L_k, and its Euler–Lagrange expressions with respect to nL, a and b, at nL = 1, are −2 a³b⁴ E^(k)0_0, −6 a²b⁴ E^(k)1_1 '
     'and −8 a³b³ E^(k)4_4 — exactly the field-equation components (δg^00 = −2δnL, three δg^ii = 2δa/a³, four δg^kk = 2δb/b³). For '
     'k = 4 all three vanish identically: √|g| L_4 is a total derivative.')
inp(r'''
gLapse = DiagonalMatrix[{nL[x0]^2, -a3[x0]^2, -a3[x0]^2, -a3[x0]^2, -b4[x0]^2, -b4[x0]^2, -b4[x0]^2, -b4[x0]^2}];
rudLapse = rudAssoc[curvature[gLapse, Together]["rud"]];
eulerLagrange[lag_, q_] := D[lag, q[x0]] - D[D[lag, q'[x0]], x0] + D[D[lag, q''[x0]], {x0, 2}];
miniCheck[k_] := Module[{lag = Together[nL[x0] a3[x0]^3 b4[x0]^4 lovelockL[rudLapse, k]], ek = lovelockE[rudLapse, k]},
   Together[{eulerLagrange[lag, nL] + 2 a3[x0]^3 b4[x0]^4 ek[[1, 1]], eulerLagrange[lag, a3] + 6 a3[x0]^2 b4[x0]^4 ek[[2, 2]],
      eulerLagrange[lag, b4] + 8 a3[x0]^3 b4[x0]^3 ek[[5, 5]]} /. nL -> (1 &)]];
check["variation of nL a^3 b^4 L_k: EL_nL = -2 a^3 b^4 E^(k)0_0, EL_a = -6 a^2 b^4 E^(k)1_1, EL_b = -8 a^3 b^3 E^(k)4_4 (k = 1..4)",
  And @@ Table[miniCheck[k] === {0, 0, 0}, {k, 1, 4}]]
''')

subsection('12.4 An inhomogeneous metric: the Lovelock tensors are divergence-free (exact check at a point)')
text('For an inhomogeneous metric with the explicit polynomials nP, aP, bP of x0 and x4 below, the divergence of E^(k) is evaluated '
     'EXACTLY (rational numbers) at the point x0 = 1/2, x4 = 1/3. The derivative of E^(k), a polynomial of degree k in the curvature, is '
     'obtained from the curvature at the point and its first derivatives there: ∂_λ E^(k) = d/dε E^(k)[R + ε ∂_λR] at ε = 0.')
inp(r'''
nP = 1 + x0 x4/5 + x0^2/11; aP = 1 + x0/3 + x4^2/7; bP = 2 - x0 x4/4 + x0^2/9 + x4/13;
gPoly = DiagonalMatrix[{nP^2, -aP^2, -aP^2, -aP^2, -bP^2, -bP^2, -bP^2, -bP^2}];
cvPoly = curvature[gPoly, Together];
ptP = {x0 -> 1/2, x4 -> 1/3};
rudPoly = rudAssoc[cvPoly["rud"]];
rudAt = Map[(# /. ptP) &, rudPoly];
dRud[lam_] := Map[(D[#, X[[lam + 1]]] /. ptP) &, rudPoly];
chrAt = cvPoly["chr"] /. ptP;
divLovelock[k_] := Module[{e0 = lovelockE[rudAt, k], d0, d4},
   d0 = D[lovelockE[rudAt + \[Epsilon] dRud[0], k], \[Epsilon]] /. \[Epsilon] -> 0;
   d4 = D[lovelockE[rudAt + \[Epsilon] dRud[4], k], \[Epsilon]] /. \[Epsilon] -> 0;
   Table[d0[[1, n]] + d4[[5, n]] + Sum[chrAt[[m, m, l]] e0[[l, n]] - chrAt[[l, m, n]] e0[[m, l]], {m, 8}, {l, 8}], {n, 8}]];
check["inhomogeneous metric: nabla_mu E^(k)mu_nu == 0 EXACTLY at (x0, x4) = (1/2, 1/3) for k = 1, 2, 3", And @@ Table[divLovelock[k] === ConstantArray[0, 8], {k, 1, 3}]]
check["inhomogeneous metric: E^(1) == G at the point, and E^(k)_mu nu symmetric (k = 1, 2, 3)",
  lovelockE[rudAt, 1] == (cvPoly["G"] /. ptP) && And @@ Table[With[{el = (gPoly /. ptP).lovelockE[rudAt, k]}, el == Transpose[el]], {k, 1, 3}]]
lzPoly = lanczosDiag[<|"g" -> (gPoly /. ptP), "ginv" -> (cvPoly["ginv"] /. ptP), "ric" -> (cvPoly["ric"] /. ptP), "R" -> (cvPoly["R"] /. ptP), "riem" -> (cvPoly["riem"] /. ptP)|>];
check["inhomogeneous metric: E^(2) == Lanczos tensor at the point (exact)", lzPoly[[2]] == lovelockE[rudAt, 2] && lzPoly[[1]] == lovelockL[rudAt, 2]]
{"E^(3)0_4 at the point" -> lovelockE[rudAt, 3][[1, 5]]}
''')

subsection('12.5 A plane gravitational wave: only the Einstein tensor survives')
text('The pp-wave ds² = dx0² − dx1² − … − dx7² + fPP(x0 − x1, x2, x3) (dx0 − dx1)² with an arbitrary profile fPP is not diagonal. '
     'With the null vector l_μ = (1, −1, 0, …, 0): E^(1)_μν = +(1/2)(∂2² + ∂3²) fPP l_μ l_ν (so the wave is a vacuum solution when '
     'the profile is harmonic in the transverse directions), and E^(2) = E^(3) = 0 for EVERY profile: plane waves solve the '
     'Einstein–Lovelock equations with any α_2, α_3.')
inp(r'''
gPP = DiagonalMatrix[{1, -1, -1, -1, -1, -1, -1, -1}] + fPP[x0 - x1, x2, x3] Outer[Times, {1, -1, 0, 0, 0, 0, 0, 0}, {1, -1, 0, 0, 0, 0, 0, 0}];
cvPP = curvature[gPP];
rudPP = rudAssoc[cvPP["rud"]];
lNull = {1, -1, 0, 0, 0, 0, 0, 0};
check["pp-wave: E^(1)_mu nu == +(1/2)(d2^2 + d3^2) fPP l_mu l_nu",
  Simplify[gPP.lovelockE[rudPP, 1] - (D[fPP[x0 - x1, x2, x3], {x2, 2}] + D[fPP[x0 - x1, x2, x3], {x3, 2}])/2 Outer[Times, lNull, lNull]] == ConstantArray[0, {8, 8}]]
check["pp-wave: E^(2) == E^(3) == 0 for every profile fPP", Simplify[lovelockE[rudPP, 2]] == ConstantArray[0, {8, 8}] && Simplify[lovelockE[rudPP, 3]] == ConstantArray[0, {8, 8}]]
''')

# =====================================================================================
section('13. Generalized frame matrices, the canonical spin connection, and the covariant Dirac equation of dirac16complex')
text('A frame field (vielbein) e^A_μ(X) is a square root of the metric, g_μν = η17_AB e^A_μ e^B_ν; its inverse is e_A^μ (code ei[[μ+1, A+1]]). '
     'The GENERALIZED FRAME MATRICES are the position-dependent real 16×16 matrices')
text('      Γ_μ(X) = e^A_μ(X) Γ17[A],    Γ^μ(X) = e_A^μ(X) Γ17up[A],    Γ_μ Γ_ν + Γ_ν Γ_μ = 2 g_μν(X) ID16.')
text('The CANONICAL (Levi-Civita) SPIN CONNECTION is ω_μ^A_B = e^A_ν (∂_μ e_B^ν + {ν μλ} e_B^λ), antisymmetric when its first frame index '
     'is lowered, ω_μAB = −ω_μBA; on 16-component spinors it acts through')
text('      Ω_μ = (1/4) ω_μAB Γ17up[A] Γ17up[B],   D_μ ψ = ∂_μ ψ + Ω_μ ψ,   D_μ ψbar = ∂_μ ψbar − ψbar Ω_μ   (ψbar = ψ† Γ17[0]).')
text('The TOTAL covariant derivative acts on every index of an object: Christoffel symbols on coordinate indices, ω on frame indices, '
     'and the commutator with Ω on the two spinor indices of a 16×16 matrix. The canonical spin connection is the one for which '
     'the total covariant derivative of the frame field and of the frame matrices vanishes:')
text('      D_μ e^A_ν = ∂_μ e^A_ν − {λ μν} e^A_λ + ω_μ^A_B e^B_ν = 0   (vielbein postulate),')
text('      ∇_μ Γ_ν = ∂_μ Γ_ν − {λ μν} Γ_λ + [Ω_μ, Γ_ν] = 0,     ω_μ^A_B Γ17up[B] + [Ω_μ, Γ17up[A]] = 0.')
text('Its curvature F_μν = ∂_μΩ_ν − ∂_νΩ_μ + [Ω_μ, Ω_ν] = [D_μ, D_ν] is the Riemann tensor in spinor form, '
     'F_μν = (1/4) R_ABμν Γ17up[A] Γ17up[B] with R_ABμν = η17_AC e^C_ρ R^ρ_σμν e_B^σ; it satisfies the Bianchi identity '
     'D_[λ F_μν] = 0, and (1/2) Γ^μ Γ^ν F_μν = −(R/4) ID16, the term of the Lichnerowicz formula (Γ^μ D_μ)² = g^μν ∇_μ D_ν − R/4.')
text('The covariant Dirac equation of dirac16complex and its Lagrangian density are')
text('      (i Γ^μ(X) D_μ − m) ψ = 0,      √|g| ψbar (i Γ^μ D_μ − m) ψ,')
text('the spin connection Ω_μ being the coupling of gravity to the field. The current J^μ = ψbar Γ^μ ψ is conserved, '
     '∂_μ(√|g| J^μ) = 0, because ∂_μ(√|g| ψbar Γ^μ ψ) = √|g| [(D_μψbar) Γ^μ ψ + ψbar Γ^μ D_μ ψ] and the Dirac equation and its adjoint '
     'make the bracket vanish. All of this is checked below for three metrics: the cosmological one, the pp-wave (not diagonal) and '
     'the inhomogeneous one with arbitrary nI, aI, bI.')

subsection('13.1 The frame data of a metric')
inp(r'''
frameData[cv_Association, e_, simp_: Simplify] := Module[{g = cv["g"], chr = cv["chr"], ei, gam, gamU, om, omL, Om, F},
  ei = simp[Inverse[e]];
  gam = Table[Sum[e[[A + 1, m]] \[CapitalGamma]17[A], {A, 0, 7}], {m, 8}];
  gamU = Table[Sum[ei[[m, A + 1]] \[CapitalGamma]17up[A], {A, 0, 7}], {m, 8}];
  om = Table[simp[Sum[e[[A, n]] (D[ei[[n, B]], X[[m]]] + Sum[chr[[n, m, l]] ei[[l, B]], {l, 8}]), {n, 8}]], {m, 8}, {A, 8}, {B, 8}];
  omL = Table[\[Eta]17.om[[m]], {m, 8}];
  Om = Table[simp[(Zero16 + Sum[If[omL[[m, A + 1, B + 1]] === 0, 0, omL[[m, A + 1, B + 1]] \[CapitalGamma]17up[A].\[CapitalGamma]17up[B]], {A, 0, 7}, {B, 0, 7}])/4], {m, 8}];   (* Zero16 + ...: a 16x16 matrix even when every term vanishes *)
  F = Table[simp[D[Om[[n]], X[[m]]] - D[Om[[m]], X[[n]]] + Om[[m]].Om[[n]] - Om[[n]].Om[[m]]], {m, 8}, {n, 8}];
  <|"e" -> e, "ei" -> ei, "Gam" -> gam, "GamUp" -> gamU, "om" -> om, "omL" -> omL, "Om" -> Om, "F" -> F|>];
frameChecks[name_String, cv_Association, fd_Association, simp_: Simplify] := Module[{g = cv["g"], chr = cv["chr"], e = fd["e"], ei = fd["ei"],
   gam = fd["Gam"], gamU = fd["GamUp"], om = fd["om"], omL = fd["omL"], Om = fd["Om"], F = fd["F"], z16 = ConstantArray[0, {16, 16}], rABmn},
  check[name <> ": eta17 e e == g and {Gamma_mu, Gamma_nu} == 2 g_mu nu ID16 (generalized frame matrices, real)",
   simp[Transpose[e].\[Eta]17.e - g] == ConstantArray[0, {8, 8}] &&
    And @@ Flatten[Table[simp[gam[[m]].gam[[n]] + gam[[n]].gam[[m]] - 2 g[[m, n]] ID16] == z16, {m, 8}, {n, m, 8}]] && FreeQ[gam, Complex]];
  check[name <> ": canonical spin connection antisymmetric, omega_mu AB == -omega_mu BA", And @@ Table[simp[omL[[m]] + Transpose[omL[[m]]]] == ConstantArray[0, {8, 8}], {m, 8}]];
  check[name <> ": vielbein postulate d_mu e^A_nu - {l mu nu} e^A_l + omega_mu^A_B e^B_nu == 0",
   And @@ Flatten[Table[simp[D[e[[A, n]], X[[m]]] - Sum[chr[[l, m, n]] e[[A, l]], {l, 8}] + Sum[om[[m, A, B]] e[[B, n]], {B, 8}]] === 0, {m, 8}, {n, 8}, {A, 8}]]];
  check[name <> ": TOTAL covariant derivative of the frame matrices vanishes, d_mu Gamma_nu - {l mu nu} Gamma_l + [Omega_mu, Gamma_nu] == 0",
   And @@ Flatten[Table[simp[D[gam[[n]], X[[m]]] - Sum[chr[[l, m, n]] gam[[l]], {l, 8}] + Om[[m]].gam[[n]] - gam[[n]].Om[[m]]] == z16, {m, 8}, {n, 8}]]];
  check[name <> ": constant frame matrices are covariantly constant, omega_mu^A_B Gamma17up[B] + [Omega_mu, Gamma17up[A]] == 0",
   And @@ Flatten[Table[simp[Sum[om[[m, A + 1, B + 1]] \[CapitalGamma]17up[B], {B, 0, 7}] + Om[[m]].\[CapitalGamma]17up[A] - \[CapitalGamma]17up[A].Om[[m]]] == z16, {m, 8}, {A, 0, 7}]]];
  check[name <> ": Omega_mu real and Gamma17[0].Transpose[Omega_mu].Gamma17[0] == -Omega_mu (so D_mu psibar = d_mu psibar - psibar Omega_mu)",
   FreeQ[Om, Complex] && And @@ Table[simp[\[CapitalGamma]17[0].Transpose[Om[[m]]].\[CapitalGamma]17[0] + Om[[m]]] == z16, {m, 8}]];
  rABmn[m_, n_] := Table[simp[Sum[\[Eta]17[[A, A]] e[[A, r]] cv["riem"][[r, s, m, n]] ei[[s, B]], {r, 8}, {s, 8}]], {A, 8}, {B, 8}];
  check[name <> ": spin curvature F_mu nu == (1/4) R_AB mu nu Gamma17up[A].Gamma17up[B]",
   And @@ Flatten[Table[With[{rr = rABmn[m, n]}, simp[F[[m, n]] - Sum[rr[[A + 1, B + 1]] \[CapitalGamma]17up[A].\[CapitalGamma]17up[B], {A, 0, 7}, {B, 0, 7}]/4] == z16], {m, 8}, {n, m + 1, 8}]]];
  check[name <> ": Bianchi identity d_l F_mn + [Omega_l, F_mn] + cyclic == 0",
   And @@ Flatten[Table[simp[Sum[D[F[[t[[2]], t[[3]]]], X[[t[[1]]]]] + Om[[t[[1]]]].F[[t[[2]], t[[3]]]] - F[[t[[2]], t[[3]]]].Om[[t[[1]]]],
          {t, {{l, m, n}, {m, n, l}, {n, l, m}}}]] == z16, {l, 8}, {m, l + 1, 8}, {n, m + 1, 8}]]];
  check[name <> ": (1/2) Gamma^mu Gamma^nu F_mu nu == -(R/4) ID16", simp[Sum[gamU[[m]].gamU[[n]].F[[m, n]], {m, 8}, {n, 8}]/2 + cv["R"] ID16/4] == z16];
  ];
''')

subsection('13.2 The cosmological metric: frame, spin connection, Dirac equation')
text('Frame e^0 = dx0, e^i = a dx^i (i = 1, 2, 3), e^k = b dx^k (k = 4 … 7). The spin connection is Ω_0 = 0, '
     'Ω_i = (a′/2) Γ17[i].Γ17[0], Ω_k = (b′/2) Γ17[k].Γ17[0], NOT zero whenever the universe expands, and Γ^μ Ω_μ = (3Ha/2 + 2Hb) Γ17[0]. '
     'The rescaled field χ = (a³b⁴)^(1/2) ψ obeys an equation without connection; for a plane wave χ = e^{i(p⃗·x⃗ + q⃗·y⃗)} c(x0) '
     '(p⃗ = (p1, p2, p3), q⃗ = (q4, …, q7), x⃗ = (x1, x2, x3), y⃗ = (x4, …, x7)) the Dirac equation becomes the Schrödinger-type equation')
text('      i dc/dx0 = h(x0) c,     h(x0) = αD·(p⃗/a, q⃗/b) + βD m = Σ_i αD[i] p_i/a + Σ_k αD[k] q_k/b + βD m,')
text('the flat-space Hamiltonian matrix of section 7.3 with the physical (redshifted) momenta p⃗/a and q⃗/b.')
inp(r'''
eCos = DiagonalMatrix[{1, a3[x0], a3[x0], a3[x0], b4[x0], b4[x0], b4[x0], b4[x0]}];
fdCos = frameData[cvCos, eCos];
frameChecks["cosmological metric", cvCos, fdCos];
check["cosmological metric: Omega_0 == 0, Omega_i == (a'/2) Gamma17[i].Gamma17[0], Omega_k == (b'/2) Gamma17[k].Gamma17[0]",
  fdCos["Om"][[1]] == ConstantArray[0, {16, 16}] && And @@ Table[Simplify[fdCos["Om"][[i + 1]] - a3'[x0]/2 \[CapitalGamma]17[i].\[CapitalGamma]17[0]] == ConstantArray[0, {16, 16}], {i, 1, 3}] &&
  And @@ Table[Simplify[fdCos["Om"][[k + 1]] - b4'[x0]/2 \[CapitalGamma]17[k].\[CapitalGamma]17[0]] == ConstantArray[0, {16, 16}], {k, 4, 7}]]
check["cosmological metric: Gamma^mu Omega_mu == (3a'/(2a) + 2b'/b) Gamma17[0]",
  Simplify[Sum[fdCos["GamUp"][[m]].fdCos["Om"][[m]], {m, 8}] - (3 a3'[x0]/(2 a3[x0]) + 2 b4'[x0]/b4[x0]) \[CapitalGamma]17[0]] == ConstantArray[0, {16, 16}]]
''')
inp(r'''
psiF = Table[fpsi[j][x0, x1, x5], {j, 16}]; psibF = Table[hpsi[j][x0, x1, x5], {j, 16}];
DpsiCos = Table[D[psiF, X[[m]]] + fdCos["Om"][[m]].psiF, {m, 8}];
DpsibCos = Table[D[psibF, X[[m]]] - psibF.fdCos["Om"][[m]], {m, 8}];
sqgCos = a3[x0]^3 b4[x0]^4;
check["current identity: d_mu(sqrt|g| psibar Gamma^mu psi) == sqrt|g| [(D_mu psibar) Gamma^mu psi + psibar Gamma^mu D_mu psi] (psi, psibar arbitrary)",
  Together[Sum[D[sqgCos psibF.fdCos["GamUp"][[m]].psiF, X[[m]]], {m, 8}] - sqgCos Sum[DpsibCos[[m]].fdCos["GamUp"][[m]].psiF + psibF.fdCos["GamUp"][[m]].DpsiCos[[m]], {m, 8}]] === 0]
DslCos[f_] := Sum[fdCos["GamUp"][[m]].(D[f, X[[m]]] + fdCos["Om"][[m]].f), {m, 8}];
boxCos = Sum[cvCos["ginv"][[m, m]] (D[DpsiCos[[m]], X[[m]]] + fdCos["Om"][[m]].DpsiCos[[m]] - Sum[cvCos["chr"][[l, m, m]] DpsiCos[[l]], {l, 8}]), {m, 8}];
check["Lichnerowicz formula: (Gamma^mu D_mu)^2 psi == (g^mu nu nabla_mu D_nu - R/4) psi", Together[DslCos[DslCos[psiF]] - boxCos + cvCos["R"] psiF/4] === ConstantArray[0, 16]]
chiF = Table[cchi[j][x0, x1, x5], {j, 16}];
check["rescaling psi = (a^3 b^4)^(-1/2) chi removes the spin connection: Gamma^mu D_mu psi == (a^3 b^4)^(-1/2) Gamma^mu d_mu chi",
  Simplify[DslCos[chiF/Sqrt[sqgCos]] - Sum[fdCos["GamUp"][[m]].D[chiF, X[[m]]], {m, 8}]/Sqrt[sqgCos]] === ConstantArray[0, 16]]
cModes = Table[cm[j][x0], {j, 16}];
hCos = Sum[\[Alpha]D[i] pmom[i]/a3[x0], {i, 1, 3}] + Sum[\[Alpha]D[k] qmom[k]/b4[x0], {k, 4, 7}] + \[Beta]D mD;
phaseCos = Exp[I (Sum[pmom[i] X[[i + 1]], {i, 1, 3}] + Sum[qmom[k] X[[k + 1]], {k, 4, 7}])];
check["mode equation: Gamma17[0].(i Gamma^mu D_mu - m) psi == (a^3 b^4)^(-1/2) e^{i(p.x+q.y)} (i dc/dx0 - h(x0) c) for psi = (a^3 b^4)^(-1/2) e^{i(p.x+q.y)} c(x0)",
  Simplify[\[CapitalGamma]17[0].(I DslCos[phaseCos cModes/Sqrt[sqgCos]] - mD phaseCos cModes/Sqrt[sqgCos]) Sqrt[sqgCos]/phaseCos - (I D[cModes, x0] - hCos.cModes)] === ConstantArray[0, 16]]
''')

subsection('13.3 The pp-wave and the inhomogeneous metric: frames that are not diagonal and not homogeneous')
text('pp-wave frame: e^0 = (1 + fPP/2) dx0 − (fPP/2) dx1, e^1 = (fPP/2) dx0 + (1 − fPP/2) dx1, e^A = dx^A otherwise. '
     'Inhomogeneous frame: e^0 = nI dx0, e^i = aI dx^i, e^k = bI dx^k. Every identity of section 13.1 holds for both.')
inp(r'''
ePP = IdentityMatrix[8]; ePP[[1, 1]] = 1 + fPP[x0 - x1, x2, x3]/2; ePP[[1, 2]] = -fPP[x0 - x1, x2, x3]/2;
ePP[[2, 1]] = fPP[x0 - x1, x2, x3]/2; ePP[[2, 2]] = 1 - fPP[x0 - x1, x2, x3]/2;
frameChecks["pp-wave", cvPP, frameData[cvPP, ePP, Together], Together];
eInh = DiagonalMatrix[{nI[x0, x4], aI[x0, x4], aI[x0, x4], aI[x0, x4], bI[x0, x4], bI[x0, x4], bI[x0, x4], bI[x0, x4]}];
frameChecks["inhomogeneous metric", cvInh, frameData[cvInh, eInh, Together], Together];
''')

subsection('13.4 Conformal weight 7/2 of the massless field')
text('For a conformally flat metric g = C(x0)² η17 (x0 is then conformal time) the massless covariant Dirac operator satisfies '
     'Γ^μ D_μ (C^(−7/2) f) = C^(−9/2) Γ17up[A] ∂_A f for every f: a massless dirac16complex in a conformally flat 7+1 spacetime is '
     'a rescaled flat-space solution, and no pairs are created. (The cosmological metric is conformally flat only when a = b.)')
inp(r'''
gConf = cConf[x0]^2 \[Eta]17;
cvConf = curvature[gConf];
fdConf = frameData[cvConf, cConf[x0] IdentityMatrix[8]];
kF = Table[kflat[j][x0, x1, x5], {j, 16}];
check["conformal weight 7/2: Gamma^mu D_mu (C^(-7/2) f) == C^(-9/2) Gamma17up[A] d_A f for g = C(x0)^2 eta17",
  Simplify[Sum[fdConf["GamUp"][[m]].(D[cConf[x0]^(-7/2) kF, X[[m]]] + fdConf["Om"][[m]].(cConf[x0]^(-7/2) kF)), {m, 8}] -
     cConf[x0]^(-9/2) Sum[\[CapitalGamma]17up[A].D[kF, X[[A + 1]]], {A, 0, 7}]] === ConstantArray[0, 16]]
''')

# =====================================================================================
section('14. The three epochs: primordial (near x0 = 0), modern, far future')
text('For the cosmological metric with a homogeneous source T^μ_ν = diag(ε, −pa, −pa, −pa, −pb, −pb, −pb, −pb) (energy density ε, '
     'pressure pa along ordinary space and pb along the auxiliary space) the Einstein–Lovelock equations of section 12 are three '
     'ordinary differential equations,')
text('      Σ_k α_k E^(k)0_0 = κ ε,    Σ_k α_k E^(k)1_1 = −κ pa,    Σ_k α_k E^(k)4_4 = −κ pb,')
text('with the closed forms of section 12.3, and the conservation law ε′ + 3Ha (ε + pa) + 4Hb (ε + pb) = 0 follows from them. '
     'The three epochs of this notebook are exact solutions or limits of these equations.')

subsection('14.1 Primordial epoch near x0 = 0: Kasner solutions')
text('Near x0 = 0 a power law a = x0^p, b = x0^q is the natural ansatz. Every term of E^(k) then scales as x0^(−2k), so for x0 → 0 '
     'the highest Lovelock order with a non-zero coupling dominates (the cubic one if α_3 ≠ 0), and matter with ε growing more slowly '
     'than x0^(−2k) can be neglected. The vacuum solutions of each order separately are:')
item('Einstein (k = 1): 3p + 4q = 1 and 3p² + 4q² = 1, i.e. p = (1 ± 2√2)/7, q = (2 ∓ 3√2)/14. On the branch p = 0.547, q = −0.160 '
     'ordinary space expands while the auxiliary space contracts: dynamical compactification of x4 … x7.')
item('Gauss–Bonnet (k = 2): 3p + 4q = 3 and a cubic equation for p with three real roots.')
item('Cubic Lovelock (k = 3): the non-degenerate solution (p, q) = (−15/7, 20/7) with 3p + 4q = 5, and the degenerate families '
     'q = 0 (static auxiliary space) and p = 0, for which E^(3) vanishes identically.')
inp(r'''
(* a = x0^pK, b = x0^qK: Ha = pK/x0, Hb = qK/x0, Ha' = -pK/x0^2, Hb' = -qK/x0^2 in the closed forms of section 12.3 *)
EKas = Table[Factor[Expand[x0^(2 k) toH[ECos[[k + 1, c, c]]] /. {Ha -> pK/x0, Hb -> qK/x0, dHa -> -pK/x0^2, dHb -> -qK/x0^2}]], {k, 1, 3}, {c, {1, 2, 5}}];
gKas = DiagonalMatrix[{1, -x0^(2 pK), -x0^(2 pK), -x0^(2 pK), -x0^(2 qK), -x0^(2 qK), -x0^(2 qK), -x0^(2 qK)}];
check["Kasner metric: the Einstein tensor computed directly from the metric equals the closed form (x0^2 G^mu_mu, mu = 0, 1, 4)",
  Simplify[x0^2 Diagonal[curvature[gKas, PowerExpand[Together[#]] &]["G"]][[{1, 2, 5}]] - EKas[[1]]] === {0, 0, 0}]
check["Einstein Kasner: (p, q) = ((1+2Sqrt[2])/7, (2-3Sqrt[2])/14) and ((1-2Sqrt[2])/7, (2+3Sqrt[2])/14) solve E^(1) == 0",
  And @@ Flatten[Table[Simplify[EKas[[1]] /. {pK -> s[[1]], qK -> s[[2]]}] == {0, 0, 0}, {s, {{(1 + 2 Sqrt[2])/7, (2 - 3 Sqrt[2])/14}, {(1 - 2 Sqrt[2])/7, (2 + 3 Sqrt[2])/14}}}]]]
gbCubic = Expand[64 (4 pK^3 + 18 pK^2 qK + 12 pK qK^2 + qK^3) /. qK -> (3 - 3 pK)/4];
gbRoots = Sort[pK /. NSolve[gbCubic == 0, pK, Reals, WorkingPrecision -> 30]];
check["Gauss-Bonnet Kasner: 3 real roots of the cubic; (p, (3-3p)/4) solve E^(2) == 0",
  Length[gbRoots] == 3 && And @@ Table[Max[Abs[N[EKas[[2]] /. {pK -> r, qK -> (3 - 3 r)/4}, 25]]] < 10^-18, {r, gbRoots}]]
check["cubic Lovelock Kasner: (p, q) = (-15/7, 20/7) solves E^(3) == 0; E^(3) == 0 identically for q = 0 and for p = 0",
  Simplify[EKas[[3]] /. {pK -> -15/7, qK -> 20/7}] == {0, 0, 0} && Simplify[EKas[[3]] /. qK -> 0] == {0, 0, 0} && Simplify[EKas[[3]] /. pK -> 0] == {0, 0, 0}]
{"x0^2 E^(1) (00, 11, 44)" -> EKas[[1]], "Gauss-Bonnet cubic" -> gbCubic, "Gauss-Bonnet Kasner p" -> N[gbRoots, 15], "Einstein Kasner p, q" -> N[{(1 + 2 Sqrt[2])/7, (2 - 3 Sqrt[2])/14}, 15]}
''')

subsection('14.2 The mode equation in two dimensions')
text('With unit vectors n⃗a (in ordinary space) and n⃗b (in the auxiliary space) along p⃗ and q⃗, h(x0) = |p⃗|/a A + |q⃗|/b B + m βD with the '
     'CONSTANT real matrices A = αD·n⃗a, B = αD·n⃗b. A, B and βD square to ID16 and anticommute, and Zc = A.B.βD squares to −ID16, '
     'is antisymmetric and commutes with all three and with h(x0), so its eigenvalue ±i is conserved. On the eigenspace Zc = +i '
     '(eight dimensions) A, B, βD act like the Pauli matrices σ1, σ2, σ3 (σ1σ2σ3 = +i), on Zc = −i like σ1, σ2, −σ3, so the '
     '16-component mode equation is four copies of each of')
text('      i dc/dx0 = (|p⃗|/a σ1 + |q⃗|/b σ2 ± m σ3) c,     c ∈ ℂ².')
text('For m = 0 the two signs give the same evolution; for m ≠ 0 they do not, and the anisotropic expansion creates the two classes '
     'of states at different rates.')
text('The direction of the physical momentum (|p⃗|/a, |q⃗|/b) turns as the two spaces expand at different rates. Pairs are created when '
     'the turning is not slow compared with the frequency ω = √(|p⃗|²/a² + |q⃗|²/b² + m²): the Bogoliubov coefficient |β|² of a mode '
     'is the probability that a particle (positive-energy state) at the initial time ends in a negative-energy state at the final time.')
inp(r'''
naVec = {1, 2, 2}/3; nbVec = {1, -1, 1, 1}/2;
Aop = Sum[naVec[[i]] \[Alpha]D[i], {i, 1, 3}]; Bop2 = Sum[nbVec[[k - 3]] \[Alpha]D[k], {k, 4, 7}]; Zc = Aop.Bop2.\[Beta]D;
check["2x2 reduction: A, B, betaD square to ID16 and anticommute; Zc = A.B.betaD: Zc.Zc == -ID16, Transpose[Zc] == -Zc, Zc commutes with A, B, betaD",
  Aop.Aop == ID16 && Bop2.Bop2 == ID16 && Aop.Bop2 + Bop2.Aop == Zero16 && Aop.\[Beta]D + \[Beta]D.Aop == Zero16 && Bop2.\[Beta]D + \[Beta]D.Bop2 == Zero16 &&
  Zc.Zc == -ID16 && Transpose[Zc] == -Zc && Zc.Aop == Aop.Zc && Zc.Bop2 == Bop2.Zc && Zc.\[Beta]D == \[Beta]D.Zc]
pauli = {{{0, 1}, {1, 0}}, {{0, -I}, {I, 0}}, {{1, 0}, {0, -1}}};
bogoliubov[hf_, t0_, t1_, n_: 2] := Module[{e0, v0, e1, v1, half = n/2, sol, cc, tt, b2 = 0, a2 = 0, ord0, ord1},
  {e0, v0} = Eigensystem[N[hf[t0]]]; ord0 = Ordering[Re[e0]]; v0 = Normalize /@ v0[[ord0]];
  {e1, v1} = Eigensystem[N[hf[t1]]]; ord1 = Ordering[Re[e1]]; v1 = Normalize /@ v1[[ord1]];
  Do[sol = NDSolveValue[{cc'[tt] == -I hf[tt].cc[tt], cc[t0] == v0[[j]]}, cc[t1], {tt, t0, t1}, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 8},
     PrecisionGoal -> 12, AccuracyGoal -> 14, MaxSteps -> Infinity];
   b2 += Total[Abs[Conjugate[v1[[1 ;; half]]].sol]^2]; a2 += Total[Abs[Conjugate[v1[[half + 1 ;; n]]].sol]^2], {j, half + 1, n}];
  {b2, a2}/half];
''')

subsection('14.3 Pair creation in the primordial epoch (Einstein–Kasner branch p = 0.547, q = −0.160)')
text('For a massless mode, x0 → λ x0 maps the problem with momenta (|p⃗|, |q⃗|) onto one with the single scaling variable '
     'K = |p⃗| λ^(1−p), λ = (|q⃗|/|p⃗|)^(1/(q−p)), so |β|² depends on K only. Creation happens around the time at which |p⃗|/a ≈ |q⃗|/b; '
     'the evolution is adiabatic before and after, so the integration runs from 10⁻¹⁵ λ to 300 λ (or wider). An isotropic expansion '
     '(a = b) creates no massless pairs (section 13.4); the anisotropy of the primordial epoch does.')
inp(r'''
pKE = N[(1 + 2 Sqrt[2])/7, 20]; qKE = N[(2 - 3 Sqrt[2])/14, 20];
hKas[ka_, kb_][t_] := ka t^-pKE pauli[[1]] + kb t^-qKE pauli[[2]];
kasData = Table[With[{ka = kk[[1]], kb = kk[[2]], lam = (kk[[2]]/kk[[1]])^(1/(qKE - pKE))},
    Join[{ka, kb, ka lam^(1 - pKE)}, bogoliubov[hKas[ka, kb], 10^-15 Min[1, lam], 300 Max[1, lam]]]], {kk, {{1., 1.}, {2., 1.}, {1., 3.}}}];
check["Kasner pair creation: |alpha|^2 + |beta|^2 == 1 (to 10^-8) for (|p|, |q|) = (1,1), (2,1), (1,3)", And @@ Table[Abs[r[[4]] + r[[5]] - 1] < 10^-8, {r, kasData}]]
exK = (1 - pKE)/(qKE - pKE); kaSc = (1/2.5^exK)^(1/(1 - exK)); lamSc = (2.5/kaSc)^(1/(qKE - pKE));
check["Kasner scaling: (|p|, |q|) = (1.4302.., 2.5) has K = 1 and the same |beta|^2 as (1, 1) (to 10^-6)",
  Abs[bogoliubov[hKas[kaSc, 2.5], 10^-15 Min[1, lamSc], 300 Max[1, lamSc]][[1]] - kasData[[1, 4]]] < 10^-6]
mK16 = 0.3; pv16 = {0.6, 0., 0.8}; qv16 = 1.5 {0.5, 0.5, -0.5, 0.5};
h16Kas[t_] := Sum[pv16[[i]] t^-pKE \[Alpha]D[i], {i, 1, 3}] + Sum[qv16[[k - 3]] t^-qKE \[Alpha]D[k], {k, 4, 7}] + mK16 \[Beta]D;
h2Kas[sg_][t_] := Norm[pv16] t^-pKE pauli[[1]] + Norm[qv16] t^-qKE pauli[[2]] + sg mK16 pauli[[3]];
b16Kas = bogoliubov[h16Kas, 0.05, 40., 16]; b2Kas = bogoliubov[h2Kas[1], 0.05, 40.]; b2KasM = bogoliubov[h2Kas[-1], 0.05, 40.];
check["16x16 evolution == 4 copies of the 2x2 evolution with +m sigma3 and 4 with -m sigma3 (Kasner, m = 0.3, |p| = 1, |q| = 1.5, x0 from 0.05 to 40): average |beta|^2 agrees to 10^-8; unitarity; the two classes differ",
  Abs[b16Kas[[1]] - (b2Kas[[1]] + b2KasM[[1]])/2] < 10^-8 && Abs[Total[b16Kas] - 1] < 10^-8 && Abs[b2Kas[[1]] - b2KasM[[1]]] > 10^-4]
check["isotropic expansion, massless: no pair creation (|beta|^2 < 10^-20)", bogoliubov[Function[t, t^-0.4 (0.6 pauli[[1]] + 0.8 pauli[[2]])], 10^-6, 100.][[1]] < 10^-20]
{"Kasner beta2 (1,1)" -> kasData[[1, 4]], "Kasner beta2 (2,1)" -> kasData[[2, 4]], "Kasner beta2 (1,3)" -> kasData[[3, 4]],
 "Kasner K (1,1)" -> kasData[[1, 3]], "Kasner K (2,1)" -> kasData[[2, 3]], "Kasner K (1,3)" -> kasData[[3, 3]], "beta2 16x16" -> b16Kas[[1]], "beta2 2x2 +m" -> b2Kas[[1]], "beta2 2x2 -m" -> b2KasM[[1]]}
''')

subsection('14.4 Modern epoch: a static auxiliary space')
text('When the auxiliary space has stopped evolving (b = b0, Hb = 0) the Lovelock tensors reduce to E^(1) = diag(3Ha², 2Ha′ + 3Ha², '
     '3Ha′ + 6Ha²), E^(2) = diag(0, 0, −12 Ha²(Ha′ + Ha²)), E^(3) = 0 (entries 00, 11, 44). The 00 and 11 equations are the Friedmann '
     'equations of ordinary cosmology, 3Ha² = κ ε + Λ, 2Ha′ + 3Ha² = −κ pa + Λ, with the four-dimensional Newton constant '
     'G4 = G8/V⊥, where V⊥ = b0⁴ Vc is the proper volume of the auxiliary space (Vc its coordinate volume; ε8 = ε4/V⊥ for matter '
     'confined to, or smeared over, ordinary space). The 44 equation then FIXES the stress that must hold b constant: '
     'κ pb = Λ − 3Ha′ − 6Ha² + 12 α_2 Ha²(Ha′ + Ha²). A field of momentum q⃗ along the auxiliary space behaves in ordinary space as a '
     'particle of mass √(m² + |q⃗|²/b0²) (a Kaluza–Klein tower when the auxiliary space is a torus of coordinate period Lc, q⃗ = 2π n⃗/Lc). '
     'The present expansion rate H0 ≈ 1.4 × 10⁻⁴² GeV is so small that pair creation today is suppressed by exp(−2πm/H0) (section 14.5).')
inp(r'''
modH = Table[Simplify[toH[ECos[[k + 1, c, c]]] /. {Hb -> 0, dHb -> 0}], {k, 1, 3}, {c, {1, 2, 5}}];
check["static auxiliary space (Hb = 0): E^(1) = (3Ha^2, 2Ha' + 3Ha^2, 3Ha' + 6Ha^2), E^(2) = (0, 0, -12 Ha^2 (Ha' + Ha^2)), E^(3) = 0",
  Expand[modH - {{3 Ha^2, 2 dHa + 3 Ha^2, 3 dHa + 6 Ha^2}, {0, 0, -12 Ha^2 (dHa + Ha^2)}, {0, 0, 0}}] == ConstantArray[0, {3, 3}]]
''')

subsection('14.5 Far future: de Sitter expansion and thermal pair creation')
text('Two de Sitter solutions: (a) eight-dimensional de Sitter a = b = e^{H x0}, an exact VACUUM solution of the Einstein–Lovelock '
     'equations when 21 α_1 H² − 420 α_2 H⁴ + 2520 α_3 H⁶ = Λ (all E^(k) are then proportional to δ^μ_ν); (b) de Sitter expansion of '
     'ordinary space over a static auxiliary space, a = e^{H x0}, b = b0, with E^(1) = diag(3, 3, 6) H², E^(2) = diag(0, 0, −12) H⁴, '
     'E^(3) = 0. In both, a mode of dirac16complex obeys i dc/dx0 = (k e^{−H x0} σ1 + M σ3) c with M = m in (a) and the Kaluza–Klein mass '
     'M = √(m² + |q⃗|²/b0²) in (b). Starting from the adiabatic vacuum at early times, the number of created pairs per mode is EXACTLY '
     'the Fermi–Dirac occupation at the Gibbons–Hawking temperature H/(2π):')
text('      |β|² = 1/(exp(2π M/H) + 1),')
text('independent of the momentum k. The check integrates the mode equation numerically (H = 1).')
inp(r'''
dS8H = Table[Simplify[toH[ECos[[k + 1, c, c]]] /. {Ha -> hdS, Hb -> hdS, dHa -> 0, dHb -> 0}], {k, 1, 3}, {c, {1, 2, 5}}];
check["8D de Sitter: E^(1) = 21 H^2, E^(2) = -420 H^4, E^(3) = 2520 H^6 (all diagonal entries)", dS8H == {{21, 21, 21} hdS^2, {-420, -420, -420} hdS^4, {2520, 2520, 2520} hdS^6}]
dS4H = Table[Simplify[toH[ECos[[k + 1, c, c]]] /. {Ha -> hdS, Hb -> 0, dHa -> 0, dHb -> 0}], {k, 1, 3}, {c, {1, 2, 5}}];
check["de Sitter over a static auxiliary space: E^(1) = (3, 3, 6) H^2, E^(2) = (0, 0, -12) H^4, E^(3) = 0", dS4H == {{3, 3, 6} hdS^2, {0, 0, -12} hdS^4, {0, 0, 0}}]
dSbeta = Table[bogoliubov[Function[t, Exp[-t] pauli[[1]] + mH pauli[[3]]], -Log[2000.], Log[10.^6/mH^2]][[1]], {mH, {0.25, 0.5, 1.}}];
check["de Sitter pair creation: |beta|^2 == 1/(Exp[2 Pi M/H] + 1) to 10^-7, M/H = 0.25, 0.5, 1",
  And @@ MapThread[Abs[#1 - 1/(Exp[2 Pi #2] + 1)] < 10^-7 &, {dSbeta, {0.25, 0.5, 1.}}]]
lzBeta = Table[bogoliubov[Function[t, Sqrt[k2] pauli[[1]] + t pauli[[3]]], -150., 150.][[1]], {k2, {0.1, 0.5, 1.}}];
check["integrator benchmark (Landau-Zener): |beta|^2 == Exp[-Pi kappa^2] to 10^-6", And @@ MapThread[Abs[#1 - Exp[-Pi #2]] < 10^-6 &, {lzBeta, {0.1, 0.5, 1.}}]]
{"dS beta2 M/H=0.25" -> dSbeta[[1]], "dS beta2 M/H=0.5" -> dSbeta[[2]], "dS beta2 M/H=1" -> dSbeta[[3]], "FD 0.25" -> N[1/(Exp[2 Pi 0.25] + 1)],
 "Landau-Zener" -> lzBeta}
''')
