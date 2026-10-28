# Section 25 of 7+1.nb: the auxiliary torus and the world brane -- their precise definitions, and worked examples that check them.
# Appended on 28 October 2026.  Executed by tools/make_7p1_nb.py, which provides section(), subsection(), text(), item(), inp() and B.
# The names of section 25 begin with tb.

section('25. The auxiliary torus and the world brane: definitions and worked examples')
text([('About section 25. ', B),
      'Section 25 was appended on 28 October 2026. The notebook as it was before (sections 0–24, with the summary of all checks numbered 25; '
      'evaluated) is kept in the workspace as backups/7+1_2026-10-28_before_torus_brane.nb (excluded from the repository by .gitignore). '
      'Sections 14 and 17–24 use two geometric objects, the AUXILIARY TORUS and the WORLD BRANE. This section defines both precisely and checks '
      'every part of the definitions with worked examples; nothing in sections 0–24 is changed. Coordinates as everywhere: X = {x0, …, x7}, '
      'x0 is time (the direction of the octonion unit e[0] = 1), x1, x2, x3 ordinary space, y = (x4, x5, x6, x7) the four auxiliary '
      'space-like coordinates. η17 = diag(+1, −1, …, −1) is simply the squares of the octonion basis (e0² = +1, e_i² = −1); xAct declares the '
      'spacetime with the signature list {7, 1, 0}: seven space-like dimensions and one time-like dimension.'])
text([('The auxiliary torus. ', B),
      'The four auxiliary coordinates are periodic: x_k ≃ x_k + Lc n_k (k = 4 … 7, n_k ∈ ℤ, Lc > 0 the coordinate period). The auxiliary '
      'space is therefore the four-dimensional torus T⁴ = ℝ⁴/(Lc ℤ)⁴ — a compact, flat space without boundary, whose fundamental domain is the '
      'hypercube [0, Lc)⁴ (equivalently [−Lc/2, Lc/2)⁴) with opposite faces identified — and the spacetime of the model, the BULK, is '
      'ℝ (time x0) × ℝ³ (ordinary space x1, x2, x3) × T⁴. The SIZE b of the torus is defined by b⁸ = det(−g_kl) (k, l = 4 … 7), so that '
      'the proper volume element of the torus is b⁴ d⁴y. For every metric of sections 11–24, −g_kl = b² δ_kl: each of the four circles has the '
      'proper circumference b·Lc, and the proper volume is V⊥ = b⁴ Lc⁴; today b = b0 and ℓ⊥ = V⊥^(1/4) = b0·Lc. Only the product b·Lc is '
      'physical (section 25.1).'])
text([('The world brane. ', B),
      'W = {x4 = x5 = x6 = x7 = 0} is the 3+1-dimensional submanifold of the bulk through the point y = 0 of the torus. It is spanned by the '
      'time x0 and ordinary space x1, x2, x3; its embedding is X(x) = (x0, x1, x2, x3, 0, 0, 0, 0), and its metric is the induced metric '
      'h_μν = g_AB ∂_μX^A ∂_νX^B = g_μν(x, y = 0) (μ, ν = 0 … 3). The Standard Model lives on W; gravity and dirac16complex live in the bulk '
      'and meet the Standard Model at W. In the counting of xAct the bulk has the signature list {7, 1, 0} and W the list {3, 1, 0}; the flat '
      'metric diag(+1, −1, −1, −1) of W is the squares of the basis e[0] = 1, e[1], e[2], e[3] of the quaternion subalgebra (section 17). In '
      'this model W is a fixed submanifold (its position is not a dynamical field).'])

# ------------------------------------------------------------------------------------------------------------------------
subsection('25.1 The auxiliary torus: size, volume and symmetries')
text('For the cosmological metric of section 11.2 the size is b(x0), and the proper volume of the fundamental domain is b⁴ Lc⁴. The change '
     'of coordinates y = λ y′ describes the SAME torus with the period Lc/λ and the size λ b, so only b·Lc (and V⊥) can be measured. The '
     'isometries of the torus are the translations y → y + c and the 384 = 2⁴·4! signed permutations of the four axes (the symmetries of the '
     'lattice Lc ℤ⁴; an orthogonal matrix with integer entries is a signed permutation). The signed permutations fix the point y = 0, and with '
     'it the world brane. They permute the vectors n ∈ ℤ⁴ of each shell |n|² = k; vectors with the same multiset of |n_i| form one orbit: the '
     'shell k = 4 consists of the 8 vectors ±2e_i and the 16 vectors (±1, ±1, ±1, ±1), the shell k = 9 of the 8 vectors ±3e_i and the 96 '
     'arrangements of (±2, ±2, ±1, 0).')
inp(r'''
tbG4 = Table[gCos[[k, l]], {k, 5, 8}, {l, 5, 8}]; tbDy = {tbd4, tbd5, tbd6, tbd7};
check["auxiliary torus: for the cosmological metric det(-g_kl) == b^8 (k, l = 4..7), and the proper volume of the fundamental domain [0, Lc)^4, Integrate[Sqrt[det(-g_kl)], {y, 0, Lc}^4], is b^4 Lc^4 = V_perp",
  Simplify[Det[-tbG4] - b4[x0]^8] === 0 &&
  Simplify[Integrate[Sqrt[Det[-tbG4]], {x4, 0, tbLc}, {x5, 0, tbLc}, {x6, 0, tbLc}, {x7, 0, tbLc}] - b4[x0]^4 tbLc^4, b4[x0] > 0 && tbLc > 0] === 0]
check["only the product b Lc of the auxiliary torus is physical: y = lambda y' maps the torus of period Lc and size b to the one of period Lc/lambda and size lambda b; the line element b^2 dy.dy, the proper period b Lc and V_perp = b^4 Lc^4 are unchanged",
  Simplify[tbB^2 (tbLam tbDy).(tbLam tbDy) - (tbLam tbB)^2 tbDy.tbDy] === 0 && Simplify[(tbLc/tbLam) (tbLam tbB) - tbLc tbB] === 0 &&
  Simplify[(tbLc/tbLam)^4 (tbLam tbB)^4 - tbLc^4 tbB^4] === 0]
tbSP = Flatten[Table[DiagonalMatrix[s].IdentityMatrix[4][[p]], {p, Permutations[Range[4]]}, {s, Tuples[{1, -1}, 4]}], 1];
tbShell[k_] := Select[Tuples[Range[-3, 3], 4], #.# == k &];
tbOrbits[k_] := GatherBy[tbShell[k], Sort[Abs[#]] &];
tbOrbitQ[k_] := AllTrue[tbOrbits[k], Function[cl, Sort[Union[(#.First[cl]) & /@ tbSP]] == Sort[cl]]];
check["symmetries of the auxiliary torus that fix y = 0 (and the brane): the 384 = 2^4 4! signed permutation matrices are distinct orthogonal integer matrices (they map Z^4 onto itself and preserve |n|^2); the shells |n|^2 = 1, 2, 3, 4, 9 split into orbits of sizes {8}, {24}, {32}, {8, 16}, {8, 96}",
  Length[tbSP] == 384 && Length[Union[tbSP]] == 384 && AllTrue[tbSP, #.Transpose[#] == IdentityMatrix[4] &] &&
  AllTrue[{1, 2, 3, 4, 9}, tbOrbitQ] && (Sort[Length /@ tbOrbits[#]] & /@ {1, 2, 3, 4, 9}) == {{8}, {24}, {32}, {8, 16}, {8, 96}}]
{"orbit sizes of the shells k = 1..12" -> Table[Sort[Length /@ tbOrbits[k]], {k, 1, 12}]}
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('25.2 Fields on the torus: Fourier modes, the Laplacian and Kaluza–Klein masses')
text('A field on the torus is periodic in each y_k and is therefore a Fourier series in the modes φ_n(y) = e^{i q_n·y}/√V⊥, '
     'q_n = 2π n/Lc, n ∈ ℤ⁴ (q_n is the wave vector in the coordinate y; the physical momentum along the torus is q_n/b). The modes are '
     'orthonormal for the proper volume element, ∫ b⁴ d⁴y φ_n* φ_m = δ_nm, and complete (Parseval: ∫ b⁴ d⁴y |f|² = Σ_n |c_n|², '
     'c_n = ∫ b⁴ d⁴y φ_n* f). They are eigenfunctions of the Laplacian of the torus, −b⁻² Σ_k ∂²φ_n/∂y_k² = (|q_n|²/b²) φ_n. A bulk field '
     'Φ = χ(x) φ_n(y) of mass m is, for a static torus (b = b0), a particle of mass M_n = (m² + μ_n²)^(1/2) in ordinary spacetime, with the '
     'Kaluza–Klein mass μ_n = |q_n|/b0 = 2π|n|/ℓ⊥ (for dirac16complex the same tower, four brane Dirac fermions per n: section 18.1).')
inp(r'''
tbY = {x4, x5, x6, x7}; tbV = tbB^4 tbLc^4;
tbQ[n_] := 2 Pi n/tbLc;
tbPhi[n_] := Exp[I tbQ[n].tbY]/Sqrt[tbV]; tbPhiC[n_] := Exp[-I tbQ[n].tbY]/Sqrt[tbV];
tbInt[f_] := Simplify[Integrate[tbB^4 f, {x4, 0, tbLc}, {x5, 0, tbLc}, {x6, 0, tbLc}, {x7, 0, tbLc}], tbLc > 0 && tbB > 0];
check["Fourier modes of the auxiliary torus: phi_n = Exp[I q_n.y]/Sqrt[V_perp] with q_n = 2 Pi n/Lc (n in Z^4) are orthonormal for the proper volume element b^4 d^4y (exact integrals: n = m = (1, 0, -2, 3); n = (1, 0, -2, 3), m = (1, 1, -2, 3); n = 0, m = (0, 2, 0, 0))",
  {tbInt[tbPhiC[{1, 0, -2, 3}] tbPhi[{1, 0, -2, 3}]], tbInt[tbPhiC[{1, 0, -2, 3}] tbPhi[{1, 1, -2, 3}]], tbInt[tbPhiC[{0, 0, 0, 0}] tbPhi[{0, 2, 0, 0}]]} == {1, 0, 0}]
tbF = Product[Cos[Pi tbY[[k]]/tbLc]^2, {k, 4}];
tbW1[j_] := If[j == 0, 1/2, 1/4];
tbCoef[n_] := Sqrt[tbV] Times @@ (tbW1 /@ n);
check["Parseval on the auxiliary torus: f(y) = Prod_k cos^2(Pi y_k/Lc) = Prod_k (1/2 + (1/4) e^(2 Pi I y_k/Lc) + (1/4) e^(-2 Pi I y_k/Lc)) has the 81 non-zero coefficients c_n = Sqrt[V_perp] Prod_k w(n_k), w(0) = 1/2, w(+-1) = 1/4 (three of them checked by integration), and Sum_n |c_n|^2 == Int b^4 d^4y f^2 == (3/8)^4 V_perp",
  AllTrue[{{0, 0, 0, 0}, {1, 0, -1, 0}, {1, 1, 1, -1}}, Simplify[tbInt[tbPhiC[#] tbF] - tbCoef[#], tbLc > 0 && tbB > 0] === 0 &] &&
  Simplify[Total[tbCoef[#]^2 & /@ Tuples[{-1, 0, 1}, 4]] - tbInt[tbF^2], tbLc > 0 && tbB > 0] === 0 && Simplify[tbInt[tbF^2] - (3/8)^4 tbV, tbLc > 0 && tbB > 0] === 0]
check["the Laplacian of the auxiliary torus: -(1/b^2) Sum_k d^2 phi_n/dy_k^2 == (|q_n|^2/b^2) phi_n (n = (2, -1, 0, 3))",
  Simplify[-Sum[D[tbPhi[{2, -1, 0, 3}], tbY[[k]], tbY[[k]]], {k, 4}]/tbB^2 - (tbQ[{2, -1, 0, 3}].tbQ[{2, -1, 0, 3}]/tbB^2) tbPhi[{2, -1, 0, 3}]] === 0]
tbGS = DiagonalMatrix[{1, -1, -1, -1, -tbB^2, -tbB^2, -tbB^2, -tbB^2}];          (* a static torus of size b0 = tbB *)
tbBox[f_] := With[{gi = Inverse[tbGS], sg = Sqrt[-Det[tbGS]]}, Sum[D[sg gi[[A, A]] D[f, X[[A]]], X[[A]]], {A, 8}]/sg];
tbWave = Exp[-I tbE x0 + I {tbk1, tbk2, tbk3}.{x1, x2, x3} + I tbQ[{1, -2, 0, 2}].tbY];
check["scalar Kaluza-Klein masses: for the static metric diag(1, -1, -1, -1, -b0^2, ..., -b0^2) the Klein-Gordon operator Box + m^2 (Box f = |g|^(-1/2) d_A(|g|^(1/2) g^AB d_B f)) on Exp[-I E x0 + I k.x + I q.y] gives (-E^2 + |k|^2 + |q|^2/b0^2 + m^2) times the wave: in ordinary spacetime a particle of mass Sqrt[m^2 + |q|^2/b0^2]",
  Simplify[(tbBox[tbWave] + tbM^2 tbWave)/tbWave - (-tbE^2 + tbk1^2 + tbk2^2 + tbk3^2 + tbQ[{1, -2, 0, 2}].tbQ[{1, -2, 0, 2}]/tbB^2 + tbM^2)] === 0]
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('25.3 Counting the levels: Jacobi’s r_4(k), Weyl’s law, and the torus today')
text('The Kaluza–Klein masses depend on k = |n|² only. The number of n ∈ ℤ⁴ with |n|² = k is Jacobi’s r_4(k) = 8σ(k) − 32σ(k/4) (σ the '
     'divisor sum, σ(k/4) = 0 unless 4 divides k; section 22.2 computes it as rIIr4 for k ≤ rIIkmax). The number N(R) of n with |n| ≤ R '
     'approaches the volume π²R⁴/2 of the four-dimensional ball (Weyl’s law), so the number of levels with μ_n ≤ μ is N(μ) ≈ V⊥ μ⁴/(32π²) '
     'and their density dN = V⊥ μ³ dμ/(8π²). With M8 = 10 TeV (κ = M8⁻⁶) and the reduced Planck mass M_Pl = 2.435 × 10¹⁸ GeV, '
     'V⊥ = M_Pl²/M8⁶ (sections 17.5 and 22.2). For comparison (CODATA 2018): the proton charge radius 0.8414 fm and the Bohr radius '
     '5.29177210903 × 10⁻¹¹ m; R = ℓ⊥/(2π) is the radius of each circle in the convention V⊥ = (2πR)⁴ of the Particle Data Group, in which '
     'the model has M_D = (2π)^(2/3) M8 (section 24.1).')
inp(r'''
check["Jacobi's four-square theorem: the number of n in Z^4 with |n|^2 = k is r_4(k) = 8 sigma(k) - 32 sigma(k/4) for k = 1..100 (direct count over {-10, ..., 10}^4)",
  With[{nn = Total[Tuples[Range[-10, 10], 4]^2, {2}]}, Table[Count[nn, kk], {kk, 100}]] == rIIr4[[1 ;; 100]]]
tbNR[R_] := 1 + Total[rIIr4[[1 ;; R^2]]];
tbRs = {10, 30, 100, 300, 1000}; tbNRs = tbNR /@ tbRs; tbDev = N[tbNRs/(Pi^2 tbRs^4/2) - 1];
check["Weyl's law: the number N(R) of n in Z^4 with |n| <= R approaches the volume Pi^2 R^4/2 of the four-ball: N(10) = 49689, N(100) = 493490641, and |N(R)/(Pi^2 R^4/2) - 1| < 10^-2, 10^-3, 3 10^-5, 5 10^-6, 10^-6 for R = 10, 30, 100, 300, 1000",
  tbNRs[[1]] == 49689 && tbNRs[[3]] == 493490641 && And @@ MapThread[Abs[#1] < #2 &, {tbDev, {10^-2, 10^-3, 3*^-5, 5*^-6, 10^-6}}]]
{"r_4(k), k = 0..16" -> Join[{1}, rIIr4[[1 ;; 16]]], "N(R) for R = 10, 30, 100, 300, 1000" -> tbNRs, "N(R)/(Pi^2 R^4/2) - 1 for R = 10, 30, 100, 300, 1000" -> tbDev}
''')
inp(r'''
tbLp = rIIVp^(1/4); tbHbarCm = 1.973269804*^-16;                  (* hbar c in GeV m *)
tbNmu[mu_] := rIIVp mu^4/(32 Pi^2);
check["the auxiliary torus today (M8 = 10 TeV, M_Pl = 2.435 10^18 GeV): V_perp = M_Pl^2/M8^6 = 5.929 10^12 GeV^-4, l_perp = V_perp^(1/4) = b0 Lc = 1560.4 GeV^-1 = 3.079 10^-13 m, lowest Kaluza-Klein mass 2 Pi/l_perp = 4.027 MeV; N(mu) = V_perp mu^4/(32 Pi^2) = 1.17 10^13 levels below m_r/2 = 5 GeV (the lattice count of section 22.2 agrees to 10^-5) and N(M8) = M_Pl^2/(32 Pi^2 M8^2)",
  Abs[rIIVp/5.929*^12 - 1] < 10^-3 && Abs[tbLp/1560.4 - 1] < 10^-4 && Abs[tbLp tbHbarCm/3.079*^-13 - 1] < 10^-3 && Abs[2 Pi/tbLp 10^3/4.027 - 1] < 10^-3 &&
  Abs[(1 + Total[rIIr4])/tbNmu[5.] - 1] < 10^-5 && Abs[tbNmu[5.]/1.17*^13 - 1] < 0.01 && Abs[tbNmu[rIIM8]/(rIIMPl^2/(32 Pi^2 rIIM8^2)) - 1] < 10^-12]
tbLevels = Table[{k, rIIr4[[k]], 4 rIIr4[[k]], 2 Pi Sqrt[k]/tbLp 10^3}, {k, 1, 8}];
{"V_perp (GeV^-4), l_perp (GeV^-1), l_perp (m), l_perp (fm), 2 Pi/l_perp (MeV)" -> {rIIVp, tbLp, tbLp tbHbarCm, tbLp tbHbarCm 10^15, 2 Pi/tbLp 10^3},
 "l_perp/proton radius, Bohr radius/l_perp, R = l_perp/(2 Pi) (GeV^-1, m), M_D = (2 Pi)^(2/3) M8 (GeV)" -> {tbLp tbHbarCm/0.8414*^-15, 5.29177210903*^-11/(tbLp tbHbarCm), tbLp/(2 Pi), tbLp tbHbarCm/(2 Pi), (2 Pi)^(2/3) rIIM8},
 "N(mu) = V_perp mu^4/(32 Pi^2) for mu = 5 GeV, m_h, M8; lattice count below 5 GeV" -> {tbNmu[5.], tbNmu[mHb], tbNmu[rIIM8], 1 + Total[rIIr4]},
 "first levels: k, r_4(k), Dirac fermions 4 r_4(k), mass 2 Pi Sqrt[k]/l_perp (MeV), k = 1..8" -> tbLevels}
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('25.4 Newton’s law on the torus: images and Kaluza–Klein modes')
text('In linearized eight-dimensional gravity (no stabilizing potential) a static mass M has the potential Φ = −5κM/(32π³ r⁵) at the '
     'distance r (κ = 8πG8): the Green’s function 1/(5 Ω_6 r⁵) of the seven-dimensional Laplacian (Ω_6 = 16π³/15 the area of the unit '
     'six-sphere) times 8πG8 (D − 3)/(D − 2) = (5/6) κ. On the torus the mass has images at the lattice points ℓ⊥ n (n ∈ ℤ⁴), and a second '
     'mass on the brane at the distance r along the brane feels Φ(r) = −(5κM/(32π³)) Σ_n (r² + ℓ⊥² |n|²)^(−5/2). Poisson summation turns the '
     'sum over images into a sum over the Kaluza–Klein modes, Σ_n (r² + |n|²)^(−5/2) = (4π²/(3r)) Σ_m e^{−2π r |m|} (lengths in units of '
     'ℓ⊥), because the four-dimensional Fourier transform of (r² + |y|²)^(−5/2) is (4π²/(3r)) e^{−p r}. Hence '
     'Φ(r) = −(5/3) G4 M Σ_m e^{−μ_m r}/r with G4 = G8/V⊥: every Kaluza–Klein level (mass μ_m = 2π|m|/ℓ⊥) adds a Yukawa potential, and the '
     'zero mode gives (5/3) G4 M/r. Of the factor 5/3, 1 is the four-dimensional graviton and 2/3 the radion (massless without the '
     'stabilizer), whose coupling (φ/(√3 M_Pl)) T^μ_μ (section 20.3) gives −(2/3) G4 M M′ e^{−m_r r}/r; with the stabilizing potential the '
     'radion has the mass m_r, and at distances r ≫ ℓ⊥, 1/m_r Newton’s law −G4 M M′/r holds. Near the source (r ≪ ℓ⊥) the '
     'eight-dimensional law returns: r⁵ S(r) → 1 + r⁵ Z(5/2), with Z(5/2) = Σ_(n≠0) |n|⁻⁵ = 8(1 − 4^(−3/2)) ζ(5/2) ζ(3/2) = 7 ζ(5/2) ζ(3/2).')
inp(r'''
check["Newton's law in eight dimensions: Omega_6 = 2 Pi^(7/2)/Gamma(7/2) = 16 Pi^3/15, r^-5 is harmonic away from r = 0 in seven space dimensions ((1/r^6) d/dr (r^6 d/dr r^-5) = 0), and grad^2 Phi = 8 Pi G8 ((D - 3)/(D - 2)) M delta^7 (D = 8) gives Phi = -(5/6) kappa M/(5 Omega_6 r^5) = -5 kappa M/(32 Pi^3 r^5); the same formula in D = 4 gives -G M/r",
  Simplify[2 Pi^(7/2)/Gamma[7/2] - 16 Pi^3/15] === 0 && Simplify[D[tbr^6 D[tbr^-5, tbr], tbr]] === 0 &&
  Simplify[tbKap (5/6)/(5 (16 Pi^3/15)) - 5 tbKap/(32 Pi^3)] === 0 && Simplify[8 Pi tbGN (1/2)/(1 (2 Pi^(3/2)/Gamma[3/2])) - tbGN] === 0]
tbS[r_?NumericQ] := Module[{kk = N[Range[Length[rIIr4]]], u0 = r^2 + Length[rIIr4] + 1/2},
   r^-5 + Total[rIIr4 (r^2 + kk)^(-5/2)] + Pi^2 (2/Sqrt[u0] - 2 r^2/(3 u0^(3/2)))];      (* images: Sum over Z^4 of (r^2 + |n|^2)^(-5/2) *)
tbSkk[r_?NumericQ] := 4 Pi^2/(3 r) (1 + Total[rIIr4[[1 ;; 20000]] Exp[-2 Pi r Sqrt[N[Range[20000]]]]]);   (* Kaluza-Klein modes *)
tbPois = Table[With[{si = tbS[r], sk = tbSkk[r]}, {si, sk, si/sk - 1}], {r, {1/4, 1/2, 1, 2}}];
check["images = Kaluza-Klein modes (Poisson summation; lengths in units of l_perp): Sum_n (r^2 + |n|^2)^(-5/2) == (4 Pi^2/(3 r)) Sum_m Exp[-2 Pi r |m|] to 10^-8 at r = 1/4, 1/2, 1, 2 (left: shells with r_4(k) to k = 1.5 10^6 and the continuum tail Pi^2 Int k (r^2 + k)^(-5/2) dk; right: shells to k = 20000)",
  AllTrue[tbPois[[All, 3]], Abs[#] < 10^-8 &]]
tbZ52 = N[7 Zeta[5/2] Zeta[3/2], 20];
tbZ52s = Total[rIIr4 N[Range[Length[rIIr4]]]^(-5/2)] + 2 Pi^2/Sqrt[Length[rIIr4] + 1/2];
tbNear = Table[(r^5 tbS[r] - 1)/(r^5 tbZ52), {r, {1/20, 1/10}}];
tbFar = Table[(3 r tbSkk[r]/(4 Pi^2) - 1)/(8 Exp[-2 Pi r]), {r, {1, 2, 3}}];
check["near and far field: Z(5/2) = Sum_(n != 0) |n|^-5 = Sum_k r_4(k) k^(-5/2) = 7 zeta(5/2) zeta(3/2) = 24.531 (shell sum with tail, to 10^-8); near the source r^5 S(r) - 1 -> r^5 Z(5/2) (the ratio lies between 0.98 and 1 at r = 1/10 and 1/20, closer to 1 at 1/20): the eight-dimensional law; far from it 3 r S/(4 Pi^2) - 1 -> 8 Exp[-2 Pi r], the 8 lightest Kaluza-Klein modes (ratio 1 to 2 10^-2 at r = 2 and 3)",
  Abs[tbZ52s/tbZ52 - 1] < 10^-8 && Abs[tbZ52 - 24.531] < 10^-3 && AllTrue[tbNear, 0.98 < # < 1 &] && tbNear[[1]] > tbNear[[2]] &&
  AllTrue[tbFar[[2 ;; 3]], Abs[# - 1] < 2*^-2 &]]
check["graviton + radion: far from the source the images add up to (5 kappa M/(32 Pi^3)) (4 Pi^2/(3 V_perp r)) == (5/3) G4 M/r with G4 = kappa/(8 Pi V_perp); the radion of section 20.3, coupled to matter by (phi/(Sqrt[3] M_Pl)) T, gives the potential -(1/(Sqrt[3] M_Pl))^2 M M' Exp[-m_r r]/(4 Pi r), and (1/(Sqrt[3] M_Pl))^2/(4 Pi) == (2/3) G4 (8 Pi G4 = 1/M_Pl^2): graviton 1 + massless radion 2/3 = 5/3",
  Simplify[(5 tbKap/(32 Pi^3)) (4 Pi^2/(3 tbVs)) - (5/3) tbKap/(8 Pi tbVs)] === 0 && Simplify[(1/(Sqrt[3] tbMP))^2/(4 Pi) - (2/3)/(8 Pi tbMP^2), tbMP > 0] === 0]
{"images S(r), Kaluza-Klein (4 Pi^2/(3 r)) Sum Exp[-2 Pi r |m|], relative difference, at r = 1/4, 1/2, 1, 2" -> tbPois, "Z(5/2) = 7 zeta(5/2) zeta(3/2), shell sum with tail" -> {tbZ52, tbZ52s},
 "near field (r^5 S - 1)/(r^5 Z(5/2)) at r = 1/20, 1/10; far field (3 r S/(4 Pi^2) - 1)/(8 Exp[-2 Pi r]) at r = 1, 2, 3" -> Join[tbNear, tbFar]}
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('25.5 Kaluza–Klein momentum: conserved in the bulk, not at the brane')
text('The translations y → y + c are isometries of the torus; their conserved charge is the Kaluza–Klein momentum n. A local interaction in '
     'the bulk conserves it: the overlap of three modes, ∫ b⁴ d⁴y φ_n φ_m φ_l*, is 1/√V⊥ if l = n + m and 0 otherwise. The world brane sits '
     'at the point y = 0 and is not invariant under the translations: there every mode has the same value 1/√V⊥, so an interaction on the '
     'brane couples all modes with the same strength (the portals of sections 18.2 and 24.4). The delta function of the torus, '
     'δ⁴(y) = Lc⁻⁴ Σ_n e^{i q_n·y}, contains every mode with the same weight; its average over the torus (its zero mode) is 1/Lc⁴, which is '
     'why a brane source smeared over the torus has the eight-dimensional density ρ_4/V⊥ (section 20.2).')
inp(r'''
tbVtx[n_, m_, l_] := tbInt[tbPhi[n] tbPhi[m] tbPhiC[l]];
check["Kaluza-Klein momentum is conserved in the bulk and not at the brane: Int b^4 d^4y phi_n phi_m phi_l^* == 1/Sqrt[V_perp] for l = n + m and 0 otherwise (exact: n = (1, 0, 0, 2), m = (0, 1, -1, 0), l = (1, 1, -1, 2) and (1, 1, -1, 1)); at y = 0, phi_n phi_m phi_l^* == V_perp^(-3/2) for every n, m, l",
  Simplify[tbVtx[{1, 0, 0, 2}, {0, 1, -1, 0}, {1, 1, -1, 2}] - 1/Sqrt[tbV], tbLc > 0 && tbB > 0] === 0 && tbVtx[{1, 0, 0, 2}, {0, 1, -1, 0}, {1, 1, -1, 1}] === 0 &&
  Simplify[((tbPhi[{1, 0, 0, 2}] tbPhi[{0, 1, -1, 0}] tbPhiC[{1, 1, -1, 1}]) /. Thread[tbY -> 0]) - tbV^(-3/2), tbLc > 0 && tbB > 0] === 0]
check["the delta function of the torus: Int_(-Lc/2)^(Lc/2) Exp[-I q_n y] delta(y) dy == 1 for every integer n, so delta^4(y) = Lc^-4 Sum_n Exp[I q_n.y] (every mode with the weight 1/Lc^4); a brane source rho4 delta^4(y)/b^4 has Int b^4 d^4y (rho4 delta^4(y)/b^4) == rho4, i.e. the average rho4/V_perp over the torus",
  Simplify[Integrate[Exp[-I 2 Pi tbn y1D/tbLc] DiracDelta[y1D], {y1D, -tbLc/2, tbLc/2}, Assumptions -> tbLc > 0 && Element[tbn, Integers]]] === 1 &&
  Simplify[Integrate[tbB^4 tbRho DiracDelta[x4] DiracDelta[x5] DiracDelta[x6] DiracDelta[x7]/tbB^4, {x4, -tbLc/2, tbLc/2}, {x5, -tbLc/2, tbLc/2},
     {x6, -tbLc/2, tbLc/2}, {x7, -tbLc/2, tbLc/2}, Assumptions -> tbLc > 0]] === tbRho]
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('25.6 The world brane: embedding, induced metric, volume element and extrinsic curvature')
text('The embedding X(x) = (x0, x1, x2, x3, 0, 0, 0, 0) of W has the Jacobian J^A_μ = ∂X^A/∂x^μ (an 8×4 matrix of zeros and ones), and the '
     'induced metric is h = Jᵀ g J. For the cosmological metric of section 11.2 the brane sees the Friedmann–Robertson–Walker metric '
     'diag(1, −a², −a², −a²): it does not see b. For the inhomogeneous test metric of section 11.3 it sees diag(n², −a², −a², −a²) evaluated '
     'at x4 = 0. Because g is block diagonal at the brane, √|g| = √|h| b⁴ there, and a brane integral is an eight-dimensional integral with '
     'the factor δ⁴(y)/b⁴. How the brane sits in the bulk is measured by its extrinsic curvature K^(k)_μν = (1/(2b)) ∂_k g_μν at W (for '
     'g_μk = 0 and −g_kl = b² δ_kl at W). The Christoffel symbols {k μν} = (1/(2b²)) ∂_k g_μν at W are the accelerations transverse to the '
     'brane of a free particle that moves along it: if they vanish, the brane is totally geodesic, and a particle that moves along the brane '
     'stays on it.')
inp(r'''
tbJ = Table[If[A == mu, 1, 0], {A, 8}, {mu, 4}];                       (* dX^A/dx^mu for X = (x0, x1, x2, x3, 0, 0, 0, 0) *)
tbOnW = {x4 -> 0, x5 -> 0, x6 -> 0, x7 -> 0};
tbInduced[g_] := Simplify[Transpose[tbJ].(g /. tbOnW).tbJ];
check["the world brane W = {x4 = x5 = x6 = x7 = 0}: the induced metric h = J^T.g.J of the embedding X = (x0, x1, x2, x3, 0, 0, 0, 0) is, for the cosmological metric of section 11.2, the Friedmann-Robertson-Walker metric diag(1, -a^2, -a^2, -a^2) (independent of b); for the inhomogeneous metric of section 11.3 it is diag(n^2, -a^2, -a^2, -a^2) at x4 = 0 with n = 1 + x0^2/11, a = 1 + x0/3",
  tbInduced[gCos] == DiagonalMatrix[{1, -a3[x0]^2, -a3[x0]^2, -a3[x0]^2}] &&
  Simplify[tbInduced[gPoly] - DiagonalMatrix[{(1 + x0^2/11)^2, -(1 + x0/3)^2, -(1 + x0/3)^2, -(1 + x0/3)^2}]] === ConstantArray[0, {4, 4}]]
tbPt = {x0 -> 1/2, x4 -> 0};
check["brane integrals in eight-dimensional form: at the brane Sqrt|det g| == Sqrt|det h| b^4 with b^4 = Sqrt[det(-g_kl)] (cosmological metric; inhomogeneous metric exactly at x0 = 1/2), so Int_W d^4x Sqrt|h| L == Int d^8X Sqrt|g| (delta^4(y)/b^4) L; test: Int d^4y b^4 (delta^4(y)/b^4) F(y) == F(0) = 1 for F = Exp[x4 + 2 x5] Cos[x6 - x7]",
  Simplify[Sqrt[-Det[gCos]] - Sqrt[-Det[tbInduced[gCos]]] Sqrt[Det[-tbG4]], a3[x0] > 0 && b4[x0] > 0] === 0 &&
  (Sqrt[-Det[gPoly /. tbPt]] - Sqrt[-Det[tbInduced[gPoly] /. tbPt]] Sqrt[Det[-(gPoly[[5 ;; 8, 5 ;; 8]] /. tbPt)]]) === 0 &&
  Simplify[Integrate[tbB^4 (DiracDelta[x4] DiracDelta[x5] DiracDelta[x6] DiracDelta[x7]/tbB^4) Exp[x4 + 2 x5] Cos[x6 - x7], {x4, -1, 1}, {x5, -1, 1},
     {x6, -1, 1}, {x7, -1, 1}]] === 1]
tbChrW[cv_, pt_] := Table[cv["chr"][[k, mu, nu]] /. pt, {k, 5, 8}, {mu, 1, 4}, {nu, 1, 4}];      (* {k mu nu}, k = 4..7, mu, nu = 0..3, at the brane *)
check["how the brane sits in the bulk: the Christoffel symbols {k mu nu} at the brane (k = 4..7, mu, nu = 0..3; = (1/(2 b^2)) d_k g_mu nu for these metrics, b times the extrinsic curvature) all vanish for the cosmological metric: the brane is totally geodesic, and a free particle moving along it stays on it. For the inhomogeneous metric at x0 = 1/2 the only non-zero one is {4 00} = n (d_4 n)/b^2 = 1458/58619: a free particle at rest on the brane is accelerated off it",
  AllTrue[Flatten[tbChrW[cvCos, tbOnW]], # === 0 &] &&
  tbChrW[cvPoly, tbPt] == ReplacePart[ConstantArray[0, {4, 4, 4}], {1, 1, 1} -> 1458/58619] &&
  ((nP D[nP, x4]/bP^2) /. tbPt) == 1458/58619]
{"Christoffel {4 00} at the brane (inhomogeneous metric, x0 = 1/2)" -> N[1458/58619, 16]}
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('25.7 dirac16complex at the brane: brane spinor ⊗ transverse spinor')
text('At the brane the eight Dirac matrices of dirac16complex split into four along W and four transverse to it (section 18): Γ17[A].U = U.R_A '
     'with R_μ = γ_μ ⊗ 1_4, R_k = γ5 ⊗ ρ_k, so ψ = U.(Σ_a ψ_a ⊗ e_a) is four brane Dirac spinors ψ_a, and the index a is a spinor index of the '
     'four transverse directions. Products of two Dirac matrices generate the rotations: Γ_μΓ_ν (the Lorentz group of the brane) act on the '
     'brane spinor index only, Γ_kΓ_l (the rotations about the brane) on the transverse index only, and Γ_μΓ_k, which would turn the brane '
     'into the transverse directions, mix the two: the brane breaks the Lorentz group of the bulk to SO(1,3) × SO(4) (on the torus only the '
     '384 signed permutations of section 25.1 remain of SO(4)). On the transverse index the rotations form Spin(4) = SU(2) × SU(2): '
     'ρ4567 = ρ4ρ5ρ6ρ7 (as a 4×4 matrix equal to −γ5) commutes with every ρ_kρ_l, and on its eigenspace s = ±1 the combinations '
     'ρ_kρ_l + s ρ_mρ_n ((k, l, m, n) an even permutation of (4, 5, 6, 7)) vanish, so only one SU(2) acts there: the four values of a form the '
     'representation (2, 1) ⊕ (1, 2). Every Γ17[A] is a signed permutation matrix, and 2√2 U has only the entries 0, ±1, ±i; both are '
     'printed below.')
inp(r'''
tbSPermQ[m_] := SubsetQ[{-1, 0, 1}, Union[Flatten[m]]] && AllTrue[m, Count[#, 0] == 15 &] && AllTrue[Transpose[m], Count[#, 0] == 15 &];
tbU8 = Simplify[2 Sqrt[2] uInt];
check["explicit forms: every Gamma17[A] is a signed permutation matrix (one entry +-1 in each row and column), and 2 Sqrt[2] U has only the entries 0, +-1, +-I, eight non-zero ones in each row and in each column",
  AllTrue[Table[\[CapitalGamma]17[A], {A, 0, 7}], tbSPermQ] && SubsetQ[{0, 1, -1, I, -I}, Union[Flatten[tbU8]]] &&
  AllTrue[tbU8, Count[#, 0] == 8 &] && AllTrue[Transpose[tbU8], Count[#, 0] == 8 &]]
tbUd = ConjugateTranspose[uInt]; tbKp = KroneckerProduct; tbZ16 = ConstantArray[0, {16, 16}];
check["brane spinor (x) transverse spinor: Gamma17[mu].Gamma17[nu] == U.(gamma_mu.gamma_nu (x) 1_4).U^dag, Gamma17[k].Gamma17[l] == U.(1_4 (x) rho_k.rho_l).U^dag and Gamma17[mu].Gamma17[k] == U.(gamma_mu.gamma5 (x) rho_k).U^dag (mu, nu = 0..3, k, l = 4..7): the Lorentz group of the brane acts on the brane spinor index, the rotations about the brane on the transverse index, and the generators that mix them are broken by the brane",
  And @@ Flatten[Table[Simplify[\[CapitalGamma]17[mu].\[CapitalGamma]17[nu] - uInt.tbKp[g4D[[mu + 1]].g4D[[nu + 1]], ID4].tbUd] == tbZ16, {mu, 0, 3}, {nu, 0, 3}]] &&
  And @@ Flatten[Table[Simplify[\[CapitalGamma]17[k + 3].\[CapitalGamma]17[l + 3] - uInt.tbKp[ID4, rhoK[[k]].rhoK[[l]]].tbUd] == tbZ16, {k, 4}, {l, 4}]] &&
  And @@ Flatten[Table[Simplify[\[CapitalGamma]17[mu].\[CapitalGamma]17[k + 3] - uInt.tbKp[g4D[[mu + 1]].g5D, rhoK[[k]]].tbUd] == tbZ16, {mu, 0, 3}, {k, 4}]]]
tbPs[s_] := (ID4 + s rho4567)/2; tbEven = {{1, 2, 3, 4}, {1, 3, 4, 2}, {1, 4, 2, 3}};       (* (4,5,6,7), (4,6,7,5), (4,7,5,6) *)
check["the transverse index carries Spin(4) = SU(2) x SU(2): rho4567 == -gamma5 (4 x 4), it commutes with every rho_k.rho_l and anticommutes with every rho_k; for each even permutation (k, l, m, n) of (4, 5, 6, 7) the combination rho_k.rho_l + s rho_m.rho_n vanishes on the eigenspace rho4567 = s, and rho_k.rho_l - s rho_m.rho_n acts there with rank 2: the four values of a form (2, 1) + (1, 2)",
  rho4567 == -g5D && And @@ Flatten[Table[rho4567.rhoK[[k]].rhoK[[l]] == rhoK[[k]].rhoK[[l]].rho4567, {k, 4}, {l, 4}]] &&
  And @@ Table[rho4567.rhoK[[k]] == -rhoK[[k]].rho4567, {k, 4}] &&
  And @@ Flatten[Table[(rhoK[[p[[1]]]].rhoK[[p[[2]]]] + s rhoK[[p[[3]]]].rhoK[[p[[4]]]]).tbPs[s] == ConstantArray[0, {4, 4}], {p, tbEven}, {s, {1, -1}}]] &&
  And @@ Flatten[Table[MatrixRank[(rhoK[[p[[1]]]].rhoK[[p[[2]]]] - s rhoK[[p[[3]]]].rhoK[[p[[4]]]]).tbPs[s]] == 2, {p, tbEven}, {s, {1, -1}}]]]
check["chirality at the brane: PL17 == U.(P_L (x) P_+ + P_R (x) P_-).U^dag with P_L,R = (1 -+ gamma5)/2 and P_+- = (1 +- rho4567)/2: the left-handed half psi_L (chi17 = -1, 8 components) is two brane-left-handed Weyl spinors with rho4567 = +1 and two brane-right-handed ones with rho4567 = -1 (each part of rank 4)",
  Simplify[PL17 - uInt.(tbKp[(ID4 - g5D)/2, tbPs[1]] + tbKp[(ID4 + g5D)/2, tbPs[-1]]).tbUd] == tbZ16 &&
  MatrixRank[tbKp[(ID4 - g5D)/2, tbPs[1]]] == 4 && MatrixRank[tbKp[(ID4 + g5D)/2, tbPs[-1]]] == 4]
tbu = {1, 2, 0, -1}; tbPsi = uInt.Flatten[Outer[Times, tbu, UnitVector[4, 1]]];
tbBil = Table[Simplify[Conjugate[tbPsi].\[CapitalGamma]17[0].\[CapitalGamma]17up[mu].tbPsi], {mu, 0, 3}];
check["an explicit spinor: psi = U.(u (x) e_1) with the brane spinor u = (1, 2, 0, -1): psibar psi == ubar u == 4 and psibar Gamma^mu psi == ubar gamma^mu u (mu = 0..3)",
  Simplify[Conjugate[tbPsi].\[CapitalGamma]17[0].tbPsi] == 4 == Conjugate[tbu].g4D[[1]].tbu &&
  tbBil == Table[Conjugate[tbu].g4D[[1]].(\[Eta]17[[mu + 1, mu + 1]] g4D[[mu + 1]]).tbu, {mu, 0, 3}]]
{"Gamma17[A] as signed permutations: column j has the entry sign(c) in row |c|, c listed for A = 0..7" -> Table[Total[\[CapitalGamma]17[A][[All, j]] Range[16]], {A, 0, 7}, {j, 16}],
 "Re(2 Sqrt[2] U)" -> Re[tbU8], "Im(2 Sqrt[2] U)" -> Im[tbU8],
 "example psi = U.(u (x) e_1), u = (1, 2, 0, -1): psibar psi, psibar Gamma^mu psi (mu = 0..3)" -> Join[{Simplify[Conjugate[tbPsi].\[CapitalGamma]17[0].tbPsi]}, tbBil]}
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('25.8 Matter on the brane, and the brane in four-dimensional physics')
text('A brane action S_W = ∫_W d⁴x √|h| L_W gives the eight-dimensional energy–momentum tensor T^AB = (δ⁴(y)/b⁴) T_W^μν δ^A_μ δ^B_ν with '
     'T_W μν = (2/√|h|) δS_W/δh^μν: it has no components along the torus. The tension (L_W = −f_W) has T_W μν = f_W h_μν: the energy density '
     'f_W and the pressure −f_W along the brane (w = −1). Smeared over the torus (δ⁴(y) → Lc⁻⁴, section 20.2) the tension gives ε = f_W/V⊥, '
     'pa = −f_W/V⊥, pb = 0, and brane radiation pa = ε/3, pb = 0, while radiation in the bulk, isotropic in the seven space directions, has '
     'pa = pb = ε/7. In the static modern regime the tension is a part of the four-dimensional cosmological constant, Λ4 = Λ + κU0 + κf_W/V⊥ '
     '(section 20.3): a four-dimensional vacuum energy f_W, for f_W = M8⁴ about 4 × 10⁶² times the observed dark-energy density. The brane '
     'plasma loses energy to the bulk by emitting Kaluza–Klein modes, at the relative rate T⁷/M8⁶ = (number of levels V⊥T⁴) × (rate per level '
     'T³/M_Pl²) (estimate, numerical factors omitted), which equals the expansion rate T²/M_Pl at the normalcy temperature '
     'T* = (M8⁶/M_Pl)^(1/5). Finally, a thin brane is an idealization: the linearized field of its tension falls like κf_W/(4π²|y|²) (the '
     'Green’s function −1/(4π²|y|²) of the four transverse directions) and is of order one within |y| ≈ (κf_W)^(1/2)/(2π) = 1/(2π M8) of the '
     'brane for f_W = M8⁴ (estimate); there the thin-brane description fails, and the smeared description of sections 20–23 holds far from it.')
inp(r'''
check["brane tension: for S_W = -Int d^4x Sqrt|h| f_W and h = diag(1/u0, -1/u1, -1/u2, -1/u3) (h^00 = u0, h^ii = -u_i), T_mu nu = (2/Sqrt|h|) dS_W/dh^mu nu == f_W h_mu nu, T^0_0 == f_W (energy density) and T^1_1 == f_W (pressure -f_W: w = -1 along the brane)",
  Module[{sq = 1/Sqrt[tbu0 tbu1 tbu2 tbu3], t00, t11},
   t00 = (2/sq) D[-tbfW sq, tbu0]; t11 = (2/sq) (-D[-tbfW sq, tbu1]);
   Simplify[{t00 - tbfW/tbu0, t11 + tbfW/tbu1, tbu0 t00, -tbu1 t11}, tbu0 > 0 && tbu1 > 0 && tbu2 > 0 && tbu3 > 0] === {0, 0, tbfW, tbfW}]]
check["sources smeared over the torus (W(a, b) the comoving energy, section 20.2): the brane tension W = a^3 f_W/Lc^4 gives eps = f_W/V_perp, pa = -eps, pb = 0; brane radiation W = a^3 r0 a^-4/Lc^4 gives pa = eps/3, pb = 0; radiation of the bulk, isotropic in the seven space directions, has <cos^2 theta> = 3/7 and <sin^2 theta> = 4/7 for the weight cos^2 theta sin^3 theta of the momentum sphere, so pa = (1/3)(3/7) eps = eps/7 and pb = (1/4)(4/7) eps = eps/7",
  Simplify[feStress[Function[{tbAa, tbBb}, tbAa^3 tbfW/tbLc^4]] - {tbfW/(bS^4 tbLc^4), -tbfW/(bS^4 tbLc^4), 0}] === {0, 0, 0} &&
  With[{st = feStress[Function[{tbAa, tbBb}, tbAa^3 tbr0 tbAa^-4/tbLc^4]]}, Simplify[{st[[2]] - st[[1]]/3, st[[3]]}] === {0, 0}] &&
  With[{c2 = Integrate[Cos[tbth]^4 Sin[tbth]^3, {tbth, 0, Pi/2}]/Integrate[Cos[tbth]^2 Sin[tbth]^3, {tbth, 0, Pi/2}]}, c2 == 3/7 && c2/3 == 1/7 && (1 - c2)/4 == 1/7]]
tbTstar = (rIIM8^6/rIIMPl)^(1/5);
check["the brane in four-dimensional physics: (i) the tension is vacuum energy, kappa f_W/V_perp == f_W/M_Pl^2 (V_perp = kappa M_Pl^2), and f_W = M8^4 is between 3 10^62 and 5 10^62 times rho_Lambda = 3 H0^2 M_Pl^2 Omega_Lambda (section 24.5); (ii) normalcy temperature: (T^7/M8^6)/(T^2/M_Pl) == 1 at T* = (M8^6/M_Pl)^(1/5) = 13.3 GeV, and V_perp T^4 T^3/M_Pl^2 == T^7/M8^6; (iii) thin brane: (1/y^3) d/dy (y^3 d/dy y^-2) == 0 and 2 Pi^2 y^3 d/dy y^-2 == -4 Pi^2 (the Green's function -1/(4 Pi^2 |y|^2) of the four transverse directions), and Sqrt[kappa f_W] == 1/M8 for f_W = M8^4",
  Simplify[(tbKap tbfW/tbVs /. tbVs -> tbKap tbMP^2) - tbfW/tbMP^2] === 0 && 3*^62 < rIIM8^4/dkRhoL < 5*^62 &&
  Simplify[((tbT^7/tbM8^6)/(tbT^2/tbMP)) /. tbT -> (tbM8^6/tbMP)^(1/5), tbM8 > 0 && tbMP > 0] === 1 && Abs[tbTstar - 13.3] < 0.05 &&
  Simplify[(tbMP^2/tbM8^6) tbT^4 tbT^3/tbMP^2 - tbT^7/tbM8^6] === 0 &&
  Simplify[D[tby^3 D[tby^-2, tby], tby]] === 0 && Simplify[2 Pi^2 tby^3 D[tby^-2, tby]] === -4 Pi^2 && Simplify[Sqrt[tbM8^-6 tbM8^4] - 1/tbM8, tbM8 > 0] === 0]
{"T* = (M8^6/M_Pl)^(1/5) (GeV), f_W/rho_Lambda for f_W = M8^4, 1/(2 Pi M8) (GeV^-1)" -> {tbTstar, rIIM8^4/dkRhoL, 1/(2 Pi rIIM8)}}
''')
