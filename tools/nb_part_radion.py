# Section 26 of 7+1.nb: the radion -- the size of the four extra dimensions -- defined precisely and checked, with the four extra
# dimensions either a 4-torus (the space of 7+1.nb) or a 4-sphere (for comparison).  Appended on 28 October 2026.  Executed by
# tools/make_7p1_nb.py, which provides section(), subsection(), text(), item(), inp() and B.  The names of section 26 begin with rd.

section('26. The radion: the size of the extra space, on the 4-torus and on the 4-sphere')
text([('About section 26. ', B),
      'Section 26 was appended on 28 October 2026. The notebook as it was before (sections 0–25, with the summary of all checks numbered 26; '
      'evaluated) is kept in the workspace as backups/7+1_2026-10-28_before_radion.nb (excluded from the repository by .gitignore); nothing '
      'in sections 0–25 is changed. The section checks the ideas of a note of the author on the radion (a Kaluza–Klein ansatz with the '
      'factors e^(2αφ) and e^(2βφ), the Einstein frame, the canonical normalization, the physical volume, the radion as a modulus, the fifth '
      'force, stabilization, the warped alternative of Randall and Sundrum with the stabilization of Goldberger and Wise, and the mixing of '
      'the radion with the Higgs boson on the brane), states the corrections that the checks require, and works '
      'every idea out for the two cases in which the four extra dimensions are a 4-TORUS T⁴ (the auxiliary space of 7+1.nb, section 25) and a '
      '4-SPHERE S⁴ (for comparison). Coordinates as everywhere: X = {x0, …, x7}, x0 is time (the direction of the octonion unit e[0] = 1), '
      'x1, x2, x3 ordinary space, x4 … x7 the four extra space-like coordinates — on T⁴ the periodic coordinates of section 25, on S⁴ the '
      'angles x4, x5, x6 ∈ [0, π], x7 ∈ [0, 2π), with the north pole at x4 = 0. η17 = diag(+1, −1, …, −1) is simply the squares of the '
      'octonion basis; xAct declares the spacetime with the signature list {7, 1, 0}: seven space-like and one time-like dimension.'])
text([('The radion. ', B),
      'Write the metric as ds² = e^(−4β(x)) g_E,μν(x) dx^μ dx^ν − e^(2β(x)) γ_ij(y) dy^i dy^j (x = (x0, …, x3), y = (x4, …, x7)), with a fixed '
      'reference metric γ of the extra space: γ = b0² δ_ij on T⁴ (period Lc) and γ = r0² γ̂ on S⁴ (γ̂ the metric of the unit sphere, '
      'γ̂ = diag(1, sin²x4, sin²x4 sin²x5, sin²x4 sin²x5 sin²x6)). The RADION is the field β(x): the size is b = b0 e^β on T⁴, the radius '
      'r = r0 e^β on S⁴, and the proper volume V = V0 e^(4β) with V0 = b0⁴ Lc⁴ (T⁴) or V0 = (8π²/3) r0⁴ (S⁴). β is the same field as '
      'β = ln(b/b0) of section 20.3. The factor e^(−4β) in front of g_E puts the four-dimensional action into the EINSTEIN FRAME: the '
      'coefficient of the four-dimensional curvature, proportional to e^(2·(−2β)) e^(4β), is then constant. In the notation of the note, '
      'ds² = e^(2αφ) g + e^(2βφ) γ with 2α + nβ = 0 (n = 4 extra dimensions).'])

# ------------------------------------------------------------------------------------------------------------------------
subsection('26.1 The definition: Einstein frame, canonical normalization, physical volume')
text('A small engine for diagonal metrics in any dimension (rdChr, rdRic, rdScal, rdEin: Christoffel symbols, Ricci tensor, scalar '
     'curvature, mixed Einstein tensor G^i_i, with the curvature conventions of section 11) is checked against curvature[] of section 11 '
     'on the cosmological metric. The claim to be checked is the REDUCTION of the eight-dimensional Einstein–Hilbert action (with the '
     'bulk cosmological constant Λ and the stabilizing potential U of section 20): with M_Pl² = V0/κ,')
text('      −(1/(2κ)) ∫ d⁴y √|G| (R + 2Λ) − ∫ d⁴y √|G| U  =  √|g_E| [ −(M_Pl²/2) R_E + 6 M_Pl² g_E^μν ∂_μβ ∂_νβ − V(β) ]  + (a total derivative),')
text('      V(β) = e^(−4β) [ M_Pl² Λ + V0 U − (M_Pl²/2) R_int ],   R_int = 0 (T⁴),   R_int = 12/(r0² e^(2β)) = 12/r² (S⁴).')
text('The kinetic term 6 M_Pl² (∂β)² = ½ (∂φ)² defines the CANONICAL RADION φ = √12 M_Pl β (section 20.3); the curvature of the sphere '
     'adds the negative potential −6 M_Pl² e^(−6β)/r0². The check is done on the homogeneous metrics g_E = diag(N², −a², −a², −a²) with '
     'N(x0), a(x0), β(x0) arbitrary: the Euler–Lagrange expressions of the difference of the two sides vanish identically, so the difference '
     'is a total derivative. For n extra dimensions (a torus Tⁿ, D = 4 + n) the coefficient 6 = 12/2 becomes n(n+2)/4, i.e. '
     '−(M_Pl²/2) R_E + (M_Pl²/2)(n(n+2)/2)(∂β)². The note writes β = β_n φ with β_n = √(2/(n(n+2))) and states the reduced action '
     '(M_Pl²/2)[R − ½(∂φ)²] (in its signature); the check shows that this β_n gives (M_Pl²/2)[R − (∂φ)²] — the form with ½ needs '
     'β_n = 1/√(n(n+2)). For n = 4 the note’s values α = −1/√3, β = 1/√12 satisfy 2α + 4β = 0, and its φ is the canonical radion in units '
     'of M_Pl: M_Pl φ_note = √12 M_Pl ln(b/b0) = φ; the physical volume is V0 e^(4β φ_note) = V0 e^((2/√3) φ_note) = V0 e^(2φ/(√3 M_Pl)).')
inp(r'''
rdChr[gd_List, xs_List] := Module[{n = Length[xs]}, Table[Which[a == b == c, D[gd[[a]], xs[[a]]]/(2 gd[[a]]), a == b, D[gd[[a]], xs[[c]]]/(2 gd[[a]]),
     a == c, D[gd[[a]], xs[[b]]]/(2 gd[[a]]), b == c, -D[gd[[b]], xs[[a]]]/(2 gd[[a]]), True, 0], {a, n}, {b, n}, {c, n}]];
rdRic[gd_List, xs_List, simp_: Simplify] := Module[{n = Length[xs], G = rdChr[gd, xs]},
   Table[simp[Sum[D[G[[r, m, s]], xs[[r]]] - D[G[[r, r, s]], xs[[m]]] + Sum[G[[r, r, l]] G[[l, m, s]] - G[[r, m, l]] G[[l, r, s]], {l, n}], {r, n}]], {s, n}, {m, n}]];
rdScal[gd_List, xs_List, simp_: Simplify] := Module[{ric = rdRic[gd, xs, simp]}, simp[Sum[ric[[i, i]]/gd[[i]], {i, Length[xs]}]]];
rdEin[gd_List, xs_List, simp_: Simplify] := Module[{ric = rdRic[gd, xs, simp], rs}, rs = simp[Sum[ric[[i, i]]/gd[[i]], {i, Length[xs]}]];
   Table[simp[ric[[i, i]]/gd[[i]] - rs/2], {i, Length[xs]}]];
rdEL[L_, q_] := Simplify[D[L, q[x0]] - D[D[L, q'[x0]], x0] + D[D[L, q''[x0]], {x0, 2}]];
rdGam4 = {1, Sin[x4]^2, Sin[x4]^2 Sin[x5]^2, Sin[x4]^2 Sin[x5]^2 Sin[x6]^2};
check["section 26 engine for diagonal metrics: rdScal and rdEin reproduce curvature[] of section 11 for the cosmological metric (R and every G^A_A)",
  Simplify[rdScal[Diagonal[gCos], X] - cvCos["R"]] === 0 && Simplify[rdEin[Diagonal[gCos], X] - Diagonal[cvCos["G"]]] === ConstantArray[0, 8]]
rdRicS4 = rdRic[Join[{1, -1, -1, -1}, -rdR0^2 rdGam4], X];
check["the 4-sphere of radius r0: Ricci_ij = 3 gamma_ij (gamma the unit-sphere metric), Ricci_ij = 0 for i != j, so R(S^4) = 12/r0^2 (in the signature of 7+1.nb the space-like block gives R = -12/r0^2, also with curvature[] of section 11); the volume of the unit S^4 is 8 Pi^2/3",
  And @@ Table[Simplify[rdRicS4[[i, i]] - 3 rdGam4[[i - 4]]] === 0, {i, 5, 8}] && And @@ Flatten[Table[If[i != j, rdRicS4[[i, j]] === 0, True], {i, 8}, {j, 8}]] &&
  Simplify[curvature[DiagonalMatrix[Join[{1, -1, -1, -1}, -rdR0^2 rdGam4]]]["R"] + 12/rdR0^2] === 0 &&
  Integrate[Sin[x4]^3 Sin[x5]^2 Sin[x6], {x4, 0, Pi}, {x5, 0, Pi}, {x6, 0, Pi}, {x7, 0, 2 Pi}] == 8 Pi^2/3]
rdMini[gam_, rint_, v0_] := Module[{gd, rs, lj, le, re, vp, mp2},
   gd = Join[{rdLp[x0]^2 E^(-4 rdBt[x0])}, -rdSa[x0]^2 E^(-4 rdBt[x0]) {1, 1, 1}, -E^(2 rdBt[x0]) gam];
   rs = rdScal[gd, X]; mp2 = v0/rdKap;
   lj = -(1/(2 rdKap)) rdLp[x0] rdSa[x0]^3 E^(-4 rdBt[x0]) v0 (rs + 2 rdLam) - rdLp[x0] rdSa[x0]^3 E^(-4 rdBt[x0]) v0 rdU[rdBt[x0]];
   re = rdScal[{rdLp[x0]^2, -rdSa[x0]^2, -rdSa[x0]^2, -rdSa[x0]^2}, Take[X, 4]];
   vp = E^(-4 rdBt[x0]) (mp2 rdLam + v0 rdU[rdBt[x0]]) - mp2/2 E^(-4 rdBt[x0]) rint;
   le = -(mp2/2) rdLp[x0] rdSa[x0]^3 re + 6 mp2 rdSa[x0]^3 rdBt'[x0]^2/rdLp[x0] - rdLp[x0] rdSa[x0]^3 vp;
   {FreeQ[rs, x4 | x5 | x6 | x7], rdEL[lj - le, #] & /@ {rdLp, rdSa, rdBt}}];
check["radion on T^4: for ds^2 = e^(-4 beta) (N^2 dx0^2 - a^2 dx.dx) - e^(2 beta) b0^2 dy.dy the reduced action -(1/(2 kappa)) Int(R + 2 Lambda) - Int U equals -(M_Pl^2/2) R_E + 6 M_Pl^2 (d beta)^2 - e^(-4 beta)(M_Pl^2 Lambda + V0 U) up to a total derivative (Euler-Lagrange expressions for N, a, beta vanish), M_Pl^2 = V0/kappa, V0 = b0^4 Lc^4",
  rdMini[rdB0^2 {1, 1, 1, 1}, 0, rdB0^4 rdLc^4] === {True, {0, 0, 0}}]
check["radion on S^4: the same reduction for ds^2 = e^(-4 beta) g_E - e^(2 beta) r0^2 gamma (unit-sphere gamma) holds with the additional potential -(M_Pl^2/2) e^(-4 beta) 12/(r0^2 e^(2 beta)) = -6 M_Pl^2 e^(-6 beta)/r0^2 from the curvature of the sphere, V0 = (8 Pi^2/3) r0^4; R is independent of the angles",
  rdMini[rdR0^2 rdGam4, 12/(rdR0^2 E^(2 rdBt[x0])), rdR0^4 8 Pi^2/3] === {True, {0, 0, 0}}]
rdKinOK[n_, beta_, coef_] := Module[{xs = Join[{x0}, Table[ToExpression["rdY" <> ToString[i]], {i, 3 + n}]], gd, lj, le, re},
   gd = Join[{rdLp[x0]^2 E^(-n beta rdPh[x0])}, -rdSa[x0]^2 E^(-n beta rdPh[x0]) {1, 1, 1}, Table[-E^(2 beta rdPh[x0]), {n}]];
   lj = -(1/(2 rdKap)) rdLp[x0] rdSa[x0]^3 E^(-n beta rdPh[x0]) rdScal[gd, xs];
   re = rdScal[{rdLp[x0]^2, -rdSa[x0]^2, -rdSa[x0]^2, -rdSa[x0]^2}, Take[xs, 4]];
   le = -(1/(2 rdKap)) rdLp[x0] rdSa[x0]^3 re + (1/(2 rdKap)) coef rdSa[x0]^3 rdPh'[x0]^2/rdLp[x0];
   AllTrue[rdEL[lj - le, #] & /@ {rdLp, rdSa, rdPh}, # === 0 &]];
rdKinC[n_] := Module[{xs = Join[{x0}, Table[ToExpression["rdY" <> ToString[i]], {i, 3 + n}]], gd, lj, le, re, cc},
   gd = Join[{rdLp[x0]^2 E^(-n rdPh[x0])}, -rdSa[x0]^2 E^(-n rdPh[x0]) {1, 1, 1}, Table[-E^(2 rdPh[x0]), {n}]];
   lj = -(1/(2 rdKap)) rdLp[x0] rdSa[x0]^3 E^(-n rdPh[x0]) rdScal[gd, xs];
   re = rdScal[{rdLp[x0]^2, -rdSa[x0]^2, -rdSa[x0]^2, -rdSa[x0]^2}, Take[xs, 4]];
   le = -(1/(2 rdKap)) rdLp[x0] rdSa[x0]^3 re + (1/(2 rdKap)) cc rdSa[x0]^3 rdPh'[x0]^2/rdLp[x0];
   Simplify[cc /. First[Solve[rdEL[lj - le, rdLp] == 0, cc]]]];
rdKinCs = Table[rdKinC[n], {n, 4}];
check["n extra dimensions (torus T^n, D = 4 + n, n = 1..4, Einstein frame e^(-n B) g_E): the kinetic term is (M_Pl^2/2)(n(n+2)/2)(dB)^2; with the note's B = beta_n phi, beta_n = Sqrt[2/(n(n+2))], it is (M_Pl^2/2)(d phi)^2, not (M_Pl^2/2)(1/2)(d phi)^2; the form with 1/2 needs beta_n = 1/Sqrt[n(n+2)]",
  rdKinCs === Table[n (n + 2)/2, {n, 4}] &&
  And @@ Table[rdKinOK[n, 1, n (n + 2)/2] && rdKinOK[n, Sqrt[2/(n (n + 2))], 1] && ! rdKinOK[n, Sqrt[2/(n (n + 2))], 1/2] && rdKinOK[n, 1/Sqrt[n (n + 2)], 1/2], {n, 4}]]
check["the note's eight-dimensional values alpha = -1/Sqrt[3], beta = 1/Sqrt[12]: 2 alpha + 4 beta == 0 (Einstein frame), 4 beta == 2/Sqrt[3] (V = V0 e^(4 beta phi)), and Sqrt[12] beta == 1, so M_Pl phi_note = Sqrt[12] M_Pl ln(b/b0) is the canonical radion of section 20.3",
  2 (-1/Sqrt[3]) + 4/Sqrt[12] == 0 && 4/Sqrt[12] == 2/Sqrt[3] && Sqrt[12]/Sqrt[12] == 1]
{"general n: kinetic coefficient n(n+2)/2 for n = 1..4" -> N[rdKinCs], "the note's beta_n^2 = 2/(n(n+2)) and the canonical 1/(n(n+2)), n = 4" -> N[{1/rdKinCs[[4]], 1/(2 rdKinCs[[4]])}]}
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('26.2 The potential of the radion on T⁴ and on S⁴, and its mass')
text('Each term of the action scales with its own power of e^β in the Einstein frame: the bulk cosmological constant and any energy density '
     'U(b) of the bulk as e^(−4β), the curvature of the sphere and any gradient energy along the extra space as e^(−6β), a magnetic 4-form '
     'flux through S⁴ (a field that 7+1.nb does not contain, used here only as an illustration) as e^(−12β). With U(b) = U0 + U1 β + '
     '(K_U/2) β² (β = ln(b/b0); section 20.3 has U1 = 0) the radion mass is m_r² = V″(φ) = V″(β)/(12 M_Pl²). Four cases are checked, each '
     'in two independent ways — from the Einstein-frame potential V(β), and from the eight-dimensional field equations: on S⁴ the Einstein '
     'tensor gains (6, 6, 6, 6, 3, 3, 3, 3)/b² over the torus (for the cosmological metric with the sphere of radius b), R gains −12/b², '
     'and the radion equation of section 20.3 becomes β″ + 3Ha β′ + 4β′² = κ(ε/6 − pa/2 + pb/3) + Λ/3 − 3/b².')
item('T⁴ with U: a static flat four-dimensional vacuum at b = b0 needs Λ + κU0 = 0 (with f_W = 0) and U1 = 0; m_r² = κ K_U/12, as in '
     'section 20.3. Without U the torus is a MODULUS (V = 0 if Λ = 0), or runs away (Λ ≠ 0).')
item('S⁴ with U: the curvature must be balanced: Λ + κU0 = 6/r0² and κU1 = −12/r0², and the mass is m_r² = κ K_U/12 − 2/r0²: a sphere '
     'needs a stiffer stabilizer than a torus.')
item('S⁴ with Λ only: dS4 × S⁴ with r0² = 9/Λ and H² = Λ/9 = 1/r0² solves the eight-dimensional equations, but it sits at a MAXIMUM of '
     'V: m² = −6H² (tachyonic) — the size of a sphere held by Λ alone is unstable.')
item('S⁴ with a 4-form flux (Freund–Rubin; illustration): Minkowski × S⁴ needs Λ = 9/(2r0²) and κρ_F = 3/(2r0²) (ρ_F the energy density '
     'of the flux, with pa = −ρ_F, pb = +ρ_F), and the radion is stable, m² = 6/r0².')
text('The note’s eight-dimensional Goldberger–Wise potential V(φ) = ∫d⁴y √γ e^(4βφ)[½ e^(−2βφ)(∇Φ)² + ½ m8² Φ²] lacks the Einstein-frame '
     'factor e^(4αφ) = e^(−8βφ) of √|G|: for a static profile Φ(y) the energy is e^(−4β)∫d⁴y √γ [½ e^(−2β)|∇Φ|²_γ + ½ m² Φ²] (e^(−4βφ) in '
     'the note’s notation) — the gradient part scales as e^(−6β), like the curvature, the mass part as e^(−4β), like U.')
inp(r'''
rdES = rdEin[Join[{1, -a3[x0]^2, -a3[x0]^2, -a3[x0]^2}, -b4[x0]^2 rdGam4], X]; rdET = rdEin[Diagonal[gCos], X];
check["the cosmological metric with the sphere of radius b(x0) (ds^2 = dx0^2 - a^2 dx.dx - b^2 dOmega_4^2): G^A_A = (torus value) + (6, 6, 6, 6, 3, 3, 3, 3)/b^2 and R = (torus value) - 12/b^2",
  toH[rdES - rdET] == {6, 6, 6, 6, 3, 3, 3, 3}/b4[x0]^2 && Simplify[toH[rdScal[Join[{1, -a3[x0]^2, -a3[x0]^2, -a3[x0]^2}, -b4[x0]^2 rdGam4], X] - cvCos["R"]] + 12/b4[x0]^2] === 0]
rdE00 = toH[rdES[[1]]] /. b4[x0] -> rdb; rdE11 = toH[rdES[[2]]] /. b4[x0] -> rdb; rdE44 = toH[rdES[[5]]] /. b4[x0] -> rdb;
rdSolS = First[Solve[{rdE11 - rdLam == -rdKap rdPa, rdE44 - rdLam == -rdKap rdPb}, {dHa, dHb}]];
rdRadS = Expand[(dHb /. rdSolS) /. Ha^2 -> (rdKap rdEps + rdLam - 12 Ha Hb - 6 Hb^2 - 6/rdb^2)/3];
check["radion equation on S^4 (Einstein gravity; 11 and 44 equations with the 00 constraint): Hb' + 3 Ha Hb + 4 Hb^2 == kappa (eps/6 - pa/2 + pb/3) + Lambda/3 - 3/b^2",
  Simplify[rdE00 - (3 (Ha^2 + 4 Ha Hb + 2 Hb^2) + 6/rdb^2)] === 0 && Simplify[rdRadS - (-3 Ha Hb - 4 Hb^2 + rdKap (rdEps/6 - rdPa/2 + rdPb/3) + rdLam/3 - 3/rdb^2)] === 0]
rdUq = rdU0 + rdU1 rdBe + rdKU/2 rdBe^2;
rdV[curv_, flux_] := rdMP2 (E^(-4 rdBe) (rdLam + rdKap rdUq) - curv 6 E^(-6 rdBe)/rdR0^2 + flux E^(-12 rdBe));
rdRHS[eps_, pa_, pb_, lam_] := rdKap (eps/6 - pa/2 + pb/3) + lam/3 - 3/(rdR0^2 E^(2 rdBe));
rdSolT = Solve[{(rdV[0, 0] /. rdBe -> 0) == 0, (D[rdV[0, 0], rdBe] /. rdBe -> 0) == 0}, {rdLam, rdU1}];
check["radion on T^4 with U = U0 + U1 beta + (K_U/2) beta^2: a flat vacuum at b = b0 needs Lambda + kappa U0 = 0 and U1 = 0, and m_r^2 = V''(beta)/(12 M_Pl^2) = kappa K_U/12 (section 20.3)",
  Length[rdSolT] == 1 && Simplify[(rdLam /. First[rdSolT]) + rdKap rdU0] === 0 && (rdU1 /. First[rdSolT]) === 0 &&
  Simplify[(D[rdV[0, 0], {rdBe, 2}] /. rdBe -> 0)/(12 rdMP2) /. First[rdSolT]] === rdKap rdKU/12]
rdSolS2 = First[Solve[{(rdV[1, 0] /. rdBe -> 0) == 0, (D[rdV[1, 0], rdBe] /. rdBe -> 0) == 0}, {rdLam, rdU1}]];
rdUb = rdUq /. rdSolS2; rdR2 = rdRHS[rdUb, -rdUb, -rdUb - D[rdUb, rdBe]/4, rdLam /. rdSolS2];
check["radion on S^4 with U: Lambda + kappa U0 == 6/r0^2, kappa U1 == -12/r0^2 and m_r^2 == kappa K_U/12 - 2/r0^2 from V(beta); the eight-dimensional radion equation gives the same static solution and the same mass",
  Simplify[(rdLam /. rdSolS2) + rdKap rdU0 - 6/rdR0^2] === 0 && Simplify[rdKap (rdU1 /. rdSolS2) + 12/rdR0^2] === 0 &&
  Simplify[(D[rdV[1, 0], {rdBe, 2}] /. rdBe -> 0)/(12 rdMP2) /. rdSolS2] === Simplify[rdKap rdKU/12 - 2/rdR0^2] &&
  Simplify[rdR2 /. rdBe -> 0] === 0 && Simplify[-(D[rdR2, rdBe] /. rdBe -> 0) - (rdKap rdKU/12 - 2/rdR0^2)] === 0]
rdVL = rdV[1, 0] /. {rdU0 -> 0, rdU1 -> 0, rdKU -> 0}; rdLamL = rdLam /. First[Solve[(D[rdVL, rdBe] /. rdBe -> 0) == 0, rdLam]];
rdHE2 = Simplify[(rdVL /. rdBe -> 0 /. rdLam -> rdLamL)/(3 rdMP2)]; rdM2L = Simplify[(D[rdVL, {rdBe, 2}] /. rdBe -> 0 /. rdLam -> rdLamL)/(12 rdMP2)];
rdDS = {rdE00, rdE11, rdE44} /. {Ha -> 1/rdR0, Hb -> 0, dHa -> 0, dHb -> 0, rdb -> rdR0};
check["S^4 with Lambda only: dS4 x S^4 with r0^2 = 9/Lambda, H^2 = 1/r0^2 = Lambda/9 solves the 00, 11 and 44 equations; it is a MAXIMUM of V: m^2 = -6 H^2 from V(beta), and m^2 = -6/r0^2 = -6 H^2 from the eight-dimensional radion equation",
  Simplify[rdLamL rdR0^2] === 9 && Simplify[rdHE2 rdR0^2] === 1 && Simplify[rdM2L/rdHE2] === -6 &&
  Simplify[rdDS - 9/rdR0^2] === {0, 0, 0} && Simplify[rdRHS[0, 0, 0, 9/rdR0^2] /. rdBe -> 0] === 0 &&
  Simplify[-(D[rdRHS[0, 0, 0, 9/rdR0^2], rdBe] /. rdBe -> 0) rdR0^2] === -6]
rdVF = rdV[1, rdFx] /. {rdU0 -> 0, rdU1 -> 0, rdKU -> 0}; rdSolF = First[Solve[{(rdVF /. rdBe -> 0) == 0, (D[rdVF, rdBe] /. rdBe -> 0) == 0}, {rdLam, rdFx}]];
rdStF = feStress[Function[{aa, rIIIbb}, aa^3 rIIIbb^4 rdQ2/(2 rIIIbb^8)]];
rdRhoF = (rdFx /. rdSolF)/rdKap E^(-8 rdBe);
check["S^4 with a 4-form flux (illustration; W = a^3 b^4 q^2/(2 b^8) gives eps = rho_F, pa = -rho_F, pb = +rho_F): Minkowski x S^4 needs Lambda = 9/(2 r0^2) and kappa rho_F = 3/(2 r0^2); m^2 = 6/r0^2 from V(beta) and from the eight-dimensional radion equation",
  Simplify[rdStF - rdQ2/(2 bS^8) {1, -1, 1}] === {0, 0, 0} && Simplify[(rdLam /. rdSolF) rdR0^2] === 9/2 && Simplify[(rdFx /. rdSolF) rdR0^2] === 3/2 &&
  Simplify[(D[rdVF, {rdBe, 2}] /. rdBe -> 0 /. rdSolF) rdR0^2/(12 rdMP2)] === 6 &&
  Simplify[rdRHS[rdRhoF, -rdRhoF, rdRhoF, rdLam /. rdSolF] /. rdBe -> 0] === 0 && Simplify[-(D[rdRHS[rdRhoF, -rdRhoF, rdRhoF, rdLam /. rdSolF], rdBe] /. rdBe -> 0) rdR0^2] === 6]
rdGdS = Join[{rdLp[x0]^2 E^(-4 rdBt[x0])}, -rdSa[x0]^2 E^(-4 rdBt[x0]) {1, 1, 1}, -rdR0^2 E^(2 rdBt[x0]) rdGam4];
rdSqG = rdLp[x0] rdSa[x0]^3 E^(-4 rdBt[x0]) rdR0^4 Sin[x4]^3 Sin[x5]^2 Sin[x6];
check["the note's Goldberger-Wise potential in eight dimensions: for a static profile Phi(y) (here Phi(x4, x5)) Sqrt|G| (1/2 G^AB dPhi dPhi - 1/2 m^2 Phi^2) == -N a^3 e^(-4 beta) Sqrt[gamma] (1/2 e^(-2 beta) |grad Phi|^2_gamma + 1/2 m^2 Phi^2): the Einstein-frame factor is e^(-4 beta) (e^(-4 beta phi) in the note's notation), not e^(+4 beta phi)",
  Simplify[rdSqG (1/2 Sum[D[rdPhi[x4, x5], X[[i]]]^2/rdGdS[[i]], {i, 5, 8}] - 1/2 rdMm^2 rdPhi[x4, x5]^2) +
     rdLp[x0] rdSa[x0]^3 E^(-4 rdBt[x0]) rdR0^4 Sin[x4]^3 Sin[x5]^2 Sin[x6] (1/2 E^(-2 rdBt[x0]) Sum[D[rdPhi[x4, x5], X[[i]]]^2/(rdR0^2 rdGam4[[i - 4]]), {i, 5, 8}] + 1/2 rdMm^2 rdPhi[x4, x5]^2)] === 0]
{"sphere with U: Lambda + kappa U0, kappa U1, m_r^2 - kappa K_U/12 (units 1/r0^2)" -> Simplify[{((rdLam /. rdSolS2) + rdKap rdU0) rdR0^2, rdKap (rdU1 /. rdSolS2) rdR0^2,
    (((D[rdV[1, 0], {rdBe, 2}] /. rdBe -> 0)/(12 rdMP2) /. rdSolS2) - rdKap rdKU/12) rdR0^2}],
 "sphere with Lambda only: Lambda r0^2, H^2 r0^2, m^2/H^2" -> {rdLamL rdR0^2, rdHE2 rdR0^2, rdM2L/rdHE2} // Simplify,
 "sphere with a 4-form flux: Lambda r0^2, kappa rho_F r0^2, m^2 r0^2" -> N[Simplify[{(rdLam /. rdSolF) rdR0^2, (rdFx /. rdSolF) rdR0^2, (D[rdVF, {rdBe, 2}] /. rdBe -> 0 /. rdSolF) rdR0^2/(12 rdMP2)}]]}
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('26.3 The radion meets the world brane: coupling, fifth force, Eöt-Wash')
text('At the brane the induced metric is h = e^(−4β) g_E, on T⁴ and on S⁴ alike (the factor depends only on n = 4). A particle of mass M at '
     'rest on the brane has the action −M ∫ dτ √(h_μν ẋ^μ ẋ^ν) = −M ∫ dx0 e^(−2β) = −M ∫ dx0 (1 − 2β + …): with β = φ/(√12 M_Pl) the '
     'coupling is (φ/(√3 M_Pl)) T^μ_μ (T^μ_μ = M for the particle), the coupling of section 20.3. For n extra dimensions h = e^(−nβ) g_E and '
     'the coupling is √(n/(2(n+2)))/M_Pl. A massless radion then adds to Newton’s potential a force of the same range with the strength '
     'α = 2 (coupling × M_Pl)² = n/(n+2) = 2/3 — the FIFTH FORCE of the note (sections 25.4 and the factor 5/3 = 1 + 2/3 of the image sum). '
     'The Eöt-Wash torsion balance (J. G. Lee et al., Phys. Rev. Lett. 124 (2020) 101101) excludes gravitational-strength (|α| = 1) Yukawa '
     'forces with a range λ > 38.6 μm (95% CL); a radion of strength 2/3 must therefore be heavier than about ħc/38.6 μm ≈ 5 meV. The 10 GeV '
     'radion of sections 22–24 has the range ħc/(10 GeV) = 2.0 × 10⁻¹⁷ m. In warped models the radion couples with 1/⟨μ⟩ (TeV) and mixes with '
     'the Higgs boson through the parameter γ = v/⟨μ⟩ (Giudice, Rattazzi and Wells); in 7+1 the scale is √3 M_Pl, so γ = v/(√3 M_Pl).')
inp(r'''
rdPP = Series[-rdMs E^(-rdNn rdBe/2), {rdBe, 0, 1}] // Normal;
check["the radion couples to a mass M on the brane as (phi/(Sqrt[3] M_Pl)) M: -M e^(-n beta/2) = -M + (n/2) M beta + ..., with beta = phi/(M_Pl Sqrt[n (n + 2)/2]) the coupling is Sqrt[n/(2 (n + 2))]/M_Pl = 1/(Sqrt[3] M_Pl) for n = 4, and the fifth-force strength 2 (coupling M_Pl)^2 = n/(n + 2) = 2/3",
  Simplify[rdPP - (-rdMs + rdNn rdMs rdBe/2)] === 0 && Simplify[((rdNn/2)/Sqrt[rdNn (rdNn + 2)/2])^2 - rdNn/(2 (rdNn + 2)), rdNn > 0] === 0 &&
  ((4/2)/Sqrt[4 (4 + 2)/2])^2 == 1/3 && Simplify[2 rdNn/(2 (rdNn + 2)) - rdNn/(rdNn + 2)] === 0 && 4/(4 + 2) == 2/3]
rdEW = {tbHbarCm/38.6*^-6 10^9, tbHbarCm/10., vEW/(Sqrt[3] rIIMPl)};
check["Eot-Wash (|alpha| = 1 excluded for lambda > 38.6 micrometres): a radion of strength 2/3 must be heavier than about hbar c/38.6 micrometres = 5.1 meV; the 10 GeV radion has the range 2.0 10^-17 m; Higgs-radion mixing in 7+1 is governed by gamma = v/(Sqrt[3] M_Pl) = 5.8 10^-17",
  Abs[rdEW[[1]]/5.112*^-3 - 1] < 10^-3 && Abs[rdEW[[2]]/1.9733*^-17 - 1] < 10^-4 && Abs[rdEW[[3]]/5.838*^-17 - 1] < 10^-3]
{"Eot-Wash bound on the radion mass (eV), range of a 10 GeV radion (m), gamma = v/(Sqrt[3] M_Pl)" -> rdEW}
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('26.4 The 4-sphere and the 4-torus side by side: levels, the brane, dirac16complex')
text('Kaluza–Klein levels of a scalar. T⁴: μ = 2π|n|/ℓ⊥ with the multiplicity r_4(k) (section 25.3). S⁴ of radius r: the eigenvalues of '
     '−Δ are l(l+3)/r², l = 0, 1, 2, …, with the multiplicity d_l = (2l+3)(l+1)(l+2)/6 (the harmonic polynomials of degree l in five '
     'variables: C(l+4, 4) − C(l+2, 4)). Both obey Weyl’s law N(μ) ≈ V μ⁴/(32π²); for S⁴ this is (μr)⁴/12. The brane at the north pole x4 = 0 '
     'sees, at each level l, only the ZONAL harmonic Z_l = √(d_l/V) C_l^(3/2)(cos x4)/C_l^(3/2)(1) (C the Gegenbauer polynomial, '
     'C_l^(3/2)(1) = (l+1)(l+2)/2); all other harmonics of the level vanish there. So the brane couples to ONE state per level, with the '
     'strength √(d_l/V) — on T⁴ the brane couples to the combination Σ_{|n|²=k} φ_n of each level, with the strength √(r_4(k)/V). The '
     'symmetry: S⁴ has the 10 Killing vectors of SO(5); the brane keeps the 6 of SO(4) (the rotations that fix the north pole) and breaks 4 — on '
     'T⁴ the brane keeps the 384 signed permutations and breaks the 4 translations.')
text('dirac16complex. On T⁴ with periodic spinors the Dirac operator of the extra space has zero modes: four massless brane Dirac fermions for '
     'm = 0 (section 25). On S⁴ it has none: Killing spinors ε± (∇_i ε± = ±(i/(2r)) γ̃_i ε±, γ̃_a = i ρ_(a+3) Hermitian, γ̃_aγ̃_b + γ̃_bγ̃_a = 2δ_ab) '
     'have D̸ε± = ±(2i/r) ε±, the lowest eigenvalue allowed by Friedrich’s bound λ² ≥ n R/(4(n−1)) = 4/r²; the spectrum is ±i(l+2)/r with '
     'the degeneracy 4·C(l+3, 3) for each sign (R. Camporesi and A. Higuchi, J. Geom. Phys. 20 (1996) 1). The Kaluza–Klein masses of '
     'dirac16complex are √(m² + λ²) on both spaces; the lightest level on S⁴ has λ = 2/r: eight Dirac fermions of mass √(m² + 4/r²). With '
     'the same volume V⊥ = M_Pl²/M8⁶ (the same four-dimensional Newton constant) the sphere has the radius r = (3V⊥/(8π²))^(1/4).')
inp(r'''
rdDl[l_] := (2 l + 3) (l + 1) (l + 2)/6;
rdZon[l_] := GegenbauerC[l, 3/2, Cos[x4]];
check["zonal harmonics on S^4 (l = 0..6, exact): (1/sin^3 x4) d/dx4 (sin^3 x4 d/dx4 C_l^(3/2)(cos x4)) == -l (l + 3) C_l^(3/2), C_l^(3/2)(1) == (l + 1)(l + 2)/2, and Int over S^4 of (C_l/C_l(1))^2 == V/d_l (unit sphere), so the normalized zonal harmonic is Sqrt[d_l/V] at the north pole; d_l == Binomial[l + 4, 4] - Binomial[l + 2, 4] for l = 0..20",
  And @@ Table[Simplify[D[Sin[x4]^3 D[rdZon[l], x4], x4]/Sin[x4]^3 + l (l + 3) rdZon[l]] === 0 && GegenbauerC[l, 3/2, 1] == (l + 1) (l + 2)/2 &&
     2 Pi^2 Integrate[Expand[(rdZon[l]/GegenbauerC[l, 3/2, 1])^2 Sin[x4]^3], {x4, 0, Pi}] == (8 Pi^2/3)/rdDl[l], {l, 0, 6}] &&
  And @@ Table[rdDl[l] == Binomial[l + 4, 4] - Binomial[l + 2, 4], {l, 0, 20}]]
rdNS4[x_] := Total[rdDl /@ Select[Range[0, Ceiling[x] + 2], # (# + 3) <= x^2 &]];
rdWeyl = Table[N[rdNS4[x]/(x^4/12) - 1], {x, {10, 100, 1000}}];
check["Weyl's law on S^4: the number of levels with l(l + 3) <= (mu r)^2 approaches (mu r)^4/12 = V mu^4/(32 Pi^2) (V = 8 Pi^2 r^4/3): relative deviation -1.0 10^-2, -1.0 10^-4, -1.0 10^-6 at mu r = 10, 100, 1000 (like -1/(mu r)^2)",
  Abs[rdWeyl[[1]] + 1/100] < 10^-6 && Abs[rdWeyl[[2]] + 1/10^4] < 10^-8 && Abs[rdWeyl[[3]] + 1/10^6] < 10^-10 && Simplify[(8 Pi^2/3) rdR0^4/(32 Pi^2) - rdR0^4/12] === 0]
rdEmb = rdR0 {Sin[x4] Sin[x5] Sin[x6] Cos[x7], Sin[x4] Sin[x5] Sin[x6] Sin[x7], Sin[x4] Sin[x5] Cos[x6], Sin[x4] Cos[x5], Cos[x4]};
rdJe = Table[D[rdEmb[[a]], X[[i]]], {a, 5}, {i, 5, 8}];
check["the brane at the north pole of S^4: the embedding (r0 sin x4 sin x5 sin x6 cos x7, ..., r0 cos x4) of S^4 in R^5 induces r0^2 gamma; the north pole x4 = 0 is (0, 0, 0, 0, r0); of the 10 rotations L_ab of R^5 (the Killing vectors of S^4, SO(5)) the 6 with a, b <= 4 vanish there (SO(4), kept by the brane) and the 4 L_a5 move it",
  Simplify[Transpose[rdJe].rdJe - rdR0^2 DiagonalMatrix[rdGam4]] === ConstantArray[0, {4, 4}] && (rdEmb /. x4 -> 0) === {0, 0, 0, 0, rdR0} &&
  Count[Subsets[Range[5], {2}], {a_, b_} /; (UnitVector[5, 5][[a]] UnitVector[5, b] - UnitVector[5, 5][[b]] UnitVector[5, a]) == {0, 0, 0, 0, 0}] == 6 &&
  Length[Subsets[Range[5], {2}]] == 10]
rdGt = I rhoK; rdHf = rdR0 {1, Sin[x4], Sin[x4] Sin[x5], Sin[x4] Sin[x5] Sin[x6]};
rdOm[mu_, i_, j_] := If[mu == i, D[rdHf[[i]], X[[4 + j]]]/rdHf[[j]], 0] - If[mu == j, D[rdHf[[j]], X[[4 + i]]]/rdHf[[i]], 0];
rdNab[M_, mu_] := D[M, X[[4 + mu]]] + Sum[If[i != j, rdOm[mu, i, j] rdGt[[i]].rdGt[[j]]/4, 0], {i, 4}, {j, 4}].M;
rdKS[s_] := (Cos[x4/2] ID4 + s I Sin[x4/2] rdGt[[1]]).(Cos[x5/2] ID4 + Sin[x5/2] rdGt[[1]].rdGt[[2]]).(Cos[x6/2] ID4 + Sin[x6/2] rdGt[[2]].rdGt[[3]]).
   (Cos[x7/2] ID4 + Sin[x7/2] rdGt[[3]].rdGt[[4]]);
rdZeroQ[m_] := Simplify[TrigExpand[m /. {x4 -> 2 rdv4, x5 -> 2 rdv5, x6 -> 2 rdv6, x7 -> 2 rdv7}]] === ConstantArray[0, Dimensions[m]];      (* half angles *)
check["Killing spinors of S^4 (radius r0, frame e^a = (r0, r0 sin x4, ...) dx^(a+3), gt_a = I rho_(a+3)): eps_s = exp(s (i/2) x4 gt_1) exp((1/2) x5 gt_1 gt_2) exp((1/2) x6 gt_2 gt_3) exp((1/2) x7 gt_3 gt_4) eps0 satisfies nabla_mu eps_s == s (i/(2 r0)) e_mu^a gt_a eps_s (s = +-1, all eps0), so the Dirac operator gives (2 s i/r0) eps_s: |lambda| = 2/r0, Friedrich's bound n R/(4(n - 1)) = 4/r0^2",
  And @@ Table[rdGt[[a]] == ConjugateTranspose[rdGt[[a]]], {a, 4}] && And @@ Flatten[Table[rdGt[[a]].rdGt[[b]] + rdGt[[b]].rdGt[[a]] == 2 KroneckerDelta[a, b] ID4, {a, 4}, {b, 4}]] &&
  And @@ Flatten[Table[rdZeroQ[rdNab[rdKS[s], mu] - s I/(2 rdR0) rdHf[[mu]] rdGt[[mu]].rdKS[s]], {s, {1, -1}}, {mu, 4}]] &&
  And @@ Table[rdZeroQ[Sum[rdGt[[a]].rdNab[rdKS[s], a]/rdHf[[a]], {a, 4}] - s 2 I/rdR0 rdKS[s]], {s, {1, -1}}] &&
  4 (12/rdR0^2)/(4 (4 - 1)) == (2/rdR0)^2]
rdSpW = N[Total[Table[2 4 Binomial[l + 3, 3], {l, 0, 998}]]/(1000^4/3) - 1];
check["Dirac spectrum on S^4: Weyl's law for spinors (4 components) N(|lambda| <= L/r) ~ 4 V (L/r)^4/(32 Pi^2) = L^4/3 agrees with the degeneracies 2 x 4 Binomial[l + 3, 3] of the eigenvalues +-i(l + 2)/r (to 2.0 10^-3 at L = 1000)",
  Abs[rdSpW] < 5*^-3 && Simplify[4 (8 Pi^2/3) rdR0^4/(32 Pi^2) - rdR0^4/3] === 0]
rdRS4 = (3 rIIVp/(8 Pi^2))^(1/4);
rdToday = {rdRS4, rdRS4 tbHbarCm, 2 Pi/tbLp 10^3, 2/rdRS4 10^3, 2/rdRS4 10^3, 1 + Total[rIIr4], rdNS4[5 rdRS4], tbNmu[5.]};
check["the sphere with the volume of the torus today (V = M_Pl^2/M8^6): r = (3 V/(8 Pi^2))^(1/4) = 688.9 GeV^-1 = 1.359 10^-13 m; lowest scalar levels 4.03 MeV (torus, 8 states) and 2.90 MeV (sphere, l = 1, 5 states); lowest dirac16complex level on the sphere 2/r = 2.90 MeV (8 Dirac fermions; the torus has 4 massless ones for m = 0); levels below 5 GeV: torus 1.17335 10^13, sphere 1.17375 10^13, Weyl 1.17335 10^13",
  Abs[rdToday[[1]]/688.94 - 1] < 10^-4 && Abs[rdToday[[4]]/2.903 - 1] < 10^-3 && Abs[rdToday[[7]]/rdToday[[8]] - 1] < 10^-3 && Abs[rdToday[[6]]/rdToday[[8]] - 1] < 10^-4 &&
  rdDl[1] == 5 && 2 4 Binomial[3, 3] == 8]
{"S^4 multiplicities d_l, l = 0..8" -> Table[rdDl[l], {l, 0, 8}], "S^4 Weyl: N(x)/(x^4/12) - 1 at x = mu r = 10, 100, 1000" -> rdWeyl,
 "Dirac on S^4: degeneracy 4 Binomial[l + 3, 3] of +i(l + 2)/r, l = 0..5; Weyl ratio - 1 at L = 1000" -> Join[Table[4 Binomial[l + 3, 3], {l, 0, 5}], {rdSpW}],
 "today: r_S4 (GeV^-1, m), lowest masses (MeV): torus 2 Pi/l_perp, sphere scalar 2/r, sphere dirac16complex 2/r; levels below 5 GeV: torus lattice, sphere exact, Weyl" -> rdToday}
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('26.5 Flat versus warped: the note’s Randall–Sundrum example, checked')
text('The note’s warped radion is the five-dimensional model of Randall and Sundrum (RS1): ds² = e^(−2A(y)) η_μν dx^μ dx^ν − dy², '
     'y ∈ [−πR, πR] with y ≃ −y, branes at y = 0 (UV) and y = πR (IR). With the conventions of 7+1.nb (G^A_B − Λ δ^A_B = κ T^A_B) the '
     'Einstein tensor is G^μ_ν = 3A″ − 6A′², G^y_y = −6A′². In the bulk A = k|y| gives Λ = −6k² (anti-de Sitter); at the branes the slope A′ '
     'jumps by +2k (UV: from −k to +k) and by −2k (IR), so A″ = ±2k δ, '
     'and a brane tension f (T^μ_ν = f δ(y) δ^μ_ν) requires f_UV = +6k/κ and f_IR = −6k/κ — the UV brane has POSITIVE tension, as the note '
     'states first (λ_A = −λ_B = 6k/κ8²); its later Israel-junction computation reverses both signs. The four-dimensional Planck mass is '
     'M_Pl² = (1 − e^(−2πkR))/(κk). Promoting R to a field T(x) (ds² = e^(−2kT|y|) g − T² dy²), the bulk action gives the kinetic term '
     '(3/(κk)) (∂ e^(−kπT))², so the canonical radion is μ = f e^(−kπT) with f = √(6/(κk)). Here 1/(2κ) is the coefficient of √|G| R in the '
     'bulk action, M_Pl²/2 that of √|g| R₄ and ½ that of (∂μ)². Randall and Sundrum and Goldberger and Wise write the bulk action as '
     '2M³∫√|G| R, i.e. κ = 1/(4M³): then f = √(24M³/k) and the tensions ±24M³k. The note writes S5 = M5³∫√|G| R5 (κ = 1/(2M5³)) together '
     'with the formulas f = √(24M5³/k) and M_Pl² = (M5³/k)(1 − e^(−2πk⟨T⟩)) of those papers; in its own normalization they read '
     'f = √(12M5³/k) and M_Pl² = (2M5³/k)(1 − e^(−2πk⟨T⟩)). Its eight-dimensional version has the right M_Pl² but f = √(12/(κk)) instead '
     'of √(6/(κk)). On the IR brane h = (μ/f)² g: the radion couples to brane matter with 1/⟨μ⟩, TeV-strong, instead of 1/(√3 M_Pl). In '
     'EIGHT dimensions the note’s equations do not hold: for ds² = e^(−2A(y)) η − dy² − dz⃗² (warping along one of four extra directions) '
     'G^y_y = −6A′² but G^z_z = 4A″ − 10A′², and with only Λ there is no warped solution (A′ = 0). 7+1.nb is unwarped: its bulk is a product.')
text('The note’s warped stabilization (Goldberger and Wise, in five dimensions): a bulk scalar Φ(y) of mass m with the boundary values '
     'Φ(0) = v_A, Φ(πT) = v_B obeys Φ″ − 4kΦ′ − m²Φ = 0 between the branes (from the action ∫√|G| (½ G^AB ∂Φ∂Φ − ½ m²Φ²), √|G| = '
     'e^(−4ky)), solved by Φ = e^(2ky) (C1 e^(νky) + C2 e^(−νky)) with ν = √(4 + m²/k²) ≈ 2 + ε, ε = m²/(4k²). Goldberger and Wise find '
     'V ≈ kεv_A² + 4k e^(−4kπT) (v_B − v_A e^(−εkπT))², i.e. V(μ) = 4k (μ/f)⁴ [v_B − v_A (μ/f)^ε]² + const. Its minimum is at '
     '(⟨μ⟩/f)^ε = v_B/v_A, i.e. kπ⟨T⟩ = (4k²/m²) ln(v_A/v_B), and its curvature there is V″(⟨μ⟩) = 8kε²v_B²⟨μ⟩²/f⁴ (the note’s algebra is '
     'right). With f² = 6/(κk) the radion mass is m_r² = (4/3) κ k² ε² v_B² (⟨μ⟩/f)² = k² v_B² ε² e^(−2πk⟨T⟩)/(3M³) (κ = 1/(4M³): the '
     'result of Goldberger and Wise, Phys. Lett. B 475 (2000) 275, whose footnote adds that a more careful treatment of the ε terms changes '
     'ε² into ε^(3/2)); the note’s factor 2/3 instead of 4/3 comes from its f² = 12/(κk).')
inp(r'''
rdXs5 = {rdt, rdx1, rdx2, rdx3, rdy}; rdG5 = rdEin[{E^(-2 rdA[rdy]), -E^(-2 rdA[rdy]), -E^(-2 rdA[rdy]), -E^(-2 rdA[rdy]), -1}, rdXs5];
rdc5 = Coefficient[rdG5[[1]], rdA''[rdy]]; rdJmp = D[rdk rdy, rdy] - D[-rdk rdy, rdy];
rdfUV = rdf /. First[Solve[rdc5 rdJmp == rdKap rdf, rdf]]; rdfIR = rdf /. First[Solve[rdc5 (-rdJmp) == rdKap rdf, rdf]];
check["RS1 in five dimensions (mostly minus): G^mu_nu = 3 A'' - 6 A'^2, G^y_y = -6 A'^2; bulk A = k y needs Lambda = -6 k^2; the delta functions A'' = +-2 k delta at the branes need the tensions f_UV = +6 k/kappa, f_IR = -6 k/kappa (UV positive)",
  Simplify[rdG5 - Join[Table[3 rdA''[rdy] - 6 rdA'[rdy]^2, {4}], {-6 rdA'[rdy]^2}]] === ConstantArray[0, 5] &&
  Simplify[(rdG5 /. {rdA'[rdy] -> rdk, rdA''[rdy] -> 0}) - (-6 rdk^2)] === ConstantArray[0, 5] &&
  rdc5 === 3 && rdJmp === 2 rdk && rdfUV === 6 rdk/rdKap && rdfIR === -6 rdk/rdKap]
rdKmix = Table[(1/2) D[{E^(-2 rdA[rdy]), -E^(-2 rdA[rdy]), -E^(-2 rdA[rdy]), -E^(-2 rdA[rdy])}[[i]], rdy]/{E^(-2 rdA[rdy]), -E^(-2 rdA[rdy]), -E^(-2 rdA[rdy]), -E^(-2 rdA[rdy])}[[i]], {i, 4}];
rdKjump = Simplify[(rdKmix /. rdA'[rdy] -> rdk) - (rdKmix /. rdA'[rdy] -> -rdk)];
check["Israel junction at the UV brane of RS1: with K^mu_nu = (1/2) g^mu alpha d_y g_alpha nu = -A' delta^mu_nu (normal d_y), [K^mu_nu] - delta^mu_nu [K] = 6 k = (coefficient of A'' in G^mu_nu) x (jump of A'), the delta-function part of G^mu_nu; in the signature (-,+,+,+,+) of the note G^A_B changes sign, so there [K_mu nu - g_mu nu K] = -kappa S_mu nu (the note uses +kappa S_mu nu, hence its reversed tensions)",
  Simplify[rdKmix + rdA'[rdy]] === {0, 0, 0, 0} && Simplify[rdKjump - Total[rdKjump]] === Table[rdc5 rdJmp, {4}] && rdc5 rdJmp === 6 rdk]
check["RS1: M_Pl^2 = (1/kappa) Int_(-Pi R)^(Pi R) e^(-2 k |y|) dy = (1 - e^(-2 Pi k R))/(kappa k)",
  Simplify[Integrate[E^(-2 rdk Abs[rdy]), {rdy, -Pi rdRR, Pi rdRR}, Assumptions -> rdk > 0 && rdRR > 0]/rdKap - (1 - E^(-2 Pi rdk rdRR))/(rdKap rdk)] === 0]
rdR5T = rdScal[{E^(-2 rdk rdT[rdt] rdy), -E^(-2 rdk rdT[rdt] rdy), -E^(-2 rdk rdT[rdt] rdy), -E^(-2 rdk rdT[rdt] rdy), -rdT[rdt]^2}, rdXs5];
rdLb = 2 Integrate[Expand[(-(1/(2 rdKap)) E^(-4 rdk rdTs rdy) rdTs rdR5T) /. {rdT''[rdt] -> rdT2, rdT'[rdt] -> rdT1, rdT[rdt] -> rdTs}], {rdy, 0, Pi},
    Assumptions -> rdk > 0 && rdTs > 0];
rdC2 = Coefficient[Expand[rdLb], rdT2]; rdL1 = Expand[rdLb - rdC2 rdT2 - D[rdC2, rdTs] rdT1^2];
check["RS1 radion: with R promoted to T(x), the bulk action is (3/(kappa k)) (d e^(-k Pi T))^2 + V(T) + (a total derivative), so the canonical radion is mu = Sqrt[6/(kappa k)] e^(-k Pi T) (the note: Sqrt[12/(kappa k)])",
  FreeQ[rdC2, rdT1 | rdT2] && Simplify[Coefficient[rdL1, rdT1, 2] - (3/(rdKap rdk)) (rdk Pi E^(-rdk Pi rdTs))^2] === 0 && Coefficient[rdL1, rdT1, 1] === 0]
rdG8w = rdEin[{E^(-2 rdA[rdy]), -E^(-2 rdA[rdy]), -E^(-2 rdA[rdy]), -E^(-2 rdA[rdy]), -1, -1, -1, -1}, {rdt, rdx1, rdx2, rdx3, rdy, rdz1, rdz2, rdz3}];
check["eight dimensions, warped along one extra direction (ds^2 = e^(-2 A(y)) eta - dy^2 - dz.dz): G^mu_nu = 3 A'' - 6 A'^2, G^y_y = -6 A'^2, G^z_z = 4 A'' - 10 A'^2; with only Lambda (all G^A_A equal) A'' = 0 and A' = 0: the note's five-dimensional equations do not carry over",
  Simplify[rdG8w - Join[Table[3 rdA''[rdy] - 6 rdA'[rdy]^2, {4}], {-6 rdA'[rdy]^2}, Table[4 rdA''[rdy] - 10 rdA'[rdy]^2, {3}]]] === ConstantArray[0, 8] &&
  DeleteDuplicates[Solve[{3 rdA2 - 6 rdA1^2 == -6 rdA1^2, -6 rdA1^2 == 4 rdA2 - 10 rdA1^2}, {rdA1, rdA2}]] === {{rdA1 -> 0, rdA2 -> 0}}]
rdGWL = E^(-4 rdk rdy) (-(1/2) rdPhy'[rdy]^2 - (1/2) rdMg^2 rdPhy[rdy]^2);
rdGWsol = E^(2 rdk rdy) (rdC1 E^(rdNu rdk rdy) + rdC2 E^(-rdNu rdk rdy)) /. rdNu -> Sqrt[4 + rdMg^2/rdk^2];
check["Goldberger-Wise in RS1: the static bulk scalar (Sqrt|G| = e^(-4 k y), G^yy = -1) obeys Phi'' - 4 k Phi' - m^2 Phi = 0, solved by e^(2 k y)(C1 e^(nu k y) + C2 e^(-nu k y)), nu = Sqrt[4 + m^2/k^2] = 2 + eps + O(m^4), eps = m^2/(4 k^2)",
  Simplify[D[rdGWL, rdPhy[rdy]] - D[D[rdGWL, rdPhy'[rdy]], rdy] - E^(-4 rdk rdy) (rdPhy''[rdy] - 4 rdk rdPhy'[rdy] - rdMg^2 rdPhy[rdy])] === 0 &&
  Simplify[D[rdGWsol, {rdy, 2}] - 4 rdk D[rdGWsol, rdy] - rdMg^2 rdGWsol, rdk > 0] === 0 &&
  Normal[Series[Sqrt[4 + rdMg^2/rdk^2], {rdMg, 0, 2}]] === 2 + rdMg^2/(4 rdk^2)]
rdVgw = 4 rdk (rdmu/rdfm)^4 (rdvB - rdvA (rdmu/rdfm)^rdeps)^2; rdmu0 = rdfm (rdvB/rdvA)^(1/rdeps);      (* all of k, f, v_A, v_B, eps, kappa positive: PowerExpand is exact here *)
rdV2gw = Simplify[PowerExpand[D[rdVgw, {rdmu, 2}] /. rdmu -> rdmu0]];
rdMgw = Simplify[PowerExpand[rdV2gw/(rdk^2 rdeps^2 rdvB^2 (rdmu0/rdfm)^2) /. rdfm -> Sqrt[6/(rdKap rdk)]]/rdKap];
rdMgwNote = Simplify[PowerExpand[rdV2gw/(rdk^2 rdeps^2 rdvB^2 (rdmu0/rdfm)^2) /. rdfm -> Sqrt[12/(rdKap rdk)]]/rdKap];
check["Goldberger-Wise radion mass: V(mu) = 4 k (mu/f)^4 (v_B - v_A (mu/f)^eps)^2 has V'(<mu>) = 0 at (<mu>/f)^eps = v_B/v_A (k Pi <T> = (4 k^2/m^2) Log[v_A/v_B]) and V''(<mu>) = 8 k eps^2 v_B^2 <mu>^2/f^4; with f^2 = 6/(kappa k): m_r^2 = (4/3) kappa k^2 eps^2 v_B^2 (<mu>/f)^2 (the note's f^2 = 12/(kappa k) gives 2/3)",
  Simplify[PowerExpand[D[rdVgw, rdmu] /. rdmu -> rdmu0]] === 0 &&
  Simplify[PowerExpand[rdV2gw - 8 rdk rdeps^2 rdvB^2 rdmu0^2/rdfm^4]] === 0 &&
  Simplify[PowerExpand[-Log[rdmu0/rdfm] - (4 rdk^2/rdMg^2) (Log[rdvA] - Log[rdvB])] /. rdeps -> rdMg^2/(4 rdk^2)] === 0 &&
  rdMgw === 4/3 && rdMgwNote === 2/3]
{"RS1: coefficient of A'' in G^0_0, f_UV kappa/k, f_IR kappa/k, f^2 kappa k (canonical radion)" -> Simplify[{rdc5, rdfUV rdKap/rdk, rdfIR rdKap/rdk,
    2 Coefficient[rdL1, rdT1, 2] rdKap rdk/(rdk Pi E^(-rdk Pi rdTs))^2}],
 "strength of the radion coupling: warped with <mu> = 3 TeV over 7+1, i.e. Sqrt[3] M_Pl/(3 TeV)" -> Sqrt[3] rIIMPl/3000.,
 "Goldberger-Wise: m_r^2/(kappa k^2 eps^2 v_B^2 (<mu>/f)^2) with f^2 = 6/(kappa k) and with the note's f^2 = 12/(kappa k)" -> N[{rdMgw, rdMgwNote}]}
''')

# ------------------------------------------------------------------------------------------------------------------------
subsection('26.6 The Higgs boson on the brane and the radion: mixing (the note’s formulas, checked), and its size in 7+1')
text('The Higgs doublet H lives on the world brane (section 25.4) and feels the induced metric h_μν = Ω² g_E,μν. In 7+1 Ω = e^(−2β) = '
     'e^(−φ/(√3 M_Pl)) on T⁴ and on S⁴ alike (section 26.3); in RS1 Ω = μ/⟨μ⟩ after the brane fields are rescaled so that v is the '
     'physical vacuum value. Write Ω = 1 + r/Λ + O(r²) with the bare radion r (Λ = ⟨μ⟩ in RS1, Λ = −√3 M_Pl in 7+1; the sign of Λ can '
     'be absorbed into r → −r) and γ = v/Λ. In unitary gauge H†H = φ_H²/2 with φ_H = v + h. The note adds the curvature term '
     '−ξ ∫√|h| R(h) H†H of Giudice, Rattazzi and Wells (Nucl. Phys. B 595 (2001) 250), written in the signature (−,+,+,+); in the signature '
     '(+,−,−,−) of 7+1.nb the scalar curvature changes sign and the brane Lagrangian reads')
text('      L = √|h| [ ½ h^μν ∂_μφ_H ∂_νφ_H − (λ/4)(φ_H² − v²)² + ½ ξ R(h) φ_H² ] + ½ (∂r)² − ½ m_r0² r²   (ξ = 1/6: conformal coupling),')
text('with the canonical radion kinetic term of section 26.1 and the stabilized radion mass m_r0 (section 26.2). The checks: (i) for '
     'h = Ω²η, √|h| R(h) = −6 Ω □Ω with □ = ∂0² − ∇² (four dimensions); (ii) to second order in the fields and up to total derivatives '
     'L = ½(∂h)² + ½(1 + 6ξγ²)(∂r)² + 6ξγ ∂h·∂r − ½m_h0² h² − ½m_r0² r² (m_h0² = 2λv²), and the terms linear in the fields are total '
     'derivatives (the vacuum is a solution) — the note’s kinetic Lagrangian lacks the 6ξγ² of the radion, which its Z contains; (iii) the '
     'shift h = h′ − (6ξγ/Z) r′, r = r′/Z makes the kinetic terms canonical if and only if Z² = 1 + 6ξγ²(1 − 6ξ) (with the note’s literal '
     'kinetic terms it would be Z² = 1 − 36ξ²γ²); (iv) the mass matrix of (h′, r′) is diagonal in the fields (h_m, r_m) given by '
     'h′ = cos θ h_m − sin θ r_m, r′ = sin θ h_m + cos θ r_m with tan 2θ = 12γξZ m_h0²/(m_r0² − m_h0²(Z² − 36ξ²γ²)) (the note writes the same '
     'matrix with the primed and unprimed fields interchanged, which would reverse the sign of θ), and then h = C_h h_m + …, r = C_r h_m + … '
     'with C_h = cos θ − (6ξγ/Z) sin θ, C_r = sin θ/Z, the coefficients of the note. The note’s rates Γ(h → γγ)/Γ_SM = |C_h + C_r γ b_EM/F_SM|² '
     'and Γ(h → gg)/Γ_SM = |C_h + C_r γ b_QCD/F_heavy|² (its values b_EM = −11/3, F_SM = −6.5, b_QCD = 7, F_heavy = −4/3) are evaluated '
     'below for 7+1. There γ = v/(√3 M_Pl) = 5.8 × 10⁻¹⁷: for ξ = 1 and m_r0 = 10 GeV the mixing angle is 3.5 × 10⁻¹⁶ and the rates change '
     'by 10⁻³¹ — the flat radion of 7+1 (T⁴ or S⁴) cannot be seen through the Higgs boson; the note’s LHC bounds (⟨μ⟩ > 2.5–3 TeV, '
     '−0.1 < ξ < 0.4) concern the warped radion (γ = v/(3 TeV) = 0.082).')
inp(r'''
rdEL2[L_, q_] := Simplify[D[L, q[x0, x1]] - D[D[L, Derivative[1, 0][q][x0, x1]], x0] - D[D[L, Derivative[0, 1][q][x0, x1]], x1] +
    D[D[L, Derivative[2, 0][q][x0, x1]], {x0, 2}] + D[D[L, Derivative[1, 1][q][x0, x1]], x0, x1] + D[D[L, Derivative[0, 2][q][x0, x1]], {x1, 2}]];
rdc4 = Simplify[rdOmg[x0, x1]^4 rdScal[rdOmg[x0, x1]^2 {1, -1, -1, -1}, {x0, x1, x2, x3}]/(rdOmg[x0, x1] (D[rdOmg[x0, x1], {x0, 2}] - D[rdOmg[x0, x1], {x1, 2}]))];
check["conformally flat brane metric h = Omega^2 eta (four dimensions, Omega(x0, x1)): Sqrt|h| R(h) == -6 Omega box Omega with box = d0^2 - d1^2 (curvature conventions of 7+1.nb)",
  rdc4 === -6]
rdOmE = 1 + rdO1 rdEp rdRr[x0, x1] + rdO2 rdEp^2 rdRr[x0, x1]^2; rdPhE = rdVv + rdEp rdHh[x0, x1];
rdRh = rdScal[rdOmE^2 {1, -1, -1, -1}, {x0, x1, x2, x3}];
rdLbr = rdOmE^4 ((1/2) rdOmE^-2 (D[rdPhE, x0]^2 - D[rdPhE, x1]^2) - rdLamH/4 (rdPhE^2 - rdVv^2)^2 + (1/2) rdXi rdRh rdPhE^2) +
   rdEp^2 ((1/2) (D[rdRr[x0, x1], x0]^2 - D[rdRr[x0, x1], x1]^2) - (1/2) rdMr0^2 rdRr[x0, x1]^2);
rdHL1 = Simplify[D[rdLbr, rdEp] /. rdEp -> 0]; rdHL2 = Simplify[D[rdLbr, {rdEp, 2}]/2 /. rdEp -> 0];
rdGam = rdO1 rdVv;
rdLq = (1/2) (D[rdHh[x0, x1], x0]^2 - D[rdHh[x0, x1], x1]^2) + (1/2) (1 + 6 rdXi rdGam^2) (D[rdRr[x0, x1], x0]^2 - D[rdRr[x0, x1], x1]^2) +
   6 rdXi rdGam (D[rdHh[x0, x1], x0] D[rdRr[x0, x1], x0] - D[rdHh[x0, x1], x1] D[rdRr[x0, x1], x1]) - (1/2) (2 rdLamH rdVv^2) rdHh[x0, x1]^2 - (1/2) rdMr0^2 rdRr[x0, x1]^2;
rdKhr = -Coefficient[Expand[rdEL2[rdHL2, rdRr]], Derivative[2, 0][rdHh][x0, x1]]; rdKrr = -Coefficient[Expand[rdEL2[rdHL2, rdRr]], Derivative[2, 0][rdRr][x0, x1]];
check["the Higgs boson on the brane (h = Omega^2 eta, Omega = 1 + r/Lambda + a2 r^2, phi_H = v + h, gamma = v/Lambda): the terms linear in (h, r) are total derivatives, and to second order L == (1/2)(dh)^2 + (1/2)(1 + 6 xi gamma^2)(dr)^2 + 6 xi gamma dh.dr - (1/2) m_h0^2 h^2 - (1/2) m_r0^2 r^2 up to total derivatives (Euler-Lagrange expressions in x0 and x1); a2 drops out",
  rdEL2[rdHL1, rdRr] === 0 && rdEL2[rdHL1, rdHh] === 0 && rdEL2[rdHL2 - rdLq, rdRr] === 0 && rdEL2[rdHL2 - rdLq, rdHh] === 0 &&
  Simplify[rdKhr - 6 rdXi rdGam] === 0 && Simplify[rdKrr - (1 + 6 rdXi rdGam^2)] === 0 && FreeQ[rdEL2[rdHL2, rdRr], rdO2] && FreeQ[rdEL2[rdHL2, rdHh], rdO2]]
rdKm = {{1, rdKhr}, {rdKhr, rdKrr}}; rdKnote = {{1, 6 rdXi rdGam}, {6 rdXi rdGam, 1}};
rdSm = {{1, -6 rdXi rdGam/rdZ}, {0, 1/rdZ}}; rdZ2 = 1 + 6 rdXi rdGam^2 (1 - 6 rdXi);
check["canonical kinetic terms: Transpose[S].K.S == 1 for S = {{1, -6 xi gamma/Z}, {0, 1/Z}} if and only if Z^2 = 1 + 6 xi gamma^2 (1 - 6 xi) (the note's Z); with the note's literal kinetic terms (no 6 xi gamma^2 for the radion) it would be Z^2 = 1 - 36 xi^2 gamma^2",
  Simplify[Transpose[rdSm].rdKm.rdSm /. rdZ -> Sqrt[rdZ2]] === IdentityMatrix[2] && Simplify[(Transpose[rdSm].rdKm.rdSm)[[2, 2]] rdZ^2 - rdZ2] === 0 &&
  Simplify[(Transpose[rdSm].rdKnote.rdSm)[[2, 2]] rdZ^2 - (1 - 36 rdXi^2 rdGam^2)] === 0 && Simplify[(Transpose[rdSm].rdKnote.rdSm)[[1, 2]]] === 0]
rdMp = Transpose[rdSm].DiagonalMatrix[{rdMh0^2, rdMr0^2}].rdSm;
rdQ[th_] := {{Cos[th], -Sin[th]}, {Sin[th], Cos[th]}};
check["mass mixing: with (h', r') = Q(theta).(h_m, r_m), Q = {{cos, -sin}, {sin, cos}}, the off-diagonal element of Transpose[Q].M'.Q is -(1/2)(M'_11 - M'_22) sin 2 theta + M'_12 cos 2 theta, zero for tan 2 theta = 2 M'_12/(M'_11 - M'_22) = 12 gamma xi Z m_h0^2/(m_r0^2 - m_h0^2 (Z^2 - 36 xi^2 gamma^2)) (the note's formula); with the inverse matrix (the note's display) it is +(1/2)(M'_11 - M'_22) sin 2 theta + M'_12 cos 2 theta (theta -> -theta); (h, r) = S.Q.(h_m, r_m) gives C_h = cos theta - (6 xi gamma/Z) sin theta, C_r = sin theta/Z",
  Simplify[TrigReduce[(Transpose[rdQ[rdth]].rdMp.rdQ[rdth])[[1, 2]]] + (1/2) (rdMp[[1, 1]] - rdMp[[2, 2]]) Sin[2 rdth] - rdMp[[1, 2]] Cos[2 rdth]] === 0 &&
  Simplify[2 rdMp[[1, 2]]/(rdMp[[1, 1]] - rdMp[[2, 2]]) - 12 rdGam rdXi rdZ rdMh0^2/(rdMr0^2 - rdMh0^2 (rdZ^2 - 36 rdXi^2 rdGam^2))] === 0 &&
  Simplify[TrigReduce[(rdQ[rdth].rdMp.Transpose[rdQ[rdth]])[[1, 2]]] - (1/2) (rdMp[[1, 1]] - rdMp[[2, 2]]) Sin[2 rdth] - rdMp[[1, 2]] Cos[2 rdth]] === 0 &&
  Simplify[(rdSm.rdQ[rdth])[[All, 1]] - {Cos[rdth] - 6 rdXi rdGam Sin[rdth]/rdZ, Sin[rdth]/rdZ}] === {0, 0}]
rdHRnum[gam_, xi_, mr_, mh_] := Module[{z, t, ch, cr}, z = Sqrt[1 + 6 xi gam^2 (1 - 6 xi)];
   t = ArcTan[12 gam xi z mh^2/(mr^2 - mh^2 (z^2 - 36 xi^2 gam^2))]/2; ch = Cos[t] - 6 xi gam Sin[t]/z; cr = Sin[t]/z;
   {z - 1, t, (ch + cr gam (-11/3)/(-13/2))^2 - 1, (ch + cr gam 7/(-4/3))^2 - 1}];
rdHR71 = rdHRnum[SetPrecision[rdEW[[3]], 60], 1, 10, SetPrecision[mHb, 60]];
check["Higgs-radion mixing in 7+1 (gamma = v/(Sqrt[3] M_Pl), xi = 1, m_r0 = 10 GeV, m_h0 = 125.20 GeV; 60-digit arithmetic): |Z - 1| < 10^-30, |theta| = 3.5 10^-16, the rates of h -> gamma gamma and h -> g g change by less than 10^-30; in RS1 with <mu> = 3 TeV gamma = v/(3 TeV) = 0.082",
  Abs[rdHR71[[1]]] < 10^-30 && Abs[Abs[rdHR71[[2]]]/3.53*^-16 - 1] < 10^-2 && Abs[rdHR71[[3]]] < 10^-30 && Abs[rdHR71[[4]]] < 10^-30 && Abs[vEW/3000/0.0821 - 1] < 10^-2]
{"Higgs-radion mixing: c in Sqrt|h| R(h) = c Omega box Omega, K_hr/(xi gamma), (K_rr - 1)/(xi gamma^2)" -> N[Simplify[{rdc4, rdKhr/(rdXi rdGam), (rdKrr - 1)/(rdXi rdGam^2)}]],
 "Higgs-radion mixing in 7+1 (gamma = v/(Sqrt[3] M_Pl), xi = 1, m_r0 = 10 GeV): Z - 1, theta, Gamma(h -> gamma gamma)/SM - 1, Gamma(h -> g g)/SM - 1; warped gamma = v/(3 TeV)" ->
   Join[SetPrecision[rdHR71, MachinePrecision], {vEW/3000.}]}
''')
