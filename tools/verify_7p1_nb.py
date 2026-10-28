#!/usr/bin/env python3
"""Independent cross-check of the EVALUATED notebook 7+1.nb against the Python reference (ref71.py).

Usage (from the repository root):  python3 tools/verify_7p1_nb.py [7+1.nb]
Checks: no message cells; every check[...] output is True; Length[checks], failedChecks, allChecksPassed;
the matrices tau[A], tauO[A], Gamma17[A] printed by Mathematica equal the Python ones; gridx44 and gridO tables;
Lambda matrices; the eight real on-shell spinors; the classification table; $Version.
Sections 7-9 (dirac16complex, Kohn-Sham, chiral Kohn-Sham): every check of each section is True; the printed real spinors
u_s(p) are recomputed; the printed closed forms of the uniform-matter integrals are differentiated with sympy; every printed
number (uniform matter in the Hartree and Hartree-Fock approximations, chiral form, NJL vacuum) is compared with the
independent Python computation of tools/ref_dirac16.py.
Sections 11-19 (gravity, Kohn-Sham in gravitational fields, Standard Model): every check of each section is True; the
summary table of the last section lists Length[checks] checks, all True; sections 0-6 and 0-10 equal the backups taken before
they were extended, and the whole notebook equals the backup taken before the signature explanation of 27 October 2026, except
for the cells of that explanation (each changed or added cell states that diag(+1,-1,...,-1) is simply the squares of the
octonion basis, or contains the xAct signature list {7, 1, 0}); no cell of a backup is missing; the printed numbers of
sections 14-19 are compared with tools/ref_gravity.py, tools/ref_ks_gravity.py and tools/ref_sm.py (their outputs are
computed by running them, or read from build/docs/runs.json when it exists).
Sections 20-23 (the coupled field equations and their solutions in the three gravity regimes, appended on 27 October 2026): every
check is True; with sections 20-23 removed and the title of the summary section restored, the notebook equals the backup taken
before they were appended (cells and outputs; the summary counts, the private context name and measured times excepted), and the
older comparisons above are made on that reduced notebook; the printed numbers of sections 20-23 are compared with
tools/ref_coupled.py.
Section 24 (dirac16complex as dark matter and dark energy, appended on 27 October 2026): every check is True; with section 24 removed and
the title of the summary section restored, the notebook equals the backup taken before it was appended (cells and outputs, with the
same exceptions and one documented correction: in the table of section 23.3 the cold-dark-matter density reads "Omega_c = 0.264", where
the backup had "0.265"; 0.120/0.674^2 = 0.26416); the printed numbers of section 24 are compared with tools/ref_dark.py.
"""
import os, re, sys
sys.setrecursionlimit(100000)
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import nbparse as P
import sympy as sp
from ref71 import *

path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(here), '7+1.nb')
nb = P.parse_notebook(path)
cells = P.iter_cells(nb)
outs = []; msgs = []; inputs = 0; check_calls = 0; outsec = []; section = ''
for d, c in cells:
    st, content, opts = P.cell_info(c)
    if st == 'Section': section = P.box_to_text(content).split('.')[0].strip('"')
    if st == 'Input':
        inputs += 1; check_calls += P.box_to_text(content).count('check["')
    if st == 'Output': outs.append(P.box_to_text(content)); outsec.append(section)
    if st == 'Message': msgs.append(P.box_to_text(content))
R = {}
R['input cells'] = inputs
R['output cells'] = len(outs)
R['message cells (must be 0)'] = len(msgs)
chk = [o for o in outs if re.fullmatch(r'"[^"]*(?:"[^"]*)*"->(True|False)', o.strip(), flags=re.S)]
R['check outputs'] = len(chk)
R['check outputs that are True'] = sum(o.strip().endswith('->True') for o in chk)
R['check["..."] calls in the input cells'] = check_calls
def find_out(pred):
    return [o for o in outs if pred(o)]
R['$Version output'] = [o for o in outs if '15.0.1' in o and 'Windows' in o][:1]
R['private notebook context'] = any('Notebook$$' in o for o in outs)
# summary outputs (the last four outputs of section 6 precede the final Grid)
n_checks = [o for o in outs if re.fullmatch(r'\d+', o.strip())]
R['Length[checks]'] = n_checks[-1] if n_checks else None
R['failedChecks == {}'] = '{}' in [o.strip() for o in outs]
R['allChecksPassed output True'] = outs[-2].strip() == 'True' if len(outs) > 1 else None

def balanced(s, i):
    """return substring starting at s[i] == '{' up to the matching '}'"""
    depth = 0
    for j in range(i, len(s)):
        if s[j] == '{': depth += 1
        elif s[j] == '}':
            depth -= 1
            if depth == 0: return s[i:j+1]
    raise ValueError('unbalanced')
def parse_int_matrix(txt):
    rows = re.findall(r'\{([^{}]*)\}', txt)
    return sp.Matrix([[sp.Integer(x) for x in r.split(',')] for r in rows])
def labelled_matrices(label):
    found = {}
    for o in outs:
        for m in re.finditer(r'"' + re.escape(label) + r'\[(\d)\]"', o):
            k = int(m.group(1))
            i = o.index('{{', m.end())
            found[k] = parse_int_matrix(balanced(o, i))
    return found
mt = labelled_matrices('τ'); mo = labelled_matrices('τO'); mg = labelled_matrices('Γ17')
R['tau[A] printed == 4+4 tau[A] (python)'] = len(mt) == 8 and all(mt[A] == tau[A] for A in range(8))
R['tauO[A] printed == python tauO[A]'] = len(mo) == 8 and all(mo[A] == tauO[A] for A in range(8))
R['Gamma17[A] printed == python G17[A]'] = len(mg) == 8 and all(mg[A] == G17[A] for A in range(8))
# tables
def table_text(C, name, base):
    rows = []
    hdr = ['x'] + ['%s[%d]' % (name, k + base) for k in range(8)]
    rows.append(hdr)
    for Bi in range(8):
        row = ['%s[%d]' % (name, Bi + base)]
        for k in range(8):
            p = mul(C, basis(Bi), basis(k)); d = [i for i in range(8) if p[i]][0]
            row.append(('-' if p[d] < 0 else '') + '%s[%d]' % (name, d + base))
        rows.append(row)
    return '{' + ', '.join('{' + ', '.join(r) + '}' for r in rows) + '}'
R['gridO printed == python octonion table'] = table_text(CO, 'e', 0) in [o.strip() for o in outs]
R['gridx44 printed == python split table (4+4.nb)'] = table_text(CS, 'eA', 1) in [o.strip() for o in outs]
# Lambda matrices
lam = [o for o in outs if 'Cosh[β]' in o and 'Cos[θ]' in o]
R['Lambda boost/rotation printed'] = bool(lam) and '{Cosh[β], -Sinh[β], 0, 0, 0, 0, 0, 0}, {-Sinh[β], Cosh[β], 0' in lam[0] and '{0, 0, 0, 0, Cos[θ], -Sin[θ], 0, 0}, {0, 0, 0, 0, Sin[θ], Cos[θ], 0, 0}' in lam[0]
# real on-shell amplitudes (section 5.4): the MatrixForm output with 8 rows of 16 integers
maj = [o for o in outs if o.startswith('(\\[NoBreak]{{') and o.count('},') == 7 and len(re.findall(r'\{([^{}]*)\}', o)[0].split(',')) == 16]
ok = False
if maj:
    W = parse_int_matrix(maj[0])
    p = [3, 2, 2, 0, 0, 0, 0, 0]
    ps = sum((p[A] * G17[A] for A in range(8)), sp.zeros(16))
    ok = W.shape == (8, 16) and W.rank() == 8 and all(((ps - sp.eye(16)) * W.row(i).T) == sp.zeros(16, 1) for i in range(8))
R['8 real on-shell spinors printed solve (pslash - m) w = 0 (python check)'] = ok
cl = [o for o in outs if 'WRONG in [1B]' in o]
R['classification table printed'] = bool(cl) and cl[0].count('WRONG in [1B]') == 3 and cl[0].count('agrees') == 6
# ---------------- sections 7-9: dirac16complex, Kohn-Sham, chiral Kohn-Sham ----------------
chk_re = re.compile(r'"[^"]*(?:"[^"]*)*"->(True|False)', flags=re.S)
for sec in ('7', '8', '9'):
    so = [o for o, q in zip(outs, outsec) if q == sec and chk_re.fullmatch(o.strip())]
    R['section %s: check outputs (all True)' % sec] = '%d (%s)' % (len(so), all(o.strip().endswith('->True') for o in so))
# the eight real spinors u_s(p), p = (2,2,4,0,0,0,0), printed as Sqrt[5] u0Matrix (16 x 8)
um = [o for o, q in zip(outs, outsec) if q == '7' and o.startswith('(\\[NoBreak]{{')]
ok = False
if um:
    U5 = parse_int_matrix(um[0])
    p0 = [2, 2, 4, 0, 0, 0, 0]
    M0 = sp.eye(8) + sum((p0[i-1]*tauO[i] for i in range(1, 8)), sp.zeros(8))
    U5py = sp.Matrix.hstack(*[sp.Matrix.vstack(5*sp.eye(8)[:, s_], M0.T*sp.eye(8)[:, s_]) for s_ in range(8)])
    ps = 5*G17[0] + sum((p0[i-1]*G17[i] for i in range(1, 8)), sp.zeros(16))
    ok = U5 == U5py and (ps - sp.eye(16))*U5 == sp.zeros(16, 8)
R['section 7: printed Sqrt[5] u_s(p) == python (5 e_s; M^T e_s) and (pslash - m) u_s == 0'] = ok
# closed forms of the uniform-matter integrals
from sympy.parsing.mathematica import parse_mathematica
cf = [o for o, q in zip(outs, outsec) if q == '8' and 'ArcTanh' in o and o.startswith('{')]
ok = False
if cf:
    a, b = cf[0][1:-1].split('),((1)/(384))')
    eS = parse_mathematica(a + ')'); eE = parse_mathematica('((1)/(384))' + b)
    ms_, kF_ = sp.Symbol('ms', positive=True), sp.Symbol('kF', positive=True)
    eS = eS.subs({sp.Symbol('ms'): ms_, sp.Symbol('kF'): kF_}); eE = eE.subs({sp.Symbol('ms'): ms_, sp.Symbol('kF'): kF_})
    ok = (sp.simplify(sp.diff(eS, kF_) - kF_**6*ms_/sp.sqrt(kF_**2 + ms_**2)) == 0 and sp.simplify(sp.diff(eE, kF_) - kF_**6*sp.sqrt(kF_**2 + ms_**2)) == 0
          and sp.limit(eS, kF_, 0) == 0 and sp.limit(eE, kF_, 0) == 0)
R['section 8: printed closed forms of int k^6 m*/E and int k^6 E: d/dkF == integrands, zero at kF = 0 (sympy)'] = ok
# numbers versus tools/ref_dirac16.py
num_re = re.compile(r'"([^"]+)"->(-?[0-9]+\.?[0-9]*)(?:`[0-9.]*)?(?:\*\^(-?[0-9]+))?')
def numbers(sec, key):
    for o, q in zip(outs, outsec):
        if q == sec and key in o:
            return {m.group(1): sp.Float(m.group(2) + ('e' + m.group(3) if m.group(3) else ''), 30) for m in num_re.finditer(o)}
    return {}
import io, contextlib
ref = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(os.path.join(here, 'ref_dirac16.py'), encoding='utf-8').read().replace("if __name__ == '__main__':", 'if False:'),
                 'ref_dirac16.py', 'exec'), ref)
RR = ref['R']
def close(x, y, rel=1e-17):
    x, y = sp.Float(x, 30), sp.Float(y, 30)
    return abs(x - y) <= rel*max(abs(y), sp.Float('1e-300'))
hs = numbers('8', '"kF0"'); h1 = numbers('8', '"epsF"->1.98'); hf = numbers('8', '"HF mstar"')
cu = numbers('9', '"kappa = rhoS/2"'); nj = numbers('9', '"g_c"')
pairs_num = [(hs.get('kF0'), RR['NUM Hartree saturation kF0']), (hs.get('rho0'), RR['NUM Hartree rho0']), (hs.get('mstar'), RR['NUM Hartree m*']),
             (hs.get('rhoS'), RR['NUM Hartree rhoS']), (hs.get('V0'), RR['NUM Hartree V0']), (hs.get('E/A - m'), RR['NUM Hartree E/A - m']),
             (hs.get('epsF'), RR['NUM Hartree eps_F']), (h1.get('mstar'), RR['NUM Hartree kF=1: m*']), (h1.get('E/A - m'), RR['NUM Hartree kF=1: E/A - m']),
             (h1.get('epsF'), RR['NUM Hartree kF=1: eps_F']), (hf.get('HF mstar'), RR['NUM HF kF=3/4: m*']), (hf.get('HF V0'), RR['NUM HF kF=3/4: V0']),
             (hf.get('HF E/A - m'), RR['NUM HF kF=3/4: E/A - m']), (hf.get('Hartree mstar'), RR['NUM Hartree kF=3/4: m*']),
             (hf.get('Hartree E/A - m'), RR['NUM Hartree kF=3/4: E/A - m']), (cu.get('kappa = rhoS/2'), RR['NUM Hartree rhoS']/2),
             (cu.get('M = m*'), RR['NUM Hartree m*']), (nj.get('g_c'), RR['NUM NJL g_c']), (nj.get('g'), RR['NUM NJL g = 6/5 g_c']),
             (nj.get('|M| (m = 0)'), RR['NUM NJL |M| (m = 0)']), (nj.get('kappa (m = 0, M real)'), RR['NUM NJL kappa = -4 M J(|M|) (M real)']),
             (nj.get('E(M) - E(0)'), RR['NUM NJL E(|M|) - E(0)']), (nj.get('M (m = 1/100)'), RR['NUM NJL M (m = 1/100)'])]
R['sections 8-9: 23 printed numbers agree with tools/ref_dirac16.py to 17 digits'] = (all(a is not None for a, b in pairs_num)
    and all(close(a, b) for a, b in pairs_num))
R['sections 8-9: G_S == -98750, G_V == 20000 printed'] = hf.get('G_S') is not None and float(hf['G_S']) == -98750.0 and float(hf['G_V']) == 20000.0
R['tools/ref_dirac16.py: all its boolean checks True'] = all(v for v in RR.values() if isinstance(v, bool)) and sum(isinstance(v, bool) for v in RR.values()) >= 80

# sections 0-6 must be exactly the notebook of the first stage: compare every non-output cell with the backup that was
# taken before sections 7-10 were appended (the backup is kept in backups/, which .gitignore excludes)
def _cells(pth):
    out = []
    for d0, c0 in P.iter_cells(P.parse_notebook(pth)):
        st0, ct0, op0 = P.cell_info(c0)
        if st0 not in ('Output', 'Message', 'Print'): out.append((st0, P.box_to_text(ct0)))
    return out
# sections 20-23 (stage 4, 27 October 2026) are removed before the comparisons with the older backups, and the title of the summary
# section ("24. Summary of all checks (sections 1-23)") is given the title it had before ("20. Summary of all checks (sections 1-19)")
STAGE4_FIRST = '20. The coupled field equations'
STAGE4_SECS = (STAGE4_FIRST, '21. Regime', '22. Regime', '23. Regime')
# section 24 (dirac16complex as dark matter and dark energy, 27 October 2026): it lies between section 23 and the summary, so _strip4
# removes it too; the comparisons of outputs leave it out as well.  The one cell of sections 0-23 that changed with it:
STAGE5_FIRST = '24. dirac16complex'
STAGE5_SECS = (STAGE5_FIRST,)
CORRECTIONS5 = [('(Omega_c = 0.264 if it is the dark matter)', '(0.265 if it is the dark matter)')]   # (new text, text of the backup)
def _is_summary(t): return 'Summary of all checks' in t
def _strip4(cl, old_title):
    """cl: list of (style, text); drop the cells of sections 20-23 and give the summary section that follows them (the final one) its
    old title; the earlier summary sections (6 and 10) are left alone"""
    out = []; skip = False
    for st0, t0 in cl:
        if st0 == 'Section':
            tt0 = t0.strip('"')
            if tt0.startswith(STAGE4_FIRST): skip = True
            elif skip and _is_summary(tt0):
                skip = False
                if old_title is not None: t0 = old_title
        if not skip: out.append((st0, t0))
    return out
# the signature explanation of 27 October 2026 (the author's request): every cell it changed or added says that
# diag(+1,-1,...,-1) is simply the squares of the octonion basis, or contains the xAct signature list {7, 1, 0}
SIG_MARKERS = ('squares of the octonion basis', '{7, 1, 0}')
def _compare_sig(old, new):
    """(ok, equal, changed, added): every cell of old is present in new in the same order, unchanged or changed by the signature
    explanation; every changed or added cell of new contains one of SIG_MARKERS"""
    import difflib
    sm = difflib.SequenceMatcher(a=old, b=new, autojunk=False)
    eq = ch = ad = 0; ok = True
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == 'equal': eq += i2 - i1; continue
        from collections import Counter
        if tag == 'delete' or (i2 - i1) > (j2 - j1) or (Counter(c[0] for c in old[i1:i2]) - Counter(c[0] for c in new[j1:j2])):
            ok = False                     # a cell of the backup is missing, or a cell changed its style
        ch += i2 - i1; ad += (j2 - j1) - (i2 - i1)
        if not all(any(m in new[j][1] for m in SIG_MARKERS) for j in range(j1, j2)): ok = False
    return ok, eq, ch, ad
_b = _cells(path)
bk = os.path.join(os.path.dirname(here), 'backups', '7+1_2026-10-26_before_dirac16complex.nb')
if os.path.exists(bk):
    _a = _cells(bk)
    _k = next((i for i, (s0, t0) in enumerate(_b) if s0 == 'Section' and t0.strip('"').startswith('7.')), len(_b))
    _ok, _eq, _ch, _ad = _compare_sig(_a, _b[:_k])
    R['sections 0-6: the %d input and text cells of the backup are all present, %d unchanged, %d changed and %d added by the '
      'signature explanation of 27 October 2026' % (len(_a), _eq, _ch, _ad)] = _ok
else:
    R['sections 0-6: backup not present, comparison skipped'] = 'skipped'

# ---------------- sections 11-19: gravity, Kohn-Sham in gravitational fields, Standard Model ----------------
stage3 = any(q in ('11', '12', '13', '14', '15', '16', '17', '18', '19') for q in outsec)
if stage3:
    for sec in ('11', '12', '13', '14', '15', '16', '17', '18', '19'):
        so = [o for o, q in zip(outs, outsec) if q == sec and chk_re.fullmatch(o.strip())]
        R['section %s: check outputs (all True)' % sec] = '%d (%s)' % (len(so), all(o.strip().endswith('->True') for o in so))
    grid = outs[-1] if outs else ''
    allc = re.findall(r'"((?:[^"\\]|\\.)+)",\s*(True|False)', grid)
    R['summary table: entries == Length[checks], all True'] = (str(len(allc)) == str(R['Length[checks]']) and all(r == 'True' for n_, r in allc))
    names = set(n_ for n_, r in allc)
    printed = [re.match(r'"(.*)"->(True|False)$', o.strip(), flags=re.S).group(1) for o in chk]
    R['summary table: every printed check is listed'] = all(n_ in names for n_ in printed)
    bk3 = os.path.join(os.path.dirname(here), 'backups', '7+1_2026-10-26_before_gravity.nb')
    if os.path.exists(bk3):
        _a3 = _cells(bk3)
        _k3 = next((i for i, (s0, t0) in enumerate(_b) if s0 == 'Section' and t0.strip('"').startswith('11.')), len(_b))
        _ok, _eq, _ch, _ad = _compare_sig(_a3, _b[:_k3])
        R['sections 0-10: the %d input and text cells of the backup taken before sections 11-20 are all present, %d unchanged, %d '
          'changed and %d added by the signature explanation of 27 October 2026' % (len(_a3), _eq, _ch, _ad)] = _ok
    else:
        R['sections 0-10: backup not present, comparison skipped'] = 'skipped'
    bk4 = os.path.join(os.path.dirname(here), 'backups', '7+1_2026-10-27_before_signature.nb')
    if os.path.exists(bk4):
        _a4 = _cells(bk4)
        _t4 = next((t0 for st0, t0 in reversed(_a4) if st0 == 'Section'), None)
        _ok, _eq, _ch, _ad = _compare_sig(_a4, _strip4(_b, _t4))
        R['sections 0-20: the %d input and text cells of the backup taken before the signature explanation are all present, %d '
          'unchanged, %d changed and %d added by that explanation' % (len(_a4), _eq, _ch, _ad)] = _ok
    else:
        R['sections 0-20: backup not present, comparison skipped'] = 'skipped'
    if os.path.exists(bk4):
        # the OUTPUTS: the explanation of 27 October 2026 must not change any result.  Allowed differences: the new check outputs
        # (each with a marker), the check counts and tables of the summary sections, the declared list {7,1,0} itself, the name of
        # the private notebook context ($Context, section 0) and measured times (bare real numbers).
        def _outs(pth):
            out = []; sec = ''
            for d0, c0 in P.iter_cells(P.parse_notebook(pth)):
                st0, ct0, op0 = P.cell_info(c0)
                if st0 == 'Section': sec = P.box_to_text(ct0).strip('"')
                if st0 in ('Output', 'Message', 'Print'): out.append((sec, P.box_to_text(ct0)))
            return out
        import difflib
        _o4 = _outs(bk4); _n4 = [x_ for x_ in _outs(path) if not x_[0].startswith(STAGE4_SECS + STAGE5_SECS)]
        _sm = difflib.SequenceMatcher(a=[t for q, t in _o4], b=[t for q, t in _n4], autojunk=False)
        _same = 0; _new = 0; _ok = True; _bad = []
        _num = re.compile(r'\d+\.\d*`?')
        def _allowed(o, n):
            _ctx = lambda t: re.sub(r'Notebook\$\$[0-9$]+`', 'Notebook`', t)
            if n[0].startswith('0.') and o[0].startswith('0.') and 'Notebook$$' in n[1] and _ctx(o[1]) == _ctx(n[1]): return True
            if _num.fullmatch(o[1].strip()) and _num.fullmatch(n[1].strip()): return True
            if 'Summary of all checks' in n[0] and 'Summary of all checks' in o[0]: return True
            if '"SignatureOfMetric[g71]"->{7,1,0}' in n[1] and o[1].replace('{1,7,0}', '{7,1,0}') == n[1]: return True
            return False
        for tag, i1, i2, j1, j2 in _sm.get_opcodes():
            if tag == 'equal': _same += i2 - i1; continue
            olds = _o4[i1:i2]; news = _n4[j1:j2]
            # pair the changed outputs in order; the surplus of new outputs must be new check outputs with a marker
            for k in range(len(olds)):
                if k >= len(news) or not _allowed(olds[k], news[k]): _ok = False; _bad.append((olds[k][0][:30], olds[k][1][:80]))
            for n_ in news[len(olds):]:
                if re.fullmatch(r'"[^"]*(?:"[^"]*)*"->True', n_[1].strip(), flags=re.S) and any(m in n_[1] for m in SIG_MARKERS): _new += 1
                else: _ok = False; _bad.append((n_[0][:30], n_[1][:80]))
        R['sections 0-20, outputs: of the %d outputs of the backup taken before the signature explanation %d are identical; the others '
          'differ only in the check counts and tables of the summary sections, the declared list {7,1,0}, the private context name and '
          'measured times; %d new check outputs, all True' % (len(_o4), _same, _new)] = _ok
        if _bad: R['  outputs that differ otherwise'] = _bad[:10]
    # printed numbers versus the independent Python programs
    import json, subprocess
    rfile = os.path.join(os.path.dirname(here), 'build', 'docs', 'runs.json')
    pyout = {}
    if os.path.exists(rfile):
        for t_, o_, c_ in json.load(open(rfile, encoding='utf-8')):
            pyout[t_.split(' ')[0]] = o_
    for f_ in ('ref_gravity.py', 'ref_ks_gravity.py', 'ref_sm.py'):
        if 'tools/' + f_ not in pyout:
            p_ = subprocess.run([sys.executable, '-B', os.path.join(here, f_)], capture_output=True, text=True, cwd=os.path.dirname(here))
            pyout['tools/' + f_] = p_.stdout
    from doc_stage3_common import gather3, parse_ref
    runs_ = [['tools/%s (x)' % f_, pyout['tools/' + f_], 0] for f_ in ('ref_gravity.py', 'ref_ks_gravity.py', 'ref_sm.py')]
    D3 = gather3(runs_)
    import doc_ks_gravity, doc_sm
    worst = []
    for lab, a_, b_ in doc_ks_gravity.compare(D3):
        r_ = doc_ks_gravity.max_rel(a_, b_)
        worst.append(('KS ' + lab, r_, 1e-6))
    for lab, a_, b_, note in doc_sm.compare(D3):
        r_ = doc_sm.max_rel(a_, b_)
        # rounded prints: 1e-4; ODE integration through the Sun and the singular values of a matrix with entries from 1e-2 eV
        # to 1e9 eV (machine precision times the condition number): 1e-5; everything else 1e-6
        tol = 1e-4 if note else (1e-5 if ('Sun' in lab or 'portal' in lab or 'N^' in lab) else 1e-6)
        worst.append(('SM ' + lab, r_, tol))
    g_ = D3['nums']['grav']; refg = D3['ref'].get('grav', ('', '', {}))[2]
    from doc_stage3_common import nums_in, mnum
    kb_nb = [mnum(g_.get(k_, '') or '') for k_ in ('Kasner beta2 (1,1)', 'Kasner beta2 (2,1)', 'Kasner beta2 (1,3)')]
    kb_py = nums_in(refg.get('Kasner (Einstein branch), massless: |beta|^2 for (|p|, |q|) = (1,1), (2,1), (1,3)', ''))
    worst.append(('gravity: Kasner |beta|^2', doc_ks_gravity.max_rel(kb_nb, kb_py), 1e-5))
    bad = [(l_, r_) for l_, r_, tol in worst if r_ is None or r_ > tol]
    R['sections 14-19: %d groups of printed numbers agree with the Python programs (relative 1e-6; rounded prints 1e-4; ODE and ill-conditioned results 1e-5)' % len(worst)] = not bad
    if bad:
        R['  disagreeing groups'] = bad

# ---------------- sections 20-23: the coupled field equations and their solutions in the three gravity regimes ----------------
stage4 = any(q in ('21', '22', '23') for q in outsec)      # before stage 4 the last section (the summary) was section 20
if stage4:
    for sec in ('20', '21', '22', '23'):
        so = [o for o, q in zip(outs, outsec) if q == sec and chk_re.fullmatch(o.strip())]
        R['section %s: check outputs (all True)' % sec] = '%d (%s)' % (len(so), all(o.strip().endswith('->True') for o in so))
    bk5 = os.path.join(os.path.dirname(here), 'backups', '7+1_2026-10-27_before_coupled.nb')
    if os.path.exists(bk5):
        _a5 = _cells(bk5)
        _t5 = next((t0 for st0, t0 in reversed(_a5) if st0 == 'Section'), None)
        _b5 = _strip4(_b, _t5)
        R['stage 4: without sections 20-23 (and with the old title of the summary section) the %d input and text cells equal those of '
          'the backup taken before sections 20-23 were appended' % len(_a5)] = (_a5 == _b5)
        def _outs5(pth):
            out = []; sec = ''
            for d0, c0 in P.iter_cells(P.parse_notebook(pth)):
                st0, ct0, op0 = P.cell_info(c0)
                if st0 == 'Section': sec = P.box_to_text(ct0).strip('"')
                if st0 in ('Output', 'Message', 'Print'): out.append((sec, P.box_to_text(ct0)))
            return out
        _o5 = _outs5(bk5); _n5 = [x_ for x_ in _outs5(path) if not x_[0].startswith(STAGE4_SECS + STAGE5_SECS)]
        _ctx5 = lambda t: re.sub(r'Notebook\$\$[0-9$]+`', 'Notebook`', t)
        _numre5 = re.compile(r'\d+\.\d*`?')
        _bad5 = []
        if len(_o5) != len(_n5): _bad5.append(('count', len(_o5), len(_n5)))
        else:
            for (qa, ta), (qb, tb) in zip(_o5, _n5):
                if ta == tb: continue
                if _ctx5(ta) == _ctx5(tb): continue
                if _numre5.fullmatch(ta.strip()) and _numre5.fullmatch(tb.strip()): continue
                if _is_summary(qa) and _is_summary(qb): continue
                _bad5.append((qa[:30], ta[:80]))
        R['stage 4: the %d outputs of sections 0-19 and of the summary equal those of the backup (summary counts and tables, the private '
          'context name and measured times excepted)' % len(_o5)] = not _bad5
        if _bad5: R['  stage-4 outputs that differ'] = _bad5[:10]
    else:
        R['stage 4: backup not present, comparison skipped'] = 'skipped'
    import json, subprocess
    rfile4 = os.path.join(os.path.dirname(here), 'build', 'docs', 'runs.json')
    pyc = None
    if os.path.exists(rfile4):
        for t_, o_, c_ in json.load(open(rfile4, encoding='utf-8')):
            if t_.startswith('tools/ref_coupled.py'): pyc = o_
    if pyc is None:
        pyc = subprocess.run([sys.executable, '-B', os.path.join(here, 'ref_coupled.py')], capture_output=True, text=True, cwd=os.path.dirname(here)).stdout
    from doc_stage3_common import gather3
    D4 = gather3([['tools/ref_coupled.py (x)', pyc, 0]])
    import doc_coupled
    cmp4 = doc_coupled.compare(D4)
    bad4 = [(lab, r_, tol) for lab, a_, b_, tol, r_ in cmp4 if r_ is None or r_ > tol]
    R['sections 20-23: %d groups of printed numbers agree with tools/ref_coupled.py (tolerances per group, from 1e-14 for exact constants to 2e-3 for numbers read off numerical solutions)' % len(cmp4)] = not bad4
    if bad4: R['  disagreeing groups (sections 20-23)'] = bad4

# ---------------------------------------------------------------- section 24: dirac16complex as dark matter and dark energy ----------------
stage5 = any(s0 == 'Section' and t0.strip('"').startswith(STAGE5_FIRST) for s0, t0 in _b)
if stage5:
    so = [o for o, q in zip(outs, outsec) if q == '24' and chk_re.fullmatch(o.strip())]
    R['section 24: check outputs (all True)'] = '%d (%s)' % (len(so), all(o.strip().endswith('->True') for o in so))
    def _strip5(cl, old_title):
        """drop the cells of section 24 and give the summary section that follows it its old title"""
        out = []; skip = False
        for st0, t0 in cl:
            if st0 == 'Section':
                tt0 = t0.strip('"')
                if tt0.startswith(STAGE5_FIRST): skip = True
                elif skip and _is_summary(tt0):
                    skip = False
                    if old_title is not None: t0 = old_title
            if not skip: out.append((st0, t0))
        return out
    def _undo5(t):
        n = 0
        for new_, old_ in CORRECTIONS5:
            if new_ in t: t = t.replace(new_, old_); n += 1
        return t, n
    bk6 = os.path.join(os.path.dirname(here), 'backups', '7+1_2026-10-27_before_dark.nb')
    if os.path.exists(bk6):
        _a6 = _cells(bk6)
        _t6 = next((t0 for st0, t0 in reversed(_a6) if st0 == 'Section'), None)
        _b6 = []; _nc6 = 0
        for st0, t0 in _strip5(_b, _t6):
            t1, n1 = _undo5(t0); _b6.append((st0, t1)); _nc6 += n1
        R['stage 5: without section 24 (and with the old title of the summary section) the %d input and text cells equal those of the '
          'backup taken before section 24 was appended, after undoing the %d documented correction(s) of section 23.3' % (len(_a6), _nc6)] = (
            _a6 == _b6 and _nc6 == len(CORRECTIONS5))
        def _outs6(pth):
            out = []; sec = ''
            for d0, c0 in P.iter_cells(P.parse_notebook(pth)):
                st0, ct0, op0 = P.cell_info(c0)
                if st0 == 'Section': sec = P.box_to_text(ct0).strip('"')
                if st0 in ('Output', 'Message', 'Print'): out.append((sec, P.box_to_text(ct0)))
            return out
        _o6 = _outs6(bk6); _n6 = [x_ for x_ in _outs6(path) if not x_[0].startswith(STAGE5_SECS)]
        _ctx6 = lambda t: re.sub(r'Notebook\$\$[0-9$]+`', 'Notebook`', t)
        _numre6 = re.compile(r'\d+\.\d*`?')
        _bad6 = []; _oc6 = 0
        if len(_o6) != len(_n6): _bad6.append(('count', len(_o6), len(_n6)))
        else:
            for (qa, ta), (qb, tb) in zip(_o6, _n6):
                if ta == tb: continue
                tb1, n1 = _undo5(tb)
                if n1 and tb1 == ta: _oc6 += 1; continue
                if _ctx6(ta) == _ctx6(tb): continue
                if _numre6.fullmatch(ta.strip()) and _numre6.fullmatch(tb.strip()): continue
                if _is_summary(qa) and _is_summary(qb): continue
                _bad6.append((qa[:30], ta[:80]))
        R['stage 5: the %d outputs of sections 0-23 and of the summary equal those of the backup (summary counts and tables, the private '
          'context name and measured times excepted; %d output(s) differ only by the documented correction)' % (len(_o6), _oc6)] = not _bad6
        if _bad6: R['  stage-5 outputs that differ'] = _bad6[:10]
    else:
        R['stage 5: backup not present, comparison skipped'] = 'skipped'
    import json, subprocess
    rfile5 = os.path.join(os.path.dirname(here), 'build', 'docs', 'runs.json')
    pyd = None
    if os.path.exists(rfile5):
        for t_, o_, c_ in json.load(open(rfile5, encoding='utf-8')):
            if t_.startswith('tools/ref_dark.py'): pyd = o_
    if pyd is None:
        pyd = subprocess.run([sys.executable, '-B', os.path.join(here, 'ref_dark.py')], capture_output=True, text=True, cwd=os.path.dirname(here)).stdout
    from doc_stage3_common import gather3
    D5 = gather3([['tools/ref_dark.py (x)', pyd, 0]])
    import doc_dark
    cmp5 = doc_dark.compare(D5)
    bad5 = [(lab, r_, tol) for lab, a_, b_, tol, r_ in cmp5 if r_ is None or r_ > tol]
    R['section 24 numbers: %d groups of printed numbers agree with tools/ref_dark.py (tolerances per group, from 1e-14 for closed forms to 1e-5 for a number read off a numerical solution)' % len(cmp5)] = not bad5
    if bad5: R['  disagreeing groups (section 24)'] = bad5

width = max(len(k) for k in R)
for k, v in R.items(): print(f'{k:{width}s} : {v}')
good = (R['message cells (must be 0)'] == 0 and R['check outputs'] == R['check outputs that are True'] and R['check outputs'] > 0
        and (stage3 or (R['check["..."] calls in the input cells'] == R['check outputs'] and str(R['Length[checks]']) == str(R['check outputs'])))
        and R['failedChecks == {}'] and R['allChecksPassed output True']
        and all(R[k] for k in R if k.endswith('(python)') or 'printed' in k or k.startswith('private'))
        and all(str(R[k]).endswith('(True)') for k in R if k.startswith('section ') and 'check outputs' in k)
        and all(R[k] for k in R if k.startswith('section 7:') or k.startswith('section 8:') or k.startswith('sections 8-9') or k.startswith('tools/ref_dirac16'))
        and all(R[k] is True or R[k] == 'skipped' for k in R if k.startswith('sections 0-6'))
        and all(R[k] is True or R[k] == 'skipped' for k in R if k.startswith('sections 0-10') or k.startswith('sections 0-20') or k.startswith('summary table')
                or k.startswith('sections 14-19') or k.startswith('stage 4') or k.startswith('sections 20-23') or k.startswith('stage 5')
                or k.startswith('section 24 numbers')))
print('\nVERDICT:', 'ALL CROSS-CHECKS PASSED' if good else 'SOMETHING FAILED')
sys.exit(0 if good else 1)
