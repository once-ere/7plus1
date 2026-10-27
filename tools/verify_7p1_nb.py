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
bk = os.path.join(os.path.dirname(here), 'backups', '7+1_2026-10-26_before_dirac16complex.nb')
if os.path.exists(bk):
    def _cells(pth):
        out = []
        for d0, c0 in P.iter_cells(P.parse_notebook(pth)):
            st0, ct0, op0 = P.cell_info(c0)
            if st0 not in ('Output', 'Message', 'Print'): out.append((st0, P.box_to_text(ct0)))
        return out
    _a = _cells(bk); _b = _cells(path)
    _k = next((i for i, (s0, t0) in enumerate(_b) if s0 == 'Section' and t0.strip('"').startswith('7.')), len(_b))
    R['sections 0-6: the %d input and text cells equal those of the backup' % len(_a)] = (_a == _b[:_k])
else:
    R['sections 0-6: backup not present, comparison skipped'] = 'skipped'

width = max(len(k) for k in R)
for k, v in R.items(): print(f'{k:{width}s} : {v}')
good = (R['message cells (must be 0)'] == 0 and R['check outputs'] == R['check outputs that are True'] and R['check outputs'] > 0
        and R['check["..."] calls in the input cells'] == R['check outputs']
        and str(R['Length[checks]']) == str(R['check outputs']) and R['failedChecks == {}'] and R['allChecksPassed output True']
        and all(R[k] for k in R if k.endswith('(python)') or 'printed' in k or k.startswith('private'))
        and all(str(R[k]).endswith('(True)') for k in R if k.startswith('section ') and 'check outputs' in k)
        and all(R[k] for k in R if k.startswith('section 7:') or k.startswith('section 8:') or k.startswith('sections 8-9') or k.startswith('tools/ref_dirac16'))
        and all(R[k] is True or R[k] == 'skipped' for k in R if k.startswith('sections 0-6')))
print('\nVERDICT:', 'ALL CROSS-CHECKS PASSED' if good else 'SOMETHING FAILED')
sys.exit(0 if good else 1)
