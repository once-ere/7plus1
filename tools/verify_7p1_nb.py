#!/usr/bin/env python3
"""Independent cross-check of the EVALUATED notebook 7+1.nb against the Python reference (ref71.py).

Usage (from the repository root):  python3 tools/verify_7p1_nb.py [7+1.nb]
Checks: no message cells; every check[...] output is True; Length[checks], failedChecks, allChecksPassed;
the matrices tau[A], tauO[A], Gamma17[A] printed by Mathematica equal the Python ones; gridx44 and gridO tables;
Lambda matrices; the eight real on-shell spinors; the classification table; $Version.
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
outs = []; msgs = []; inputs = 0; check_calls = 0
for d, c in cells:
    st, content, opts = P.cell_info(c)
    if st == 'Input':
        inputs += 1; check_calls += P.box_to_text(content).count('check["')
    if st == 'Output': outs.append(P.box_to_text(content))
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
    W = parse_int_matrix(maj[-1])
    p = [3, 2, 2, 0, 0, 0, 0, 0]
    ps = sum((p[A] * G17[A] for A in range(8)), sp.zeros(16))
    ok = W.shape == (8, 16) and W.rank() == 8 and all(((ps - sp.eye(16)) * W.row(i).T) == sp.zeros(16, 1) for i in range(8))
R['8 real on-shell spinors printed solve (pslash - m) w = 0 (python check)'] = ok
cl = [o for o in outs if 'WRONG in [1B]' in o]
R['classification table printed'] = bool(cl) and cl[0].count('WRONG in [1B]') == 3 and cl[0].count('agrees') == 6
width = max(len(k) for k in R)
for k, v in R.items(): print(f'{k:{width}s} : {v}')
good = (R['message cells (must be 0)'] == 0 and R['check outputs'] == R['check outputs that are True'] and R['check outputs'] > 0
        and R['check["..."] calls in the input cells'] == R['check outputs']
        and str(R['Length[checks]']) == str(R['check outputs']) and R['failedChecks == {}'] and R['allChecksPassed output True']
        and all(R[k] for k in R if k.endswith('(python)') or 'printed' in k or k.startswith('private')))
print('\nVERDICT:', 'ALL CROSS-CHECKS PASSED' if good else 'SOMETHING FAILED')
sys.exit(0 if good else 1)
