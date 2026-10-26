"""Data shared by the documentation builders: every matrix and table printed in the documents is
computed here from the Python reference implementation (ref44.py, ref71.py) or parsed from the
evaluated notebook 7+1.nb, never typed by hand."""
import os, re, sys, hashlib, datetime
sys.setrecursionlimit(100000)
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import sympy as sp
import nbparse as P
from ref71 import *          # ref44 objects and the 7+1 constructions
from docgen import pmat

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''): h.update(chunk)
    return h.hexdigest()

def table_cells(C, sym, base, lead):
    """multiplication table as rows of LaTeX inline math strings"""
    rows = []
    for Bi in range(8):
        row = ['$%s_{%d}$' % (sym, Bi + base)]
        for k in range(8):
            p = mul(C, basis(Bi), basis(k)); d = [i for i in range(8) if p[i]][0]
            row.append('$%s%s_{%d}$' % ('-' if p[d] < 0 else '', sym, d + base))
        rows.append(row)
    header = [lead] + ['$%s_{%d}$' % (sym, k + base) for k in range(8)]
    return header, rows

def eight_square():
    a = sp.symbols('a0:8'); b = sp.symbols('b0:8')
    xy = mul(CO, list(a), list(b))
    return [sp.expand(c) for c in xy], a, b

def parse_evaluated_notebook(path=None):
    path = path or os.path.join(ROOT, '7+1.nb')
    nb = P.parse_notebook(path)
    outs = []; msgs = 0; inputs = 0; checks = []; section = ''; subsection = ''
    for d, c in P.iter_cells(nb):
        st, content, opts = P.cell_info(c)
        if st == 'Section': section = P.box_to_text(content); subsection = ''
        elif st == 'Subsection': subsection = P.box_to_text(content)
        elif st == 'Input': inputs += 1
        elif st == 'Output':
            o = P.box_to_text(content); outs.append(o)
            m = re.fullmatch(r'"(.*)"->(True|False)', o.strip(), flags=re.S)
            if m: checks.append((m.group(1), m.group(2), section, subsection))
        elif st == 'Message': msgs += 1
    version = next((o for o in outs if '15.0.1' in o and 'Windows' in o), '')
    ver = re.search(r'"([^"]*15\.0\.1[^"]*)"', version).group(1) if version else ''
    ctx = re.search(r'"(Notebook\$\$[^"]*)"', version).group(1).rstrip('`') if version else ''
    comm_time = next((o.strip().rstrip('`') for o in outs if re.fullmatch(r'\d+\.\d+`?', o.strip())), '')
    # the 8 real on-shell spinors printed by Mathematica
    maj = [o for o in outs if o.startswith('(\\[NoBreak]{{') and o.count('},') == 7 and len(re.findall(r'\{([^{}]*)\}', o)[0].split(',')) == 16]
    W = None
    if maj:
        rows = re.findall(r'\{([^{}]*)\}', maj[-1])
        W = sp.Matrix([[sp.Integer(x) for x in r.split(',')] for r in rows])
    return dict(inputs=inputs, outputs=len(outs), messages=msgs, checks=checks, version=ver, context=ctx,
                commutant_time=comm_time, W=W, path=path)

def eval_log():
    p = os.path.join(ROOT, 'build', 'eval_log.txt')
    d = {}
    if os.path.exists(p):
        for line in open(p, encoding='utf-8', errors='replace'):
            if ':' in line:
                k, v = line.split(':', 1); d[k.strip()] = v.strip()
    return d

def gather():
    D = {}
    D['s4'] = s4; D['t4'] = t4; D['tau'] = tau; D['tauO'] = tauO; D['EO'] = EAaO; D['E44'] = EAa(8)
    D['oct_table'] = table_cells(CO, 'e', 0, '$\\cdot$')
    D['split_table'] = table_cells(CS, '\\varepsilon', 1, '$\\cdot$')
    D['xy'], D['a'], D['b'] = eight_square()
    D['nb'] = parse_evaluated_notebook()
    D['log'] = eval_log()
    D['sha'] = {name: sha256(os.path.join(ROOT, name)) for name in ['7+1.nb', '4+4.nb']}
    D['sha_tools'] = {name: sha256(os.path.join(HERE, name)) for name in sorted(os.listdir(HERE))
                      if os.path.isfile(os.path.join(HERE, name)) and (name.endswith('.py') or name.endswith('.wl'))}
    D['size44'] = os.path.getsize(os.path.join(ROOT, '4+4.nb'))
    D['size71'] = os.path.getsize(os.path.join(ROOT, '7+1.nb'))
    D['today'] = datetime.date.today().isoformat()
    # numerical octonion example of section 3.4 of the notebook
    def ecomb(v):
        terms = []
        for k, c in enumerate(v):
            c = sp.nsimplify(c)
            if c == 0: continue
            s = ('+' if c > 0 else '-')
            mag = abs(c)
            coef = '' if mag == 1 else sp.latex(mag)
            terms.append((s, coef + 'e_{%d}' % k))
        if not terms: return '0'
        out = ''.join((s if i or s == '-' else '') + t for i, (s, t) in enumerate(terms))
        return out
    xN = [1, 2, 0, 0, -1, 0, 0, 3]; yN = [-2, 0, 1, 0, 0, 1, 0, 0]; zN = [0, 0, 0, 1, 0, 0, 1, 0]
    xy = mul(CO, xN, yN); yx = mul(CO, yN, xN)
    Nx = sum(t * t for t in xN); Ny = sum(t * t for t in yN); Nxy = sum(t * t for t in xy)
    xinv = [sp.Rational(c, Nx) for c in conj(xN)]
    assoc_xyz = sub(mul(CO, mul(CO, xN, yN), zN), mul(CO, xN, mul(CO, yN, zN)))
    D['num'] = dict(x=ecomb(xN), y=ecomb(yN), z=ecomb(zN), xy=ecomb(xy), yx=ecomb(yx), Nx=Nx, Ny=Ny, Nxy=Nxy,
                    xinv=ecomb(xinv), assoc=ecomb(assoc_xyz))
    D['ecomb'] = ecomb
    return D

if __name__ == '__main__':
    D = gather()
    print('checks:', len(D['nb']['checks']), 'all True:', all(c[1] == 'True' for c in D['nb']['checks']))
    print('version:', D['nb']['version'], ' context:', D['nb']['context'], ' commutant time:', D['nb']['commutant_time'])
    print('messages:', D['nb']['messages'], ' W shape:', D['nb']['W'].shape if D['nb']['W'] is not None else None)
    print('log:', D['log'])
    print('sha:', D['sha'])
    print('xy[0] =', D['xy'][0])
