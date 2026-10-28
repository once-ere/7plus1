"""Data for the stage-3, stage-4, dark-sector, torus, brane and radion documents (GRAVITY*, KOHN_SHAM_GRAVITY*, STANDARD_MODEL*,
COUPLED_SOLUTIONS*, DARK_SECTOR*, AUXILIARY_TORUS*, WORLD_BRANE*, RADION*): parsed from the evaluated 7+1.nb (sections 11-26 and the final
summary of all checks) and from the outputs of the independent Python programs tools/ref_gravity.py, tools/ref_ks_gravity.py, tools/ref_sm.py,
tools/ref_coupled.py, tools/ref_dark.py, tools/ref_torus_brane.py and tools/ref_radion.py.  Nothing in the documents that can be computed is
typed by hand."""
import os, re, sys, platform
from collections import OrderedDict
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import nbparse as P
from doc_common import sha256, eval_log

STAGE3_SECTIONS = {'grav': ('11.', '12.', '13.', '14.'), 'ks': ('15.', '16.'), 'sm': ('17.', '18.', '19.'), 'cp': ('20.', '21.', '22.', '23.'),
                   'dk': ('24. dirac16complex',), 'tb': ('25. The auxiliary torus',), 'rd': ('26. The radion',)}
REF = {'grav': 'ref_gravity.py', 'ks': 'ref_ks_gravity.py', 'sm': 'ref_sm.py', 'cp': 'ref_coupled.py', 'dk': 'ref_dark.py',
       'tb': 'ref_torus_brane.py', 'rd': 'ref_radion.py'}
FRAME_PREFIXES = ('cosmological metric:', 'pp-wave:', 'inhomogeneous metric', 'static radial metric:')

def notebook_rows(path=None):
    """(style, section, subsection, text) for every Output and Message cell of the evaluated notebook"""
    path = path or os.path.join(ROOT, '7+1.nb')
    nb = P.parse_notebook(path)
    rows = []; section = ''; subsection = ''; inputs = 0; sections = []
    for d, c in P.iter_cells(nb):
        st, content, opts = P.cell_info(c)
        if st == 'Section':
            section = P.box_to_text(content).strip('"'); subsection = ''; sections.append(section)
        elif st == 'Subsection': subsection = P.box_to_text(content).strip('"')
        elif st == 'Input': inputs += 1
        elif st in ('Output', 'Message', 'Print'):
            rows.append((st, section, subsection, P.box_to_text(content)))
    return rows, inputs, sections

def split_top(s, sep=','):
    """split s at top-level separators (outside brackets, association brackets <| |> and strings); the > of -> is not a bracket"""
    out, depth, cur, instr = [], 0, '', False
    for i, ch in enumerate(s):
        if ch == '"' and not cur.endswith('\\'): instr = not instr
        if not instr:
            if ch in '{[(' or (ch == '<' and s[i + 1:i + 2] == '|'): depth += 1
            elif ch in '}])' or (ch == '>' and cur.endswith('|')): depth -= 1
            elif ch == sep and depth == 0:
                out.append(cur); cur = ''; continue
        cur += ch
    out.append(cur)
    return [x.strip() for x in out]

def rules_in(text):
    """the rules "name" -> value of an output that is a list of rules (or a single rule)"""
    t = text.strip()
    if t.startswith('{') and t.endswith('}'): t = t[1:-1]
    res = OrderedDict()
    for part in split_top(t):
        m = re.match(r'"((?:[^"\\]|\\.)*)"\s*->\s*(.*)$', part, flags=re.S)
        if m: res[m.group(1)] = m.group(2).strip()
    return res

def mnum(v):
    """a Mathematica number (with precision marks and *^ exponents) as a float, or None"""
    v = v.strip()
    m = re.fullmatch(r'(-?\d+\.?\d*)(?:`[\d.]*)?(?:\*\^(-?\d+))?', v)
    if not m: return None
    return float(m.group(1) + ('e' + m.group(2) if m.group(2) else ''))

def mlist(v):
    """a flat Mathematica list of numbers as floats"""
    v = v.strip()
    if not (v.startswith('{') and v.endswith('}')): return None
    items = split_top(v[1:-1])
    vals = [mnum(x) for x in items]
    return vals if all(x is not None for x in vals) else None

def parse_ref(out):
    d = OrderedDict()
    for line in out.splitlines():
        if ' : ' in line:
            k, v = line.split(' : ', 1); d[k.strip()] = v.strip()
    return d

def gather3(runs=None):
    rows, inputs, sections = notebook_rows()
    D = {'rows': rows, 'inputs': inputs, 'sections': sections}
    secnums = [int(t.split('.')[0]) for t in sections if t[:1].isdigit() and t.split('.')[0].isdigit()]
    last_is_summary = bool(sections) and 'Summary of all checks' in sections[-1]
    D['n_last'] = (max(secnums) - 1 if last_is_summary else max(secnums)) if secnums else 19
    outs = [(sec, sub, t) for st, sec, sub, t in rows if st == 'Output']
    D['messages'] = [(sec, sub, t) for st, sec, sub, t in rows if st == 'Message']
    checks = OrderedDict()
    for key, pre in STAGE3_SECTIONS.items():
        lst = []
        for sec, sub, t in outs:
            if not sec.startswith(pre): continue
            m = re.fullmatch(r'"(.*)"->(True|False)', t.strip(), flags=re.S)
            if m: lst.append((m.group(1), m.group(2), sec, sub))
        checks[key] = lst
    # the final summary (last section) lists every check, including those made inside frameChecks
    last = sections[-1] if sections else ''
    summary_outs = [t for sec, sub, t in outs if sec == last]
    D['summary_section'] = last
    D['summary_outs'] = summary_outs
    grid = next((t for t in summary_outs if 'Grid' in t or 'check' in t and 'result' in t), '')
    allc = re.findall(r'"((?:[^"\\]|\\.)+)",\s*(True|False)', grid)
    D['all_checks_grid'] = allc
    known = set(n for v in checks.values() for n, *_ in v)
    frame = [(n, r, '13. Generalized frame matrices, the canonical spin connection, and the covariant Dirac equation of dirac16complex', '13.1')
             for n, r in allc if n.startswith(FRAME_PREFIXES) and n not in known]
    checks['grav'] = checks['grav'] + frame
    D['checks'] = checks
    D['n_checks_total'] = next((mnum(t) for t in summary_outs if mnum(t) is not None), None)
    if D['n_checks_total'] is not None and float(D['n_checks_total']).is_integer():
        D['n_checks_total'] = int(D['n_checks_total'])        # a count: 346, not 346.0
    D['all_passed'] = next((t.strip() for t in summary_outs if t.strip() in ('True', 'False')), '')
    D['failed'] = next((t.strip() for t in summary_outs if t.strip().startswith('{') and 'Grid' not in t and len(t) < 4000), '')
    nums = OrderedDict()
    for key, pre in STAGE3_SECTIONS.items():
        d = OrderedDict()
        for sec, sub, t in outs:
            if sec.startswith(pre) and t.strip().startswith('{') and '->' in t:
                for k, v in rules_in(t).items():
                    if v not in ('True', 'False'): d[k] = v
        nums[key] = d
    D['nums'] = nums
    D['log'] = eval_log()
    try:
        D['closed_forms_tex'], D['kasner_tex'], D['closed_forms'] = lovelock_tex()
    except Exception as exc:
        D['closed_forms_tex'], D['kasner_tex'], D['closed_forms'] = {}, {}, {}
    # xAct: package versions and the Lovelock tensors with abstract indices (section 12.2-12.3)
    g = nums.get('grav', {})
    ver = OrderedDict()
    for pk in ('xPerm', 'xTensor', 'xCoba', 'xTras'):
        m = re.search(r'"([0-9][0-9.]*)"', g.get(pk, '') or '')
        if m: ver[pk] = m.group(1)
    D['xact_version_map'] = ver
    D['xact_versions'] = ', '.join('%s %s' % kv for kv in ver.items()) if ver else ''
    D['xact_lovelock_inputform'], D['xact_lovelock_tex'] = {}, {}
    for k in (1, 2, 3):
        v = g.get('Lovelock%d (InputForm)' % k)
        if not v: continue
        v = v.strip()
        if v.startswith('"') and v.endswith('"'): v = v[1:-1]
        v = v.replace('\\"', '"')
        D['xact_lovelock_inputform'][k] = v
        try:
            D['xact_lovelock_tex'][k] = xact_to_latex(v)[0]
        except Exception as exc:
            D.setdefault('xact_tex_errors', []).append('Lovelock%d: %r' % (k, exc))
    try:
        D['expo'] = exponential_solutions()
    except Exception as exc:
        D['expo'] = {}
        D.setdefault('xact_tex_errors', []).append('exponential solutions: %r' % (exc,))
    ref_files = [os.path.join(ROOT, f) for f in ('Einstein-Lovelock-7+1-Nash.nb', 'Einstein-Lovelock-7+1-Nash.pdf')]
    D['sha_reference'] = tuple(sha256(f) for f in ref_files) if all(os.path.exists(f) for f in ref_files) else None
    D['ref'] = {}
    if runs:
        for key, f in REF.items():
            for title, out, *rest in runs:
                if title.startswith('tools/' + f): D['ref'][key] = (title, out, parse_ref(out))
    files = ['7+1.nb', 'tools/make_7p1_nb.py', 'tools/nb_part_gravity.py', 'tools/nb_part_ks_gravity.py', 'tools/nb_part_sm.py',
             'tools/nb_part_coupled.py', 'tools/nb_part_dark.py', 'tools/nb_part_torus_brane.py', 'tools/nb_part_radion.py', 'tools/lovelock71.py',
             'tools/ref_gravity.py', 'tools/ref_ks_gravity.py', 'tools/ref_sm.py', 'tools/ref_coupled.py', 'tools/ref_dark.py', 'tools/ref_torus_brane.py',
             'tools/ref_radion.py', 'tools/ref71.py', 'tools/ref44.py']
    D['sha'] = OrderedDict((f, sha256(os.path.join(ROOT, f))) for f in files if os.path.exists(os.path.join(ROOT, f)))
    bk = os.path.join(ROOT, 'backups', '7+1_2026-10-26_before_gravity.nb')
    D['backup'] = (sha256(bk), os.path.getsize(bk)) if os.path.exists(bk) else None
    D['versions'] = dict(python=platform.python_version())
    for mod in ('sympy', 'numpy', 'scipy', 'mpmath'):
        try: D['versions'][mod] = __import__(mod).__version__
        except Exception: D['versions'][mod] = '?'
    return D

def lovelock_tex():
    """LaTeX of the closed forms E^(k)c_c (c = 0, 1, 4) of the cosmological metric and of x0^(2k) E^(k) for a = x0^p, b = x0^q"""
    import sympy as sp
    from lovelock71 import cosmological_closed_forms
    CF = cosmological_closed_forms()
    Ha, Hb, dHa, dHb = sp.symbols('H_a H_b dH_a dH_b', real=True)
    pa, pb = sp.Symbol("H'_a"), sp.Symbol("H'_b")
    p, q, x = sp.symbols('p q x_0', positive=True)
    cf, ks = {}, {}
    for (k, c), v in CF.items():
        cf[(k, c)] = sp.latex(sp.factor(v).subs({dHa: pa, dHb: pb}))
        w = sp.factor(sp.expand(x**(2*k)*v.subs({Ha: p/x, Hb: q/x, dHa: -p/x**2, dHb: -q/x**2})))
        ks[(k, c)] = sp.latex(w)
    return cf, ks, CF

def fmt(x, n=6):
    """a float as a short decimal string (LaTeX exponent for small or large numbers)"""
    if x is None: return '?'
    if x == 0: return '0'
    import math
    e = int(math.floor(math.log10(abs(x))))
    if -4 <= e <= 5:
        s = ('%.' + str(max(n - 1 - e, 0)) + 'f') % x
        return s
    m = x/10**e
    return (('%.' + str(n - 1) + 'f') % m) + r'\times10^{%d}' % e

if __name__ == '__main__':
    D = gather3()
    for k, v in D['checks'].items():
        print(k, len(v), 'all True:', all(r == 'True' for n, r, *_ in v))
    print('summary section:', D['summary_section'], 'grid checks:', len(D['all_checks_grid']))
    print('messages:', len(D['messages']))
    for k, v in D['nums'].items():
        print(k, list(v.items())[:6])

def break_sum(lhs, rhs, width=95, end=','):
    """align lines 'lhs &= rhs' with the LaTeX sum rhs broken at top-level + and - signs into lines of about width characters"""
    terms, depth, cur = [], 0, ''
    i = 0
    while i < len(rhs):
        ch = rhs[i]
        if ch in '{([': depth += 1
        elif ch in '})]': depth -= 1
        if depth == 0 and ch in '+-' and cur.strip() and not cur.rstrip().endswith(('^', '_', '(', '{', '\\left(')):
            terms.append(cur); cur = ch
        else:
            cur += ch
        i += 1
    terms.append(cur)
    lines, line = [], ''
    for t in terms:
        if len(line) + len(t) > width and line.strip():
            lines.append(line); line = t
        else:
            line += t
    lines.append(line)
    out = ['%s&=%s' % (lhs, lines[0].strip())]
    for l in lines[1:]:
        out.append('&\\quad %s' % l.strip())
    out[-1] += end
    return out


# ----------------------------------------------------------------------------- xAct InputForm -> LaTeX
_XHEAD = {'delta': (r'\delta', 'D0'), 'RicciScalarCD': ('R', 'D1'), 'RicciCD': ('R', 'D2'), 'RiemannCD': ('R', 'D3'), 'EinsteinCD': ('G', 'D4')}
_XTOK = re.compile(r'(RiemannCD|RicciCD|RicciScalarCD|delta|EinsteinCD)\[([^\]]*)\]')

def _xtensor_tex(head, args):
    idx = [a.strip() for a in args.split(',') if a.strip()]
    base = _XHEAD[head][0]
    if not idx:
        return base
    groups, cur_pos, cur = [], None, []
    for a in idx:
        pos = 'd' if a.startswith('-') else 'u'
        if pos != cur_pos and cur:
            groups.append((cur_pos, cur)); cur = []
        cur_pos = pos; cur.append(a.lstrip('-'))
    groups.append((cur_pos, cur))
    out = base
    for k, (pos, names) in enumerate(groups):
        out += ('' if k == 0 else '{}') + ('_{%s}' if pos == 'd' else '^{%s}') % ''.join(names)
    return out

def _strip_heads(s, head):
    """replace every head[...] (bracket-matched) by (...)"""
    out, i, key = '', 0, head + '['
    while True:
        j = s.find(key, i)
        if j < 0:
            return out + s[i:]
        out += s[i:j] + '('
        depth, k = 1, j + len(key)
        while depth:
            if s[k] == '[': depth += 1
            elif s[k] == ']': depth -= 1
            if depth: out += s[k]
            k += 1
        out += ')'
        i = k

def xact_to_latex(inputform):
    """LaTeX of an xAct expression in InputForm built from delta, RicciScalarCD, RicciCD, RiemannCD (abstract indices,
    a minus sign marks a lower index).  Products are ordered delta, R (scalar), Ricci, Riemann."""
    import sympy as sp
    names = {}
    def rep(m):
        tex = _xtensor_tex(m.group(1), m.group(2))
        for k, v in names.items():
            if v == tex and k.startswith(_XHEAD[m.group(1)][1]):
                return k
        key = '%s_%d' % (_XHEAD[m.group(1)][1], len(names))
        names[key] = tex
        return key
    s = _strip_heads(inputform, 'Scalar')
    s = _XTOK.sub(rep, s).replace('^', '**')
    syms = {k: sp.Symbol(k) for k in names}
    expr = sp.expand(sp.sympify(s, locals=syms))
    return sp.latex(expr, symbol_names={syms[k]: v for k, v in names.items()}), expr, names


# ----------------------------------------------------------------------------- exponential vacuum solutions (section 12.9)
def exponential_solutions():
    """For a = e^{H x0}, b = e^{qH H x0}: the three distinct vacuum equations (00, 11, 44) of
    0 = -Lambda delta + H^-2 w1 E1 + H^-4 w2 E2 + H^-6 w3 E3, the rank identity, the general solution for (w2, w3), and the
    case qH = -1.  Computed with sympy from the closed forms of tools/lovelock71.py; returns LaTeX strings and sympy objects."""
    import sympy as sp
    from lovelock71 import cosmological_closed_forms
    CF = cosmological_closed_forms()
    Ha, Hb, dHa, dHb = sp.symbols('H_a H_b dH_a dH_b', real=True)
    H, qH, Lam, w1, w2, w3 = sp.symbols('H q_H Lambda w_1 w_2 w_3')
    sub = {Ha: H, Hb: qH*H, dHa: 0, dHb: 0}
    eq = {c: sp.expand(-Lam + w1*CF[(1, c)].subs(sub)/H**2 + w2*CF[(2, c)].subs(sub)/H**4 + w3*CF[(3, c)].subs(sub)/H**6)
          for c in (0, 1, 4)}
    rank_identity = sp.expand(3*H*(eq[0] - eq[1]) + 4*qH*H*(eq[0] - eq[4])) == 0
    sol = sp.solve([eq[1], eq[4]], [w2, w3], dict=True)[0]
    sol = {k: sp.factor(v) for k, v in sol.items()}
    all_zero = all(sp.simplify(eq[c].subs(sol)) == 0 for c in (0, 1, 4))
    m1 = {c: sp.expand(eq[c].subs(qH, -1)) for c in (0, 1, 4)}
    s1 = {k: sp.expand(v.subs(qH, -1)) for k, v in sol.items()}
    return {'eq': eq, 'rank_identity': rank_identity, 'sol': sol, 'all_zero': all_zero, 'eq_qm1': m1, 'sol_qm1': s1,
            'tex_eq_qm1': {c: sp.latex(m1[c]) for c in (0, 1, 4)},
            'tex_sol': {str(k): sp.latex(v) for k, v in sol.items()},
            'tex_sol_qm1': {str(k): sp.latex(v) for k, v in s1.items()},
            'alpha_qm1': (sp.simplify(s1[w2].subs(Lam, 0)/w1), sp.simplify(s1[w3].subs(Lam, 0)/w1))}


# ----------------------------------------------------------------------------- numbers in printed values
_NUM = re.compile(r'-?\d+(?:\.\d*)?(?:[eE][-+]?\d+|\*\^-?\d+)?')

def nums_in(s):
    """all numbers in a printed value (a Mathematica output or a Python repr), as floats, in order of appearance;
    Mathematica precision marks (`) and numpy wrappers np.float64(...) are removed first"""
    s = re.sub(r'np\.(?:float|int|complex)\d*\(([^()]*)\)', r'\1', s or '')
    s = re.sub(r'`[\d.]*', '', s)
    s = re.sub(r'\b[A-Za-z_][A-Za-z_0-9]*\d', ' ', s)       # identifiers that end in a digit (np.float64, gg1, ...)
    return [float(x.replace('*^', 'e')) for x in _NUM.findall(s)]

def fierz_matrix():
    """the 9x9 Fierz matrix F_rs = (1/16) Sum_l (-1)^(rs-l) C(s,l) C(8-s,r-l) as Fractions, and G = (1 - F^T) g as
    LaTeX strings for G_0 and G_1"""
    from fractions import Fraction
    from math import comb
    F = [[Fraction(sum((-1)**(r*s - l)*comb(s, l)*comb(8 - s, r - l) for l in range(0, min(r, s) + 1) if r - l <= 8 - s), 16)
          for s in range(9)] for r in range(9)]
    def gtex(s):
        terms = []
        for r in range(9):
            c = (1 if r == s else 0) - F[r][s]
            if c == 0: continue
            sign = '-' if c < 0 else '+'
            a = abs(c)
            coef = '' if a == 1 else (r'\tfrac{%d}{%d}' % (a.numerator, a.denominator) if a.denominator != 1 else '%d' % a.numerator)
            terms.append('%s%sg_%d' % (sign, coef, r))
        t = ''.join(terms)
        return t[1:] if t.startswith('+') else t
    return F, {0: gtex(0), 1: gtex(1)}
