"""Data for the dirac16complex documents (DIRAC16COMPLEX*.md/.tex/.pdf): parsed from the evaluated 7+1.nb
(sections 7-10: checks, numbers, printed matrices) and recomputed independently by tools/ref_dirac16.py.
Nothing in the documents is typed by hand that can be computed."""
import os, re, sys, io, contextlib, platform
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import sympy as sp
import nbparse as P
from doc_common import gather as gather_common, sha256

NUM_RE = re.compile(r'"([^"]+)"->(-?[0-9]+\.?[0-9]*)(?:`[0-9.]*)?(?:\*\^(-?[0-9]+))?')

def notebook_parts(path=None):
    """outputs of the evaluated notebook with their section/subsection; checks of sections 7-9"""
    path = path or os.path.join(ROOT, '7+1.nb')
    nb = P.parse_notebook(path)
    rows = []; section = ''; subsection = ''
    for d, c in P.iter_cells(nb):
        st, content, opts = P.cell_info(c)
        if st == 'Section': section = P.box_to_text(content).strip('"'); subsection = ''
        elif st == 'Subsection': subsection = P.box_to_text(content).strip('"')
        elif st in ('Output', 'Message'):
            rows.append((st, section, subsection, P.box_to_text(content)))
    return rows

def numbers_from(text):
    return {m.group(1): sp.Float(m.group(2) + ('e' + m.group(3) if m.group(3) else ''), 25) for m in NUM_RE.finditer(text)}

def run_ref():
    """execute tools/ref_dirac16.py (without its __main__ block) and return its result dictionary"""
    ns = {}
    src = open(os.path.join(HERE, 'ref_dirac16.py'), encoding='utf-8').read().replace("if __name__ == '__main__':", 'if False:')
    cwd = os.getcwd(); os.chdir(HERE)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(src, 'ref_dirac16.py', 'exec'), ns)
    finally:
        os.chdir(cwd)
    return ns

def fmt(x, n=15):
    """a sympy Float as a decimal string with n significant digits, scientific notation as LaTeX"""
    s = sp.Float(x, 30)
    if s == 0: return '0'
    e = int(sp.floor(sp.log(abs(s), 10)))
    if -4 <= e <= 5:
        return sp.N(s, n).__format__('') if False else ('%.*g' % (n, float(s)) if n <= 15 else str(sp.N(s, n)))
    mant = s / sp.Integer(10)**e
    return r'%s\times10^{%d}' % (str(sp.N(mant, n)).rstrip('0').rstrip('.'), e)

def gather():
    D = gather_common()
    rows = notebook_parts()
    D['rows'] = rows
    outs = [(sec, sub, t) for st, sec, sub, t in rows if st == 'Output']
    checks = []
    for sec, sub, t in outs:
        m = re.fullmatch(r'"(.*)"->(True|False)', t.strip(), flags=re.S)
        if m and sec[:2] in ('7.', '8.', '9.'): checks.append((m.group(1), m.group(2), sec, sub))
    D['checks79'] = checks
    D['messages_all'] = sum(1 for r in rows if r[0] == 'Message')
    def find(sec_prefix, key):
        for sec, sub, t in outs:
            if sec.startswith(sec_prefix) and key in t: return numbers_from(t)
        return {}
    D['hs'] = find('8.', '"kF0"'); D['h1'] = find('8.', '"epsF"->1.98'); D['hf'] = find('8.', '"HF mstar"')
    D['cu'] = find('9.', '"kappa = rhoS/2"'); D['nj'] = find('9.', '"g_c"')
    # the printed matrix Sqrt[5] u_s(p), p = (2,2,4,0,0,0,0)
    um = [t for sec, sub, t in outs if sec.startswith('7.') and t.startswith('(\\[NoBreak]{{')]
    if um:
        rws = re.findall(r'\{([^{}]*)\}', um[0])
        D['U5'] = sp.Matrix([[sp.Integer(x) for x in r.split(',')] for r in rws])
    else:
        D['U5'] = None
    # closed forms printed in section 8.7
    cf = [t for sec, sub, t in outs if sec.startswith('8.') and 'ArcTanh' in t and t.startswith('{')]
    D['closed_forms_text'] = cf[0] if cf else ''
    # summary numbers of section 10
    tail = [t for st, sec, sub, t in rows if st == 'Output' and sec.startswith('10.')]
    D['summary10'] = tail[:3]
    # Python reference
    ref = run_ref()
    D['ref'] = ref['R']
    D['ref_bool'] = [(k, v) for k, v in ref['R'].items() if isinstance(v, bool) or type(v).__name__ == 'bool_']
    D['versions'] = dict(python=platform.python_version(), sympy=sp.__version__, mpmath=__import__('mpmath').__version__,
                         numpy=__import__('numpy').__version__)
    D['sha_backup'] = None
    bk = os.path.join(ROOT, 'backups', '7+1_2026-10-26_before_dirac16complex.nb')
    if os.path.exists(bk): D['sha_backup'] = sha256(bk); D['size_backup'] = os.path.getsize(bk)
    return D

if __name__ == '__main__':
    D = gather()
    print('checks 7-9:', len(D['checks79']), 'all True:', all(c[1] == 'True' for c in D['checks79']))
    print('messages:', D['messages_all'])
    print('hs', D['hs']); print('h1', D['h1']); print('hf', D['hf']); print('cu', D['cu']); print('nj', D['nj'])
    print('U5', None if D['U5'] is None else D['U5'].shape)
    print('closed forms', D['closed_forms_text'][:120])
    print('summary10', D['summary10'])
    print('ref bool checks', len(D['ref_bool']), all(v for k, v in D['ref_bool']))
    print('versions', D['versions'], 'backup sha', D['sha_backup'])
