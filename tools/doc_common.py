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

EVAL_LOG = os.path.join(ROOT, 'build', 'eval_log.txt')          # written by tools/run_eval.wl; build/ is excluded by .gitignore
EVAL_LOG_COPY = os.path.join(HERE, 'eval_log_7p1.txt')           # the copy that is part of the repository

def _parse_log(txt):
    d = {}
    for line in txt.splitlines():
        if ':' in line:
            k, v = line.split(':', 1); d[k.strip()] = v.strip()
    return d

def eval_log(nbinfo=None):
    """The record of the last evaluation of 7+1.nb (Wolfram version, kernel, duration, message and output cells, finish time).
    tools/run_eval.wl writes it to build/eval_log.txt, which .gitignore excludes; a copy, with the SHA-256 of the evaluated notebook
    added, is kept in tools/eval_log_7p1.txt, which is part of the repository, so that the documents can be rebuilt from a fresh
    clone.  A record is used only if it belongs to the notebook on disk: build/eval_log.txt must give the notebook's numbers of
    output and message cells (checked when nbinfo, the result of parse_evaluated_notebook, is given), the copy its SHA-256."""
    nb = os.path.join(ROOT, '7+1.nb')
    sha = sha256(nb) if os.path.exists(nb) else ''
    if os.path.exists(EVAL_LOG):
        txt = open(EVAL_LOG, encoding='utf-8', errors='replace').read()
        d = _parse_log(txt)
        if nbinfo is not None and (d.get('output cells') != str(nbinfo['outputs']) or d.get('message cells') != str(nbinfo['messages'])):
            print('warning: build/eval_log.txt (%s output cells, %s message cells) does not belong to 7+1.nb (%d, %d); not used'
                  % (d.get('output cells'), d.get('message cells'), nbinfo['outputs'], nbinfo['messages']))
            return {}
        copy = txt.rstrip('\n') + '\nnotebook SHA-256: ' + sha + '\n'
        old = open(EVAL_LOG_COPY, encoding='utf-8', errors='replace').read() if os.path.exists(EVAL_LOG_COPY) else None
        if old != copy:
            with open(EVAL_LOG_COPY, 'w', encoding='utf-8', newline='\n') as f: f.write(copy)
        d['notebook SHA-256'] = sha
        return d
    if os.path.exists(EVAL_LOG_COPY):
        d = _parse_log(open(EVAL_LOG_COPY, encoding='utf-8', errors='replace').read())
        if d.get('notebook SHA-256') == sha:
            return d
        print('warning: tools/eval_log_7p1.txt belongs to another 7+1.nb (SHA-256 %s); not used' % d.get('notebook SHA-256'))
    return {}

def all_checks(D):
    """Every check of 7+1.nb, in the order of the summary table of its last section: the cells check["name", test] that print
    their result where they stand (D['nb']['checks']), and the checks that the function frameChecks of section 13.1 makes for
    four metrics (nine each), which only that summary table prints (taken from the third-stage data D['s3'])."""
    printed = list(D['nb']['checks'])
    s3 = D.get('s3') or {}
    grid = s3.get('all_checks_grid') or []
    if not grid:
        return printed
    byname = {c[0]: c for c in printed}
    # the subsection in which frameChecks runs for each metric (its printed title where the subsection has printed checks)
    titles = {'13.3': '13.3 The pp-wave and the inhomogeneous metric: frames that are not diagonal and not homogeneous'}
    for c in printed:
        if c[3]: titles.setdefault(c[3].split(' ')[0], c[3])
    where = (('cosmological metric:', '13.2'), ('pp-wave:', '13.3'), ('inhomogeneous metric', '13.3'), ('static radial metric:', '13.5'))
    extra = {}
    for c in s3.get('checks', {}).get('grav', []):
        if c[0] in byname: continue
        num = next((w for p, w in where if c[0].startswith(p)), None)
        extra[c[0]] = (c[0], c[1], c[2], titles.get(num, c[3]) if num else c[3])
    return [byname.get(n) or extra.get(n) or (n, r, '?', '?') for n, r in grid]

# ------------------------------------------------------------------------------------------------------------------------
# The signature of 7+1 spacetime, written into EVERY document, symbol file and provenance file (the author's explanation of
# 27 October 2026): the xAct declaration {7,1,0} declares seven space-like dimensions and one time-like dimension (none
# degenerate), while the matrix diag(+1,-1,...,-1) is simply the squares of the octonion basis.
def sig_eq(eta=r'\eta'):
    """eta_AB = Re(e_A e_B), the squares of the octonion basis, and the quadratic form of the octonions (two centred lines)"""
    return (eta + r'_{AB}=\mathrm{Re}(e_Ae_B),\qquad e_0^2=+1,\qquad e_i^2=-1\ \ (i=1,\dots,7),\\ '
            r'\mathrm{Re}(x\,x)=x_0^2-x_1^2-\cdots-x_7^2=' + eta + r'_{AB}\,x^Ax^B\qquad (x=x^Ae_A\in\mathbb O).')

def signature_block(d, eta=r'\eta', eta71=r'\eta_{71}', gam=r'\Gamma', name17='η17', name71='η71', lead=None):
    """the paragraph, the displayed equation and the closing paragraph that explain {7,1,0} and diag(+1,-1,...,-1)"""
    lead = '**The signature: $\\{7,1,0\\}$ and $\\mathrm{diag}(+1,-1,\\dots,-1)$.** ' if lead is None else lead
    d.p(lead + 'The 7+1 spacetime has one time-like dimension ($x_0$, the direction of the octonion unit $e_0=1$) and seven '
        'space-like dimensions ($x_1,x_2,x_3$ ordinary space, $x_4,\\dots,x_7$ auxiliary space). Two different objects describe '
        'this, and neither may be replaced by the other. (i) The tensor package xAct, used in section 12 of 7+1.nb exactly as in '
        'the author’s notebook Einstein-Lovelock-7+1-Nash.nb, declares the metric with `DefMetric[{7, 1, 0}, g71[-a, -b], CD, …]`: '
        'the list $\\{7,1,0\\}$ declares seven space-like dimensions, one time-like dimension and no degenerate dimension — the '
        'numbers of entries $+1$ and $-1$ of $%s=\\mathrm{diag}(-1,+1,\\dots,+1)$ (`%s`), the same spacetime written in the sign '
        'convention of 4+4.nb (space-like $=+1$, time-like $=-1$). (ii) The matrix $%s=-%s=\\mathrm{diag}(+1,-1,\\dots,-1)$ (`%s`) '
        'is simply the squares of the octonion basis: $e_0^2=+1$ for the unit $e_0=1$ (time) and $e_i^2=-1$ for the seven '
        'imaginary units $e_1,\\dots,e_7$ (space); more generally it is the real part of the product of two basis units, and it '
        'gives the quadratic form of the octonions:' % (eta71, name71, eta, eta71, name17))
    d.eq(sig_eq(eta))
    d.p('$%s$ gives the signs of the Clifford relation $%s_A%s_B+%s_B%s_A=2\\,%s_{AB}\\,1_{16}$ of the real 16×16 Dirac '
        'matrices: time squares to $+1_{16}$ and each of the seven space directions to $-1_{16}$, exactly like the octonion '
        'units (in the notation $Cl(p,q)$, with $p$ generators squaring to $+1$ and $q$ to $-1$, the algebra '
        '$Cl(1,7)\\cong\\mathbb R(16)$). It also gives the signs of the metric components: $g_{00}>0$ along the time-like '
        'direction, negative squares along the space-like directions. Every component computed by 7+1.nb follows from the '
        'component matrix of the metric; the declared list fixes for xAct the sign of the determinant of the metric, '
        '`SignDetOfMetric[g71]` $=-1$, which is the same in both sign conventions, since $\\det %s=\\det %s=-1$. 7+1.nb checks '
        '`SignatureOfMetric[g71]` $=\\{7,1,0\\}$, `SignDetOfMetric[g71]` $=-1$, '
        '$\\mathrm{Re}(e_Ae_B)=%s_{AB}$ and $\\mathrm{Re}(x\\,x)=%s_{AB}x^Ax^B$.' % (eta, gam, gam, gam, gam, eta, eta, eta71, eta, eta))

def _sig_backup_text():
    p = os.path.join(ROOT, 'backups', '7+1_2026-10-27_before_signature.nb')
    return (', %d bytes, SHA-256 `%s`' % (os.path.getsize(p), sha256(p))) if os.path.exists(p) else ''

def signature_provenance(d, eta=r'\eta', eta71=r'\eta_{71}', gam=r'\Gamma', name17='η17', name71='η71'):
    """the provenance record of the author's explanation of 27 October 2026, followed by the signature block itself"""
    d.p('**Source.** On 27 October 2026 the author of the project explained that the xAct list $\\{7,1,0\\}$ declares seven '
        'space-like dimensions and one time-like dimension, while $\\mathrm{diag}(+1,-1,\\dots,-1)$ is simply the squares of the '
        'octonion basis, and asked that this be made clear in all documentation, examples and provenance records. '
        '**What it corrected.** The version of 7+1.nb pushed '
        'in the commits 52d8729 to dc2f178 of github.com/once-ere/7plus1 declared the xAct metric of section 12 with '
        '`DefMetric[{1, 7, 0}, …]`, reading the list as the signs of $%s$; since then section 12 declares '
        '`DefMetric[{7, 1, 0}, g71[-a, -b], CD, …]`, exactly as the author’s notebook Einstein-Lovelock-7+1-Nash.nb does. **What '
        'was added.** Checks in 7+1.nb — section 4.7: $%s$ is the squares of the octonion basis, $e_Ae_A=%s_{AA}e_0$ and '
        '$\\mathrm{Re}(e_Ae_B)=%s_{AB}$; section 12.2: `SignatureOfMetric[g71]` $=\\{7,1,0\\}$ with `SignDetOfMetric[g71]` $=-1$, '
        '$\\mathrm{Re}(e_Ae_B)=%s_{AB}$, and $\\mathrm{Re}(x\\,x)=%s_{AB}x^Ax^B$; section 17: the flat metric '
        '$\\mathrm{diag}(+1,-1,-1,-1)$ of the world brane is the squares of the basis $e_0,e_1,e_2,e_3$ of the quaternion '
        'subalgebra — and the same statements in the Python programs `tools/ref71.py`, `tools/ref_gravity.py` and '
        '`tools/ref_sm.py`; the explanation itself in the text of 7+1.nb (sections 1.2, 1.4, 4.2, 4.7, 7, 11, 12.2, 13, 15.1, 17) '
        'and in every document, symbol file and provenance file of the repository. **What did not change.** No result. Compared '
        'output by output with the backup below, every output of 7+1.nb is the same, except the declared list itself, the new checks '
        'and the check counts (and, as in every run, the name of the private context and the measured times); in particular every '
        'component of every Lovelock tensor is unchanged. The components follow from the component matrix of the metric, and the sign '
        'of the determinant, `SignDetOfMetric`, is $-1$ for both lists. **Backup.** The evaluated notebook as it was before is kept in '
        'the workspace as `backups/7+1_2026-10-27_before_signature.nb` (the folder backups/ is excluded from the repository by '
        '.gitignore)%s; `tools/verify_7p1_nb.py` compares every input, text and section cell of 7+1.nb with it and finds every cell '
        'present, each one unchanged or changed by this explanation, and it compares every output in the same way.'
        % (eta, eta, eta, eta, eta, eta, _sig_backup_text()))
    signature_block(d, eta=eta, eta71=eta71, gam=gam, name17=name17, name71=name71,
                    lead='**The explanation, in full.** ')

def gather():
    D = {}
    D['s4'] = s4; D['t4'] = t4; D['tau'] = tau; D['tauO'] = tauO; D['EO'] = EAaO; D['E44'] = EAa(8)
    D['oct_table'] = table_cells(CO, 'e', 0, '$\\cdot$')
    D['split_table'] = table_cells(CS, '\\varepsilon', 1, '$\\cdot$')
    D['xy'], D['a'], D['b'] = eight_square()
    D['nb'] = parse_evaluated_notebook()
    D['log'] = eval_log(D['nb'])
    D['sha'] = {name: sha256(os.path.join(ROOT, name)) for name in ['7+1.nb', '4+4.nb']}
    D['sha_tools'] = {name: sha256(os.path.join(HERE, name)) for name in sorted(os.listdir(HERE))
                      if os.path.isfile(os.path.join(HERE, name)) and (name.endswith('.py') or name.endswith('.wl'))}
    D['size44'] = os.path.getsize(os.path.join(ROOT, '4+4.nb'))
    D['size71'] = os.path.getsize(os.path.join(ROOT, '7+1.nb'))
    try:
        from zoneinfo import ZoneInfo
        D['today'] = datetime.datetime.now(ZoneInfo('America/Los_Angeles')).date().isoformat()   # the user's time zone
    except Exception:
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
