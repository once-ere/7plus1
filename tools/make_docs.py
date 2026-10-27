#!/usr/bin/env python3
"""Build the documentation of 7+1.nb.

    python3 tools/make_docs.py

writes, in the repository root,
    README.md                      README.tex                      README.pdf
    README_PROVENANCE.md           README_PROVENANCE.tex           README_PROVENANCE.pdf
    README_SYMBOLS.md              README_SYMBOLS.tex              README_SYMBOLS.pdf
    DIRAC16COMPLEX.md              DIRAC16COMPLEX.tex              DIRAC16COMPLEX.pdf
    DIRAC16COMPLEX_PROVENANCE.md   DIRAC16COMPLEX_PROVENANCE.tex   DIRAC16COMPLEX_PROVENANCE.pdf
    DIRAC16COMPLEX_SYMBOLS.md      DIRAC16COMPLEX_SYMBOLS.tex      DIRAC16COMPLEX_SYMBOLS.pdf
The README documents cover the whole notebook; the DIRAC16COMPLEX documents cover its sections 7-10 (the Dirac equation
and the Kohn-Sham theories of the field dirac16complex).  The LaTeX files are compiled with latexmk -xelatex inside
build/docs (excluded by .gitignore); only the PDFs are copied back.  Every matrix, table, check list and number is
computed by doc_common.py / doc_dirac16_common.py or parsed from the evaluated 7+1.nb.
"""
import os, sys, subprocess, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True
from doc_dirac16_common import gather
import doc_readme, doc_symbols, doc_provenance
import doc_dirac16, doc_dirac16_symbols, doc_dirac16_provenance

def run(cmd):
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    return (p.stdout + p.stderr).strip(), p.returncode

def all_runs():
    """run the verification programs; their outputs go into the provenance documents"""
    out_v, rc_v = run([sys.executable, '-B', 'tools/verify_7p1_nb.py'])
    verify = ['tools/verify_7p1_nb.py (parses the evaluated 7+1.nb; exit code %d)' % rc_v, out_v, rc_v]
    runs = [verify]
    for f, what in [('ref44.py', 'rebuilds 4+4.nb'), ('ref71.py', 'rebuilds 7+1.nb objects'), ('extra71.py', 'further identities'),
                    ('ref_dirac16.py', 'independent checks of sections 7-9')]:
        out, rc = run([sys.executable, '-B', 'tools/' + f])
        runs.append(['tools/%s (%s; exit code %d)' % (f, what, rc), out, rc])
    return runs

def main(argv):
    """python3 tools/make_docs.py [--reuse-runs] [NAME ...]: --reuse-runs takes the program outputs saved by the previous run
    (build/docs/runs.json) instead of running the programs again; NAMEs restrict the build to some documents."""
    import json
    reuse = '--reuse-runs' in argv
    only = [a for a in argv if not a.startswith('--')]
    D = gather()
    bdir = os.path.join(ROOT, 'build', 'docs'); os.makedirs(bdir, exist_ok=True)
    rfile = os.path.join(bdir, 'runs.json')
    if reuse and os.path.exists(rfile):
        runs = json.load(open(rfile, encoding='utf-8'))
    else:
        runs = all_runs()
        json.dump(runs, open(rfile, 'w', encoding='utf-8'))
    rc_v, rc_r = runs[0][2], runs[-1][2]
    runs_readme = [(t, o) for t, o, rc in runs]
    runs16 = [(runs[0][0], runs[0][1]), (runs[-1][0], runs[-1][1])]
    builders = {'README': lambda: doc_readme.build(D),
                'README_SYMBOLS': lambda: doc_symbols.build(D),
                'README_PROVENANCE': lambda: doc_provenance.build(D, runs_readme),
                'DIRAC16COMPLEX': lambda: doc_dirac16.build(D),
                'DIRAC16COMPLEX_SYMBOLS': lambda: doc_dirac16_symbols.build(D),
                'DIRAC16COMPLEX_PROVENANCE': lambda: doc_dirac16_provenance.build(D, runs16)}
    docs = {k: f() for k, f in builders.items() if not only or k in only}
    report = []
    for name, d in docs.items():
        md, tex = d.to_md(), d.to_tex()
        with open(os.path.join(ROOT, name + '.md'), 'w', encoding='utf-8') as f: f.write(md)
        with open(os.path.join(ROOT, name + '.tex'), 'w', encoding='utf-8') as f: f.write(tex)
        shutil.copyfile(os.path.join(ROOT, name + '.tex'), os.path.join(bdir, name + '.tex'))
        p = subprocess.run(['latexmk', '-xelatex', '-interaction=nonstopmode', '-halt-on-error', name + '.tex'],
                           cwd=bdir, capture_output=True, text=True)
        pdf = os.path.join(bdir, name + '.pdf')
        ok = p.returncode == 0 and os.path.exists(pdf)
        if ok: shutil.copyfile(pdf, os.path.join(ROOT, name + '.pdf'))
        logf = os.path.join(bdir, name + '.log')
        log = open(logf, encoding='utf-8', errors='replace').read() if os.path.exists(logf) else ''
        report.append((name, len(md), len(tex), p.returncode, ok, log.count('Overfull \\hbox'), log.count('Underfull \\hbox'),
                       log.count('\n! ')))
    for r in report:
        print('%-26s md %6d chars  tex %6d chars  latexmk exit %d  pdf %s  overfull %d  underfull %d  errors %d' % r)
    print('verify_7p1_nb.py exit code %d, ref_dirac16.py exit code %d' % (rc_v, rc_r))
    return all(r[4] for r in report) and rc_v == 0 and rc_r == 0

if __name__ == '__main__':
    sys.exit(0 if main(sys.argv[1:]) else 1)
