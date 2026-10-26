#!/usr/bin/env python3
"""Build the documentation of 7+1.nb.

    python3 tools/make_docs.py

writes, in the repository root,
    README.md             README.tex             README.pdf
    README_PROVENANCE.md  README_PROVENANCE.tex  README_PROVENANCE.pdf
    README_SYMBOLS.md     README_SYMBOLS.tex     README_SYMBOLS.pdf
The LaTeX files are compiled with latexmk -xelatex inside build/docs (excluded by .gitignore); only the PDFs are
copied back.  Every matrix, table and check list is computed by doc_common.py or parsed from the evaluated 7+1.nb.
"""
import os, sys, subprocess, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
from doc_common import gather
import doc_readme, doc_symbols, doc_provenance

def run(cmd):
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    return (p.stdout + p.stderr).strip(), p.returncode

def main():
    D = gather()
    runs = []
    out, rc = run([sys.executable, 'tools/verify_7p1_nb.py'])
    runs.append(('tools/verify_7p1_nb.py (parses the evaluated 7+1.nb; exit code %d)' % rc, out))
    out, rc = run([sys.executable, 'tools/ref44.py'])
    runs.append(('tools/ref44.py (rebuilds 4+4.nb; exit code %d)' % rc, out))
    out, rc = run([sys.executable, 'tools/ref71.py'])
    runs.append(('tools/ref71.py (rebuilds 7+1.nb objects; exit code %d)' % rc, out))
    out, rc = run([sys.executable, 'tools/extra71.py'])
    runs.append(('tools/extra71.py (further identities; exit code %d)' % rc, out))
    docs = {'README': doc_readme.build(D),
            'README_SYMBOLS': doc_symbols.build(D),
            'README_PROVENANCE': doc_provenance.build(D, runs)}
    bdir = os.path.join(ROOT, 'build', 'docs'); os.makedirs(bdir, exist_ok=True)
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
        report.append((name, len(md), len(tex), p.returncode, ok))
    for r in report: print('%-18s md %6d chars  tex %6d chars  latexmk exit %d  pdf %s' % r)
    return all(r[4] for r in report)

if __name__ == '__main__':
    sys.exit(0 if main() else 1)
